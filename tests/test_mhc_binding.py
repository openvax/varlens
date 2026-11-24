from io import StringIO

import pandas
import pytest
import mhctools
from mhctools.base_predictor import BasePredictor
import topiary
from varcode import Variant

from varlens.commands import variants as variants_cmd
from varlens.mhc_binding import (
    binding_affinities,
    fragment_predictions,
    mhc_predictor_from_args,
    parse_hla_alleles,
)
from varlens.variant_includes import MHCBindingAffinity


class FixedBindingPredictor(BasePredictor):
    """Exercise MHCTools 4's native results without an external executable."""

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        return [
            mhctools.PeptideResult(
                preds=tuple(
                    mhctools.Prediction(
                        peptide=peptide,
                        allele=allele,
                        kind=kind,
                        score=score,
                        value=value,
                        predictor_name="fixed-test",
                    )
                    for allele in self.alleles
                    for kind, value, score in [
                        ("pMHC_affinity", 50.25, 0.8),
                        ("pMHC_presentation", None, 0.01),
                    ]
                )
            )
            for peptide in peptides
        ]


def test_topiary_mhctools4_dataframe_contract():
    model = FixedBindingPredictor(
        alleles=["HLA-A*02:01"], default_peptide_lengths=[9]
    )
    predictions = topiary.TopiaryPredictor(models=model).predict_from_named_sequences(
        {"fixture": "ACDEFGHIKLMN"}
    )
    assert set(predictions["kind"]) == {"pMHC_affinity", "pMHC_presentation"}
    affinities = predictions.loc[predictions["kind"].eq("pMHC_affinity")]
    assert affinities["affinity"].tolist() == [50.25] * 4
    assert affinities["allele"].tolist() == ["HLA-A*02:01"] * 4
    assert predictions.loc[
        predictions["kind"].eq("pMHC_presentation"), "affinity"
    ].isna().all()


@pytest.fixture
def variants():
    return [
        Variant("22", 46931060, "A", "C", genome="GRCh37"),
        Variant("22", 46931062, "G", "A", genome="GRCh37"),
    ]


def install_predictions(monkeypatch, predictions):
    calls = []

    def model(**kwargs):
        calls.append(kwargs)
        return kwargs

    class Predictor:
        def __init__(self, models, only_novel_epitopes):
            assert only_novel_epitopes is True
            assert models == calls[-1]

        def predict_from_variants(self, variants):
            assert len(variants) > 0
            return predictions

    monkeypatch.setattr(mhctools, "NetMHCpan", model)
    monkeypatch.setattr(topiary, "TopiaryPredictor", Predictor)
    return calls


def test_affinity_summary_excludes_other_scores_and_missing_alleles(
    monkeypatch, variants
):
    first, second = variants
    predictions = pandas.DataFrame(
        [
            [first.short_description, "HLA-A*02:01", "pMHC_affinity", float("nan")],
            [first.short_description, "HLA-A*02:02", "pMHC_affinity", 200.0],
            [first.short_description, "HLA-A*02:02", "pMHC_affinity", 50.0],
            [first.short_description, "HLA-A*02:01", "pMHC_presentation", 0.01],
            [second.short_description, "HLA-A*02:01", "pMHC_affinity", float("nan")],
        ],
        columns=["variant", "allele", "kind", "affinity"],
    )
    calls = install_predictions(monkeypatch, predictions)
    result = binding_affinities(
        iter([second, first, first]), ["HLA-A*02:01", "HLA-A*02:02"]
    )
    assert result["variant"].tolist() == [second, first, first]
    assert pandas.isna(result.loc[0, "binding_affinity"])
    assert pandas.isna(result.loc[0, "binding_allele"])
    assert result["binding_affinity"].iloc[1:].tolist() == [50, 50]
    assert result["binding_allele"].iloc[1:].tolist() == ["HLA-A*02:02"] * 2
    binding_affinities(variants, ["HLA-A*02:02"], epitope_lengths=[9])
    assert calls[-1]["default_peptide_lengths"] == (9,)


@pytest.mark.parametrize(
    "predictions",
    [
        pandas.DataFrame(),
        pandas.DataFrame(columns=["variant", "allele", "kind", "affinity"]),
    ],
)
def test_empty_predictions(monkeypatch, variants, predictions):
    install_predictions(monkeypatch, predictions)
    result = binding_affinities(variants, ["HLA-A*02:01"])
    assert result["binding_affinity"].isna().all()
    assert result["binding_allele"].isna().all()
    assert list(binding_affinities([], ["HLA-A*02:01"]).columns) == list(result.columns)


def test_presentation_only_predictions_are_missing(monkeypatch, variants):
    predictions = pandas.DataFrame(
        [
            [variants[0].short_description, "HLA-A*02:01", "pMHC_presentation", 0.01],
        ],
        columns=["variant", "allele", "kind", "affinity"],
    )
    install_predictions(monkeypatch, predictions)
    result = binding_affinities(variants, ["HLA-A*02:01"])
    assert result["binding_affinity"].isna().all()
    assert result["binding_allele"].isna().all()


def test_reject_mixed_reference_genomes(variants):
    other = Variant("22", 46931060, "A", "C", genome="GRCh38")
    with pytest.raises(ValueError, match="reference genomes"):
        binding_affinities([variants[0], other], ["HLA-A*02:01"])


def test_chunked_mhc_annotation(monkeypatch, variants):
    predictions = pandas.DataFrame(
        [
            [variants[0].short_description, "HLA-A*02:01", "pMHC_affinity", 50.0],
        ],
        columns=["variant", "allele", "kind", "affinity"],
    )
    install_predictions(monkeypatch, predictions)
    frame = pandas.DataFrame({"variant": variants}, index=[4, 9])
    annotation = MHCBindingAffinity(hla=["HLA-A*02:01"])
    assert list(annotation.compute(frame, chunk_rows=1)) == [1, 1]
    assert frame.loc[4, "binding_allele"] == "HLA-A*02:01"
    assert frame.loc[4, "binding_affinity"] == 50
    assert pandas.isna(frame.loc[9, "binding_allele"])


@pytest.mark.parametrize("chunk_rows", [None, 1])
@pytest.mark.parametrize(
    "affinity", [50.25, float("nan")], ids=["fractional", "missing"]
)
def test_reannotate_integer_affinity_from_csv(
    monkeypatch, variants, chunk_rows, affinity
):
    frame = pandas.read_csv(
        StringIO("binding_affinity,binding_allele\n50,\n100,HLA-A*02:01\n")
    )
    assert frame["binding_affinity"].dtype == "int64"
    frame.index = [4, 9]
    frame["variant"] = variants
    predictions = pandas.DataFrame(
        [
            [variants[0].short_description, "HLA-A*02:01", "pMHC_affinity", affinity],
        ],
        columns=["variant", "allele", "kind", "affinity"],
    )
    install_predictions(monkeypatch, predictions)

    annotation = MHCBindingAffinity(hla=["HLA-A*02:01"])
    assert list(annotation.compute(frame, chunk_rows=chunk_rows)) == [1]
    if pandas.isna(affinity):
        assert pandas.isna(frame.loc[4, "binding_affinity"])
        assert pandas.isna(frame.loc[4, "binding_allele"])
    else:
        assert frame.loc[4, "binding_affinity"] == affinity
        assert frame.loc[4, "binding_allele"] == "HLA-A*02:01"
    # Fully annotated rows must retain their existing values.
    assert frame.loc[9, "binding_affinity"] == 100
    assert frame.loc[9, "binding_allele"] == "HLA-A*02:01"


@pytest.mark.ensembl
@pytest.mark.parametrize(
    "predictor_class", [mhctools.RandomBindingPredictor, FixedBindingPredictor]
)
def test_real_topiary_variant_pipeline(monkeypatch, variants, predictor_class):
    # Replace only the external executable; exercise current Topiary's real
    # annotation, peptide generation and DataFrame result contract.
    monkeypatch.setattr(mhctools, "NetMHCpan", predictor_class)
    result = binding_affinities(variants, ["HLA-A*02:01"])
    assert result["variant"].tolist() == variants
    assert result["binding_affinity"].notna().all()
    assert result["binding_allele"].tolist() == ["HLA-A*02:01"] * 2
    if predictor_class is FixedBindingPredictor:
        # The presentation score must not enter this affinity summary.
        assert result["binding_affinity"].tolist() == [50.25] * 2


@pytest.mark.ensembl
def test_intronic_primary_uses_coding_candidate(monkeypatch):
    from pyensembl import EnsemblRelease

    # Varcode 11 selects an intronic primary for this CFTR splice variant but
    # retains an exon-skip deletion; Topiary predicts from that alternative.
    variant = Variant("7", 117171171, "A", "C", EnsemblRelease(75))
    assert variant.effects().top_priority_effect().short_description == "intronic"
    monkeypatch.setattr(mhctools, "NetMHCpan", FixedBindingPredictor)
    result = binding_affinities([variant], ["HLA-A*02:01"])
    assert result["binding_affinity"].tolist() == [50.25]
    assert result["binding_allele"].tolist() == ["HLA-A*02:01"]


@pytest.mark.ensembl
def test_binding_affinities_accepts_predictor_and_mhc_predictor(variants):
    model = FixedBindingPredictor(
        alleles=["HLA-A*02:01"], default_peptide_lengths=[9]
    )
    from_model = binding_affinities(variants, predictor=model)
    assert from_model["binding_affinity"].tolist() == [50.25, 50.25]
    assert from_model["binding_allele"].tolist() == ["HLA-A*02:01", "HLA-A*02:01"]

    topiary_predictor = topiary.TopiaryPredictor(
        models=model, only_novel_epitopes=True
    )
    from_topiary = binding_affinities(variants, predictor=topiary_predictor)
    pandas.testing.assert_frame_equal(from_model, from_topiary)

    from_name = binding_affinities(
        variants, ["HLA-A*02:01"], epitope_lengths=[9], mhc_predictor="random"
    )
    assert from_name["binding_affinity"].notna().all()
    assert from_name["binding_allele"].tolist() == ["HLA-A*02:01", "HLA-A*02:01"]

    from_class = binding_affinities(
        variants,
        ["HLA-A*02:01"],
        epitope_lengths=[9],
        mhc_predictor=FixedBindingPredictor,
    )
    assert from_class["binding_affinity"].tolist() == [50.25, 50.25]


def test_binding_affinities_argument_validation(variants):
    model = FixedBindingPredictor(
        alleles=["HLA-A*02:01"], default_peptide_lengths=[9]
    )
    with pytest.raises(ValueError, match="not both"):
        binding_affinities(variants, ["HLA-A*02:01"], predictor=model)
    with pytest.raises(ValueError, match="not both"):
        binding_affinities(variants, predictor=model, mhc_predictor="random")
    with pytest.raises(ValueError, match="Unknown MHC predictor"):
        binding_affinities(variants, ["HLA-A*02:01"], mhc_predictor="no-such-model")
    with pytest.raises(TypeError, match="mhc_predictor"):
        binding_affinities(variants, ["HLA-A*02:01"], mhc_predictor=123)
    with pytest.raises(ValueError, match="positive integers"):
        binding_affinities(variants, ["HLA-A*02:01"], epitope_lengths=[True])
    with pytest.raises(ValueError, match="Specify at least one MHC allele"):
        binding_affinities(variants, [])


def test_fragment_predictions_accepts_mhc_predictor_and_mhctools_model():
    fragment = topiary.ProteinFragment(
        fragment_id="frag-1",
        source_type="sv:fusion",
        sequence="ACDEFGHIKLMN",
        target_intervals=[(4, 6)],
    )
    by_name = fragment_predictions(
        [fragment], ["HLA-A*02:01"], epitope_lengths=[9], mhc_predictor="random"
    )
    assert not by_name.empty
    assert set(by_name["allele"]) == {"HLA-A*02:01"}

    model = FixedBindingPredictor(
        alleles=["HLA-A*02:01"], default_peptide_lengths=[9]
    )
    by_model = fragment_predictions([fragment], predictor=model)
    assert not by_model.empty
    affinities = by_model.loc[by_model["kind"].eq("pMHC_affinity"), "affinity"]
    assert (affinities == 50.25).all()


def test_parse_hla_alleles_unifies_string_and_sequence_inputs():
    assert parse_hla_alleles("'A:02:01 A:02:02'") == ["A:02:01", "A:02:02"]
    assert parse_hla_alleles(["HLA-A*02:01 HLA-A*02:02"]) == [
        "HLA-A*02:01",
        "HLA-A*02:02",
    ]
    assert parse_hla_alleles(["HLA-A*02:01", "HLA-A*02:02", "HLA-A*02:01"]) == [
        "HLA-A*02:01",
        "HLA-A*02:02",
    ]
    assert parse_hla_alleles(["HLA-A*02:01,HLA-B*07:02"]) == [
        "HLA-A*02:01",
        "HLA-B*07:02",
    ]


def test_variants_mhc_cli_args_support_mhctools_and_deprecated_hla(tmp_path):
    alleles_file = tmp_path / "alleles.txt"
    alleles_file.write_text("HLA-A*02:01\nHLA-B*07:02\n")

    args = variants_cmd.parser.parse_args(
        [
            "--include-mhc-binding",
            "--mhc-predictor",
            "random",
            "--mhc-alleles",
            "HLA-A*02:01",
            "HLA-B*07:02",
            "--mhc-peptide-lengths",
            "8-9",
        ]
    )
    includeable = MHCBindingAffinity.from_args(args)
    assert isinstance(includeable.predictor, mhctools.RandomBindingPredictor)
    assert set(includeable.predictor.alleles) == {"HLA-A*02:01", "HLA-B*07:02"}
    assert includeable.predictor.default_peptide_lengths == [8, 9]

    file_args = variants_cmd.parser.parse_args(
        [
            "--include-mhc-binding",
            "--mhc-predictor",
            "random",
            "--mhc-alleles-file",
            str(alleles_file),
        ]
    )
    file_includeable = MHCBindingAffinity.from_args(file_args)
    assert set(file_includeable.predictor.alleles) == {"HLA-A*02:01", "HLA-B*07:02"}

    # Deprecated --hla accepts both quoted space-separated and multi-arg forms.
    for hla_argv in (
        ["--hla", "HLA-A*02:01 HLA-B*07:02"],
        ["--hla", "HLA-A*02:01", "HLA-B*07:02"],
    ):
        legacy_args = variants_cmd.parser.parse_args(
            ["--include-mhc-binding", "--mhc-predictor", "random"] + hla_argv
        )
        with pytest.warns(DeprecationWarning, match="--hla is deprecated"):
            predictor = mhc_predictor_from_args(legacy_args)
        assert isinstance(predictor, mhctools.RandomBindingPredictor)
        assert set(predictor.alleles) == {"HLA-A*02:01", "HLA-B*07:02"}

    # --mhc-alleles without --mhc-predictor requires an explicit predictor.
    missing_predictor = variants_cmd.parser.parse_args(
        ["--include-mhc-binding", "--mhc-alleles", "HLA-A*02:01"]
    )
    with pytest.raises(ValueError, match="--mhc-predictor"):
        MHCBindingAffinity.from_args(missing_predictor)

    # Combining --hla and --mhc-alleles is rejected.
    conflict_args = variants_cmd.parser.parse_args(
        [
            "--include-mhc-binding",
            "--mhc-predictor",
            "random",
            "--mhc-alleles",
            "HLA-A*02:01",
            "--hla",
            "HLA-B*07:02",
        ]
    )
    with pytest.raises(ValueError, match="Cannot combine deprecated --hla"):
        MHCBindingAffinity.from_args(conflict_args)


@pytest.mark.ensembl
def test_variants_mhc_cli_per_donor_hla_file_with_random_predictor(tmp_path, variants):
    hla_csv = tmp_path / "donors.csv"
    hla_csv.write_text("donor,hla\nD1,HLA-A*02:01\nD2,HLA-B*07:02\n")

    args = variants_cmd.parser.parse_args(
        [
            "--include-mhc-binding",
            "--mhc-predictor",
            "random",
            "--hla-file",
            str(hla_csv),
            "--mhc-peptide-lengths",
            "9",
        ]
    )
    includeable = MHCBindingAffinity.from_args(args)
    frame = pandas.DataFrame({"variant": variants, "donor": ["D1", "D2"]})
    list(includeable.compute(frame))
    assert frame["binding_allele"].tolist() == ["HLA-A*02:01", "HLA-B*07:02"]
    assert frame["binding_affinity"].notna().all()

    # --hla-file without --mhc-predictor is rejected.
    missing_pred = variants_cmd.parser.parse_args(
        ["--include-mhc-binding", "--hla-file", str(hla_csv)]
    )
    with pytest.raises(ValueError, match="--mhc-predictor"):
        MHCBindingAffinity.from_args(missing_pred)

    # --hla-file combined with --mhc-alleles is rejected.
    both_args = variants_cmd.parser.parse_args(
        [
            "--include-mhc-binding",
            "--mhc-predictor",
            "random",
            "--hla-file",
            str(hla_csv),
            "--mhc-alleles",
            "HLA-A*02:01",
        ]
    )
    with pytest.raises(ValueError, match="--hla-file"):
        MHCBindingAffinity.from_args(both_args)
