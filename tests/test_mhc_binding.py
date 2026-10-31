from io import StringIO

import pandas
import pytest
import mhctools
from mhctools.base_predictor import BasePredictor
import topiary
from varcode import Variant

from varlens.mhc_binding import binding_affinities
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
