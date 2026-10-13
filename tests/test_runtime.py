import pandas as pd
import pytest
from varcode import Variant

from varlens.locus import Locus, chromosome_name_map
from varlens.mhc_binding import binding_affinities, CACHED_BINDING_AFFINITIES
from varlens.support import EXPECTED_COLUMNS, variant_support
from varlens.variant_includes import Includeable


def test_support_dataframe():
    variant = Variant("22", 46931060, "A", "C", genome="b37")
    counts = pd.DataFrame([
        ["tumor", "chr22", "46931059", "46931060", "A", 7],
        ["tumor", "22", "46931059", "46931060", "C", 2],
        ["tumor", "22", "46931059", "46931060", "C", 1],
        ["tumor", "22", "46931059", "46931060", "G", 2],
        ["normal", "22", "46931059", "46931060", "A", 10],
    ], columns=EXPECTED_COLUMNS)
    original = counts.copy(deep=True)
    result = variant_support(iter([variant]), counts)
    assert result.columns.names == ["metric", "source"]
    assert result.loc[variant, ("num_alt", "tumor")] == 3
    assert result.loc[variant, ("num_ref", "tumor")] == 7
    assert result.loc[variant, ("num_other", "tumor")] == 2
    assert result.loc[variant, ("total_depth", "tumor")] == 12
    assert result.loc[variant, ("alt_fraction", "tumor")] == pytest.approx(3 / 12)
    assert result.loc[variant, ("any_alt_fraction", "tumor")] == pytest.approx(5 / 12)
    assert result["num_ref"].loc[variant, "normal"] == 10
    pd.testing.assert_frame_equal(counts, original)

    missing = Variant("22", 46931062, "G", "A", genome="b37")
    with pytest.raises(ValueError, match="No allele counts"):
        variant_support([missing], counts)
    assert (variant_support([missing], counts, ignore_missing=True) == 0).all().all()
    assert variant_support([], counts).shape == (0, 12)
    assert variant_support([variant], counts.iloc[:0]).shape == (1, 0)


def test_chromosome_aliases():
    assert Locus.parse("chr22:42") == Locus("22", 41, 42)
    assert Locus.parse("chrM:42") == Locus("MT", 41, 42)
    assert Locus.from_inclusive_coordinates("chrUn_KI270438v1", 42).contig == "chrUn_KI270438v1"
    assert chromosome_name_map(["chr22", "chrM"])["22"] == "chr22"
    assert chromosome_name_map(["22", "chr22"])["22"] == "22"


def test_resume_empty_csv_annotation_column():
    class Annotation(Includeable):
        name = "test annotation"
        columns = ["annotation"]

        def process_chunk(self, frame):
            frame["annotation"] = "example"
            return frame

    frame = pd.DataFrame({"annotation": [float("nan"), float("nan")]})
    assert list(Annotation().compute(frame, chunk_rows=1)) == [1, 1]
    assert frame.annotation.tolist() == ["example", "example"]


def test_binding_affinity_lengths_and_missing_values(binding_predictor):
    variant = Variant("22", 46931060, "A", "C", genome="b37")
    alleles = ["A:02:01", "A:02:02"]
    first = binding_affinities(iter([variant]), iter(alleles), [9])
    assert first.iloc[0].binding_affinity == 29
    assert first.iloc[0].binding_allele == "A:02:02"
    calls = len(binding_predictor.calls)
    pd.testing.assert_frame_equal(first, binding_affinities([variant], alleles, [9]))
    assert len(binding_predictor.calls) == calls
    assert binding_affinities([variant], alleles, [10]).iloc[0].binding_affinity == 30
    CACHED_BINDING_AFFINITIES[(variant, alleles[0], (9,))] = float("nan")
    assert binding_affinities([variant], alleles, [9]).iloc[0].binding_affinity == 29
    assert list(binding_affinities([], alleles).columns) == [
        "variant", "binding_affinity", "binding_allele"]
    with pytest.raises(ValueError, match="At least one"):
        binding_affinities([variant], [])


def test_binding_affinities_do_not_mix_genomes(binding_predictor, monkeypatch):
    import topiary

    variants = [Variant("22", 46931060, "A", "C", genome=genome)
                for genome in ("b37", "GRCh38")]
    assert variants[0].short_description == variants[1].short_description

    def predict(self, collection):
        return pd.DataFrame([
            {"variant": variant.short_description, "kind": "pMHC_affinity",
             "value": 100 if variant.genome.reference_name == "GRCh37" else 200}
            for variant in collection
        ])

    monkeypatch.setattr(topiary.TopiaryPredictor, "predict_from_variants", predict)
    result = binding_affinities(variants, ["A:02:01"])
    assert result.binding_affinity.tolist() == [100, 200]
