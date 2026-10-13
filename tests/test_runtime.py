import pandas as pd
import pytest
from varcode import Variant

from varlens.locus import Locus, chromosome_name_map
from varlens.mhc_binding import binding_affinities, CACHED_BINDING_AFFINITIES
from varlens.support import EXPECTED_COLUMNS, variant_support
from varlens.variant_includes import Includeable


def test_insertion_cli_and_csv_round_trip(tmp_path):
    from . import data_path
    from varlens.commands import allele_support, variants
    from varlens.loci_util import Loci
    from varlens.read_source import ReadSource
    from varlens.read_evidence import PileupCollection
    from varlens.read_evidence.pileup_collection import to_locus
    from varlens.variants_util import load_as_dataframe

    bam = data_path("gatk_mini_bundle_extract.bam")
    variant = Variant("20", 10008921, "C", "CA", genome="GRCh37")
    locus = to_locus(variant)
    assert locus == Locus("20", 10008921, 10008921)
    evidence = PileupCollection.from_bam(bam, [locus])
    assert evidence.filter(drop_duplicates=True).allele_summary(locus) == [
        ("", 2), ("A", 1), ("AA", 1)]
    assert evidence.at(locus).num_reads() == evidence.num_reads()
    loci = Loci([locus]).union(Loci([Locus("20", 10008920, 10008921)]))
    assert len(loci) == 2
    assert loci.intersects(locus)
    assert Loci([Locus("20", 10008920, 10008921)]).intersects(locus)
    source = ReadSource("test", bam)
    try:
        assert list(source.reads(Loci([locus])))
    finally:
        source.handle.close()

    common = ["--single-variant", "20:10008921", "C", "CA", "--genome", "GRCh37",
              "--reads", bam, "--not-is-duplicate"]
    counts_path = tmp_path / "alleles.csv"
    allele_support.run(common + ["--out", str(counts_path)])
    counts = pd.read_csv(counts_path, keep_default_na=False)
    assert dict(zip(counts.allele, counts["count"])) == {"": 2, "A": 1, "AA": 1}
    result_path = tmp_path / "variants.csv"
    variants.run(common + ["--include-read-evidence", "--out", str(result_path)])
    frame = load_as_dataframe(str(result_path))
    assert frame.variant.tolist() == [variant]
    assert frame.num_alt.tolist() == [1]
    assert frame.num_ref.tolist() == [2]
    assert frame.total_depth.tolist() == [4]


def test_read_tags_and_sam_output(tmp_path):
    import pysam
    from . import data_path
    from varlens.commands import reads

    bam = data_path("gatk_mini_bundle_extract.bam")
    with pysam.AlignmentFile(bam) as handle:
        read = next(handle)
        read.set_tag("XY", "example")
        assert reads.read_field(None, read, "tag:XY") == "example"
        assert reads.read_field(None, read, "tag:XZ") is None
    output = tmp_path / "reads.sam"
    reads.run([bam, "--locus", "20:10008921", "--out", str(output)])
    assert output.read_text().startswith("@HD")
    with pysam.AlignmentFile(str(output), "r") as handle:
        assert list(handle)


def test_pileup_does_not_append_terminal_soft_clips(tmp_path):
    import pysam
    from varlens.read_evidence import PileupCollection

    path = tmp_path / "clipped.bam"
    with pysam.AlignmentFile(str(path), "wb", header={
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "1", "LN": 1000}]}) as handle:
        for name, sequence, cigar in [
                ("clipped", "ACGTTT", "3M3S"),
                ("inserted", "ACGATTT", "3M1I3S"),
                ("deleted", "ACTTT", "2M1D3S")]:
            read = pysam.AlignedSegment(handle.header)
            read.query_name = name
            read.query_sequence = sequence
            read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            read.reference_id = 0
            read.reference_start = 100
            read.mapping_quality = 60
            read.cigarstring = cigar
            handle.write(read)
    pysam.index(str(path))
    locus = Locus("1", 102, 103)
    evidence = PileupCollection.from_bam(str(path), [locus])
    assert dict(evidence.allele_summary(locus)) == {"G": 1, "GA": 1, "": 1}
    insertion = Locus("1", 103, 103)
    assert dict(evidence.allele_summary(insertion)) == {"A": 1}


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
