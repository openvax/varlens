"""Input completeness and structural-event preservation, without reference data."""

import json

import pandas
import pytest
import varcode

from varlens import variants_util
from varlens.commands import variants, allele_support, reads
from varlens.errors import InputError
from varlens.mhc_binding import binding_affinities
from varlens.read_evidence.pileup_collection import to_locus
from . import data_path, run_and_parse_csv


HEADER = (
    '##fileformat=VCFv4.2\n##reference=GRCh38\n'
    '##contig=<ID=chr7,length=159345973>\n'
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Start uncertainty">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
)
SNV = 'chr7\t100\tsnv\tC\tT\t.\tPASS\t.\n'
DEL = 'chr7\t200\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300;CIPOS=-2,3\n'
BND = 'chr7\t400\tbnd\tN\tN[chr8:500[\t.\tPASS\tSVTYPE=BND\n'
SINGLE_BND = 'chr7\t600\tsingle\tN\tN.\t.\tPASS\tSVTYPE=BND\n'


def vcf(tmp_path, text, name="input.vcf"):
    filename = tmp_path / name
    filename.write_text(HEADER + text)
    return str(filename)


@pytest.mark.parametrize("rows", ["", SNV.replace("PASS", "LowQual")])
@pytest.mark.parametrize("chunk", [[], ["--chunk-rows", "1", "--include-effect"]])
def test_empty_output_has_headers(tmp_path, rows, chunk):
    result = run_and_parse_csv(variants.run, [vcf(tmp_path, rows)] + chunk)
    assert result.empty
    assert set(variants_util.STANDARD_DATAFRAME_COLUMNS + ["sources"]) <= set(result)


def test_empty_sources_and_empty_csv(tmp_path):
    empty = vcf(tmp_path, "", "empty.vcf")
    second = vcf(tmp_path, "", "second.vcf")
    result = run_and_parse_csv(variants.run, [empty, second])
    assert result.empty
    csv = tmp_path / "empty.csv"
    result.to_csv(csv, index=False)
    assert variants_util.load_as_dataframe(str(csv)).empty
    populated = vcf(tmp_path, SNV, "populated.vcf")
    combined = run_and_parse_csv(variants.run, [empty, populated])
    assert len(combined) == 1
    assert combined.iloc[0].sources == "populated"


@pytest.mark.parametrize("rows", [DEL + BND + SINGLE_BND, SNV + DEL + BND + SINGLE_BND])
def test_structural_csv_round_trip(tmp_path, rows):
    filename = vcf(tmp_path, rows)
    original = variants_util.load_as_dataframe(filename)
    result = run_and_parse_csv(variants.run, [filename])
    assert len(result) == len(rows.splitlines())
    csv = tmp_path / "roundtrip.csv"
    result.to_csv(csv, index=False)
    restored = variants_util.load_csv(csv)
    assert list(restored) == list(original.variant)
    for before, after in zip(original.variant, restored):
        if isinstance(before, varcode.StructuralVariant):
            fields = before.to_dict()
            fields["genome"] = before.reference_name
            assert after.to_dict() == fields
    structural = result[result.structural_variant.notna()]
    assert all(json.loads(value)[0]["source"] == filename for value in structural.vcf_records)
    assert [json.loads(value)[0]["metadata"]["id"] for value in structural.vcf_records] == [
        "del", "bnd", "single"]


def test_same_sv_in_two_sources_preserves_provenance(tmp_path):
    first = vcf(tmp_path, DEL, "first.vcf")
    second = vcf(tmp_path, DEL.replace("\tdel\t", "\tother_id\t"), "second.vcf")
    result = run_and_parse_csv(variants.run, [first, second])
    assert len(result) == 1
    assert result.iloc[0].sources == "first second"
    records = json.loads(result.iloc[0].vcf_records)
    assert [r["metadata"]["id"] for r in records] == ["del", "other_id"]


@pytest.mark.parametrize("limit,ids", [(None, ["del", "second"]), (1, ["del"])])
def test_duplicate_sv_records_in_one_vcf_keep_each_id(tmp_path, limit, ids):
    filename = vcf(tmp_path, DEL + DEL.replace("\tdel\t", "\tsecond\t"))
    df = variants_util.load_as_dataframe(filename, max_variants=limit)
    assert len(df) == 1
    assert [r["metadata"]["id"] for r in json.loads(df.vcf_records.iloc[0])] == ids
    csv = tmp_path / "records.csv"
    df.drop(columns="variant").to_csv(csv, index=False)
    assert variants_util.load_as_dataframe(csv).vcf_records.tolist() == df.vcf_records.tolist()


def test_sv_records_with_different_uncertainty_stay_separate(tmp_path):
    filename = vcf(tmp_path, DEL + DEL.replace("CIPOS=-2,3", "CIPOS=-5,6").replace("\tdel\t", "\tsecond\t"))
    df = variants_util.load_as_dataframe(filename)
    assert len(df) == 2
    assert sorted((row.variant.ci_start, json.loads(row.vcf_records)[0]["metadata"]["id"])
                  for _, row in df.iterrows()) == [((-5, 6), "second"), ((-2, 3), "del")]


def test_breakpoint_filter_includes_remote_mate_not_deletion_interior(tmp_path):
    filename = vcf(tmp_path, SNV + DEL + BND)
    mate = run_and_parse_csv(variants.run, [filename, "--variant-locus", "chr8:500"])
    assert mate.alt.tolist() == ["N[chr8:500["]
    interior = run_and_parse_csv(variants.run, [filename, "--variant-locus", "chr7:250"])
    assert interior.empty


@pytest.mark.parametrize("prefix", ["", SNV])
def test_unrecognized_structural_allele_fails_before_output(tmp_path, capsys, prefix):
    filename = vcf(tmp_path, prefix + BND.replace("N[chr8:500[", "N[chr8:500"))
    out = tmp_path / "out.csv"
    with pytest.raises(SystemExit) as error:
        variants.run([filename, "--out", str(out)])
    assert error.value.code == 2
    assert not out.exists()
    assert "No partial output" in capsys.readouterr().err


@pytest.mark.parametrize("flag", ["--include-effect", "--include-gene",
                                 "--include-context", "--include-read-evidence",
                                 "--include-mhc-binding"])
def test_unsupported_sv_annotation_is_explicit(tmp_path, capsys, flag):
    with pytest.raises(SystemExit) as error:
        variants.run([vcf(tmp_path, DEL), flag])
    assert error.value.code == 2
    assert "does not support structural variants" in capsys.readouterr().err


def test_structural_read_selection_and_python_apis_reject_interval_counts(tmp_path, capsys):
    filename = vcf(tmp_path, BND)
    bam = data_path("CELSR1/bams/bam_0.bam")
    for command, arguments in ((allele_support.run, ["--reads", bam]),
                               (reads.run, [bam])):
        with pytest.raises(SystemExit) as error:
            command(arguments + ["--variants", filename])
        assert error.value.code == 2
        assert "junction support" in capsys.readouterr().err
    variant = variants_util.load_as_dataframe(filename).variant.iloc[0]
    with pytest.raises(InputError, match="junction support"):
        to_locus(variant)
    with pytest.raises(InputError, match="structural variants"):
        binding_affinities([variant], ["HLA-A*02:01"])


def test_structural_csv_cannot_lose_or_contradict_payload(tmp_path):
    df = variants_util.load_as_dataframe(vcf(tmp_path, DEL))
    csv = tmp_path / "bad.csv"
    df.drop(columns=["variant", "structural_variant"]).to_csv(csv, index=False)
    with pytest.raises(InputError, match="require the structural_variant"):
        variants_util.load_as_dataframe(str(csv))
    df.loc[0, "interbase_end"] = 301
    df.drop(columns="variant").to_csv(csv, index=False)
    with pytest.raises(InputError, match="disagree"):
        variants_util.load_as_dataframe(str(csv))


def test_unsupported_format_has_actionable_diagnostic(tmp_path, capsys):
    filename = tmp_path / "fusion.json.gz"
    filename.write_bytes(b"unused")
    with pytest.raises(SystemExit) as error:
        variants.run([str(filename)])
    assert error.value.code == 2
    assert "varlens-fusions" in capsys.readouterr().err


def test_chunked_unannotated_output_and_column_selection(tmp_path):
    filename = vcf(tmp_path, SNV)
    selected = run_and_parse_csv(variants.run, [filename, "--chunk-rows", "1", "--columns", "alt"])
    pandas.testing.assert_frame_equal(selected, pandas.DataFrame({"alt": ["T"]}))


@pytest.mark.parametrize("sv_type", ["DEL", "DUP", "INV", "INS", "CNV", "BND"])
def test_python_dataframe_csv_preserves_explicit_sv_fields(tmp_path, sv_type):
    variant = varcode.StructuralVariant(
        "chr7", 200, sv_type, end=300, affected_start=201, affected_end=300,
        ci_start=(-2, 3), ci_end=(-4, 5), alt_assembly="ACGT",
        info={"caller": "test", "support": [2, 3]}, genome="GRCh38",
        **({"alt": "N[chr8:500["} if sv_type == "BND" else {}))
    df = variants_util.variants_to_dataframe([variant])
    filename = tmp_path / "python.csv"
    # Deliberately include the internal `variant` column to verify it is
    # reconstructed from the typed data, never trusted as a string object.
    df.to_csv(filename, index=False)
    restored = variants_util.load_as_dataframe(filename).variant.iloc[0]
    assert restored == variant
    assert restored.to_dict() == variant.to_dict()
