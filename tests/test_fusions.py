"""Real supplied Sid RNA: retain evidence and uncertainty through prediction."""

import gzip
from hashlib import sha256
import json
from pathlib import Path
import shutil

import pysam
import pytest

from varlens.commands import fusions, reads
from varlens.fusions import analyze_fusion
from . import data_path, run_and_parse_csv


ROOT = Path(data_path("fusions"))
CODING = ROOT / "ATP5MG--KMT2A.input.json.gz"
UNRESOLVED = ROOT / "GABBR1--SLC29A1-T1.json.gz"


def payload(filename):
    with gzip.open(filename, "rt") as stream:
        return json.load(stream)


@pytest.fixture
def random_predictor():
    pytest.importorskip("isovar.fusion")
    from mhctools import RandomBindingPredictor
    from topiary import TopiaryPredictor
    return TopiaryPredictor(models=RandomBindingPredictor(
        alleles=["HLA-A*02:01"], default_peptide_lengths=[8, 9]))


def test_fixture_checksums():
    for entry in json.loads((ROOT / "manifest.json").read_text()):
        packed = (ROOT / entry["file"]).read_bytes()
        assert sha256(packed).hexdigest() == entry["sha256"]
        assert sha256(gzip.decompress(packed)).hexdigest() == entry["uncompressed_sha256"]


def test_sid_coding_hypothesis_keeps_ambiguity_and_junction_evidence(random_predictor):
    supplied = payload(CODING)
    result = analyze_fusion(supplied, predictor=random_predictor, epitope_lengths=[8, 9])
    assert result["input"] == supplied
    reconstruction = result["reconstruction"]
    assert reconstruction["event_id"] == "ATP5MG--KMT2A"
    assert reconstruction["status"] == "ambiguous"
    assert reconstruction["reasons"]
    assert result["prediction_status"] == "completed"
    translation, = reconstruction["paths"][0]["translations"]
    fragment, = result["fragments"]
    assert fragment["sequence"] == translation["amino_acids"]
    assert fragment["n_rna_alt_reads"] == fragment["n_rna_alt_fragments"] == 2
    assert fragment["n_rna_ref_reads"] is None
    assert fragment["n_rna_alt_reads_supporting_protein_sequence"] is None
    assert {row["peptide"] for row in result["predictions"]} == {
        peptide["sequence"] for peptide in translation["candidate_peptides"]}
    for row in result["predictions"]:
        assert row["fusion_reconstruction_status"] == "ambiguous"
        assert row["fusion_provenance"] == supplied["fusion"]["provenance"]
        assert row["sample_name"] == "Sid-T0-Personalis"
        assert row["source_type"] == "sv:fusion"
        assert row["overlaps_target"] is True
        assert row["fusion_translation_observed"] is False
        assert row["wt_peptide"] is None
    json.dumps(result, allow_nan=False)


def test_sid_unresolved_result_preserves_reads_and_reasons(random_predictor):
    # Do not call a predictor at all when Isovar cannot establish a frame.
    def unexpected(_):
        pytest.fail("Unresolved RNA must not be sent to an MHC predictor")
    random_predictor.predict_from_fragments = unexpected
    supplied = payload(UNRESOLVED)
    result = analyze_fusion(supplied, predictor=random_predictor)
    reconstruction = result["reconstruction"]
    assert reconstruction["status"] == "unresolved"
    assert reconstruction["reasons"] == ["no_exact_collinear_annotated_donor"]
    assert reconstruction["paths"][0]["sequence"] == supplied["fusion"]["sequence"]
    assert reconstruction["paths"][0]["junctions"][0]["direct_support"]["fragments"] == 2
    assert result["fragments"] == result["predictions"] == []
    assert result["prediction_status"] == "no_coding_products"
    assert result["input"]["original_records"] == supplied["original_records"]


def test_missing_reads_are_insufficient_support_not_a_negative_protein(random_predictor):
    supplied = payload(CODING)
    supplied["reads"] = []
    result = analyze_fusion(supplied, predictor=random_predictor)
    assert result["reconstruction"]["status"] == "insufficient_support"
    assert result["reconstruction"]["reasons"] == ["insufficient_direct_junction_fragments"]
    assert result["predictions"] == []


def test_fusion_cli_retains_all_outcomes(tmp_path):
    pytest.importorskip("isovar.fusion")
    out = tmp_path / "report.json.gz"
    fusions.run([str(CODING), str(UNRESOLVED), "--out", str(out)])
    result = payload(out)
    assert result["schema"] == "varlens.fusions.v1"
    assert [r["reconstruction"]["status"] for r in result["events"]] == ["ambiguous", "unresolved"]
    assert all(r["prediction_status"] == "not_requested" for r in result["events"])


@pytest.mark.parametrize("bad_input", [[], {"fusion": {}}, {"unknown_key": True}])
def test_bad_fusion_input_is_actionable_and_does_not_write_partial_report(tmp_path, capsys, bad_input):
    pytest.importorskip("isovar.fusion")
    filename = tmp_path / "invalid.json"
    filename.write_text(json.dumps(bad_input))
    out = tmp_path / "result.json"
    out.write_text("previous result\n")
    with pytest.raises(SystemExit) as error:
        fusions.run([str(CODING), str(filename), "--out", str(out)])
    assert error.value.code == 2
    assert out.read_text() == "previous result\n"
    assert "invalid.json" in capsys.readouterr().err


def test_sid_raw_alignment_export_does_not_collapse_partner_records(tmp_path):
    supplied = payload(UNRESOLVED)
    header = pysam.AlignmentHeader.from_references(["chr6", "chr15"], [170805979, 101991189])
    unsorted = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as output:
        for record in supplied["original_records"]:
            for field in ("sam", "partner_sam"):
                output.write(pysam.AlignedSegment.fromstring(record[field], header))
    bam = tmp_path / "observations.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    result = run_and_parse_csv(reads.run, [str(bam), "--field", "is_supplementary", "mapping_quality"])
    assert len(result) == 4
    assert result.is_supplementary.sum() == 2


@pytest.mark.netmhcpan
def test_sid_fragment_predictions_with_real_netmhcpan():
    pytest.importorskip("isovar.fusion")
    if not shutil.which("netMHCpan"):
        pytest.skip("NetMHCpan is not installed")
    result = analyze_fusion(payload(CODING), alleles=["HLA-A*02:01"], epitope_lengths=[8, 9])
    assert result["prediction_status"] == "completed"
    assert result["predictions"]
    assert any(row["kind"] == "pMHC_affinity" and row["affinity"] > 0
               for row in result["predictions"])
    assert all(row["fusion_reconstruction_status"] == "ambiguous" for row in result["predictions"])
