"""Regression tests for the small bugs collected in issue #50."""

import logging
import re

import pytest
from pysam import Samfile

from varlens import util
from varlens.commands import variants
from varlens.locus import Locus
from varlens.read_evidence import PileupCollection

from . import data_path, run_and_parse_csv


def test_variants_limit_keeps_first_rows():
    path = data_path("CELSR1/vcfs/vcf_1.vcf")
    everything = run_and_parse_csv(variants.run, [path, "--genome", "b37"])
    assert len(everything) > 2

    limited = run_and_parse_csv(
        variants.run, [path, "--genome", "b37", "--limit", "2"])
    assert limited.to_dict("records") == everything.iloc[:2].to_dict("records")

    beyond_end = run_and_parse_csv(
        variants.run,
        [path, "--genome", "b37", "--limit", str(len(everything) + 10)])
    assert len(beyond_end) == len(everything)


def test_variants_limit_applies_before_annotation():
    result = run_and_parse_csv(
        variants.run,
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome", "b37",
            "--include-read-evidence",
            "--reads", data_path("CELSR1/bams/bam_0.bam"),
            "--limit", "1",
        ])
    assert len(result) == 1
    assert "total_depth" in result.columns


def test_variants_limit_rejects_negative(capsys):
    with pytest.raises(SystemExit):
        variants.run(
            [data_path("CELSR1/vcfs/vcf_1.vcf"), "--genome", "b37",
             "--limit", "-1"])
    assert "--limit must not be negative" in capsys.readouterr().err


@pytest.mark.parametrize("path, expected", [
    ("x.vcf", "x"),
    ("./x.vcf", "./x"),
    ("/data/run.1/x.vcf", "/data/run.1/x"),
    ("/data/run.1/x.vcf.gz", "/data/run.1/x"),
    ("a/x", "a/x"),
    (".hidden", ".hidden"),
])
def test_strip_extensions(path, expected):
    assert util.strip_extensions(path) == expected


@pytest.mark.parametrize("paths, expected", [
    (["./x.vcf"], ["x"]),
    (["/data/run.1/x.vcf"], ["x"]),
    (["a/x.vcf", "a/y.vcf"], ["x", "y"]),
    (["/data/run.1/x.vcf", "/data/run.1/y.vcf"], ["x", "y"]),
    (["./x.vcf", "./y.vcf"], ["x", "y"]),
    (["a/x.vcf.gz", "a/y.bam"], ["x", "y"]),
])
def test_drop_prefix_ignores_dots_in_directories(paths, expected):
    assert util.drop_prefix(paths) == expected


def test_drop_prefix_returns_input_on_collision():
    paths = ["a/x.vcf", "a/x.bam"]
    assert util.drop_prefix(paths) == paths


def test_locus_parse_error_names_the_input():
    with pytest.raises(ValueError) as info:
        Locus.parse("not a locus")
    message = str(info.value)
    assert "Couldn't parse locus: not a locus." in message
    assert "%s" not in message


def test_include_timing_log_is_formatted(caplog):
    caplog.set_level(logging.INFO)
    run_and_parse_csv(
        variants.run,
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome", "b37",
            "--include-read-evidence",
            "--reads", data_path("CELSR1/bams/bam_0.bam"),
        ])
    timing_messages = [
        record.getMessage() for record in caplog.records
        if record.getMessage().startswith("Processed in")
    ]
    assert timing_messages
    for message in timing_messages:
        assert re.fullmatch(r"Processed in \d+\.\d\d sec", message), message


def test_pileup_collection_filter_does_not_mutate_filters():
    locus = Locus.from_inclusive_coordinates("20", 10008951)
    evidence = PileupCollection.from_bam(
        Samfile(data_path("gatk_mini_bundle_extract.bam")), [locus])

    def keep_everything(element):
        return True

    filters = [keep_everything]
    first = evidence.filter(
        drop_duplicates=True,
        drop_improper_mate_pairs=True,
        min_mapping_quality=10,
        min_base_quality=10,
        filters=filters)
    assert filters == [keep_everything]

    # The same list can be reused: the earlier call's criteria do not carry
    # over to a call that only passes the user filters.
    assert (
        evidence.filter(filters=filters).allele_summary(locus)
        == evidence.allele_summary(locus))
    assert (
        first.allele_summary(locus)
        == evidence.filter(
            drop_duplicates=True,
            drop_improper_mate_pairs=True,
            min_mapping_quality=10,
            min_base_quality=10).allele_summary(locus))
