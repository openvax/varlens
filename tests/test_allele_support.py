# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging

import pytest
from varlens.commands import allele_support

from . import data_path, run_and_parse_csv, cols_concat


def run(args):
    logging.info("Running with args: " + " ".join(args))
    return run_and_parse_csv(allele_support.run, args)


expected_cols = [
    "contig",
    "interbase_start",
    "interbase_end",
    "allele",
    "count",
]


def test_basic():
    result = run(
        [
            "--reads",
            data_path("CELSR1/bams/bam_5.bam"),
            "--locus",
            "chr22:46929963",
            "chr22:46929964",
        ]
    )
    assert (cols_concat(result, expected_cols)) == (
        {"22-46929962-46929963-C-60", "22-46929963-46929964-A-81"}
    )

    result = run(
        [
            "--reads",
            data_path("CELSR1/bams/bam_5.bam"),
            "--locus",
            "chr22:46929963",
            "chr22:46929964",
            "--is-reverse",
        ]
    )
    assert (cols_concat(result, expected_cols)) == (
        {"22-46929962-46929963-C-37", "22-46929963-46929964-A-47"}
    )

    result = run(
        [
            "--reads",
            data_path("gatk_mini_bundle_extract.bam"),
            "--locus",
            "chr20:10008951",
            "--is-reverse",
        ]
    )
    assert (cols_concat(result, expected_cols)) == ({"20-10008950-10008951-C-1"})


def test_simple():
    result = run(
        [
            "--reads",
            data_path("CELSR1/bams/bam_0.bam"),
            "--genome",
            "b37",
            "--variants",
            data_path("CELSR1/vcfs/vcf_1.vcf"),
        ]
    )
    assert (
        cols_concat(
            result, ["contig", "interbase_start", "interbase_end", "allele", "count"]
        )
    ) == (
        {
            "22-50636217-50636218-N-0",
            "22-50875932-50875933-N-0",
            "22-21829554-21829555-N-0",
            "22-46931059-46931060-A-50",
            "22-46931061-46931062-G-51",
        }
    )


@pytest.mark.parametrize(
    "variant_filter",
    [
        ["--ref", "G"],
        ["--alt", "A"],
        ["--variant-locus", "22/46931061"],
        ["--variant-locus", "22/46931061-46931062"],
    ],
    ids=["reference", "alternate", "position", "interval"],
)
def test_variant_filter(variant_filter):
    result = run(
        [
            "--reads",
            data_path("CELSR1/bams/bam_0.bam"),
            "--genome",
            "b37",
            "--variants",
            data_path("CELSR1/vcfs/vcf_1.vcf"),
        ]
        + variant_filter
    )
    assert cols_concat(result, ["contig", "interbase_start", "interbase_end"]) == {
        "22-46931061-46931062"
    }


def test_positional_bam_and_reads_flag_consistency(capsys):
    bam_5 = data_path("CELSR1/bams/bam_5.bam")
    bam_0 = data_path("CELSR1/bams/bam_0.bam")

    positional_result = run([bam_5, "--locus", "chr22:46929963", "chr22:46929964"])
    flag_result = run(
        ["--reads", bam_5, "--locus", "chr22:46929963", "chr22:46929964"]
    )
    assert cols_concat(positional_result, expected_cols) == cols_concat(
        flag_result, expected_cols
    ) == {"22-46929962-46929963-C-60", "22-46929963-46929964-A-81"}

    # Multiple positional BAMs and positional + --reads combinations.
    multi_positional = run([bam_0, bam_5, "--locus", "chr22:46929963"])
    combined = run([bam_0, "--reads", bam_5, "--locus", "chr22:46929963"])
    repeated_flag = run(
        ["--reads", bam_0, "--reads", bam_5, "--locus", "chr22:46929963"]
    )
    assert {str(s) for s in multi_positional["source"]} == {"0", "5"}
    assert cols_concat(multi_positional, ["source"] + expected_cols) == cols_concat(
        combined, ["source"] + expected_cols
    ) == cols_concat(repeated_flag, ["source"] + expected_cols)

    # Positional BAM with --variants.
    with_variants = run(
        [
            bam_0,
            "--genome",
            "b37",
            "--variants",
            data_path("CELSR1/vcfs/vcf_1.vcf"),
        ]
    )
    assert cols_concat(with_variants, expected_cols) == {
        "22-50636217-50636218-N-0",
        "22-50875932-50875933-N-0",
        "22-21829554-21829555-N-0",
        "22-46931059-46931060-A-50",
        "22-46931061-46931062-G-51",
    }

    # Missing BAMs produces exit code 2.
    with pytest.raises(SystemExit) as error:
        allele_support.run(["--locus", "chr22:46929963"])
    assert error.value.code == 2
    assert "No read sources specified" in capsys.readouterr().err
