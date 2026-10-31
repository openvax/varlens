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

import shutil

import pytest
import logging

import pandas
import numpy
from varcode import Variant

from varlens.commands import variants

from . import data_path, run_and_parse_csv, cols_concat, temp_file


def run(args):
    logging.info("Running with args: " + " ".join(args))
    return run_and_parse_csv(variants.run, args)


expected_cols = [
    "genome",
    "contig",
    "interbase_start",
    "interbase_end",
    "ref",
    "alt",
]


def test_basic():
    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome",
            "b37",
        ]
    )
    assert (sorted(cols_concat(result, expected_cols))) == (
        sorted(
            {
                "GRCh37-22-46931059-46931060-A-C",
                "GRCh37-22-21829554-21829555-T-G",
                "GRCh37-22-46931061-46931062-G-A",
                "GRCh37-22-50636217-50636218-A-C",
                "GRCh37-22-50875932-50875933-A-C",
            }
        )
    )


@pytest.mark.ensembl
def test_genes_and_effects():
    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome",
            "b37",
            "--include-effect",
            "--include-gene",
            "--rename-column",
            "gene",
            "genez",
        ]
    )
    # Varcode owns the splice-effect description and can add outcomes between
    # releases. Verify that Varlens forwards its complete current description.
    splice_effect = Variant("22", 50875933, "A", "C", genome="GRCh37").effects()
    splice_description = splice_effect.top_priority_effect().short_description
    assert (sorted(cols_concat(result, expected_cols + ["effect", "genez"]))) == (
        sorted(
            {
                "GRCh37-22-21829554-21829555-T-G-non-coding-transcript-PI4KAP2",
                "GRCh37-22-46931059-46931060-A-C-p.S670A-CELSR1",
                "GRCh37-22-46931061-46931062-G-A-p.S669F-CELSR1",
                "GRCh37-22-50636217-50636218-A-C-intronic-TRABD",
                f"GRCh37-22-50875932-50875933-A-C-{splice_description}-PPP6R2",
            }
        )
    )


@pytest.mark.ensembl
def test_primary_intronic_annotation_retains_protein_changing_candidates():
    from pyensembl import EnsemblRelease
    from varlens.variant_includes import Effect

    variant = Variant("7", 117171171, "A", "C", EnsemblRelease(75))
    frame = Effect().process_chunk(pandas.DataFrame({"variant": [variant]}))
    assert frame.effect.iloc[0] == "intronic"
    assert "Intronic: intronic" in frame.effect_candidates.iloc[0]
    assert "Deletion:" in frame.effect_candidates.iloc[0]
    assert "Unresolved:" in frame.effect_candidates.iloc[0]


@pytest.mark.ensembl
def test_context(tmp_path):
    # pyfaidx writes an index beside the FASTA; keep fixtures read-only.
    reference_fasta = tmp_path / "chr22.fa"
    shutil.copyfile(data_path("chr22.no_line_wrap.fa"), reference_fasta)
    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome",
            "b37",
            "--include-context",
            "--context-num-bases",
            "5",
            "--reference",
            str(reference_fasta),
        ]
    )
    assert (
        sorted(
            cols_concat(
                result,
                expected_cols
                + ["context_5_prime", "context_3_prime", "context_mutation"],
            )
        )
    ) == (
        sorted(
            {
                "GRCh37-22-46931059-46931060-A-C-GCTCC-CCACC-T>G",
                "GRCh37-22-21829554-21829555-T-G-CATGA-AGTGA-T>G",
                "GRCh37-22-46931061-46931062-G-A-GAGCT-CTCCA-C>T",
                "GRCh37-22-50636217-50636218-A-C-AGGGA-GGGCA-T>G",
                "GRCh37-22-50875932-50875933-A-C-AGGCC-GGGAG-T>G",
            }
        )
    )


@pytest.mark.ensembl
@pytest.mark.netmhcpan
def test_mhc_binding_affinity():
    if not shutil.which("netMHCpan"):
        pytest.skip("NetMHCpan is not installed")

    with temp_file(".csv") as out_csv:
        run(
            [
                data_path("CELSR1/vcfs/vcf_1.vcf"),
                "--genome",
                "b37",
                "--include-mhc-binding",
                "--hla",
                "HLA-A*02:01 HLA-A*02:02",
                "--out",
                out_csv,
            ]
        )
        results = pandas.read_csv(out_csv).set_index(expected_cols)
        assert numpy.isnan(
            results.loc[("GRCh37", 22, 21829554, 21829555, "T", "G")].binding_affinity
        )
        assert numpy.isnan(
            results.loc[("GRCh37", 22, 21829554, 21829555, "T", "G")].binding_allele
        )
        for start, end, ref, alt in [
            (46931059, 46931060, "A", "C"),
            (46931061, 46931062, "G", "A"),
        ]:
            row = results.loc[("GRCh37", 22, start, end, ref, alt)]
            assert row.binding_allele in {"HLA-A*02:01", "HLA-A*02:02"}
            assert numpy.isfinite(row.binding_affinity)
            assert row.binding_affinity > 0


def test_read_evidence():
    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--include-read-evidence",
            "--reads",
            data_path("CELSR1/bams/bam_0.bam"),
            "--genome",
            "b37",
        ]
    )
    allele_groups = ["num_ref", "num_alt", "total_depth"]
    for allele_group in allele_groups:
        result[allele_group] = result[allele_group].astype(int)
    assert (cols_concat(result, ["contig", "interbase_start"] + allele_groups)) == (
        {
            "22-50636217-0-0-0",
            "22-50875932-0-0-0",
            "22-21829554-0-0-0",
            "22-46931059-50-0-50",
            "22-46931061-51-0-51",
        }
    )

    # Same thing but with chunk rows = 1
    with temp_file(".csv") as out_csv:
        run(
            [
                data_path("CELSR1/vcfs/vcf_1.vcf"),
                "--include-read-evidence",
                "--reads",
                data_path("CELSR1/bams/bam_0.bam"),
                "--genome",
                "b37",
                "--chunk-rows",
                "1",
                "--out",
                out_csv,
            ]
        )
        result = pandas.read_csv(out_csv)

        allele_groups = ["num_ref", "num_alt", "total_depth"]
        for allele_group in allele_groups:
            result[allele_group] = result[allele_group].astype(int)
        assert (cols_concat(result, ["contig", "interbase_start"] + allele_groups)) == (
            {
                "22-50636217-0-0-0",
                "22-50875932-0-0-0",
                "22-21829554-0-0-0",
                "22-46931059-50-0-50",
                "22-46931061-51-0-51",
            }
        )

    result = run(
        [
            "--include-read-evidence",
            "--reads",
            data_path("gatk_mini_bundle_extract.bam"),
            "--read-source-name",
            "foo",
            "--single-variant",
            "chr20:10008951",
            "C",
            "A",
            "--genome",
            "b37",
        ]
    )
    for allele_group in allele_groups:
        result[allele_group] = result[allele_group].astype(int)
    assert (cols_concat(result, expected_cols + allele_groups)) == (
        {"GRCh37-20-10008950-10008951-C-A-4-1-5"}
    )

    result = run(
        [
            "--include-read-evidence",
            "--reads",
            data_path("gatk_mini_bundle_extract.bam"),
            "--read-source-name",
            "foo",
            "--single-variant",
            "chr20:10008951",
            "C",
            "A",
            "--genome",
            "b37",
            "--is-reverse",
        ]
    )
    for allele_group in allele_groups:
        result[allele_group] = result[allele_group].astype(int)
    assert (cols_concat(result, expected_cols + allele_groups)) == (
        {"GRCh37-20-10008950-10008951-C-A-1-0-1"}
    )


def test_filtering():
    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome",
            "b37",
            "--ref",
            "A",
        ]
    )
    assert (sorted(cols_concat(result, expected_cols))) == (
        sorted(
            {
                "GRCh37-22-46931059-46931060-A-C",
                "GRCh37-22-50636217-50636218-A-C",
                "GRCh37-22-50875932-50875933-A-C",
            }
        )
    )

    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome",
            "b37",
            "--ref",
            "A",
            "--variant-locus",
            "22:50636218",
        ]
    )
    assert (sorted(cols_concat(result, expected_cols))) == (
        sorted(
            {
                "GRCh37-22-50636217-50636218-A-C",
            }
        )
    )

    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            data_path("CELSR1/vcfs/vcf_2.vcf"),
            "--alt",
            "C",
            "G",
            "--genome",
            "b37",
        ]
    )
    assert (sorted(cols_concat(result, expected_cols))) == (
        sorted(
            {
                "GRCh37-22-21829554-21829555-T-G",
                "GRCh37-22-45309892-45309893-T-G",
                "GRCh37-22-46931059-46931060-A-C",
                "GRCh37-22-50636217-50636218-A-C",
                "GRCh37-22-50875932-50875933-A-C",
            }
        )
    )


@pytest.mark.ensembl
def test_round_trip():
    with temp_file(".csv") as out_csv:
        variants.run(
            [
                data_path("CELSR1/vcfs/vcf_1.vcf"),
                "--out",
                out_csv,
                "--genome",
                "b37",
                "--ref",
                "A",
                "--include-gene",
            ]
        )
        result1 = pandas.read_csv(out_csv)
        assert (sorted(cols_concat(result1, expected_cols + ["gene"]))) == (
            sorted(
                {
                    "GRCh37-22-46931059-46931060-A-C-CELSR1",
                    "GRCh37-22-50636217-50636218-A-C-TRABD",
                    "GRCh37-22-50875932-50875933-A-C-PPP6R2",
                }
            )
        )

        result2 = run(
            [
                out_csv,
                "--include-gene",
            ]
        )
        assert (sorted(cols_concat(result2, expected_cols + ["gene"]))) == (
            sorted(
                {
                    "GRCh37-22-46931059-46931060-A-C-CELSR1",
                    "GRCh37-22-50636217-50636218-A-C-TRABD",
                    "GRCh37-22-50875932-50875933-A-C-PPP6R2",
                }
            )
        )


def test_distinct_variants():
    result = run(
        [
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome",
            "b37",
            "--ref",
            "A",
            "T",
            "--variant-source-name",
            "first",
            "second",
        ]
    )
    assert (sorted(cols_concat(result, expected_cols + ["sources"]))) == (
        sorted(
            {
                "GRCh37-22-21829554-21829555-T-G-first second",
                "GRCh37-22-46931059-46931060-A-C-first second",
                "GRCh37-22-50636217-50636218-A-C-first second",
                "GRCh37-22-50875932-50875933-A-C-first second",
            }
        )
    )
