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

from __future__ import absolute_import

import subprocess
import warnings
import logging

import pandas
from nose.tools import eq_

from varlens.commands import variants

from . import data_path, run_and_parse_csv, cols_concat, temp_file

def run(args):
    logging.info("Running with args: " + ' '.join(args))
    return run_and_parse_csv(variants.run, args)

reference_fasta = data_path("chr22.no_line_wrap.fa")

expected_cols = [
    "genome", "contig", "interbase_start", "interbase_end", "ref", "alt",
]

def test_basic():
    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--genome", "b37",
    ])
    eq_(sorted(cols_concat(result, expected_cols)), sorted({
        "GRCh37-22-46931059-46931060-A-C",
        "GRCh37-22-21829554-21829555-T-G",
        "GRCh37-22-46931061-46931062-G-A",
        "GRCh37-22-50636217-50636218-A-C",
        "GRCh37-22-50875932-50875933-A-C",
    }))

def test_genes_and_effects():
    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--genome", "b37",
        "--include-effect",
        "--include-gene",
        "--rename-column", "gene", "genez",
    ])
    eq_(sorted(cols_concat(result, expected_cols + ["effect", "genez"])),
        sorted({
           'GRCh37-22-21829554-21829555-T-G-non-coding-transcript-PI4KAP2',
           'GRCh37-22-46931059-46931060-A-C-p.S670A-CELSR1',
           'GRCh37-22-46931061-46931062-G-A-p.S669F-CELSR1',
           'GRCh37-22-50636217-50636218-A-C-intronic-TRABD',
           'GRCh37-22-50875932-50875933-A-C-splice-acceptor-PPP6R2',
        }))

def test_context():
    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--genome", "b37",
        "--include-context",
        "--context-num-bases", "5",
        "--reference", reference_fasta,
    ])
    eq_(sorted(cols_concat(result,
            expected_cols + [
                "context_5_prime", "context_3_prime", "context_mutation"])),
        sorted({
            "GRCh37-22-46931059-46931060-A-C-GCTCC-CCACC-T>G",
            "GRCh37-22-21829554-21829555-T-G-CATGA-AGTGA-T>G",
            "GRCh37-22-46931061-46931062-G-A-GAGCT-CTCCA-C>T",
            "GRCh37-22-50636217-50636218-A-C-AGGGA-GGGCA-T>G",
            "GRCh37-22-50875932-50875933-A-C-AGGCC-GGGAG-T>G",
        }))

def test_mhc_binding_affinity():
    # If netMHC is not installed, we skip this test
    try:
        # If this succeeds (no exception), we do nothing.
        subprocess.call(
            "netMHC", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        warnings.warn("netMHC not installed, skipping mhc binding test")
        return

    with temp_file(".csv") as out_csv:
        run([
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--genome", "b37",
            "--include-mhc-binding",
            "--hla", "A:02:01 A:02:02",
            "--out", out_csv,
        ])
        eq_(sorted(cols_concat(pandas.read_csv(out_csv),
                expected_cols + ["binding_affinity", "binding_allele"])),
            sorted({
               'GRCh37-22-21829554-21829555-T-G-nan-nan',
               'GRCh37-22-46931059-46931060-A-C-142.13-A:02:02',
               'GRCh37-22-46931061-46931062-G-A-115.77-A:02:02',
               'GRCh37-22-50636217-50636218-A-C-nan-nan',
               'GRCh37-22-50875932-50875933-A-C-nan-nan',
            }))

def test_read_evidence():
    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--include-read-evidence",
        "--reads", data_path("CELSR1/bams/bam_0.bam"),
        "--genome", "b37",
    ])
    allele_groups = ["num_ref", "num_alt", "total_depth"]
    for allele_group in allele_groups:
        result[allele_group] = result[allele_group].astype(int)
    eq_(cols_concat(
            result,
            ["contig", "interbase_start"] + allele_groups),
        {
            '22-50636217-0-0-0',
            '22-50875932-0-0-0',
            '22-21829554-0-0-0',
            "22-46931059-50-0-50",
            "22-46931061-51-0-51",
    })

    # Same thing but with chunk rows = 1
    with temp_file(".csv") as out_csv:
        run([
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--include-read-evidence",
            "--reads", data_path("CELSR1/bams/bam_0.bam"),
            "--genome", "b37",
            "--chunk-rows", "1",
            "--out", out_csv,
        ])
        result = pandas.read_csv(out_csv)
    
        allele_groups = ["num_ref", "num_alt", "total_depth"]
        for allele_group in allele_groups:
            result[allele_group] = result[allele_group].astype(int)
        eq_(cols_concat(
                result,
                ["contig", "interbase_start"] + allele_groups),
            {
                '22-50636217-0-0-0',
                '22-50875932-0-0-0',
                '22-21829554-0-0-0',
                "22-46931059-50-0-50",
                "22-46931061-51-0-51",
        })

    result = run([
        "--include-read-evidence",
        "--reads", data_path("gatk_mini_bundle_extract.bam"),
        "--read-source-name", "foo",
        "--single-variant", "chr20:10008951", "C", "A",
        "--genome", "b37",
    ])
    for allele_group in allele_groups:
        result[allele_group] = result[allele_group].astype(int)
    eq_(cols_concat(result, expected_cols + allele_groups),
        {"GRCh37-20-10008950-10008951-C-A-4-1-5"})

    result = run([
        "--include-read-evidence",
        "--reads", data_path("gatk_mini_bundle_extract.bam"),
        "--read-source-name", "foo",
        "--single-variant", "chr20:10008951", "C", "A",
        "--genome", "b37",
        "--is-reverse",
    ])
    for allele_group in allele_groups:
        result[allele_group] = result[allele_group].astype(int)
    eq_(cols_concat(result, expected_cols + allele_groups),
        {"GRCh37-20-10008950-10008951-C-A-1-0-1"})


def test_filtering():
    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--genome", "b37",
        "--ref", "A",
    ])
    eq_(sorted(cols_concat(result, expected_cols)), sorted({
        "GRCh37-22-46931059-46931060-A-C",
        "GRCh37-22-50636217-50636218-A-C",
        "GRCh37-22-50875932-50875933-A-C",
    }))

    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--genome", "b37",
        "--ref", "A",
        "--variant-locus", "22:50636218",
    ])
    eq_(sorted(cols_concat(result, expected_cols)), sorted({
        "GRCh37-22-50636217-50636218-A-C",
    }))

    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        data_path("CELSR1/vcfs/vcf_2.vcf"),
        "--alt", "C", "G",
        "--genome", "b37"
    ])
    eq_(sorted(cols_concat(result, expected_cols)), sorted({
        "GRCh37-22-21829554-21829555-T-G",
        "GRCh37-22-45309892-45309893-T-G",
        "GRCh37-22-46931059-46931060-A-C",
        "GRCh37-22-50636217-50636218-A-C",
        "GRCh37-22-50875932-50875933-A-C",
    }))

'''
def test_fields():
    result = run([
        "--field",
        "foo:ref.lower()",
        "gene_names[0]",
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--variant-filter", "ref=='A'",
        "--variant-genome", "b37"
    ])
    eq_(sorted(cols_concat(result, expected_cols + ["foo", "gene_names[0]"])),
        sorted({
            "GRCh37-22-46931059-46931060-A-C-a-CELSR1",
            "GRCh37-22-50636217-50636218-A-C-a-TRABD",
            "GRCh37-22-50875932-50875933-A-C-a-PPP6R2",
        }))
'''
def test_round_trip():
    with temp_file(".csv") as out_csv:
        variants.run([
            data_path("CELSR1/vcfs/vcf_1.vcf"),
            "--out", out_csv,
            "--genome", "b37",
            "--ref", "A",
            "--include-gene",
        ])
        result1 = pandas.read_csv(out_csv)
        eq_(sorted(cols_concat(
                result1, expected_cols + ["gene"])),
            sorted({
                "GRCh37-22-46931059-46931060-A-C-CELSR1",
                "GRCh37-22-50636217-50636218-A-C-TRABD",
                "GRCh37-22-50875932-50875933-A-C-PPP6R2",
            }))

        result2 = run([
            out_csv,
            "--include-gene",
        ])
        eq_(sorted(cols_concat(
                result2,
                expected_cols + ["gene"])),
            sorted({
                "GRCh37-22-46931059-46931060-A-C-CELSR1",
                "GRCh37-22-50636217-50636218-A-C-TRABD",
                "GRCh37-22-50875932-50875933-A-C-PPP6R2",
            }))

def test_distinct_variants():
    result = run([
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        data_path("CELSR1/vcfs/vcf_1.vcf"),
        "--genome", "b37",
        "--ref", "A", "T",
        "--variant-source-name", "first", "second",
    ])
    eq_(sorted(cols_concat(result, expected_cols + ["sources"])),
        sorted({
            "GRCh37-22-21829554-21829555-T-G-first second",
            "GRCh37-22-46931059-46931060-A-C-first second",
            "GRCh37-22-50636217-50636218-A-C-first second",
            "GRCh37-22-50875932-50875933-A-C-first second",
        }))

