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

import functools

from nose.tools import eq_, assert_raises

from varlens.commands import allele_support

from . import data_path, run_and_parse_csv, cols_concat

run = functools.partial(run_and_parse_csv, allele_support.run)

def test_simple():
    result = run([
        "--reads", data_path("CELSR1/bams/bam_0.bam"),
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
    ])
    eq_(cols_concat(result, ["contig", "interbase_start", "interbase_end"]),
        {"22-46931059-46931060", "22-46931061-46931062"})

    pick_first_variant = [
        "ref=='A'",
        "alt=='C'",
        "inclusive_start==46931060",
        "interbase_start==46931059",
        "interbase_end==46931060",
#        "filter==",
    ]
    for variant_filter in pick_first_variant:
        result = run([
            "--reads", data_path("CELSR1/bams/bam_0.bam"),
            "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
            "--variant-filter", variant_filter,
        ])
        yield (
            eq_,
            cols_concat(
                result, ["contig", "interbase_start", "interbase_end"]),
            {"22-46931059-46931060"})

