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

import os
import tempfile

import pandas

from nose.tools import eq_, assert_raises

from varlens.commands import allele_support

from . import data_path

def run_command(args):
    assert '--out' not in args
    out_file = tempfile.NamedTemporaryFile(
        suffix=".csv",
        prefix="test_varlens_support",
        delete=False).name
    print("Running variant support with arguments: %s" % str(args))
    try:
        args.extend(['--out', out_file])
        allele_support.run(args)
        return pandas.read_csv(out_file)
    finally:
        os.unlink(out_file)

def cols_concat(df, columns, delimiter="-"):
    zipped = zip(*[df[c] for c in columns])
    return set([delimiter.join(str(item) for item in row) for row in zipped])

def test_simple():
    result = run_command([
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
        result = run_command([
            "--reads", data_path("CELSR1/bams/bam_0.bam"),
            "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
            "--variant-filter", variant_filter,
        ])
        yield (
            eq_,
            cols_concat(
                result, ["contig", "interbase_start", "interbase_end"]),
            {"22-46931059-46931060"})

