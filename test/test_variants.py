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

from pandas.util.testing import assert_frame_equal
from nose.tools import eq_, assert_raises

from varlens.commands import variants

from . import data_path, run_and_parse_csv, cols_concat

run = functools.partial(run_and_parse_csv, variants.run)

def test_basic():
    expected_cols = [
        "genome", "contig", "interbase_start", "interbase_end", "ref", "alt",
    ]

    result = run([
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
    ])
    eq_(sorted(cols_concat(result, expected_cols)), sorted({
        "GRCh37-22-46931059-46931060-A-C",
        "GRCh37-22-21829554-21829555-T-G",
        "GRCh37-22-46931061-46931062-G-A",
        "GRCh37-22-50636217-50636218-A-C",
        "GRCh37-22-50875932-50875933-A-C",
    }))

    # Loading the same file twice should have no effect as we're working
    # with sets of variants.
    result_twice = run([
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
    ])
    assert_frame_equal(result, result_twice)

def test_mixing_genomes():
    # Different genomes should raise an error.
    assert_raises(ValueError, run, [
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b37"),
        "--variants", data_path("CELSR1/vcfs/vcf_1.vcf#genome=b38"),
    ])
