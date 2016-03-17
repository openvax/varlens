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

from nose.tools import eq_

from varlens.commands import reads

from . import data_path, run_and_parse_csv, cols_concat, temp_file

run = functools.partial(run_and_parse_csv, reads.run)

expected_cols = (
    "query_name,reference_start,reference_end,cigarstring").split(',')

def test_basic():
    result = run([
        data_path("CELSR1/bams/bam_0.bam"),
    ])
    eq_(result.shape, (953, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_0.bam"),
        "--is-duplicate",
    ])
    eq_(result.shape, (173, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_0.bam"),
        "--is-read1",
    ])
    eq_(result.shape, (481, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_0.bam"),
        "--is-read2",
    ])
    eq_(result.shape, (472, len(expected_cols)))

def test_loci_filtering():
    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
    ])
    eq_(result.shape, (37053, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--locus", "chr22:46930257-46930259"
    ])
    eq_(result.shape, (1795, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--locus", "chr22/46930256-46930259"
    ])
    eq_(result.shape, (1795, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--locus", "chr22:46930257-46930257"
    ])
    eq_(result.shape, (1753, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--locus", "chr22:46930257"
    ])
    eq_(result.shape, (1753, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--locus", "chr22/46930256"
    ])
    eq_(result.shape, (1753, len(expected_cols)))

def test_read_filtering():
    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--reference-start", '46932059',
    ])
    eq_(result.shape, (26, len(expected_cols)))

    result = run([
        data_path("CELSR1/bams/bam_5.bam"),
        "--reference-start", '46932059',
        "--query-name-contains", '57841',
    ])
    eq_(result.shape, (1, len(expected_cols)))

def test_round_trip():
    with temp_file(".bam") as out:
        reads.run([
            data_path("CELSR1/bams/bam_5.bam"),
            "--locus", "chr22/46930276",
            "--locus", "chr22/46930256",
            "--out", out,
        ])
        result1 = run([
            out,
        ])
        result2 = run([
            data_path("CELSR1/bams/bam_5.bam"),
            "--locus", "chr22/46930276",
            "--locus", "chr22/46930256",
        ])
        eq_(sorted(cols_concat(result1, expected_cols)),
            sorted(cols_concat(result2, expected_cols)))
