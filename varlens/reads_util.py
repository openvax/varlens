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

from .read_source import ReadSource
from . import util

def add_args(parser):
    """
    Extends a commandline argument parser with arguments for specifying
    a read source:
        --reads : One or more paths to SAM or BAM files
        --read-filter : Python expression for filtering reads
    """
    parser.add_argument("--reads", nargs="+", default=[])

    parser.add_argument(
        "--read-filter",
        action="append",
        default=[],
        nargs="+",
        help="Read filter expression, can be specified any number of times")

    parser.add_argument(
        "--read-source-name",
        nargs="+",
        help="Read source name")

def load_from_args(args):
    """
    Given parsed commandline arguments, returns a list of ReadSource objects
    """
    if not args.reads:
        return None

    if args.read_source_name:
        read_source_names = util.expand(
            args.read_source_name,
            'read_source_name',
            'read source',
            len(args.reads))
    else:
        read_source_names = util.drop_prefix(args.reads)

    read_filters = zip(*[
        util.expand(
            value, 'read_filter', 'read source', len(args.reads))
        for value in args.read_filter
    ])
    if not read_filters:
        read_filters = [[]] * len(args.reads)

    assert len(read_filters) == len(args.reads)

    return [
        load_bam(filename, name, read_filter)
        for (filename, name, read_filter)
        in zip(args.reads, read_source_names, read_filters)
    ]

def load_bam(filename, name=None, read_filters=[]):
    if not name:
        name = filename
    return ReadSource(name, filename, read_filters)

def flatten_header(header):
    for (group, rows) in header.items():
        for (index, row) in enumerate(rows):
            if not isinstance(row, dict):
                key_values = [(row, "")]
            else:
                key_values = row.items()
            for (key, value) in key_values:
                yield (str(group), index, str(key), str(value))
