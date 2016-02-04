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

import os

from . import evaluation
from .read_source import ReadSource


def add_args(parser):
    """
    Extends a commandline argument parser with arguments for specifying
    a read source:
        --reads : One or more paths to SAM or BAM files
        --read-filter : Python expression for filtering reads
    """
    parser.add_argument("--reads", action="append", default=[])

    parser.add_argument(
        "--read-filter",
        action="append",
        default=[],
        help="Read filter expression, can be specified any number of times")


def load_from_args(args):
    """
    Given parsed commandline arguments, returns a list of ReadSource objects
    """
    if not args.reads:
        return None

    default_names = drop_prefix(args.reads)
    return [
        load_bam(url, name, args.read_filter)
        for (url, name) in zip(args.reads, default_names)
    ]


def drop_prefix(strings):
    """
    Removes common prefix from a collection of strings
    """
    if len(strings) == 1:
        return [os.path.basename(strings[0])]
    prefix_len = len(os.path.commonprefix(strings))
    return [string[prefix_len:] for string in strings]


def load_bam(url, name=None, read_filters=[]):
    (url_without_fragment, fragment) = evaluation.parse_url_fragment(url)
    filters = []
    for (key, value) in fragment:
        if key == 'filter':
            filters.append(value)
        elif key == "name":
            name = value
        else:
            raise ValueError("Unsupported operation: %s" % key)

    filters.extend(read_filters)
    if not name:
        name = url
    return ReadSource(name, url_without_fragment, filters)


def flatten_header(header):
    for (group, rows) in header.items():
        for (index, row) in enumerate(rows):
            if not isinstance(row, dict):
                key_values = [(row, "")]
            else:
                key_values = row.items()
            for (key, value) in key_values:
                yield (str(group), index, str(key), str(value))
