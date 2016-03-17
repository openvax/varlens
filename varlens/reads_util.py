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

import collections
import functools

from .read_source import ReadSource
from . import util

BOOLEAN_PROPERTIES = """
is_paired is_proper_pair is_qcfail is_read1 is_read2 is_reverse is_secondary
is_unmapped mate_is_reverse mate_is_unmapped is_duplicate
""".split()

STRING_PROPERTIES = """
cigarstring query_alignment_sequence query_name
""".split()

INT_PROPERTIES = """
inferred_length mapping_quality query_alignment_length query_alignment_start
query_length reference_length reference_start template_length
""".split()

# name -> (type, help, filter function)
READ_FILTERS = collections.OrderedDict()

for prop in BOOLEAN_PROPERTIES:
    READ_FILTERS[prop] = (
        bool,
        "Only reads where %s is True" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                bool(getattr(read, field_name))),
            prop)
    )

    READ_FILTERS["not_" + prop] = (
        bool,
        "Only reads where %s is False" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                not getattr(read, field_name)),
            prop)
    )

for prop in STRING_PROPERTIES:
    READ_FILTERS["%s" % prop] = (
        str,
        "Only reads with the specified %s" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                getattr(read, field_name) == parsed_value),
            prop)
    )

    READ_FILTERS["%s_contains" % prop] = (
        str,
        "Only reads where %s contains the given string" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                parsed_value in getattr(read, field_name)),
            prop)
    )

for prop in INT_PROPERTIES:
    READ_FILTERS["%s" % prop] = (
        int,
        "Only reads with the specified %s" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                getattr(read, field_name) == parsed_value),
            prop)
    )

    READ_FILTERS["min_%s" % prop] = (
        int,
        "Only reads where %s >=N" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                getattr(read, field_name) >= parsed_value),
            prop)
    )

    READ_FILTERS["max_%s" % prop] = (
        int,
        "Only reads where %s <=N" % prop,
        functools.partial(
            (lambda field_name, parsed_value, read:
                getattr(read, field_name) <= parsed_value),
            prop)
    )

def add_args(parser, positional=False):
    """
    Extends a commandline argument parser with arguments for specifying
    read sources.
    """
    group = parser.add_argument_group("read loading")
    group.add_argument("reads" if positional else "--reads",
        nargs="+", default=[],
        help="Paths to bam files. Any number of paths may be specified.")

    group.add_argument(
        "--read-source-name",
        nargs="+",
        help="Names for each read source. The number of names specified "
        "must match the number of bam files. If not specified, filenames are "
        "used for names.")

    # Add filters
    group = parser.add_argument_group(
        "read filtering",
        "A number of read filters are available. See the pysam "
        "documentation (http://pysam.readthedocs.org/en/latest/api.html) "
        "for details on what these fields mean. When multiple filter "
        "options are specified, reads must match *all* filters.")

    for (name, (kind, message, function)) in READ_FILTERS.items():
        extra = {}
        if kind is bool:
            extra["action"] = "store_true"
            extra["default"] = None
        elif kind is int:
            extra["type"] = int
            extra["metavar"] = "N"
        elif kind is str:
            extra["metavar"] = "STRING"
        group.add_argument("--" + name.replace("_", "-"),
            help=message,
            **extra)

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

    filters = []
    for (name, info) in READ_FILTERS.items():
        value = getattr(args, name)
        if value is not None:
            filters.append(functools.partial(info[-1], value))

    return [
        load_bam(filename, name, filters)
        for (filename, name)
        in zip(args.reads, read_source_names)
    ]

def load_bam(filename, name=None, filters=[]):
    if not name:
        name = filename
    return ReadSource(name, filename, filters)

def flatten_header(header):
    for (group, rows) in header.items():
        for (index, row) in enumerate(rows):
            if not isinstance(row, dict):
                key_values = [(row, "")]
            else:
                key_values = row.items()
            for (key, value) in key_values:
                yield (str(group), index, str(key), str(value))
