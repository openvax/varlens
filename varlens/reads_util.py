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

from __future__ import annotations

import argparse
import collections
from collections.abc import Callable
import functools
from typing import Any

from .commands.arguments import (
    READ_BOOLEAN_PROPERTIES,
    READ_INT_PROPERTIES,
    READ_STRING_PROPERTIES,
    add_read_args,
)
from .read_source import ReadSource
from . import util

BOOLEAN_PROPERTIES = list(READ_BOOLEAN_PROPERTIES)

STRING_PROPERTIES = list(READ_STRING_PROPERTIES)

INT_PROPERTIES = list(READ_INT_PROPERTIES)

# name -> (type, help, filter function)
READ_FILTERS: collections.OrderedDict[str, tuple[type[Any], str, Callable[..., bool]]] = (
    collections.OrderedDict()
)

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

def field_contains(field_name, parsed_value, read):
    field_value = getattr(read, field_name)
    return field_value is not None and parsed_value in field_value

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
        functools.partial(field_contains, prop))

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

def add_args(parser: argparse.ArgumentParser, positional: bool = False) -> None:
    """
    Extends a commandline argument parser with arguments for specifying
    read sources.

    When ``positional`` is True, both optional positional BAM paths and
    ``--reads`` are accepted and combined in order.
    """
    add_read_args(parser, positional=positional)


def read_paths_from_args(args: argparse.Namespace) -> list[str]:
    """Return all BAM/SAM paths given via positional ``reads`` and/or ``--reads``."""
    positional_reads = list(getattr(args, "reads", None) or [])
    optional_reads = list(getattr(args, "optional_reads", None) or [])
    return positional_reads + optional_reads


def load_from_args(args: argparse.Namespace) -> list[ReadSource] | None:
    """
    Given parsed commandline arguments, returns a list of ReadSource objects
    """
    read_paths = read_paths_from_args(args)
    if not read_paths:
        return None

    read_source_names = util.name_sources(
        read_paths, args.read_source_name, 'read_source_name', 'read source')

    filters = []
    for (name, info) in READ_FILTERS.items():
        value = getattr(args, name)
        if value is not None:
            filters.append(functools.partial(info[-1], value))

    return [
        load_bam(filename, name, filters)
        for (filename, name)
        in zip(read_paths, read_source_names)
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
