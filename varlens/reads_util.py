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

import functools
import os

import typechecks
import pysam

from . import read_evidence
from . import evaluation

def add_args(parser):
    parser.add_argument("--reads", action="append", default=[])
    parser.add_argument("--read-filter", action="append", default=[])

def load_from_args(args):
    default_names = drop_prefix(args.reads)
    return [
        load_bam(url, name, args.read_filter)
        for (url, name) in zip(args.reads, default_names)
    ]

def drop_prefix(strings):
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
    return ReadSource(name if name else url, url_without_fragment, filters)

class ReadSource(object):
    def __init__(self, name, filename, read_filters=[]):
        self.name = name
        self.handle = pysam.Samfile(filename)
        self.read_filters = read_filters

    def pileups(self, loci):
        collection = read_evidence.PileupCollection.from_bam(self.handle, loci)
        if self.read_filters:
            for (locus, pileup) in collection.pileups.items():
                collection.pileups[locus] = pileup.filter([
                    functools.partial(
                        evaluate_pileup_element_expression,
                        expression,
                        collection,
                        pileup)
                    for expression in self.read_filters
                ])
        return collection

def evaluate_pileup_element_expression(
        expression,
        collection,
        pileup,
        element,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        bindings = eval.AttributeToKeyWrapper(
            [element, element.alignment, pileup],
            extra={
                'element': element,
                'pileup': pileup,
                'collection': collection,
            }
        )
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(pileup)   

