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
import pyensembl
import logging

from . import read_evidence
from . import evaluation

def add_args(parser):
    parser.add_argument("--reads", action="append", default=[])
    parser.add_argument("--read-filter", action="append", default=[])

def load_from_args(args):
    if not args.reads:
        return None

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
        self.filename = filename
        self.handle = pysam.Samfile(filename)
        self.read_filters = read_filters

        self.chromosome_name_map = {}
        for name in self.handle.references:
            normalized = pyensembl.locus.normalize_chromosome(name)
            self.chromosome_name_map[normalized] = name
            self.chromosome_name_map[name] = name

    def reads(self, loci=None):
        if self.handle.is_bam and not self.handle._hasIndex():
            # pysam strangely requires and index even to iterate through a bam.
            logging.info("Attempting to create BAM index for file: %s" %
                self.filename)
            pysam.index(self.filename)

            # Reopen
            self.handle.close()
            self.handle = pysam.Samfile(self.filename)

        if loci is None:
            def reads_iterator():
                return self.handle.fetch()
        else:
            def reads_iterator():
                seen = set()
                for locus in loci:
                    logging.warn(locus)
                    try:
                        chromosome = self.chromosome_name_map[locus.contig]
                    except KeyError:
                        logging.warn(
                            "No such contig in bam: %s" % locus.contig)
                        continue
                    for read in self.handle.fetch(
                            chromosome,
                            locus.start,
                            locus.end):
                        key = alignment_key(read)
                        if key not in seen:
                            yield read
                            seen.add(key)

        for alignment in reads_iterator():
            if all(evaluate_read_expression(expression, alignment)
                    for expression in self.read_filters):
                yield alignment      

    def pileups(self, loci):
        collection = read_evidence.PileupCollection.from_bam(self.handle, loci)
        if self.read_filters:
            for (locus, pileup) in collection.pileups.items():
                def evaluate(expression, element):
                    return evaluate_read_expression(
                        expression, element.alignment)
                collection.pileups[locus] = pileup.filter([
                    functools.partial(evaluate, expression)
                    for expression in self.read_filters
                ])
        return collection

def evaluate_read_expression(
        expression,
        alignment,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        bindings = evaluation.EvaluationEnvironment(
            [alignment],
            extra={})
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(alignment) 


def evaluate_pileup_element_expression(
        expression,
        collection,
        pileup,
        element,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        bindings = evaluation.EvaluationEnvironment(
            [element, element.alignment, pileup],
            extra={
                'element': element,
                'pileup': pileup,
                'collection': collection,
            })
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(pileup)   

def alignment_key(pysam_alignment_record):
    '''
    Return the identifying attributes of a `pysam.AlignedSegment` instance.
    This is necessary since these objects do not support a useful notion of
    equality (they compare on identify by default).
    '''
    return (
        read_key(pysam_alignment_record),
        pysam_alignment_record.query_alignment_start,
        pysam_alignment_record.query_alignment_end,
    )

def read_key(pysam_alignment_record):
    '''
    Given a `pysam.AlignedSegment` instance, return the attributes identifying
    the *read* it comes from (not the alignment). There may be more than one
    alignment for a read, e.g. chimeric and secondary alignments.
    '''
    return (
        pysam_alignment_record.query_name,
        pysam_alignment_record.is_duplicate,
        pysam_alignment_record.is_read1,
        pysam_alignment_record.is_read2,
    )
