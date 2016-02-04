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

import logging
import functools

import pyensembl
import pysam

from . import read_evidence
from .read_source_helpers import alignment_key, evaluate_read_expression


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

    def index_if_needed(self):
        if self.filename.endswith(".bam") and not self.handle.has_index():
            # pysam strangely requires and index even to iterate through a bam.
            logging.info(
                "Attempting to create BAM index for file: %s" % self.filename)
            samtools_output = pysam.index(self.filename)
            logging.info(
                "Done indexing" + ((": " + samtools_output) if samtools_output else ''))

            # Reopen
            self.handle.close()
            self.handle = pysam.Samfile(self.filename)

    def reads(self, loci=None):
        if loci is None:
            def reads_iterator():
                return self.handle.fetch(until_eof=True)
        elif self.filename.endswith(".sam"):
            # Inefficient.
            chromosome_intervals = {}
            for (contig, intervals) in loci.contigs.items():
                try:
                    chromosome = self.chromosome_name_map[contig]
                except KeyError:
                    logging.warn(
                        "No such contig in bam: %s" % contig)
                    continue
                chromosome_intervals[chromosome] = intervals

            def reads_iterator():
                seen = set()
                for read in self.handle.fetch(until_eof=True):
                    intervals = chromosome_intervals.get(read.reference_name)
                    if not intervals or not intervals.overlaps_range(
                            read.reference_start,
                            read.reference_end):
                        continue
                    key = alignment_key(read)
                    if key not in seen:
                        yield read
                        seen.add(key)
        else:
            self.index_if_needed()

            def reads_iterator():
                seen = set()
                for locus in loci:
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
        self.index_if_needed()
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
