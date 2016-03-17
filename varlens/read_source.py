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

import pyensembl
import pysam

from . import read_evidence

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
                "Done indexing" + (
                    (": " + samtools_output) if samtools_output else ''))

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

        return (
            read for read in reads_iterator()
            if self.read_passes_filters(read))

    def read_passes_filters(self, read):
        return all(read_filter(read) for read_filter in self.read_filters)

    def pileups(self, loci):
        self.index_if_needed()
        collection = read_evidence.PileupCollection.from_bam(self.handle, loci)
        if self.read_filters:
            for (locus, pileup) in collection.pileups.items():
                collection.pileups[locus] = pileup.filter(
                    [lambda element:
                        self.read_passes_filters(element.alignment)])
        return collection

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


