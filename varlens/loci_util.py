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
import itertools
import re

import intervaltree

from . import read_evidence
from . import variants_util
from .locus import Locus

def add_args(parser):
    # TODO:
    # - Load intervals_list files
    variants_util.add_args(parser)
    parser.add_argument('--locus', action="append", default=[],
        help="Genomic locus, like chr1:2342332 or chr1:2342-23423. "
        "Can be listed multiple times.")
    parser.add_argument("--neighbor-offsets",
        nargs="+", type=int, default=[])

def load_from_args(args):
    def loci():
        for locus in args.locus:
            match = re.match(r'(\w+):(\d+)(-(\d+))?', locus)
            if match is None:
                raise ValueError("Couldn't parse locus: %s. "
                    "Expected format is: chr5:3332 or chr5:3332-5555." % locus)
            
            (contig, start, _, maybe_end) = match.groups()
            start = int(start)
            end = int(maybe_end) if maybe_end is not None else start

            yield Locus.from_inclusive_coordinates(contig, start, end)

    def expand_with_neighbors(iterator):
        for locus in iterator:
            for offset in sorted(set(args.neighbor_offsets + [0])):
                if offset == 0:
                    yield locus
                else:
                    yield Locus(
                        locus.contig, locus.start + offset, locus.end + offset)

    variants = variants_util.load_from_args(args)
    return Loci(expand_with_neighbors(itertools.chain(
        loci(),
        (read_evidence.pileup_collection.to_locus(variant)
            for variant in variants))))

class Loci(object):
    def __init__(self, locus_iterator):
        self.contigs = collections.defaultdict(intervaltree.IntervalTree)
        for locus in locus_iterator:
            self.contigs[locus.contig].addi(locus.start, locus.end)
        
        for tree in self.contigs.values():
            tree.merge_overlaps()

    def __iter__(self):
        for contig in sorted(self.contigs):
            for interval in self.contigs[contig]:
                yield Locus(contig, interval.begin, interval.end)

    def __len__(self):
        return sum(len(tree) for tree in self.contigs.values())


