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
    """
    Return a Loci object giving the loci specified on the command line.

    If no loci-related arguments are specified, return None. This makes it
    possible to distinguish an empty set of loci, for example due to filters
    removing all loci, from the case where the user didn't specify any
    arguments.
    """
    variants = variants_util.load_from_args(args)

    if variants is None:
        if not args.locus:
            return None
        variant_loci = []
    else:
        variant_loci = (
            read_evidence.pileup_collection.to_locus(variant)
            for variant in variants)

    loci_iterator = itertools.chain(
        (Locus.parse(locus) for locus in args.locus),
        variant_loci)

    if args.neighbor_offsets:
        loci_iterator = expand_with_neighbors(
            loci_iterator, args.neighbor_offsets)

    return Loci(loci_iterator)

def expand_with_neighbors(loci_iterator, neighbor_offsets):
    offsets = sorted(set(neighbor_offsets + [0]))
    for locus in loci_iterator:
        for offset in offsets:
            if offset == 0:
                yield locus
            else:
                yield Locus(
                    locus.contig, locus.start + offset, locus.end + offset)

class Loci(object):
    def __init__(self, locus_iterator):
        self.contigs = collections.defaultdict(intervaltree.IntervalTree)
        for locus in locus_iterator:
            self.contigs[locus.contig].addi(locus.start, locus.end)

    def __iter__(self):
        for contig in sorted(self.contigs):
            for interval in self.contigs[contig]:
                yield Locus(contig, interval.begin, interval.end)

    def __len__(self):
        return sum(len(tree) for tree in self.contigs.values())


