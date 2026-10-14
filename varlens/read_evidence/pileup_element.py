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

from __future__ import absolute_import

from . import alignment_key

class PileupElement(object):
    '''
    A PileupElement represents the segment of an alignment that aligns to a
    particular base in the reference.

    Attributes
    ----------
    locus : Varcode.Locus
        The reference locus. Must be length 1, i.e. a single base.

    offset_start : int
        0-based start offset into the alignment sequence, inclusive

    offset_end : int
        0-based end offset into the alignment sequence, exclusive

    alignment : pysam.AlignedSegment
        pysam alignment instance

    alignment_key : tuple
        value computed from the alignment instance that uniquely specifies its
        properties. Used for comparisons since pysam.AlignedSegment instances
        do not support a useful notion of equality (they compare using object
        identity). See `read_evidence.alignment_key` for the implementation of
        this key.
    '''
    def __init__(self, locus, offset_start, offset_end, alignment):
        '''
        Construct a PileupElement object.
        '''
        assert offset_end >= offset_start, \
            "offset_start=%d > offset_end=%d" % (offset_start, offset_end)
        self.locus = locus
        self.offset_start = offset_start
        self.offset_end = offset_end
        self.alignment = alignment
        self.alignment_key = alignment_key(self.alignment)

    def fields(self):
        '''
        Fields that should be considered for our notion of object equality.
        '''
        return (
            self.locus, self.offset_start, self.offset_end, self.alignment_key)

    def __eq__(self, other):
        return hasattr(other, "fields") and self.fields() == other.fields()

    def __hash__(self):
        return hash(self.fields())

    @property
    def bases(self):
        '''
        The sequenced bases in the alignment that align to this locus in the
        genome, as a string.

        Empty string in the case of a deletion. String of length > 1 if there
        is an insertion here.
        '''
        sequence = self.alignment.query_sequence
        assert self.offset_end <= len(sequence), \
            "End offset=%d > sequence length=%d. CIGAR=%s. SEQUENCE=%s" % (
                self.offset_end,
                len(sequence),
                self.alignment.cigarstring,
                sequence)
        return sequence[self.offset_start:self.offset_end]

    @property
    def base_qualities(self):
        '''
        The phred-scaled base quality scores corresponding to ``self.bases``, as
        a list.
        '''
        return self.alignment.query_qualities[
            self.offset_start:self.offset_end]

    @property
    def min_base_quality(self):
        '''
        The minimum of the base qualities. In the case of a deletion, in which
        case there are no bases in this PileupElement, the minimum is taken
        over the sequenced bases immediately before and after the deletion.
        '''
        try:
            return min(self.base_qualities)
        except ValueError:
            # We are mid-deletion. We return the minimum of the adjacent bases.
            assert self.offset_start == self.offset_end
            adjacent_qualities = [
                self.alignment.query_qualities[offset]
                for offset in [self.offset_start - 1, self.offset_start]
                if 0 <= offset < len(self.alignment.query_qualities)
            ]
            return min(adjacent_qualities)

    @staticmethod
    def from_pysam_alignment(locus, pileup_read):
        '''
        Factory function to create a new PileupElement from a pysam
        ``PileupRead``.

        Parameters
        ----------
        locus : varcode.Locus
            Reference locus for which to construct a PileupElement. Must
            include exactly one base.

        pileup_read : pysam.calignmentfile.PileupRead
            pysam PileupRead instance. Its alignment must overlap the locus.

        Returns
        -------
        PileupElement

        '''
        assert not pileup_read.is_refskip, (
            "Can't create a PileupElement in a refskip (typically an intronic "
            "gap in an RNA alignment)")

        # PileupRead already resolves the CIGAR at this reference base. Walking
        # aligned_pairs again expands long introns and can incorrectly append
        # a terminal soft clip to the allele.
        if pileup_read.is_del:
            offset_start = offset_end = pileup_read.query_position_or_next
        else:
            offset_start = pileup_read.query_position
            offset_end = offset_start + 1 + max(0, pileup_read.indel)

        assert offset_end >= offset_start
        result = PileupElement(
            locus, offset_start, offset_end, pileup_read.alignment)
        return result   

