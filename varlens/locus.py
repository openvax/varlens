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

from collections import namedtuple
from collections.abc import Iterable
import re

from varcode.contigs import ContigResolver, annotation_contig_name


def normalize_contig(contig: str, reference_name: str | None = None) -> str:
    """Use Varcode's contig normalization (including UCSC accessions and GRCh37 chrM)."""
    return annotation_contig_name(contig, reference_name)


def chromosome_name_map(
    references: Iterable[str], reference_name: str | None = None
) -> dict[str, str]:
    """Resolve aliases to BAM references via Varcode's ContigResolver, giving exact names priority."""
    references = list(references)
    resolver = ContigResolver(references, reference_name=reference_name)
    result = {}
    for name in references:
        normalized = normalize_contig(name, reference_name=reference_name)
        resolution = resolver.resolve(normalized)
        if resolution.resolved_name is not None:
            result[normalized] = resolution.resolved_name
    result.update((name, name) for name in references)
    return result

class Locus(namedtuple("Locus", "contig start end")):
    '''
    A genomic interval in 0-indexed interbase coordinates.

    See this blog post for a discussion on coordinate systems:
        http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html
    '''

    @property
    def inclusive_start(self):
        return self.start + 1

    @property
    def inclusive_end(self):
        return self.end
    
    @property
    def positions(self):
        '''
        A Python range object giving the bases included in this locus.
        '''
        return range(self.start, self.end)

    @property
    def pileup_positions(self):
        """Reference bases needed to query this interval (left anchor for insertions)."""
        if self.start == self.end:
            if self.start < 1:
                raise ValueError("An insertion needs a reference base to its left")
            return range(self.start - 1, self.start)
        return self.positions

    @property
    def position(self):
        '''
        If this locus spans a single base, this property gives that position.
        Otherwise, raises a ValueError.
        '''
        if self.end != self.start + 1:
            raise ValueError("Not a single base: %s" % str(self))
        return self.start

    # Factory functions.
    @staticmethod
    def from_inclusive_coordinates(
        contig: str, start: int, end: int | None = None
    ) -> Locus:
        '''
        Given coordinates in 1-based coordinates that are inclusive on start
        and end, return a Locus instance. Locus instances are always 0-based
        "interbase" coordinates.
        '''
        if not isinstance(contig, str):
            raise TypeError(
                f"contig must be a str, got {type(contig).__name__}: {contig!r}"
            )
        if not isinstance(start, int):
            raise TypeError(
                f"start must be an int, got {type(start).__name__}: {start!r}"
            )
        if end is None:
            end = start
        if not isinstance(end, int):
            raise TypeError(
                f"end must be an int, got {type(end).__name__}: {end!r}"
            )
        contig = normalize_contig(contig)
        return Locus(contig, start - 1, end)

    @staticmethod
    def from_interbase_coordinates(
        contig: str, start: int, end: int | None = None
    ) -> Locus:
        '''
        Given coordinates in 0-based interbase coordinates, return a Locus
        instance.
        '''
        if not isinstance(contig, str):
            raise TypeError(
                f"contig must be a str, got {type(contig).__name__}: {contig!r}"
            )
        if not isinstance(start, int):
            raise TypeError(
                f"start must be an int, got {type(start).__name__}: {start!r}"
            )
        if end is None:
            end = start + 1
        if not isinstance(end, int):
            raise TypeError(
                f"end must be an int, got {type(end).__name__}: {end!r}"
            )
        contig = normalize_contig(contig)
        return Locus(contig, start, end)

    @staticmethod
    def parse(string: str) -> Locus:
        match = re.match(r'(\w+)([:/])(\d+)(-(\d+))?', string)
        if match is None:
            raise ValueError("Couldn't parse locus: %s. "
                "Expected format is: chr5:3332 or chr5:3332-5555 for "
                "inclusive 1-based coordinates and chr5/3331 or "
                "chr5/3331-5554 for half-open 0-based coordinates." % string)

        (contig, symbol, start, _, maybe_end) = match.groups()
        start = int(start)
        end = int(maybe_end) if maybe_end is not None else None

        if symbol == ":":
            # inclusive coordinatess
            return Locus.from_inclusive_coordinates(contig, start, end)
        else:
            # interbase coordinates
            assert symbol == "/"
            return Locus.from_interbase_coordinates(contig, start, end)
