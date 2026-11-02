Coordinate systems
==================

Which coordinates varlens uses in its output and which it expects on the
command line. Read this before comparing varlens positions with VCF or
samtools positions.

``varlens`` uses 0-based half-open coordinates internally. Many tools
(including samtools and VCF files) use inclusive 1-based coordinates. We try to
keep the confusion to a minimum by using the term "interbase" whenever we're
using 0-based half open coordinates and "inclusive" when we're using 1-based
inclusive coordinates.

One particularly sticky place this comes up is when specifying loci on the
commandline using e.g. ``--locus chr22:43243-43244``. To maintain consistency
with the most common other tools, when you specify a locus like
``chr22:10-20``, we interpret that as a 1-based inclusive coordinate. To
specify 0-based half-open coordinates, use a slash instead of a colon:
``chr22/9-20`` is the same 11 bases as ``chr22:10-20``.

See this `blog post <http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html>`_
for more details on coordinate systems.
