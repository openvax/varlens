.. image:: https://travis-ci.org/hammerlab/varlens.svg?branch=master
    :target: https://travis-ci.org/hammerlab/varlens

varlens
======================

A collection of Python tools for working with genomic variants and
next-generation sequencing reads. Not particularly fast for large datasets. The
emphasis is on exracting what you need from BAMs and VCFs into a CSV file for
further analysis.

Built on `varcode <https://github.com/hammerlab/varcode>`_ and `pysam <https://github.com/pysam-developers/pysam>`_.

varlens-variants
    Combine, annotate, and filter variants from VCF or CSV files. Available
    annotations include genes, variant effects, surrounding sequence context,
    counts of supporting reads from specified BAM files, and MHC I binding
    affinity prediction of mutant peptides.

varlens-reads
    Display, filter, and copy reads from a SAM/BAM file. Partial replacement for ``samtools view``.

varlens-allele-support
    Count reads supporting each allele at specified sites in BAM files.


Installation
-------------

From a git checkout:

::

    pip install .

To run the tests:

::

    nosetests .

To build the documentation (just this README plus the commandline tool help):

::

    pip install -e .
    pip install Sphinx
    cd docs
    make clean setup rst html

The docs will be written to the ``_build/html`` directory.


varlens-variants
----------------------

Given variants from one or more VCF or CSV files, apply filters, add additional
columns, and output to CSV.

Currently we can only output to CSV, not VCF.

A number of useful annotations can be added for each variant by specifying
options of the form '--include-XXX', e.g. '--include-gene'. See detailed help
(run with -h).

Examples
`````````````

Print basic info for the variants found in two VCF files. Note that variants
found in both files are listed in one row, and the 'sources' column lists
the files each variant was found in:

::

    $ varlens-variants test/data/CELSR1/vcfs/vcf_1.vcf test/data/CELSR1/vcfs/vcf_2.vcf

    genome,contig,interbase_start,interbase_end,ref,alt,sources
    GRCh37,22,21829554,21829555,T,G,1.vcf
    GRCh37,22,46931059,46931060,A,C,1.vcf
    GRCh37,22,46931061,46931062,G,A,1.vcf 2.vcf
    GRCh37,22,50636217,50636218,A,C,1.vcf
    GRCh37,22,50875932,50875933,A,C,1.vcf
    GRCh37,22,45309892,45309893,T,G,2.vcf

Same as the above but include additional columns giving varcode variant effect
annotations and the genes the variants overlap, and write to a file:

::

    $ varlens-variants test/data/CELSR1/vcfs/vcf_1.vcf test/data/CELSR1/vcfs/vcf_2.vcf \
        --include-effect \
        --include-gene \
        --out /tmp/result.csv

    Wrote: /tmp/result.csv

    $ cat /tmp/result.csv

    genome,contig,interbase_start,interbase_end,ref,alt,sources,effect,gene
    GRCh37,22,21829554,21829555,T,G,1.vcf,non-coding-transcript,PI4KAP2
    GRCh37,22,46931059,46931060,A,C,1.vcf,p.S670A,CELSR1
    GRCh37,22,46931061,46931062,G,A,1.vcf 2.vcf,p.S669F,CELSR1
    GRCh37,22,50636217,50636218,A,C,1.vcf,intronic,TRABD
    GRCh37,22,50875932,50875933,A,C,1.vcf,splice-acceptor,PPP6R2
    GRCh37,22,45309892,45309893,T,G,2.vcf,p.T214P,PHF21B

Print counts for number of reads supporting reference/variant/other alleles
from the specified BAM, counting only reads with mapping quality >= 10:

::

    $ varlens-variants test/data/CELSR1/vcfs/vcf_1.vcf \
        --include-read-evidence \
        --reads test/data/CELSR1/bams/bam_1.bam \
        --min-mapping-quality 10

    genome,contig,interbase_start,interbase_end,ref,alt,sources,num_alt,num_ref,total_depth
    GRCh37,22,21829554,21829555,T,G,vcf_1.vcf,0,0,0
    GRCh37,22,46931059,46931060,A,C,vcf_1.vcf,0,216,320
    GRCh37,22,46931061,46931062,G,A,vcf_1.vcf,0,321,321
    GRCh37,22,50636217,50636218,A,C,vcf_1.vcf,0,0,0
    GRCh37,22,50875932,50875933,A,C,vcf_1.vcf,0,0,0


varlens-reads
----------------------

Filter reads from one or more BAMs and output a CSV or a new BAM.

Loci and VCF files may be specified, in which case reads are filtered to
overlap the specified loci or variants.

Examples
`````````````

Print basic fields for the reads in a BAM:

::

    $ varlens-reads test/data/CELSR1/bams/bam_0.bam

    query_name,reference_start,reference_end,cigarstring
    HISEQ:142:C5822ANXX:3:2116:16538:101199,46929962,46930062,100M
    HISEQ:142:C5822ANXX:3:1106:18985:32932,46929964,46930064,100M
    HISEQ:142:C5822ANXX:3:2201:21091:67220,46929966,46930066,100M
    HISEQ:142:C5822ANXX:4:1304:5363:12786,46929966,46930066,100M
    HISEQ:142:C5822ANXX:4:1104:9008:85114,46929969,46930069,100M
    HISEQ:142:C5822ANXX:3:2304:9921:94828,46929970,46930070,100M
    HISEQ:142:C5822ANXX:3:2211:6266:74633,46929973,46930073,100M
    HISEQ:142:C5822ANXX:3:1305:8982:42729,46929974,46930074,100M
    HISEQ:142:C5822ANXX:4:2316:5630:7371,46929978,46930078,100M
    ...

Same as above but filter only to reads aligned on the (-) strand, write to a
file instead of stdout, and also include the mapping quality and sequenced
bases in the output:

::

    $ varlens-reads test/data/CELSR1/bams/bam_0.bam \
        --is-reverse \
        --field mapping_quality query_alignment_sequence \
        --out /tmp/result.csv

    Wrote: /tmp/result.csv

    $ head /tmp/result.csv

    query_name,reference_start,reference_end,cigarstring,mapping_quality,query_alignment_sequence
    HISEQ:142:C5822ANXX:3:2116:16538:101199,46929962,46930062,100M,60,CATGATCTGGGCATTAGGGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCG
    HISEQ:142:C5822ANXX:3:1106:18985:32932,46929964,46930064,100M,60,TGATCTGGGCATTAGGGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTC
    HISEQ:142:C5822ANXX:4:1104:9008:85114,46929969,46930069,100M,60,TGGGCATTAGGGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCT
    HISEQ:142:C5822ANXX:4:1202:18451:91174,46929979,46930079,100M,60,GGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCTCAAACATGGG
    HISEQ:142:C5822ANXX:3:1211:18522:54773,46929987,46930087,100M,60,TCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCTCAAACATGGGGGCATTGT
    HISEQ:142:C5822ANXX:3:2114:19455:45093,46929987,46930087,100M,60,TCAGGGTCGTTAGCACGAATCTTTGCCACCGCCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCTCAAACATGGGGGCATTGT
    HISEQ:142:C5822ANXX:4:2115:9153:21593,46929994,46930094,100M,60,CGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCTCAAACATGGGGGCATTGTCATTAAT
    HISEQ:142:C5822ANXX:4:1212:15644:87227,46929995,46930095,100M,60,GTTAGCACGTATGTTTGCCACCACCGACCCCACTGAGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTGCTTCTCAAACATGGGGGCAGTGTCATTAATG
    HISEQ:142:C5822ANXX:3:1103:4717:26369,46929997,46930097,100M,60,TAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCTCAAACATGGGGGCATTGTCATTAATGTC


Write a bam file consisting of reads with mapping quality >=30 and
overlapping a certain locus:

::

    $ varlens-reads test/data/CELSR1/bams/bam_0.bam \
        --min-mapping-quality 30 \
        --locus 22:46932040-46932050 \
        --out /tmp/result.bam

Write a bam file consisting of reads overlapping variants from a VCF:

::

    $ varlens-reads test/data/CELSR1/bams/bam_0.bam \
        --variants test/data/CELSR1/vcfs/vcf_1.vcf \
        --out /tmp/result.bam

Print just the header for a BAM in csv format:

::

    $ varlens-reads test/data/CELSR1/bams/bam_0.bam --header

varlens-allele-support
----------------------

Given one or more BAMs and some genomic sites to consider, write a csv file
giving counts of reads supporting each allele at each site for each BAM.

The genomic sites to consider may be specified by locus (--locus option), or via
one or more VCF files.

The positions outputted by this command are in *interbase coordinates*, i.e.
starting at 0, inclusive on first index, exclusive on second (as opposed to
the one-based inclusive coordinates used in VCF files).

Examples
`````````````

:: 

    varlens-allele-support \
        --reads test/data/CELSR1/bams/bam_1.bam \
        --locus 22:46931061 22:46931063

    source,contig,interbase_start,interbase_end,allele,count
    bam_1.bam,22,46931060,46931061,,1
    bam_1.bam,22,46931060,46931061,G,329
    bam_1.bam,22,46931062,46931063,A,327
    bam_1.bam,22,46931062,46931063,AC,1
    bam_1.bam,22,46931062,46931063,AG,2

Note on coordinate systems
-----------------------------------

``varlens`` uses 0-based half-open coordinates internally. Many tools
(including samtools and VCF files) use inclusive 1-based coordinates. We try to
keep the confusion to a minimum by using the term "interbase" whenever we're
using 0-based half open coordinates and "inclusive" when we're using 1-based
inclusive coordinates.

One particularly sticky place this comes up is when specifying loci on the
commandline using e.g. ``--locus chr22:43243-43244``. To maintain consistency
with the most common other tools, when you specify a locus like
``chr22:10-20``, we interpret that as a 1-based inclusive coordinate. To
specify 0-based half-open coordinates, use this syntax: ``chr22/11-20`` (i.e. a
slash instead of a colon).

See this `blog post <http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html>`_
for more details on coordinate systems.

.. Documentation
    -------------
    The docs are just this readme and the commandline tool help.
    They are available here: http://hammerlab.github.io/varlens/docs/html


