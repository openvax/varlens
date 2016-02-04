.. image:: https://travis-ci.org/hammerlab/varlens.svg?branch=master
    :target: https://travis-ci.org/hammerlab/varlens

varlens
======

A collection of Python tools for working with genomic variants and next-generation sequencing reads.

varlens-reads
    Display, filter, and copy reads from a SAM/BAM file (partial replacement for ``samtools view``)

varlens-variants
    Display and filter variants from VCF or CSV files
    
varlens-allele-support
    Count reads supporting different alleles at specified sites


``varlens`` is built on `varcode <https://github.com/hammerlab/varcode>`_ and `pysam <https://github.com/pysam-developers/pysam>`_.

Installation
-------------

From a git checkout:

::

    pip install .

To run the tests:

::

    nosetests -s .

On Python 3, running nosetests with the `-s` argument to disable output buffering is currently required due to a `pysam` issue.
    
Examples
-------------

varlens-reads
``````````````

Display reads in a BAM file:

::

    $ varlens-reads --reads test/data/gatk_mini_bundle_extract.bam

    query_name,query_alignment_start,query_alignment_end,cigar
    20FUKAAXX100202:6:27:4968:125377,0,101,95M2I4M
    20FUKAAXX100202:8:48:6663:137967,0,101,84M2I15M
    20FUKAAXX100202:3:3:17662:56884,0,101,51M2I48M
    20GAVAAXX100126:1:22:20478:72423,0,88,44M2I42M13S
    20GAVAAXX100126:3:8:19969:70026,0,101,35M2I64M
    20FUKAAXX100202:4:1:11664:56566,0,69,20M2I47M32S
    20FUKAAXX100202:8:28:11647:104607,0,101,19M2I80M
    20GAVAAXX100126:3:67:18892:73513,0,101,19M2I80M
    20FUKAAXX100202:8:23:18439:98544,0,101,12M2I87M
    ...

By default, the CSV output is written to stdout. Specify ``--out`` to write to a file.

With the ``--locus`` or ``--variants`` arguments we can subselect to reads
overlapping the specified sites:

::

    $ varlens-reads \
        --reads test/data/CELSR1/bams/bam_1.bam \
        --variants test/data/CELSR1/vcfs/vcf_1.vcf \
        --variant-genome b37

(The ``--variant-genome`` argument is needed if your VCF file doesn't have a header specifying the genome reference being used.)

Read filters are supported with a flexible Python-based filtering language:

::

    $ varlens-reads \
        --reads test/data/CELSR1/bams/bam_1.bam \
        --read-filter 'not is_duplicate and mapping_quality > 30'

The variables in the filtering expression are coming from a
`pysam.AlignedSegment
<http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment>`_ object.

By default, we get a few columns giving basic information on each read. We can add additional columns by specifying positional arguments:

::

    $ varlens-reads \
        --reads test/data/gatk_mini_bundle_extract.bam \
        query_sequence \
        mate_is_unmapped

    query_name,query_alignment_start,query_alignment_end,cigar,query_sequence,mate_is_unmapped
    20FUKAAXX100202:6:27:4968:125377,0,101,95M2I4M,TGTTAAAATTCAACTGACCATAGGTGTATTGGTTTATTTCTGTACTCTTAGTAGATTCCATTGACCTATATCTCTATCCTTATGCCAGTACCACACTGTTT,False
    20FUKAAXX100202:8:48:6663:137967,0,101,84M2I15M,AACTGACCATAGGTGTATTGGTTTATTTCTGTACTCTTAGTAGATTCCATTGACCTATATCTCTATCCTTATGCCAGTACCACACTGTTTTGTTTACTACA,False
    20FUKAAXX100202:3:3:17662:56884,0,101,51M2I48M,CTCTTAGTAGATTCCATTGACCTATATCTCTATCCTTATGCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTA,False
    20GAVAAXX100126:1:22:20478:72423,0,88,44M2I42M13S,TAGATTCCATTGACCTATATCTCTATCCTTATGCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTCT,False
    20GAVAAXX100126:3:8:19969:70026,0,101,35M2I64M,TTGACCTATATCTCTATCCTTATGCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTCTAACTTTGTT,False
    20FUKAAXX100202:4:1:11664:56566,0,69,20M2I47M32S,ATCCTTATGCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTGTAACTGTGTTTGTTTTTGAAGCGTG,False
    20FUKAAXX100202:8:28:11647:104607,0,101,19M2I80M,TCCTTATGCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTCTAACTTTGTTTGTTTTTCAAGAGTGT,False
    20GAVAAXX100126:3:67:18892:73513,0,101,19M2I80M,TCCTTATGCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTCTAACTTTGTTTGTTTTTCAAGAGTGT,False
    20FUKAAXX100202:8:23:18439:98544,0,101,12M2I87M,GCCAGTACCACACTGTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTCTAACTTTGTTTGTTTTTCAAGAGTGTTTTGACT,False
    ...

Here ``query_sequence`` and ``mate_is_unmapped`` are expressions that are evaluated in the same way as the read filters, and, like read filters, the names are coming from 
`pysam.AlignedSegment
<http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment>`_.

Since these are full-fledged Python expressions, we can do things like:

::

    $ varlens-reads --reads test/data/gatk_mini_bundle_extract.bam 'min(query_qualities)'

Which would include a column giving the minimum of the base qualities for each read.


varlens-variants
``````````````

Here we use the ``varlens-variants`` tool to take the union of the variants in
two VCF files and filter to only those where the reference nucleotide is 'A':

::

    $ varlens-variants \
        --variants test/data/CELSR1/vcfs/vcf_1.vcf \
        --variants test/data/CELSR1/vcfs/vcf_2.vcf \
        --variant-genome b37 \
        --variant-filter 'ref=="A"'

    genome,contig,interbase_start,interbase_end,ref,alt
    GRCh37,22,46931059,46931060,A,C
    GRCh37,22,50636217,50636218,A,C
    GRCh37,22,50875932,50875933,A,C

Similarly to ``varlens-reads``, we can use Python expressions to filter variants and extract additional properties.
The variables available to us are the attributes of a `varcode.Variant <https://github.com/hammerlab/varcode/blob/master/varcode/variant.py>`_ object.

Here we extract the names of the genes each variant overlaps:

::

    $ varlens-variants \
            --variants test/data/CELSR1/vcfs/vcf_1.vcf \
            --variants test/data/CELSR1/vcfs/vcf_2.vcf \
            --variant-genome b37 \
            gene_names

varlens-allele-support
``````````````

Here's an example of the ``varlens-allele-support`` tool:

::

    $ varlens-allele-support \
        --reads test/data/CELSR1/bams/bam_5.bam \
        --locus chr22:46930257 \
        --locus chr22:46930259-46930260

    source,contig,interbase_start,interbase_end,allele,count
    bam_5.bam,22,46930256,46930257,GCC,1
    bam_5.bam,22,46930256,46930257,G,1751
    bam_5.bam,22,46930256,46930257,N,1
    bam_5.bam,22,46930258,46930260,TG,1
    bam_5.bam,22,46930258,46930260,CG,1731
    bam_5.bam,22,46930258,46930260,AG,1
    bam_5.bam,22,46930258,46930260,NN,1
    bam_5.bam,22,46930258,46930260,CT,2

At each locus specified (you can use ``--variants`` to specify loci using a VCF
file), this tool writes out a line for each allele sequenced at that locus. By
default a simple count of the reads supporting each allele is included.
Similarly to the other tools, however, you can also include your own
expressions:

::

    $ varlens-allele-support \
        --reads test/data/CELSR1/bams/bam_5.bam \
        --locus chr22:46930257 \
        --locus chr22:46930259-46930260 \
        --count-group 'not is_duplicate'

    source,contig,interbase_start,interbase_end,allele,count,not is_duplicate
    bam_5.bam,22,46930256,46930257,GCC,1,0
    bam_5.bam,22,46930256,46930257,G,1751,1087
    bam_5.bam,22,46930256,46930257,N,1,0
    bam_5.bam,22,46930258,46930260,TG,1,1
    bam_5.bam,22,46930258,46930260,CG,1731,1081
    bam_5.bam,22,46930258,46930260,AG,1,1
    bam_5.bam,22,46930258,46930260,NN,1,0
    bam_5.bam,22,46930258,46930260,CT,2,1


Here we added a column that gives a count of the non-duplicate reads.


Coordinate systems
-------------

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

See this `blog post <http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html>`_ for more details on coordinate systems.


