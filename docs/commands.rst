Command guide
=============

Worked examples for ``varlens-variants``, ``varlens-reads`` and
``varlens-allele-support``, with the output each prints. Run them from a git
checkout: they use the test data in ``tests/data``, and the effect and gene
annotations need ``pyensembl install --release 75 --species human`` (GRCh37).
Every option is listed in the :doc:`command-line reference <cli>`.

varlens-variants
----------------

Given variants from one or more VCF or CSV files, apply filters, add additional
columns, and output to CSV.

Currently we can only output to CSV, not VCF.

A number of useful annotations can be added for each variant by specifying
options of the form '--include-XXX', e.g. '--include-gene'. See detailed help
(run with -h).

Print basic info for the variants found in two VCF files. Note that variants
found in both files are listed in one row, and the 'sources' column lists
the files each variant was found in:

::

    $ varlens-variants tests/data/CELSR1/vcfs/vcf_1.vcf tests/data/CELSR1/vcfs/vcf_2.vcf

    genome,contig,interbase_start,interbase_end,ref,alt,sources
    GRCh37,22,21829554,21829555,T,G,1
    GRCh37,22,45309892,45309893,T,G,2
    GRCh37,22,46931059,46931060,A,C,1
    GRCh37,22,46931061,46931062,G,A,1 2
    GRCh37,22,50636217,50636218,A,C,1
    GRCh37,22,50875932,50875933,A,C,1

Same as the above but include additional columns giving varcode variant effect
annotations and the genes the variants overlap, and write to a file:

::

    $ varlens-variants tests/data/CELSR1/vcfs/vcf_1.vcf tests/data/CELSR1/vcfs/vcf_2.vcf \
        --include-effect \
        --include-gene \
        --out /tmp/result.csv

    Wrote: /tmp/result.csv

    $ cat /tmp/result.csv
    genome,contig,interbase_start,interbase_end,ref,alt,sources,effect,effect_candidates,gene
    GRCh37,22,21829554,21829555,T,G,1,non-coding-transcript,NoncodingTranscript: non-coding-transcript,PI4KAP2
    GRCh37,22,45309892,45309893,T,G,2,p.T214P,Substitution: p.T214P,PHF21B
    GRCh37,22,46931059,46931060,A,C,1,p.S670A,Substitution: p.S670A,CELSR1
    GRCh37,22,46931061,46931062,G,A,1 2,p.S669F,Substitution: p.S669F,CELSR1
    GRCh37,22,50636217,50636218,A,C,1,intronic,Intronic: intronic,TRABD
    GRCh37,22,50875932,50875933,A,C,1,p.A562fs,FrameShift: p.A562fs; Unresolved: unresolved-intron-retention; Unresolved: unresolved-cryptic; Intronic: intronic,PPP6R2

``effect`` is the producer-preferred consequence on the selected transcript.
``effect_candidates`` lists that consequence and its alternatives in Varcode's
order, including unresolved outcomes. For example, an intronic primary can
retain a protein-changing exon-skip deletion. These are predicted alternatives, not observed products or
calibrated probabilities. MHC predictions use Topiary's shared candidate-aware
selection and can therefore use a coding alternative to that primary.

Print counts for number of reads supporting reference/variant/other alleles
from the specified BAM, counting only reads with mapping quality >= 10. Here
``sources`` is ``vcf_1``, derived from the input file name, and the counts are
written as floating-point numbers:

::

    $ varlens-variants tests/data/CELSR1/vcfs/vcf_1.vcf \
        --include-read-evidence \
        --reads tests/data/CELSR1/bams/bam_1.bam \
        --min-mapping-quality 10

    genome,contig,interbase_start,interbase_end,ref,alt,sources,num_alt,num_ref,total_depth
    GRCh37,22,21829554,21829555,T,G,vcf_1,0.0,0.0,0.0
    GRCh37,22,46931059,46931060,A,C,vcf_1,0.0,216.0,320.0
    GRCh37,22,46931061,46931062,G,A,vcf_1,0.0,321.0,321.0
    GRCh37,22,50636217,50636218,A,C,vcf_1,0.0,0.0,0.0
    GRCh37,22,50875932,50875933,A,C,vcf_1,0.0,0.0,0.0

varlens-reads
-------------

Filter reads from one or more BAMs and output a CSV or a new BAM.

Loci and VCF files may be specified, in which case reads are filtered to
overlap the specified loci or variants.

Print basic fields for the reads in a BAM. The ``source`` column names the BAM
each read came from:

::

    $ varlens-reads tests/data/CELSR1/bams/bam_0.bam

    source,query_name,reference_start,reference_end,cigarstring
    bam_0,HISEQ:142:C5822ANXX:3:2116:16538:101199,46929962,46930062,100M
    bam_0,HISEQ:142:C5822ANXX:3:1106:18985:32932,46929964,46930064,100M
    bam_0,HISEQ:142:C5822ANXX:3:2201:21091:67220,46929966,46930066,100M
    bam_0,HISEQ:142:C5822ANXX:4:1304:5363:12786,46929966,46930066,100M
    bam_0,HISEQ:142:C5822ANXX:4:1104:9008:85114,46929969,46930069,100M
    bam_0,HISEQ:142:C5822ANXX:3:2304:9921:94828,46929970,46930070,100M
    bam_0,HISEQ:142:C5822ANXX:3:2211:6266:74633,46929973,46930073,100M
    bam_0,HISEQ:142:C5822ANXX:3:1305:8982:42729,46929974,46930074,100M
    bam_0,HISEQ:142:C5822ANXX:4:2316:5630:7371,46929978,46930078,100M
    ...

Same as above but filter only to reads aligned on the (-) strand, write to a
file instead of stdout, and also include the mapping quality and sequenced
bases in the output:

::

    $ varlens-reads tests/data/CELSR1/bams/bam_0.bam \
        --is-reverse \
        --field mapping_quality query_alignment_sequence \
        --out /tmp/result.csv

    Wrote: /tmp/result.csv

    $ head -n 4 /tmp/result.csv

    source,query_name,reference_start,reference_end,cigarstring,mapping_quality,query_alignment_sequence
    bam_0,HISEQ:142:C5822ANXX:3:2116:16538:101199,46929962,46930062,100M,60,CATGATCTGGGCATTAGGGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCG
    bam_0,HISEQ:142:C5822ANXX:3:1106:18985:32932,46929964,46930064,100M,60,TGATCTGGGCATTAGGGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTC
    bam_0,HISEQ:142:C5822ANXX:4:1104:9008:85114,46929969,46930069,100M,60,TGGGCATTAGGGCCTTCATCAGGGTCGTTAGCACGAATCTTTGCCACCACCGACCCCACTGGGTTGTTCTCCTCAACAAACAGCTCCAGTTCGTCCTTCT

Write a bam file consisting of reads with mapping quality >=30 and
overlapping a certain locus:

::

    $ varlens-reads tests/data/CELSR1/bams/bam_0.bam \
        --min-mapping-quality 30 \
        --locus 22:46932040-46932050 \
        --out /tmp/result.bam

    Sorting read file /tmp/result.bam
    Wrote 65 reads: /tmp/result.bam

Write a bam file consisting of reads overlapping variants from a VCF:

::

    $ varlens-reads tests/data/CELSR1/bams/bam_0.bam \
        --variants tests/data/CELSR1/vcfs/vcf_1.vcf \
        --out /tmp/result.bam

    Sorting read file /tmp/result.bam
    Wrote 51 reads: /tmp/result.bam

Print just the header for a BAM in csv format:

::

    $ varlens-reads tests/data/CELSR1/bams/bam_0.bam --header

    read_source,group,index,key,value
    bam_0,HD,0,VN,
    bam_0,HD,1,GO,
    bam_0,HD,2,SO,
    bam_0,SQ,0,SN,chrM
    ...

varlens-allele-support
----------------------

Given one or more BAMs and some genomic sites to consider, write a csv file
giving counts of reads supporting each allele at each site for each BAM.

The genomic sites to consider may be specified by locus (--locus option), or via
one or more VCF files.

The positions outputted by this command are in *interbase coordinates*, i.e.
starting at 0, inclusive on first index, exclusive on second (as opposed to
the one-based inclusive coordinates used in VCF files). See
:doc:`coordinates`.

Insertions use an empty interval at the insertion boundary. For example,
VCF ``1:100 A>AT`` becomes ``interbase_start=100, interbase_end=100`` with
``ref=""`` and ``alt="T"``. The pileup reads the base immediately to the left
and counts only the inserted sequence. Version 1.0 corrects the older
single-base insertion intervals; CSV input still accepts that older encoding.

::

    $ varlens-allele-support \
        --reads tests/data/CELSR1/bams/bam_1.bam \
        --locus 22:46931061 22:46931063

    source,contig,interbase_start,interbase_end,allele,count
    bam_1,22,46931060,46931061,G,329
    bam_1,22,46931060,46931061,,1
    bam_1,22,46931062,46931063,A,327
    bam_1,22,46931062,46931063,AG,2
    bam_1,22,46931062,46931063,AC,1
