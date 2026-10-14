.. image:: https://github.com/openvax/varlens/actions/workflows/tests.yml/badge.svg
    :target: https://github.com/openvax/varlens/actions/workflows/tests.yml
    :alt: Tests

varlens
======================

A collection of Python tools for working with genomic variants and
next-generation sequencing reads. Not particularly fast for large datasets. The
emphasis is on extracting what you need from BAMs and VCFs into a CSV file for
further analysis.

Built on `varcode <https://github.com/openvax/varcode>`_ and `pysam <https://github.com/pysam-developers/pysam>`_.

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

Varlens 1.0 requires Python 3.11 or newer and pandas 3. The package uses
the current Topiary 5 prediction API.

To install from `PyPI <https://pypi.python.org/pypi/varlens>`_:

::

    pip install varlens

Or from a git checkout:

::

    pip install .

To run the tests:

::

    pip install -e '.[test]'
    ./lint.sh
    ./test.sh

Tests that need Ensembl annotation data are skipped unless human release 75
is installed. To run these tests as well:

::

    pyensembl install --release 75 --species human --custom-mirror https://github.com/openvax/ensembl-data/releases/download/GRCh37.75/
    ./test.sh

The Topiary integration test uses the real annotation and prediction pipeline
with an in-process test predictor, so it does not need NetMHCpan. To also test
the external predictor, install NetMHCpan and its model data, put ``netMHCpan``
on ``PATH``, then run ``./test.sh --run-netmhcpan``. MHC affinity results use
Topiary's normalized allele names, such as ``HLA-A*02:01``.

CI tests Python 3.11 through 3.14 and runs reference-dependent annotation
tests in a separate job.

To build the documentation (just this README plus the commandline tool help):

::

    pip install -e .
    pip install Sphinx
    cd docs
    make clean setup rst html

The docs will be written to the ``_build/html`` directory.


Python API
------------

The same read-evidence tools are available directly from Python:

.. code-block:: python

    from varlens.locus import Locus
    from varlens.read_evidence import PileupCollection

    locus = Locus.from_inclusive_coordinates("chr17", 41244936)
    evidence = PileupCollection.from_bam("sample.bam", [locus])
    counts = evidence.allele_summary(locus)
    filtered = evidence.filter(min_base_quality=20, min_mapping_quality=30)

The BAM must be indexed. ``Locus`` stores zero-based, half-open coordinates;
``from_inclusive_coordinates`` accepts one-based inclusive coordinates.

``varlens.support.allele_support_df(loci, sources)`` returns an allele-count
DataFrame from ``ReadSource`` objects. ``variant_support(variants, counts)``
now returns a DataFrame indexed by variant, with two column levels named
``measurement`` and ``source``. This replaces the removed ``pandas.Panel``
return type. For example, ``result["num_alt"]["tumor"]`` selects alternate
read counts for the source named ``tumor``.


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

Structural variants and empty inputs
-----------------------------------

``varlens-variants`` loads symbolic SV alleles and paired/single breakends from
VCF using Varcode. A mixed SNV/SV file retains both kinds of event. A malformed
structural allele that Varcode cannot parse causes an error before output is
written. Empty or all-filtered sources produce a header-only CSV, including
when ``--chunk-rows`` is used. By default, non-PASS calls are excluded; use
``--include-failing-variants`` to include them. VCF ``.`` and ``*`` ALT
placeholders are not independent variant calls. ``--max-variants-per-source``
explicitly limits the number of loaded alleles.

Structural rows add two JSON-valued CSV columns:

* ``structural_variant``: the ``varlens.structural_variant.v1`` envelope and
  Varcode constructor fields, including type, original alleles, affected span,
  partner breakpoint, orientation, uncertainty, assembly and INFO data. Genome
  objects are represented by their canonical assembly name.
* ``vcf_records``: source paths and retained Varcode VCF metadata, including
  record ID, filter, INFO and sample fields. Records from different sources are
  retained when identical events are merged.

These columns are retained without ``--include-metadata`` and are required for
a structural CSV round trip. Do not edit the displayed locus/alleles independently
of the typed payload. Custom ``--columns`` output can deliberately omit fields,
but then cannot necessarily be reloaded. Ordinary variant CSVs retain their
existing schema.

For SVs, ``--variant-locus`` selects **known breakpoint bases**, including a
remote mate, rather than every base in a deleted/duplicated interval. Confidence
intervals are preserved but do not widen this filter. Raw contig spellings are
retained in the payload; displayed contigs follow Varcode normalization.

SV export, reload, REF/ALT filtering and breakpoint filtering are supported.
Gene/effect/context annotations, the best-affinity variant summary, and
single-interval read evidence currently reject SV rows with an actionable
diagnostic. Likewise, ``varlens-reads --variants`` and
``varlens-allele-support --variants`` reject structural events. Explicit
``--locus`` read export remains available: those alignments are raw observations,
not counts of molecules supporting a junction.

varlens-fusions
---------------

Install the optional Isovar adapter and reconstruct supplied fusion RNA::

    pip install 'varlens[fusion]'
    varlens-fusions fusion.json.gz --out reconstruction.json

Optionally predict junction-overlapping peptides with Topiary and mhctools::

    varlens-fusions fusion.json.gz --hla 'HLA-A*02:01' 'HLA-B*07:02' \
        --peptide-lengths 8 9 10 11 --out predictions.json

NetMHCpan and its model data must be installed for the second command. The
reconstruction-only command does not invoke it. Multiple JSON/JSON.gz input
files can be supplied, and ``--out`` also accepts JSON.gz; omitted output goes
to stdout. A failed input or prediction prevents writing a partial report.

Each input is an `Isovar supplied-fusion input
<https://github.com/openvax/isovar/blob/v1.53.5/isovar/fusion.py>`_: a JSON object
with ``fusion``, ``references`` and ``reads`` (plus optional ``reference_names``).
``fusion`` supplies event identity, assembly, oriented donor/acceptor breakpoints,
RNA sequence, mapped blocks, half-open junction coordinates and provenance with
a sample ID. References supply annotated transcript sequence/mapping and CDS
coordinates; reads supply observed RNA sequence/mapping and sample/library/read/
fragment identity. Isovar validates these fields and transfers frames only from
compatible annotated donors. The ``note`` and ``original_records`` fields in
the Sid fixture audit envelopes are retained in the report, but are not counted
as additional RNA evidence. Unknown fields are rejected.

The output ``varlens.fusions.v1`` report contains an ``events`` list. Each
``varlens.fusion.v1`` event retains the entire input, the unmodified
``isovar.fusion_rna.v3`` reconstruction, software versions, Topiary fragments,
and full peptide prediction rows. ``prediction_status`` distinguishes
``not_requested``, ``no_coding_products`` and ``completed``; a completed run may
have no target-overlapping predictions.

Unresolved and insufficient-support outcomes keep their reasons and measured
evidence. Ambiguous or partial-CDS coding hypotheses remain separate fragments
with explicit status on every prediction; Varlens does not choose an isoform.
Junction read/fragment counts retain their supplied-observation scope. Unknown
reference counts, protein-observation evidence, UMI/cell counts and comparator
sequences are not filled with zero or inferred from alignment counts.
Junction overlap is not a claim of proteome-wide novelty or tumor specificity.
This command does not discover fusions, assemble arbitrary BAMs, translate
structural VCFs, or establish whole-sample expression from selected reads.

Python callers can use the same workflow::

    from varlens.fusions import analyze_fusion
    result = analyze_fusion(payload, alleles=["HLA-A*02:01"])

``analyze_fusion(payload)`` performs reconstruction only. A configured
``topiary.TopiaryPredictor`` can be passed as ``predictor=`` instead of alleles.
``varlens.mhc_binding.fragment_predictions`` also accepts existing
``topiary.ProteinFragment`` objects and returns target-overlapping prediction
rows with source, sample and evidence metadata intact.

The bundled Sid regressions cover ATP5MG--KMT2A (a coding hypothesis with
ambiguous donor annotation), GABBR1--SLC29A1 T1 (unresolved frame with direct
junction support), and the latter's unchanged raw alignment export. Their
upstream sources and checksums are in ``tests/data/fusions/manifest.json``.

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
    They are available here: http://openvax.github.io/varlens/docs/html


