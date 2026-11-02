.. image:: https://github.com/openvax/varlens/actions/workflows/tests.yml/badge.svg
    :target: https://github.com/openvax/varlens/actions/workflows/tests.yml
    :alt: Tests

.. image:: https://img.shields.io/pypi/v/varlens.svg
    :target: https://pypi.org/project/varlens/
    :alt: PyPI

varlens
======================

Varlens is a set of command-line tools that pull what you need out of VCF and
BAM files into CSV tables for further analysis: combine and annotate variants,
export or filter reads, and count the reads supporting each allele. It is not
particularly fast for large datasets.

It is built on `Varcode <https://github.com/openvax/varcode>`_ (variants and
effects) and `pysam <https://github.com/pysam-developers/pysam>`_ (reads), with
`Topiary <https://github.com/openvax/topiary>`_ and
`MHCtools <https://github.com/openvax/mhctools>`_ for MHC binding and
`Isovar <https://github.com/openvax/isovar>`_ for fusion RNA.

varlens-variants
    Combine, annotate, and filter variants from VCF or CSV files. Available
    annotations include genes, variant effects, surrounding sequence context,
    counts of supporting reads from specified BAM files, and MHC I binding
    affinity prediction of mutant peptides.

varlens-reads
    Display, filter, and copy reads from a SAM/BAM file. Partial replacement for ``samtools view``.

varlens-allele-support
    Count reads supporting each allele at specified sites in BAM files.

varlens-fusions
    Reconstruct a supplied fusion transcript from RNA reads with Isovar, and
    optionally predict MHC binding of junction peptides.

Install
-------

::

    pip install varlens

Varlens 1.1 requires Python 3.11 or newer, pandas 3, and
MHCTools 4.1.1 or newer (below 5). The package uses the current Topiary 5
prediction API; MHC binding predictions also need a local NetMHCpan
installation. ``pip install 'varlens[fusion]'`` adds Isovar for
``varlens-fusions``. Effect and gene annotations need a PyEnsembl release
matching your variants' genome build; the examples below use GRCh37
(``pyensembl install --release 75 --species human``).

Quick start
-----------

From a git checkout, combine the variants in two test VCFs, annotate each with
its predicted effect and the genes it overlaps, and print a CSV:

::

    $ varlens-variants tests/data/CELSR1/vcfs/vcf_1.vcf tests/data/CELSR1/vcfs/vcf_2.vcf \
        --include-effect \
        --include-gene

    genome,contig,interbase_start,interbase_end,ref,alt,sources,effect,effect_candidates,gene
    GRCh37,22,21829554,21829555,T,G,1,non-coding-transcript,NoncodingTranscript: non-coding-transcript,PI4KAP2
    GRCh37,22,45309892,45309893,T,G,2,p.T214P,Substitution: p.T214P,PHF21B
    GRCh37,22,46931059,46931060,A,C,1,p.S670A,Substitution: p.S670A,CELSR1
    GRCh37,22,46931061,46931062,G,A,1 2,p.S669F,Substitution: p.S669F,CELSR1
    GRCh37,22,50636217,50636218,A,C,1,intronic,Intronic: intronic,TRABD
    GRCh37,22,50875932,50875933,A,C,1,p.A562fs,FrameShift: p.A562fs; Unresolved: unresolved-intron-retention; Unresolved: unresolved-cryptic; Intronic: intronic,PPP6R2

A variant found in both files is listed once, and ``sources`` names the files
it came from (``1 2``). Positions are 0-based half-open (*interbase*), unlike
the 1-based positions in VCF. ``effect`` is Varcode's preferred consequence and
``effect_candidates`` lists it with its alternatives; these are predictions,
not observed products.

Command line
------------

Each command prints its options with ``--help`` and the installed Varlens
version with ``--version``. One example per command:

::

    # Count reads supporting each variant's ref and alt alleles in a BAM
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

    # Write the reads overlapping a locus (1-based, inclusive) to a new BAM
    $ varlens-reads tests/data/CELSR1/bams/bam_0.bam \
        --min-mapping-quality 30 \
        --locus 22:46932040-46932050 \
        --out /tmp/result.bam

    # Count the reads supporting each allele at two sites
    $ varlens-allele-support \
        --reads tests/data/CELSR1/bams/bam_1.bam \
        --locus 22:46931061 22:46931063

    source,contig,interbase_start,interbase_end,allele,count
    bam_1,22,46931060,46931061,G,329
    bam_1,22,46931060,46931061,,1
    bam_1,22,46931062,46931063,A,327
    bam_1,22,46931062,46931063,AG,2
    bam_1,22,46931062,46931063,AC,1

    # Reconstruct a supplied fusion transcript (needs varlens[fusion])
    $ varlens-fusions fusion.json.gz --out reconstruction.json

.. end-of-docs-front-page

Documentation
-------------

Full documentation: `openvax.github.io/varlens <https://openvax.github.io/varlens/>`_.

- `Command guide <https://openvax.github.io/varlens/docs/html/commands.html>`_:
  more examples for each command, with their output
- `Structural variants and empty inputs <https://openvax.github.io/varlens/docs/html/structural_variants.html>`_
- `Fusion RNA with varlens-fusions <https://openvax.github.io/varlens/docs/html/fusions.html>`_
- `Using varlens from Python <https://openvax.github.io/varlens/docs/html/python.html>`_
- `Coordinate systems <https://openvax.github.io/varlens/docs/html/coordinates.html>`_
- `Command-line reference <https://openvax.github.io/varlens/docs/html/cli.html>`_

Development
-----------

::

    pip install -e '.[test,fusion]'
    pyensembl install --release 75 --species human
    ./lint.sh
    ./test.sh

NetMHCpan tests are opt-in (``./test.sh --run-netmhcpan``). The
`development guide <https://openvax.github.io/varlens/docs/html/development.html>`_
covers the shared OpenVax dependency check and building the docs.

License
-------

Varlens is released under the
`Apache License 2.0 <https://github.com/openvax/varlens/blob/master/LICENSE>`_.
