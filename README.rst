.. image:: https://travis-ci.org/hammerlab/varlens.svg?branch=master
    :target: https://travis-ci.org/hammerlab/varlens

varlens
======

A collection of Python tools for working with genomic variants and next-generation sequencing reads.

varlens-variants
    Combine, annotate, and filter variants from VCF or CSV files

varlens-reads
    Display, filter, and copy reads from a SAM/BAM file (partial replacement for ``samtools view``)

varlens-allele-support
    Count reads supporting each allele at specified sites in BAM files


``varlens`` is built on `varcode <https://github.com/hammerlab/varcode>`_ and `pysam <https://github.com/pysam-developers/pysam>`_.


Documentation
-------------
Available at: http://timodonnell.github.io/varlens/docs/html

Installation
-------------

From a git checkout:

::

    pip install .

To run the tests:

::

    nosetests .

To build the documentation:

::

    pip install -e .
    pip install Sphinx
    cd docs
    make clean setup rst html

The docs will be written to the ``_build/html`` directory.

