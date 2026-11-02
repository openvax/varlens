Using varlens from Python
=========================

The read-support table and MHC binding summary that the commands build are
also available from Python. This page describes their shapes; the
:doc:`API reference <varlens>` lists every module.

Read support
------------

``varlens.support.variant_support(variants, counts)`` returns a pandas DataFrame
with variants as its row index and two column levels named ``measurement`` and
``source``. This replaces the removed pandas Panel API in Varlens 1.0::

    evidence = variant_support(variants, counts)
    alt_counts = evidence["num_alt"]  # one column per read source
    tumor_depth = evidence[("total_depth", "tumor")]

The metrics are ``num_alt``, ``num_ref``, ``num_other``, ``total_depth``,
``alt_fraction``, and ``any_alt_fraction``. Commandline CSV columns retain their
existing names. Counts include low-quality bases and overlapping mates unless
read filters are requested. Common ``chr`` aliases resolve to BAM references;
alternate and unplaced contig names are preserved.

MHC binding
-----------

MHC binding uses the current Topiary API with a local NetMHCpan installation.
It reports the lowest affinity in nM among peptides overlapping the mutation,
and the corresponding input allele. Variants without predictions have missing
values. Prediction caches distinguish peptide-length requests. Effect labels
come from the installed Varcode version and may change as its annotations improve.

For fusion transcripts, see :doc:`fusions`.
