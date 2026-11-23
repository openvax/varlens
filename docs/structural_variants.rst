Structural variants and empty inputs
====================================

How ``varlens-variants`` loads, exports and filters structural variants (SVs),
which commands reject them, and what empty or all-filtered inputs produce. Read
this if your VCFs contain symbolic alleles or breakends.

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
  record ID, filter, INFO and sample fields, plus the original VCF line. Every
  retained record is kept when identical events are merged, within or across
  sources.

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
