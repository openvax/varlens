'''
Given one or more BAMs and some genomic sites to consider, write a csv file
giving counts of reads supporting each allele at each site for each BAM.

The genomic sites to consider may be specified by locus (--locus option), or via
one or more VCF files.

The positions outputted by this command are in *interbase coordinates*, i.e.
starting at 0, inclusive on first index, exclusive on second (as opposed to
the one-based inclusive coordinates used in VCF files).

Example:

%(prog)s \\
    --reads test/data/CELSR1/bams/bam_1.bam \\
    --locus 22:46931061 22:46931063

'''
from __future__ import absolute_import

import argparse
import csv
import sys
import logging

from .. import loci_util
from .. import reads_util
from .. import variants_util

from . import configure_logging
from .. import support
from ..read_evidence.pileup_collection import to_locus

parser = argparse.ArgumentParser(usage=__doc__)
group = parser.add_argument_group("output arguments")
group.add_argument("--out")
group.add_argument("-v", "--verbose", action="store_true", default=False)
loci_util.add_args(parser.add_argument_group("loci specification"))
variants_util.add_args(parser)
reads_util.add_args(parser)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)
    configure_logging(args)

    loci = loci_util.load_from_args(args)  # may be None
    variants_df = variants_util.load_from_args_as_dataframe(args)
    if variants_df is not None:
        variant_loci = loci_util.Loci(
            to_locus(variant)
            for variant in variants_df["variant"])
        loci = variant_loci if loci is None else loci.union(variant_loci)
    
    if not loci:
        if variants_df is not None:
            parser.error("No loci: variants specified but none remained "
                "after filtering")
        else:
            parser.error("No genomic loci or variants specified.")

    logging.info("Loaded %d genomic loci." % len(loci))

    read_sources = reads_util.load_from_args(args)

    if read_sources is None:
        parser.error("No read sources (--reads argument) specified.")

    out_fd = open(args.out, "w") if args.out else sys.stdout
    writer = csv.writer(out_fd)

    rows_generator = support.allele_support_rows(loci, read_sources)
    for (i, row) in enumerate(rows_generator):
        if i == 0:
            writer.writerow(row.index.tolist())
        writer.writerow([str(x) for x in row])

    if out_fd is not sys.stdout:
        out_fd.close()
        print("Wrote: %s" % args.out)
