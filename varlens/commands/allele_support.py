'''
Count reads supporting 

Given genomic loci (e.g. from a VCF file) and one or more BAM files, write a
csv file giving read counts for each allele at each locus.

Note that the output of this command is given in *interbase coordinates*, i.e.
starting at 0, inclusive on first index, exclusive on second.

%(prog)s \
    --reads /path/to/bam1
    --reads /path/to/bam2
    --variants /path/to/first.vcf
    --variants /path/to/second.vcf

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
loci_util.add_args(parser)
reads_util.add_args(parser)
variants_util.add_args(parser)

parser.add_argument("--out")
parser.add_argument("-v", "--verbose", action="store_true", default=False)

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
