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
from . import configure_logging
from .. import support

parser = argparse.ArgumentParser(usage=__doc__)
loci_util.add_args(parser)
reads_util.add_args(parser)

parser.add_argument("--out")
support.add_args(parser)
parser.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)
    configure_logging(args)
    
    loci = loci_util.load_from_args(args)
    if not loci:
        parser.error("No genomic loci (e.g. VCF files) specified.")

    logging.info("Loaded %d genomic loci." % len(loci))

    read_sources = reads_util.load_from_args(args)

    if read_sources is None:
        parser.error("No read sources (--reads argument) specified.")

    out_fd = open(args.out, "w") if args.out else sys.stdout
    writer = csv.writer(out_fd)

    rows_generator = support.allele_support_rows(
        loci, read_sources, ["count:all"] + args.count_group)
    for (i, row) in enumerate(rows_generator):
        if i == 0:
            writer.writerow(row.index.tolist())
        writer.writerow([str(x) for x in row])

    if out_fd is not sys.stdout:
        out_fd.close()
        print("Wrote: %s" % args.out)
