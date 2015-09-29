'''
Extract reads from a BAM at specified loci.

%(prog)s \
    --reads /path/to/bam1.bam
    --variants /path/to/variants.vcf
    --out out.bam

'''
from __future__ import absolute_import

import argparse
import sys
import csv
import collections
import logging

import pysam

from .. import loci_util
from .. import reads_util
from ..evaluation import parse_labeled_expression

parser = argparse.ArgumentParser(usage=__doc__)
loci_util.add_args(parser)
reads_util.add_args(parser)

parser.add_argument("--out")
parser.add_argument("field", nargs="*")
parser.add_argument("--no-standard-fields", action="store_true", default=False)
parser.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)

    read_sources = reads_util.load_from_args(args)
    if not read_sources:
        parser.error("No read sources specified.")

    loci = loci_util.load_from_args(args)  # may be None        

    out_pysam_handle = None
    out_csv_writer = out_csv_fd = None
    if args.out and (args.out.endswith(".bam") or args.out.endswith(".sam")):
        if args.field:
            parser.error("Don't specify fields when outputting to bam or sam.")

        out_pysam_handle = pysam.AlignmentFile(
            args.out,
            "wb" if args.out.endswith(".bam") else "w",
            template=read_sources[0].handle)

    elif not args.out or args.out.endswith(".csv"):
        out_csv_fd = open(args.out, "w") if args.out else sys.stdout
        out_csv_writer = csv.writer(out_csv_fd)

        columns = collections.OrderedDict()
        if not args.no_standard_fields:
            columns["query_name"] = lambda x: x.query_name
            columns["query_alignment_start"] = (
                lambda x: x.query_alignment_start)
            columns["query_alignment_end"] = lambda x: x.query_alignment_end
            columns["cigar"] = lambda x: x.cigarstring

        for labeled_expression in args.field:
            (label, expression) = parse_labeled_expression(labeled_expression)
            columns[label] = expression

        out_csv_writer.writerow(columns.keys())
    else:
        parser.error(
            "Don't know how to write to file with output extension: %s. "
            "Supported extensions: csv, bam, sam." % args.out)

    for read_source in read_sources:
        for read in read_source.reads(loci):
            if out_pysam_handle is not None:
                out_pysam_handle.write(read)
            if out_csv_writer is not None:
                out_csv_writer.writerow([
                    str(reads_util.evaluate_read_expression(e, read))
                    for e in columns.values()
                ])

    if out_pysam_handle is not None:
        out_pysam_handle.close()
        print("Wrote: %s" % args.out)

    if out_csv_fd is not None and out_csv_fd is not sys.stdout:
        out_csv_fd.close()
        print("Wrote: %s" % args.out)

