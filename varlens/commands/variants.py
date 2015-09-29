'''
Variant manipulation.

Collect variants specified on the commandline or in VCF or CSV files, apply
filters, and write out a CSV file.

%(prog)s 
    --variants /path/to/first.vcf \
    --variants /path/to/second.vcf \
    --variant-filter "ref=='A'" \
    --out-variants result.csv

'''
from __future__ import absolute_import

import argparse
import sys
import logging
import collections

from .. import variants_util

parser = argparse.ArgumentParser(usage=__doc__)
variants_util.add_args(parser)
parser.add_argument("--field", action="append", default=[], nargs="+")
parser.add_argument("--no-standard-fields", action="store_true", default=False,
    help="Do not write standard fields (contig, genome, start, end, ref, alt)")

parser.add_argument("--out")
parser.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)

    variants = variants_util.load_from_args(args)
    if not variants:
        parser.error("No variants specified.")

    logging.info("Loaded %d variants." % len(variants))

    extra_columns = collections.OrderedDict()
    for lst in args.field:
        if len(lst) == 1:
            name = expression = lst[0]
        elif len(lst) == 2:
            (name, expression) = lst
        else:
            parser.error("Expected 1 or 2 arguments ([name], expression) for "
                "--field option, but got: %s" % lst)
        extra_columns[name] = expression

    df = variants_util.variants_to_dataframe(variants, extra_columns)
    if args.no_standard_fields:
        for column in variants.STANDARD_DATAFRAME_COLUMNS:
            del df[column]

    if args.out is None:
        # Write to stdout.
        df.to_csv(sys.stdout, index=False)
    elif args.out.endswith(".csv"):
        df.to_csv(args.out, index=False)
        print("Wrote: %s" % args.out)
    else:
        parser.error("Unsupported output file extension: %s" % args.out)

