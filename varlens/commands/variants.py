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

from .. import variants

parser = argparse.ArgumentParser(usage=__doc__)
variants.add_args(parser)
parser.add_argument("--field", action="append", default=[], nargs="+")
parser.add_argument("--no-standard-fields", action="store_true", default=False,
    help="Do not write standard fields (contig, genome, start, end, ref, alt)")

parser.add_argument("--out-variants", required=True)
parser.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)

    vc = variants.load(args)
    if not vc:
        parser.error("No variants specified.")

    logging.info("Loaded %d variants." % len(vc))

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

    df = variants.variants_to_dataframe(vc, extra_columns)
    if args.no_standard_fields:
        for column in variants.STANDARD_DATAFRAME_COLUMNS:
            del df[column]

    if args.out_csv.endswith(".csv"):
        df.write_csv(args.out_variants)
    else:
        parser.error("Unsupported output file extension: %s" % args.out_variants)

    print("Wrote: %s" % args.out_variants)
