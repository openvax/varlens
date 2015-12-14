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

from . import configure_logging
from .. import variant_includes
from .. import variants_util
from ..evaluation import parse_labeled_expression

parser = argparse.ArgumentParser(usage=__doc__)
variants_util.add_args(parser)
parser.add_argument("field", nargs="*")
parser.add_argument("--no-standard-columns",
    action="store_true", default=False,
    help="Don't write standard columns (genome, contig, start, end, ref, alt)")

parser.add_argument("--chunk-rows", metavar="N", type=int,
    help="Write out current results after processing N rows.")

parser.add_argument("--limit", metavar="N", type=int,
    help="Process only the first N variants (useful for testing)")

parser.add_argument("--columns",
    help="Column separated list of columns to output")

parser.add_argument("--out")

parser.add_argument('--include-metadata', action="store_true", default=False,
    help="Output variant metadata when loading from VCF (info column, etc).")

parser.add_argument('--include-variant-source',
    action="store_true", default=False,
    help="Output the source file each variant originates from")

for includeable in variant_includes.INCLUDEABLES:
    includeable.add_args(parser)

parser.add_argument("-v", "--verbose", action="store_true", default=False)


def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)
    configure_logging(args)

    df = variants_util.load_from_args_as_dataframe(args)
    if df is None:
        parser.error("No variants specified.")

    logging.info("Loaded %d variants." % df.shape[0])

    extra_columns = collections.OrderedDict()
    for labeled_expression in args.field:
        (label, expression) = parse_labeled_expression(labeled_expression)
        extra_columns[label] = expression

    for (column, callable_or_expression) in extra_columns.items():
        df[column] = [
            variants_util.evaluate_variant_expression(
                callable_or_expression, row.to_dict(), row.variant)
            for (i, row) in df.iterrows()
        ]

    def save():
        if args.columns:
            columns = [x.strip() for x in args.columns.split(",")]
        else:
            columns = [x for x in df.columns.tolist() if x != "variant"]
            if not args.include_metadata:
                columns = [
                    x for x in columns if not x.startswith("metadata_")
                ]
            if args.no_standard_columns:
                columns = [
                    x for x in columns
                    if x not in variants_util.STANDARD_DATAFRAME_COLUMNS
                ]
            if not args.include_variant_source:
                columns = [
                    x for x in columns
                    if x != "variant_source"
                ]

        df_save = df[columns]
        if args.out is None:
            # Write to stdout.
            df_save.to_csv(sys.stdout, index=False)
        elif args.out.endswith(".csv"):
            df_save.to_csv(args.out, index=False)
            print("Wrote: %s" % args.out)
        else:
            parser.error("Unsupported output file extension: %s" % args.out)

    for includeable in variant_includes.INCLUDEABLES:
        if includeable.requested(args):
            logging.info("Running includeable: %s" % includeable.name)
            instance = includeable.from_args(args)
            for num_rows in instance.compute(df, chunk_rows=args.chunk_rows):
                if args.chunk_rows is not None:
                    save()

    save()
 
