'''
Given variants from one or more VCF or CSV files, apply filters, add additional
columns, and output to CSV.

Currently we can only output to CSV, not VCF.

A number of useful annotations can be added for each variant by specifying
options of the form '--include-XXX', e.g. '--include-gene'. See detailed help
below.

Examples:

Print basic info for the variants found in two VCF files. Note that variants
found in both files are listed in one row, and the 'sources' column lists
the files each variant was found in:

    %(prog)s test/data/CELSR1/vcfs/vcf_1.vcf test/data/CELSR1/vcfs/vcf_2.vcf

Same as the above but include additional columns giving varcode variant effect
annotations and the genes the variants overlap, and write to a file:

    %(prog)s test/data/CELSR1/vcfs/vcf_1.vcf test/data/CELSR1/vcfs/vcf_2.vcf \\
        --include-effect \\
        --include-gene \\
        --out /tmp/result.csv

Print counts for number of reads supporting reference/variant/other alleles
from the specified BAMs, counting only reads with mapping quality >= 10:

    %(prog)s test/data/CELSR1/vcfs/vcf_1.vcf \\
        --include-read-evidence \\
        --reads test/data/CELSR1/bams/*.bam \\
        --min-mapping-quality 10

'''
from __future__ import absolute_import

import argparse
import sys
import logging

from . import configure_logging
from .. import variant_includes
from .. import variants_util

parser = argparse.ArgumentParser(usage=__doc__)
variants_util.add_args(parser, positional=True)

group = parser.add_argument_group("variant output")

group.add_argument("--no-standard-columns",
    action="store_true", default=False,
    help="Don't write standard columns (genome, contig, start, end, ref, alt)")

group.add_argument("--chunk-rows", metavar="N", type=int,
    help="Write out current results after processing N rows.")

group.add_argument("--limit", metavar="N", type=int,
    help="Process only the first N variants (useful for testing)")

group.add_argument("--columns",
    help="Column separated list of columns to output")

group.add_argument("--rename-column", nargs=2, action="append", default=[],
    metavar="COL",
    help="Rename output column first argument to second. Can be specified "
    "multiple times by repeating the --rename-column option.")

group.add_argument("--out",
    help="Output file. If not specified the CSV is written to stdout.")

group.add_argument('--include-metadata', action="store_true", default=False,
    help="Output variant metadata when loading from VCF (info column, etc).")

for includeable in variant_includes.INCLUDEABLES:
    includeable.add_args(parser)

group.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)
    configure_logging(args)

    df = variants_util.load_from_args_as_dataframe(args)
    if df is None:
        parser.error("No variants specified.")

    logging.info("Loaded %d variants." % df.shape[0])

    # We run the inverse of the column renames on the input df.
    column_renames = {}
    if args.rename_column:
        column_renames = dict(args.rename_column)
        column_renames_inverse = dict((v, k) for (k, v) in args.rename_column)
        if len(column_renames) != len(column_renames_inverse):
            raise ValueError("Column renames are not 1:1")

        df.columns = [
            column_renames_inverse.get(col, col) for col in df.columns
        ]

    def save(df):
        if column_renames:
            df = df.copy()
            df.columns = [column_renames.get(col, col) for col in df.columns]

        if args.columns:
            columns = [x.strip() for x in args.columns.split(",")]
        else:
            columns = [x for x in df.columns.tolist() if x != "variant"]
            if not args.include_metadata:
                columns = [
                    x for x in columns
                    if not x.startswith("metadata")
                ]
            if args.no_standard_columns:
                columns = [
                    x for x in columns
                    if x not in variants_util.STANDARD_DATAFRAME_COLUMNS
                ]

        df_save = df[columns].copy()
        df_save.interbase_start = df_save.interbase_start.astype(int)
        df_save.interbase_end = df_save.interbase_end.astype(int)

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
                    save(df)

    if args.chunk_rows is None:
        save(df)
 
