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

import pysam

from . import configure_logging
from .. import loci_util
from .. import reads_util
from .. import variants_util
from ..read_evidence.pileup_collection import PileupCollection, to_locus

parser = argparse.ArgumentParser(usage=__doc__)
reads_util.add_args(parser, positional=True)
loci_util.add_args(parser)
variants_util.add_args(parser)

parser.add_argument("--out")
parser.add_argument("--field", nargs="+", default=[],
    help="Possible fields include: %s" % (
        " ".join(PileupCollection._READ_ATTRIBUTE_NAMES)))

parser.add_argument("--no-standard-fields", action="store_true", default=False)
parser.add_argument("--no-sort", action="store_true", default=False)
parser.add_argument(
    "--header",
    action="store_true",
    default=False,
    help="Output BAM/SAM header only.")
parser.add_argument(
    "--header-set",
    nargs=4,
    action="append",
    help="Example --header-set RG . SM my_sample")

parser.add_argument("-v", "--verbose", action="store_true", default=False)

STANDARD_FIELDS = [
    "query_name",
    "reference_start",
    "reference_end",
    "cigarstring",
]

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)
    configure_logging(args)

    read_sources = reads_util.load_from_args(args)
    if not read_sources:
        parser.error("No read sources specified.")

    loci = loci_util.load_from_args(args)  # may be None
    variants_df = variants_util.load_from_args_as_dataframe(args)
    if variants_df is not None:
        variant_loci = loci_util.Loci(
            to_locus(variant)
            for variant in variants_df["variant"])
        loci = variant_loci if loci is None else loci.union(variant_loci)

    if args.header:
        if loci is not None:
            parser.error("If specifying --header don't specify loci.")
        if args.field:
            parser.error("If specifying --header don't specify fields.")

    out_pysam_handle = None
    out_csv_writer = out_csv_fd = None
    if args.out and (args.out.endswith(".bam") or args.out.endswith(".sam")):
        if args.field:
            parser.error("Don't specify fields when outputting to bam or sam.")

        header = update_header(args, read_sources[0].handle.header)
        out_pysam_handle = pysam.AlignmentFile(
            args.out,
            "wb" if args.out.endswith(".bam") else "w",
            header=header)

    elif not args.out or args.out.endswith(".csv"):
        out_csv_fd = open(args.out, "w") if args.out else sys.stdout
        out_csv_writer = csv.writer(out_csv_fd)

        if args.header:
            if args.field:
                parser.error("Don't specify fields when outputting header.")
            out_csv_writer.writerow([
                "read_source", "group", "index", "key", "value",
            ])
        else:
            columns = (
                ([] if args.no_standard_fields else STANDARD_FIELDS) +
                args.field)
            out_csv_writer.writerow(columns)
    else:
        parser.error(
            "Don't know how to write to file with output extension: %s. "
            "Supported extensions: csv, bam, sam." % args.out)

    num_reads = 0
    for read_source in read_sources:
        if args.header:
            header = update_header(args, read_source.handle.header)
            for (group, i, key, value) in reads_util.flatten_header(header):
                out_csv_writer.writerow(
                    [read_source.name, group, str(i), key, value])
            continue  # we don't look at reads at all.
        for read in read_source.reads(loci):
            num_reads += 1
            if out_pysam_handle is not None:
                out_pysam_handle.write(read)
            if out_csv_writer is not None:
                out_csv_writer.writerow([
                    str(read_field(read, field)) for field in columns
                ])

    if out_pysam_handle is not None:
        out_pysam_handle.close()
        if not args.no_sort:
            print("Sorting read file %s" % args.out)
            pysam.sort(
                "-o", args.out,
                "-T", "varlens_reads", args.out,
                catch_stdout=False)
        print("Wrote %d reads: %s" % (num_reads, args.out))

    if out_csv_fd is not None and out_csv_fd is not sys.stdout:
        out_csv_fd.close()
        print("Wrote: %s" % args.out)


def read_field(read, field_name):
    if field_name.startswith("tag:"):
        tag_name = field_name[len("tag:"):]
        return read.get_tags().get(tag_name)

    try:
        return getattr(read, field_name)
    except AttributeError:
        raise ValueError("Invalid read field '%s'. Valid fields include: %s"
            % (field_name, ' '.join(dir(read))))

def update_header(args, header):
    if args.header_set:
        header = dict(header)
        for (group, index_string, key, value) in args.header_set:
            if not isinstance(header[group], list):
                header[group] = [header[group]]
            if index_string == ".":
                indices = range(len(header[group]))
            else:
                indices = [int(x) for x in index_string.split(",")]
            for index in indices:
                header[group][index][key] = value
    return header
