'''
Filter reads from one or more BAMs and output a CSV or a new BAM.

Loci and VCF files may be specified, in which case reads are filtered to
overlap the specified loci or variants.

Examples:

Print basic fields for the reads in a BAM:

    %(prog)s test/data/CELSR1/bams/bam_0.bam

Same as above but filter only to reads aligned on the (-) strand, write to a 
file instead of stdout, and also include the mapping quality and sequenced
bases in the output:

    %(prog)s test/data/CELSR1/bams/bam_0.bam \\
        --is-reverse \\
        --field mapping_quality query_alignment_sequence \\
        --out /tmp/result.csv

Write a bam file consisting of reads with mapping quality >=30 and
overlapping a certain locus:

    %(prog)s test/data/CELSR1/bams/bam_0.bam \\
        --min-mapping-quality 30 \\
        --locus 22:46932040-46932050 \\
        --out /tmp/result.bam

Write a bam file consisting of reads overlapping variants from a VCF:

    %(prog)s test/data/CELSR1/bams/bam_0.bam \\
        --variants test/data/CELSR1/vcfs/vcf_1.vcf \\
        --out /tmp/result.bam

Print just the header for a BAM in csv format:

    %(prog)s test/data/CELSR1/bams/bam_0.bam --header

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

STANDARD_FIELDS = [
    "query_name",
    "reference_start",
    "reference_end",
    "cigarstring",
]

parser = argparse.ArgumentParser(usage=__doc__)
group = parser.add_argument_group("output")
group.add_argument("--out",
    help="Output file. Format is guessed from file extension: must be csv or "
    "bam. If not specified, csv is written to stdout.")
group.add_argument("--field", nargs="+", default=[],
    help="Additional read fields to output as columns in the csv. See pysam "
    "documentation (http://pysam.readthedocs.org/en/latest/api.html) for the "
    "meaning of these fields. Valid fields include: %s" % (
        " ".join(PileupCollection._READ_ATTRIBUTE_NAMES)))

group.add_argument("--no-standard-fields", action="store_true", default=False,
    help="Do not include the standard fields (%s) in csv output."
    % ', '.join(STANDARD_FIELDS))
group.add_argument("--no-sort", action="store_true", default=False,
    help="When outputting a bam, do not call samtools sort.")
group.add_argument(
    "--header",
    action="store_true",
    default=False,
    help="Output BAM/SAM header only.")
group.add_argument(
    "--header-set",
    nargs=4,
    action="append",
    help="When outputting a bam, set a particular header field to the given "
    "value. Example: --header-set RG . SM my_sample")

group.add_argument("-v", "--verbose", action="store_true", default=False)

reads_util.add_args(parser, positional=True)
loci_util.add_args(parser.add_argument_group("loci specification"))
variants_util.add_args(parser)

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
