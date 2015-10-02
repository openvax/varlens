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
import collections

from .. import loci_util
from .. import reads_util
from . import configure_logging
from ..evaluation import parse_labeled_expression
from .. import read_evidence

parser = argparse.ArgumentParser(usage=__doc__)
loci_util.add_args(parser)
reads_util.add_args(parser)

parser.add_argument("--out")
parser.add_argument("field", nargs="*")
parser.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)
    configure_logging(args)
    
    loci = loci_util.load_from_args(args)
    if not loci:
        parser.error("No genomic loci (e.g. VCF files) specified.")

    logging.info("Loaded %d genomic loci." % len(loci))

    read_sources = reads_util.load_from_args(args)

    out_fd = open(args.out, "w") if args.out else sys.stdout
    writer = csv.writer(out_fd)

    extra_columns = collections.OrderedDict()
    for labeled_expression in args.field:
        (label, expression) = parse_labeled_expression(labeled_expression)
        extra_columns[label] = expression

    writer.writerow([
        "source",
        "contig",
        "interbase_start",
        "interbase_end",
        "allele",
        "count",
    ] + extra_columns.keys())
    for source in read_sources:
        logging.info("Reading from: %s" % source.name)
        for locus in loci:
            grouped = dict(source.pileups([locus]).group_by_allele(locus))
            for (allele, group) in grouped.items():
                extra_values = []
                for expression in extra_columns.values():
                    num_reads = len(set(
                        read_evidence.read_key(element.alignment)
                        for pileup in group.pileups.values()
                        for element in pileup
                        if reads_util.evaluate_pileup_element_expression(
                            expression,
                            group,
                            pileup,
                            element)))
                    extra_values.append(str(num_reads))
                writer.writerow([
                    source.name,
                    locus.contig,
                    str(locus.start),
                    str(locus.end),
                    allele,
                    str(group.num_reads()),
                ] + extra_values)

    if out_fd is not sys.stdout:
        out_fd.close()
        print("Wrote: %s" % args.out)
