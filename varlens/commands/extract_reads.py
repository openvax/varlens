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

import pysam

from .. import loci
from .. import reads

parser = argparse.ArgumentParser(usage=__doc__)
loci.add_args(parser)
reads.add_args(parser)

parser.add_argument("--out", required=True)
parser.add_argument("-v", "--verbose", action="store_true", default=False)

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)

    sites = loci.load(args)
    if not sites:
        parser.error("No genomic loci (e.g. VCF files) specified.")

    read_sources = reads.load(args)

    if len(read_sources) != 1:
        parser.error("Exactly one read source must be specified.")

    (read_source,) = read_sources
        
    out_handle = pysam.AlignmentFile(
        args.out, "wb", template=read_source.handle)

    for locus in sites:
        pileup_collection = read_source.pileups([locus]).reads()
        for read in pileup_collection:
            out_handle.write(read)

    out_handle.close()
    print("Wrote: %s" % args.out)

