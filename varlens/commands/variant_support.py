'''
Generate a read evidence pie chart plot.

%(prog)s \
    --input-reads /path/to/bam1 "my bam 1"
    --input-reads /path/to/bam2 "my bam 2"
    --input-variants /path/tovcf.1 "mutect" "my bam 1"
    --input-variants /path/to/vcf2 "my vcf2"

'''
from __future__ import absolute_import

import argparse
import collections
import os
import logging
import csv
import sys

try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy

import varcode
import varcode.read_evidence

from .. import load_variants

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("--variant-sets", nargs="+", default=[])
parser.add_argument("--read-sets", nargs="+", default=[])
parser.add_argument("--variant-set-labels", nargs="+")
parser.add_argument("--read-set-labels", nargs="+")
parser.add_argument("--max-variants", type=int)
parser.add_argument("--associate",
    metavar=("VARIANT_SET", "READ_SET"),
    nargs=2,
    action="append",
    default=[])

parser.add_argument("--variant-filter")
parser.add_argument("--ensembl-version")

parser.add_argument("--neighboring-loci-offsets",
    nargs="+", type=int, default=[])
parser.add_argument("--extra-info-below-each-pie", nargs="+")
parser.add_argument("--sort-order",
    choices=("priority", "genomic"), default="priority")
    
parser.add_argument("--out-plot")
parser.add_argument("--out-evidence")

parser.add_argument("-v", "--verbose", action="store_true", default=False)

ReadInput = collections.namedtuple("ReadInput", "name path")
VariantInput = collections.namedtuple("VariantInput", "name path reads")

def drop_prefix(strings):
    if len(strings) == 1:
        return strings
    prefix_len = len(os.path.commonprefix(strings))
    return [string[prefix_len:] for string in strings]

def run(raw_args=sys.argv[1:]):
    args = parser.parse_args(raw_args)

    if not args.variant_set_labels:
        args.variant_set_labels = drop_prefix(args.variant_sets)

    if not args.read_set_labels:
        args.read_set_labels = drop_prefix(args.read_sets)

    read_inputs = [
        ReadInput(name, path)
        for (path, name) in zip(args.read_sets, args.read_set_labels)
    ]

    associations = collections.defaultdict(list)
    for (variant_set, read_set) in args.associate:
        associations[variant_set].extend([
            read_input for read_input in read_inputs
            if read_set in (read_input.name, read_input.path)
        ])

    variant_inputs = [
        VariantInput(
            name,
            path,
            associations.get(path, []) + associations.get(name, []))
        for (path, name) in zip(args.variant_sets, args.variant_set_labels)
    ]

    variant_to_inputs = load_variants_dict(
        variant_inputs,
        filter=args.variant_filter,
        ensembl_version=args.ensembl_version)

    print("Loaded %d variants." % len(variant_to_inputs))

    if args.max_variants and len(variant_to_inputs) > args.max_variants:
        subselected = numpy.random.choice(
            list(variant_to_inputs),
            args.max_variants,
            replace=False)
        variant_to_inputs = dict(
            (v, variant_to_inputs[v]) for v in subselected)
        print("Subselected to %d variants." % len(variant_to_inputs))

    loci = [
        varcode.read_evidence.pileup_collection.to_locus(variant)
        for variant in variant_to_inputs
    ]
    print("Collecting evidence.")
    evidence_out = collect_evidence(loci, read_inputs)

    assert evidence_out
    if args.out_evidence:
        write_evidence(args.out_evidence, evidence_out)

def write_evidence(filename, evidence):
    if filename.endswith(".pickle"):
        with open(filename, 'w') as fd:
            pickle.dump(evidence, fd, pickle.HIGHEST_PROTOCOL)
    elif filename.endswith(".csv"):
        with open(filename, 'w') as fd:
            writer = csv.writer(fd)
            writer.writerow(
                ["sample", "contig", "start", "end", "allele", "count"])
            for (i, sample) in enumerate(evidence):
                logging.info("Writing sample %d / %d" %
                    ((i + 1), len(evidence)))
                locus_to_evidence = evidence[sample]
                for locus in locus_to_evidence:
                    allele_to_evidence = locus_to_evidence[locus]
                    for (allele, count) in allele_to_evidence.items():
                        writer.writerow([
                            sample,
                            locus.contig,
                            locus.start,
                            locus.end,
                            allele,
                            str(count),
                        ])
    else:
        raise ValueError("Unsupported format: %s" % filename)
    print("Wrote: %s" % filename)

def collect_evidence(loci, read_inputs):
    '''
    Parameters
    ----------
    loci : list of varcode.Locus or varcode.Variant instances
    read_inputs : list of ReadInput
    Returns
    ----------
    dict : read label -> locus -> allele -> count
    '''
    print("Loci", loci)
    result = collections.OrderedDict()
    for (i, read_input) in enumerate(read_inputs):
        logging.info("Reading pileups for %d / %d: %s" % (
            (i + 1), len(read_inputs), read_input.path))
        evidence = read_evidence.PileupCollection.from_bam(
            read_input.path, loci)
        logging.info("Done loading evidence")
        result[read_input.name] = collections.OrderedDict()
        for (j, locus) in enumerate(loci):
            if j % 100 == 0:
                logging.info("Locus %d / %d" % ((j + 1), len(loci)))
            result[read_input.name][locus] = dict(
                evidence.allele_summary(locus))
    return result

def load_variants_dict(variant_inputs, filter=None, ensembl_version=None):
    '''

    Returns
    --------
    dict : Variant -> set(VariantInput) giving the inputs that have the variant
    '''
    result = collections.defaultdict(set)
    for (i, variant_input) in enumerate(variant_inputs):
        logging.info("Loading vcf %d / %d" % ((i + 1), len(variant_inputs)))
        vc = load_variants.load_vcf(
            variant_input.path,
            filter=filter,
            ensembl_version=ensembl_version)
        logging.info("Loaded %d variants from %s." %
            (len(vc), variant_input.path))
        for variant in vc:
            result[variant].add(variant_input.name)
    return result
