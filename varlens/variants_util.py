# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import collections

import pandas
import varcode
import varcode.reference
import logging

from . import util, loci_util
from .locus import Locus
from .read_evidence import pileup_collection

STANDARD_DATAFRAME_COLUMNS = [
    "genome",
    "contig",
    "interbase_start",
    "interbase_end",
    "ref",
    "alt",
]

def add_args(parser, positional=False):
    group = parser.add_argument_group("variant loading")
    group.add_argument("variants" if positional else "--variants",
        nargs=("*" if positional else "+"), default=[],
        help="Path to VCF file. Any number of VCF files may be specified. "
        "CSV files in the format outputted by the varlens-variants tool are "
        "also supported.")
    group.add_argument("--genome",
        help="Genome for the variants (e.g. b37). Required when the genome "
        "cannot be guessed from the metadata in the VCF.")
    group.add_argument("--include-failing-variants",
        action="store_true",
        default=False,
        help="Include variants with a non-PASS filter field.")
    group.add_argument("--variant-source-name", nargs="+",
        help="Names for variant sources. Must specify one name per variant "
        "source. If not specified, the filenames are used.")
    group.add_argument("--max-variants-per-source", type=int,
        metavar="N",
        help="Load at most N variants from each source.")
    group.add_argument("--single-variant", nargs=3, action="append",
        default=[], metavar="X",
        help="Literal variant specified as three arguments: LOCUS REF ALT. "
        "Can be specified any number of times by repeating the "
        "--single-variant option.")

    # Filters
    group = parser.add_argument_group("variant filtering",
        "If multiple filters are specified, the variants must pass *all* "
        "filters. For filtering by loci, any variants that overlap "
        "the specified loci are included.")
    group.add_argument("--ref", nargs="+",
        help="Include only variants where ref is one of the given values.")
    group.add_argument("--alt", nargs="+",
        help="Include only variants where alt is one of the given values.")
    loci_util.add_args(util.PrefixedArgumentParser(group, "variant"))

def load_from_args_as_dataframe(args):
    '''
    Given parsed variant-loading arguments, return a pandas DataFrame.

    If no variant loading arguments are specified, return None.
    '''
    if not args.variants and not args.single_variant:
        return None

    if args.variant_source_name:
        variant_source_names = util.expand(
            args.variant_source_name,
            'variant_source_name',
            'variant source',
            len(args.variants))
    else:
        variant_source_names = util.drop_prefix(args.variants)

    variant_to_sources = collections.defaultdict(list)

    dfs = []
    for i in range(len(args.variants)):
        name = variant_source_names[i]
        prefix = (
            'metadata:' if len(args.variants) == 1 else "metadata:%s:" % name)
        df = load_as_dataframe(
            args.variants[i],
            name=name,
            genome=args.genome,
            max_variants=args.max_variants_per_source,
            only_passing=not args.include_failing_variants,
            metadata_column_prefix=prefix)

        if df.shape[0] == 0:
            logging.warn("No variants loaded from: %s" % args.variants[i])
        else:
            for variant in df.variant:
                variant_to_sources[variant].append(name)
            dfs.append(df)

    if args.single_variant:
        variants = []
        extra_args = {}
        if args.genome:
            extra_args = {
                'ensembl': varcode.reference.infer_genome(args.genome)
            }
        for (locus_str, ref, alt) in args.single_variant:
            locus = Locus.parse(locus_str)
            variant = varcode.Variant(
                    locus.contig,
                    locus.inclusive_start,
                    ref,
                    alt,
                    **extra_args)
            variants.append(variant)
            variant_to_sources[variant].append("commandline")
        dfs.append(variants_to_dataframe(variants))

    df = dfs.pop(0)
    for other_df in dfs:
        df = pandas.merge(
            df,
            other_df,
            how='outer',
            on=["variant"] + STANDARD_DATAFRAME_COLUMNS)

    genomes = df["genome"].unique()
    if len(genomes) > 1:
        raise ValueError(
                "Mixing references is not supported. "
                "Reference genomes: %s" % (", ".join(genomes)))

    df["sources"] = [" ".join(variant_to_sources[v]) for v in df.variant]

    # Apply filters:
    if args.ref:
        df = df.ix[df.ref.isin(args.ref)]
    if args.alt:
        df = df.ix[df.alt.isin(args.alt)]
    loci = loci_util.load_from_args(
        util.remove_prefix_from_parsed_args(args, "variant"))
    if loci is not None:
        df = df.ix[[
            loci.intersects(pileup_collection.to_locus(v))
            for v in df.variant
        ]]
    return df

def load_as_dataframe(
        filename,
        loader=None,
        name=None,
        genome=None,
        max_variants=None,
        only_passing=True,
        metadata_column_prefix=''):

    if name is None:
        name = filename

    if loader is None:
        if (filename.endswith(".vcf") or filename.endswith(".vcf.gz")):
            # Load from VCF
            def loader(filename):
                collection = varcode.load_vcf_fast(
                    filename,
                    genome=genome,
                    max_variants=max_variants,
                    only_passing=only_passing,
                    allow_extended_nucleotides=True)
                return variants_to_dataframe(
                    collection,
                    collection.metadata,
                    metadata_column_prefix=metadata_column_prefix)

        elif (filename.endswith(".csv") or filename.endswith(".csv.gz")):
            # Load from csv
            def loader(filename):
                # Ignores only_passing
                df = pandas.read_csv(filename, nrows=max_variants)
                for column in ['ref', 'alt']:
                    df[column] = df[column].fillna('')
                df["variant"] = [
                    dataframe_row_to_variant(row) for (i, row) in df.iterrows()
                ]
                return df
        else:
            raise ValueError(
                "Unsupported input file extension for variants: %s" % filename)

    df = loader(filename)

    if 'genome' not in df:
        df["genome"] = genome

    df["variant"] = [
        dataframe_row_to_variant(row) for (i, row) in df.iterrows()
    ]
    return df

def variants_to_dataframe(
        variants, metadata=None, metadata_column_prefix=""):
    def record(variant):
        d = {
            'variant': variant,
            'genome': str(variant.reference_name),
            'contig': variant.contig,
            'interbase_start': variant.start - 1,
            'interbase_end': variant.end,
            'ref': variant.ref,
            'alt': variant.alt,
        }
        if metadata:
            for (name, value) in metadata.get(variant, {}).items():
                if name == 'info':
                    for (info_col, value) in value.items():
                        column = '%sinfo:%s' % (
                            metadata_column_prefix, info_col)
                        d[column] = value
                else:
                    d["%s%s" % (metadata_column_prefix, name.lower())] = value
        return d

    df = pandas.DataFrame.from_records([record(v) for v in variants])
    column_indices = dict(
        (column, i) for (i, column) in enumerate(STANDARD_DATAFRAME_COLUMNS))
    columns = sorted(df.columns, key=lambda col: column_indices.get(col, 100))
    return df[columns]

def dataframe_row_to_variant(row):
    return varcode.Variant(
            ensembl=row.genome,
            contig=row.contig,
            start=row.interbase_start + 1,
            ref=row.ref,
            alt=row.alt,
            allow_extended_nucleotides=True)

def dataframe_to_variants(df):
    for column in STANDARD_DATAFRAME_COLUMNS:
        if column not in df:
            raise ValueError("Missing column: %s" % column)

    extra_columns = [
        c for c in df.columns if c not in STANDARD_DATAFRAME_COLUMNS
    ]
    metadata = collections.OrderedDict()
    for (i, row) in df.iterrows():
        variant = dataframe_row_to_variant(row)
        # We ignore the interbase_end field.
        metadata[variant] = dict((c, row[c]) for c in extra_columns)

    return varcode.VariantCollection(metadata.keys(), metadata=metadata)

def load_csv(filename, genome=None):
    # Genome is ignored for now.
    df = pandas.read_csv(filename)
    return dataframe_to_variants(df)
