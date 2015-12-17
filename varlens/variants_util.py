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

import typechecks
import pandas
import varcode
import varcode.reference

from . import evaluation, Locus

STANDARD_DATAFRAME_COLUMNS = [
    "genome",
    "contig",
    "interbase_start",
    "interbase_end",
    "ref",
    "alt",
]

def add_args(parser):
    parser.add_argument("--variants", action="append", default=[],
        help="Path to VCF file. Can be specified multiple times.")
    parser.add_argument("--variant-filter")
    parser.add_argument("--variant-genome")
    parser.add_argument("--distinct-variants",
        action="store_true", default=False)
    parser.add_argument("--single-variant", nargs=3, action="append",
        default=[], metavar=("LOCUS", "REF", "ALT"),
        help="Literal variant. Can be specified any number of times.")

def load_from_args_as_dataframe(args):
    '''
    Given parsed variant-loading arguments, return a pandas DataFrame.

    If no variant loading arguments are specified, return None.
    '''
    if not args.variants and not args.single_variant:
        return None

    dfs = [
        load_as_dataframe(
            filename, filter=args.variant_filter, genome=args.variant_genome)
        for filename in args.variants
    ]
    if args.single_variant:
        variants = []
        extra_args = {}
        if args.variant_genome:
            extra_args = {
                'ensembl': varcode.reference.infer_genome(args.variant_genome)
            }
        for (locus_str, ref, alt) in args.single_variant:
            locus = Locus.parse(locus_str)
            variants.append(
                varcode.Variant(
                    locus.contig,
                    locus.inclusive_start,
                    ref,
                    alt,
                    **extra_args))
        dfs.append(variants_to_dataframe(variants))

    df = pandas.concat(dfs)
    genomes = df["genome"].unique()
    if len(genomes) > 1:
        raise ValueError(
                "Mixing references is not supported. "
                "Reference genomes: %s" % (", ".join(genomes)))

    if args.distinct_variants:
        # We return a dataframe of only unique variants. Instead of the usual
        # metadata fields, we have only one called "variant_sources", which 
        # for each variant gives a dict from source name -> the metadata for
        # the variant that was read from that source.
        sources = (
            df.groupby(["variant"] + STANDARD_DATAFRAME_COLUMNS)
            .apply(lambda group: dict(iter(group.groupby("variant_source")))))
        sources.name = "variant_sources"
        df = sources.reset_index()
    return df

ParsedVariantsURL = collections.namedtuple(
    "ParsedVariantsURL",
    "url_without_fragment filters params collection_metadata name")

def parse_variants_url(url, extra_filter=None, extra_params={}):
    (url_without_fragment, fragment) = evaluation.parse_url_fragment(url)
    filters = []
    params = {}
    collection_metadata = {}
    name = url
    for (key, value) in fragment:
        if key == 'filter':
            filters.append(value)
        elif key.startswith('metadata.'):
            metadata_key = key[len("metadata."):]
            collection_metadata[metadata_key] = value
        elif key in ("only_passing", "allow_extended_nucleotides"):
            params[key] = evaluation.string_to_boolean(value)
        elif key in ("genome", "reference_name", "reference_vcf_key"):
            params[key] = value
        elif key == "max_variants":
            params[key] = int(value)
        elif key == "name":
            name = value
        else:
            raise ValueError("Unsupported operation: %s" % key)
    
    # extra params override URL.
    for (key, value) in extra_params.items():
        if value is not None:
            params[key] = value

    # passed in filter is run after filters specified in the URL.
    if extra_filter:
        filters.append(extra_filter)

    return ParsedVariantsURL(
        url_without_fragment,
        filters,
        params,
        collection_metadata,
        name)

def load_as_dataframe(url, filter=None, loader=None, **kwargs):
    parsed = parse_variants_url(url, extra_filter=filter, extra_params=kwargs)
    
    if loader is None:
        if (parsed.url_without_fragment.endswith(".vcf") or
                parsed.url_without_fragment.endswith(".vcf.gz")):
            
            # Load from VCF
            def loader(filename, **kwargs):
                collection = varcode.load_vcf_fast(filename, **kwargs)
                return variants_to_dataframe(collection, collection.metadata)

        elif (parsed.url_without_fragment.endswith(".csv") or
                parsed.url_without_fragment.endswith(".csv.gz")):
            
            # Load from csv
            def loader(filename, genome=None, max_variants=None):
                df = pandas.read_csv(filename, nrows=max_variants)
                if genome is not None:
                    df["genome"] = genome
                df["variant"] = [
                    dataframe_row_to_variant(row) for (i, row) in df.iterrows()
                ]
                return df
        else:
            raise ValueError(
                "Unsupported input file extension for variants: %s"
                % parsed.url_without_fragment)

    df = loader(parsed.url_without_fragment, **parsed.params)
    df["variant_source"] = parsed.name

    if parsed.filters:
        df = df[[
            bool(all(
                evaluate_variant_expression(
                    expression, row.to_dict(), row.variant)
                for expression in parsed.filters))
            for (i, row) in df.iterrows()
        ]]

    return df

def evaluate_variant_expression(
        expression,
        metadata,
        variant,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        extra_bindings = {
            'inclusive_start': variant.start,
            'inclusive_end': variant.end,
            'interbase_start': variant.start - 1,
            'interbase_end': variant.end,
            'variant': variant,
            'metadata': metadata,
        }
        extra_bindings.update(metadata)
        bindings = evaluation.EvaluationEnvironment([variant], extra_bindings)
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(variant)  

def variants_to_dataframe(variants, metadata=None):
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
                d["metadata_%s" % name.lower()] = value
            if 'metadata_info' in d:
                for (info_col, value) in d['metadata_info'].items():
                    d['metadata_info_%s' % info_col] = value
                del d['metadata_info']
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
