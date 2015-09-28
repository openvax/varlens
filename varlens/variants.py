import collections

import typechecks
import pandas
import varcode

from . import util
from . import evaluation

def add_args(parser):
    parser.add_argument("--variants", action="append", default=[],
        help="Path to VCF file. Can be specified multiple times.")
    parser.add_argument("--variant-filter")
    parser.add_argument("--variant-genome")

def load(args):
    '''
    Given parsed variant-loading arguments, return a varcode.VariantCollection.
    '''
    variant_collections = [
        load_vcf(
            filename, filter=args.variant_filter, genome=args.variant_genome)
        for filename in args.variants
    ]
    if not variant_collections:
        return varcode.VariantCollection([])
    
    if len(variant_collections) == 1:
        return variant_collections[0]

    variants = set()
    metadata = {}
    for collection in variant_collections:
        variants.update(collection)
        metadata.update(collection.metadata)

    return varcode.VariantCollection(variants, metadata=metadata)

def load_vcf(url, filter=None, loader=varcode.load_vcf_fast, **kwargs):
    (url_without_fragment, fragment) = util.parse_url_fragment(url)
    filters = []
    params = {}
    collection_metadata = {}
    for (key, value) in fragment:
        if key == 'filter':
            filters.append(value)
        elif key.startswith('metadata.'):
            metadata_key = key[len("metadata."):]
            collection_metadata[metadata_key] = value
        elif key in ("only_passing", "allow_extended_nucleotides"):
            params[key] = util.string_to_boolean(value)
        elif key in ("genome", "reference_name", "reference_vcf_key"):
            params[key] = value
        elif key == "max_variants":
            params[key] = int(value)
        else:
            raise ValueError("Unsupported operation: %s" % key)
    
    # kwargs override URL.
    for (key, value) in kwargs.items():
        if value is not None:
            params[key] = value

    # passed in filter is run after filters specified in the URL.
    if filter:
        filters.append(filter)

    result = loader(url_without_fragment, **params)
    if collection_metadata:
        for variant in result:
            result.metadata[variant].update(collection_metadata)
    for expression in filters:
        result = result.filter(
            lambda variant:
                evaluate_variant_expression(expression, result, variant))
    return result

def evaluate_variant_expression(
        expression,
        collection,
        variant,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        bindings = evaluation.AttributeToKeyWrapper([variant], extra={
            'inclusive_start': variant.start,
            'inclusive_end': variant.end,
            'interbase_start': variant.start - 1,
            'interbase_end': variant.end,
            'variant': variant,
            'collection': collection,
            'metadata': collection.metadata.get(variant),
        })
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(variant)  

STANDARD_DATAFRAME_COLUMNS = [
    "genome",
    "contig",
    "interbase_start",
    "interbase_end",
    "ref",
    "alt",
]

def variants_to_dataframe(variant_collection, extra_columns={}):
    columns = collections.OrderedDict(
        (name, []) for name in STANDARD_DATAFRAME_COLUMNS)
    columns.extend(extra_columns.keys())
    for variant in variant_collection:
        columns["genome"].append(str(variant.genome))
        columns["contig"].append(variant.contig)
        columns["interbase_start"].append(variant.start - 1)
        columns["interbase_end"].append(variant.end)
        columns["ref"].append(variant.ref)
        columns["alt"].append(variant.alt)
        for (column, callable_or_expression) in extra_columns.items():
            columns[column].append(evaluate_variant_expression(
                callable_or_expression, variant_collection, variant))

    return pandas.DataFrame(columns, index=list(variant_collection))

