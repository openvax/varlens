import typechecks

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
        bindings = eval.AttributeToKeyWrapper(variant, extra={
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

