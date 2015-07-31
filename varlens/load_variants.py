import os
import sys
import re
import json
import collections
import typechecks

from future.utils import raise_

import varcode

from . import util

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
        elif key in ("ensembl_version", "reference_name", "reference_vcf_key"):
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
    #result.collection_metadata.update(collection_metadata)
    for f in filters:
        result = filter_variants(result, f)
    return result

RAISE = object()
def filter_variants(
        collection,
        expression,
        error_value=RAISE,
        extra_bindings={}):

    return collection.filter(
        lambda variant: evaluate_expression(
            expression,
            collection,
            variant,
            error_value=error_value,
            extra_bindings=extra_bindings))


STANDARD_EVALUATION_ENVIRONMENT = {
    "os": os,
    "sys": sys,
    "collections": collections,
    "re": re,
    "json": json,
}

class VariantWrapper(object):
    def __init__(self, variant):
        self.variant = variant

    def __getitem__(self, key):
        if hasattr(self.variant, key):
            return getattr(self.variant, key)
        raise KeyError("No key: %s" % key)

def evaluate_expression(
        expression,
        collection,
        variant,
        error_value=RAISE,
        extra_bindings={}):

    # Since Python 2 doesn't have a nonlocal keyword, we have to box up the
    # error_value, so we can reassign to it in the ``on_error`` function
    # below.
    error_box = [error_value] 
    try:
        if typechecks.is_string(expression):
            # Give some basic modules.
            environment = dict(STANDARD_EVALUATION_ENVIRONMENT)
            environment["variant"] = variant
            environment["collection"] = collection
            environment["metadata"] = collection.metadata.get(variant)

            # We also add our "on_error" hack.
            def on_error(value):
                error_box[0] = value
            environment["on_error"] = on_error
            environment.update(extra_bindings)

            variant_wrapper = VariantWrapper(variant)
            return eval(expression, environment, variant_wrapper)
        else:
            return expression(variant)                
    except Exception as e:
        if error_box[0] is not RAISE:
            return error_box[0]
        extra = "Error while evaluating: \n\t%s\non variant:\n%s" % (
            expression, variant)
        traceback = sys.exc_info()[2]
        raise_(ValueError, str(e) + "\n" + extra, traceback)
