import os
import sys
import collections
import re
import json
from future.utils import raise_

import typechecks

STANDARD_EVALUATION_ENVIRONMENT = {
    "os": os,
    "sys": sys,
    "collections": collections,
    "re": re,
    "json": json,
}

class AttributeToKeyWrapper(object):
    def __init__(self, wrapped_list, extra={}):
        self._wrapped_list = wrapped_list
        self._extra = extra

    def __getitem__(self, key):
        for wrapped in self._wrapped_list:
            if hasattr(wrapped, key):
                return getattr(wrapped, key)
        if key in self._extra:
            return self._extra[key]
        raise KeyError("No key: %s" % key)

RAISE = object()

def evaluate_expression(expression, bindings, error_value=RAISE):
    typechecks.require_string(expression)

    # Since Python 2 doesn't have a nonlocal keyword, we have to box up the
    # error_value, so we can reassign to it in the ``on_error`` function
    # below.
    error_box = [error_value] 
    try:
        # Give some basic modules.
        environment = dict(STANDARD_EVALUATION_ENVIRONMENT)
    
        # Add our "on_error" hack.
        def on_error(value):
            error_box[0] = value
        environment["on_error"] = on_error

        return eval(expression, environment, bindings)
    except Exception as e:
        if error_box[0] is not RAISE:
            return error_box[0]
        extra = "Error while evaluating: \n\t%s\non:\n%s" % (
            expression, bindings)
        traceback = sys.exc_info()[2]
        raise_(ValueError, str(e) + "\n" + extra, traceback)