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

import os
import sys
import collections
import re
import json
from future.utils import raise_

import typechecks
import pandas
import numpy

try:
    # Python 2
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

try:
    # Python 2
    from urlparse import parse_qsl
except ImportError:
    from urllib.parse import parse_qsl

STANDARD_EVALUATION_ENVIRONMENT = {
    "os": os,
    "sys": sys,
    "collections": collections,
    "re": re,
    "json": json,
    "pandas": pandas,
    "numpy": numpy,
}

def parse_url_fragment(url):
    parsed = urlparse(url)
    try:
        if parsed.fragment:
            # If our fragment begins with an '&' symbol, we ignore it. 
            unparsed_fragment = parsed.fragment
            if unparsed_fragment.startswith("&"):
                unparsed_fragment = unparsed_fragment[1:]
            fragment = parse_qsl(unparsed_fragment, strict_parsing=True)
        else:
            fragment = []
    except ValueError as e:
        raise ValueError("Couldn't parse fragment '%s': %s" % (
            parsed.fragment, e))

    return (parsed._replace(fragment='').geturl(), fragment)

def string_to_boolean(s):
    value = str(s).lower()
    if value in ("true", "1"):
        return True
    elif value in ("false", 0):
        return False
    raise ValueError("Not a boolean string: %s" % s)

class EvaluationEnvironment(object):
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

    def __str__(self):
        all_bindings = []
        for wrapped in self._wrapped_list:
            all_bindings.extend(dir(wrapped))
        all_bindings.extend(self._extra)
        display_bindings = [
            binding for binding in all_bindings
            if not binding.startswith("_")
        ]
        return ("<Environment with bindings: %s>"
            % " ".join(sorted(display_bindings)))

RAISE = object()

def evaluate_expression(expression, bindings, error_value=RAISE):
    typechecks.require_string(expression)

    # Since Python 2 doesn't have a nonlocal keyword, we have to box up the
    # error_value, so we can reassign to it in the ``on_error`` function
    # below.
    error_box = [error_value] 
    try:
        # Give some basic modules.
        standard_environment = dict(STANDARD_EVALUATION_ENVIRONMENT)
    
        # Add our "on_error" hack.
        def on_error(value):
            error_box[0] = value
        standard_environment["on_error"] = on_error

        return eval(expression, standard_environment, bindings)
    except Exception as e:
        if error_box[0] is not RAISE:
            return error_box[0]
        extra = "Error while evaluating: \n\t%s\non:\n%s" % (
            expression, bindings)
        traceback = sys.exc_info()[2]
        raise_(ValueError, str(e) + "\n" + extra, traceback)