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

STANDARD_EVALUATION_ENVIRONMENT = {
    "os": os,
    "sys": sys,
    "collections": collections,
    "re": re,
    "json": json,
    "pandas": pandas,
    "numpy": numpy,
}

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
            try:
                return getattr(wrapped, key)
            except AttributeError:
                pass
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
        return ("<Environment with bindings: %s>" % " ".join(
            sorted(display_bindings)))

RAISE = object()


def evaluate_expression(expression, bindings, error_value=RAISE):
    typechecks.require_string(expression)
    if not expression:
        return True

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


def parse_labeled_expression(labeled_expression):
    match = re.match(r"^([\w ]+):(.*)$", labeled_expression)
    if match is None:
        label = expression = labeled_expression
    else:
        (label, expression) = match.groups()
    return (label, expression)
