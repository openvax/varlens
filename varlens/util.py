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
import argparse

def expand(value, arg_name, input_name, length):
    if value is None or len(value) == 0:
        return [None] * length

    if len(value) == length:
        return value

    if len(value) == 1:
        return value * length

    if length == 1:
        raise ValueError(
            "With only 1 {input_name} specified, each {arg_name} argument "
            "should be length 1. If you are trying to specify multiple filters"
            " to apply consecutively, you should specify the entire argument "
            "multiple times."
            .format(
                arg_name=arg_name,
                input_name=input_name,
                length=length,
                actual=len(value)))

    else:
        raise ValueError(
            "Expected argument {arg_name} to be length 1 (i.e. apply to all "
            "{input_name} inputs) or length {length} (i.e. an individual value"
            " for each of the {length} {input_name} inputs), not {actual}."
            .format(
                arg_name=arg_name,
                input_name=input_name,
                length=length,
                actual=len(value)))


def drop_prefix(strings):
    """
    Removes common prefix from a collection of strings
    """
    if len(strings) == 1:
        return [os.path.basename(strings[0])]
    prefix_len = len(os.path.commonprefix(strings))
    return [string[prefix_len:] for string in strings]

class PrefixedArgumentParser(object):
    def __init__(self, wrapped, prefix):
        self.wrapped = wrapped
        self.prefix = prefix

    def add_argument(self, name, *args, **kwargs):
        assert name.startswith("--")
        new_name = "--" + self.prefix + "-" + name[2:]
        self.wrapped.add_argument(new_name, *args, **kwargs)


def remove_prefix_from_parsed_args(args, prefix):
    result = argparse.Namespace()
    for (arg, value) in args._get_kwargs():
        if arg.startswith(prefix + "_"):
            setattr(result, arg[len(prefix + "_"):], value)
    return result



