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
                input_name=input_name))

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
    strings_without_extensions = [
        s.split(".", 2)[0] for s in strings
    ]

    if len(strings_without_extensions) == 1:
        return [os.path.basename(strings_without_extensions[0])]
    prefix_len = len(os.path.commonprefix(strings_without_extensions))
    result = [string[prefix_len:] for string in strings_without_extensions]
    if len(set(result)) != len(strings):
        # If these operations resulted in a collision, just return the original
        # strings.
        return strings
    return result


def name_sources(paths, names, arg_name, input_name):
    """Return the name to show for each input file (e.g. each BAM or VCF).

    Used for ``--read-source-name`` and ``--variant-source-name``.

    Parameters
    ----------
    paths : list of str
        Input file paths.
    names : list of str or None
        Names given on the command line. A single name is used for every
        path; otherwise there must be one name per path. If empty or None,
        names are derived from the paths with ``drop_prefix``.
    arg_name : str
        Option name to show in the error message, e.g. ``"read_source_name"``.
    input_name : str
        Kind of input to show in the error message, e.g. ``"read source"``.

    Returns
    -------
    list of str
        One name per path, in the same order.

    Raises
    ------
    ValueError
        If the number of names is neither 1 nor ``len(paths)``.
    """
    if names:
        return expand(names, arg_name, input_name, len(paths))
    return drop_prefix(paths)


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



