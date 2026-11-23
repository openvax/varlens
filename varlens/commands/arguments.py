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

"""Command-line options that every varlens command defines the same way.

Options about variants, reads and loci are added by ``variants_util.add_args``,
``reads_util.add_args`` and ``loci_util.add_args``.
"""

from ..version import __version__


def add_version_arg(parser):
    """Add ``--version``, which prints the program name and version, then exits.

    Parameters
    ----------
    parser : argparse.ArgumentParser or argument group
        Where to add the option.

    Returns
    -------
    argparse.Action
        The added option.
    """
    return parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__)


def add_verbose_arg(parser):
    """Add ``-v/--verbose``, which turns on debug logging.

    ``varlens.commands.configure_logging`` reads the parsed ``verbose`` value.

    Parameters
    ----------
    parser : argparse.ArgumentParser or argument group
        Where to add the option.

    Returns
    -------
    argparse.Action
        The added option. Its value is False unless the flag is given.
    """
    return parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Log debug messages to stderr.")
