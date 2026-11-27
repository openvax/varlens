# PYTHON_ARGCOMPLETE_OK
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

"""Command-line manipulation of genomic variants, NGS reads, and fusion RNA."""

from __future__ import annotations

from typing import Any

from . import allele_support, fusions, reads, variants
from .arguments import add_version_arg
from .help import HelpArgumentParser, enable_argcomplete
from ..errors import InputError

SUBCOMMANDS: tuple[tuple[str, str, Any, bool], ...] = (
    (
        "variants",
        "Combine, filter, and annotate variants from VCF/CSV files",
        variants,
        True,
    ),
    (
        "reads",
        "Display, filter, and copy reads from BAM/SAM files",
        reads,
        True,
    ),
    (
        "allele-support",
        "Count reads supporting each allele at genomic sites",
        allele_support,
        True,
    ),
    (
        "fusions",
        "Reconstruct supplied fusion RNA and predict junction peptides",
        fusions,
        False,
    ),
)


def make_parser(prog: str = "varlens") -> HelpArgumentParser:
    """Build the unified ``varlens`` CLI argument parser."""
    root = HelpArgumentParser(
        prog=prog,
        description=__doc__,
        epilog="Run 'varlens COMMAND --help' for command-specific options.",
    )
    add_version_arg(root)
    subparsers = root.add_subparsers(
        title="subcommands",
        dest="subcommand",
        required=True,
        metavar="COMMAND",
    )
    for name, summary, module, doc_is_usage in SUBCOMMANDS:
        subparser_kwargs: dict[str, Any] = {"help": summary}
        if doc_is_usage:
            subparser_kwargs["usage"] = module.__doc__.lstrip()
        else:
            subparser_kwargs["description"] = module.__doc__
        subparser = subparsers.add_parser(name, **subparser_kwargs)
        module.configure_parser(subparser)
        subparser.set_defaults(_handler=module.run_parsed, _parser=subparser)
    return root


parser = make_parser()


def run(raw_args: list[str] | None = None) -> None:
    """Parse arguments and dispatch to the selected ``varlens`` subcommand."""
    enable_argcomplete(parser)
    args = parser.parse_args(raw_args)
    subparser = args._parser
    try:
        args._handler(args, subparser)
    except InputError as error:
        subparser.error(str(error))


if __name__ == "__main__":
    run()
