# PYTHON_ARGCOMPLETE_OK
"""Reconstruct supplied fusion RNA and optionally predict junction peptides."""

from __future__ import annotations

import argparse
import sys

from .arguments import add_mhc_args, add_version_arg
from .help import HelpArgumentParser, enable_argcomplete
from ..errors import InputError, input_errors


def configure_parser(parser: HelpArgumentParser) -> HelpArgumentParser:
    """Add ``fusions`` arguments to ``parser``."""
    add_version_arg(parser)
    parser.add_argument(
        "inputs", nargs="+", help="Isovar supplied-fusion JSON or JSON.gz files"
    )
    parser.add_argument("--out", help="JSON/JSON.gz report; defaults to stdout")
    parser.add_argument(
        "--min-fragments",
        type=int,
        default=2,
        help="Minimum direct junction fragments required by Isovar (default: 2)",
    )
    add_mhc_args(parser, include_peptide_lengths_alias=True)
    return parser


def make_parser(prog: str | None = None) -> HelpArgumentParser:
    """Build the command-line argument parser for ``varlens fusions``."""
    parser = HelpArgumentParser(prog=prog, description=__doc__)
    return configure_parser(parser)


parser = make_parser()


def run_parsed(
    args: argparse.Namespace, active_parser: HelpArgumentParser = parser
) -> None:
    """Execute ``fusions`` using pre-parsed arguments."""
    if args.min_fragments < 1:
        active_parser.error("--min-fragments must be positive")
    if args.out and not args.out.endswith((".json", ".json.gz")):
        active_parser.error("--out must end in .json or .json.gz")

    import gzip
    import json
    from pathlib import Path

    from ..fusions import analyze_fusion
    from ..mhc_binding import (
        DEFAULT_EPITOPE_LENGTHS,
        has_mhc_cli_args,
        mhc_predictor_from_args,
        peptide_lengths_from_args,
    )
    try:
        requested_lengths = peptide_lengths_from_args(args)
        predictor = mhc_predictor_from_args(args) if has_mhc_cli_args(args) else None
    except ValueError as error:
        raise InputError(str(error)) from error

    if requested_lengths is not None:
        epitope_lengths = requested_lengths
    elif predictor is not None and getattr(predictor, "default_peptide_lengths", None):
        epitope_lengths = tuple(predictor.default_peptide_lengths)
    else:
        epitope_lengths = DEFAULT_EPITOPE_LENGTHS

    events = []
    for filename in args.inputs:
        try:
            if not filename.endswith((".json", ".json.gz")):
                raise InputError("Expected a supplied-fusion .json or .json.gz file")
            opener = gzip.open if filename.endswith(".gz") else open
            with opener(filename, "rt") as stream:
                payload = json.load(stream)
            events.append(analyze_fusion(
                payload, epitope_lengths=epitope_lengths,
                min_fragments=args.min_fragments, predictor=predictor))
        except (ValueError, OSError) as error:
            raise InputError("%s: %s" % (filename, error)) from error
    # Finish all processing before touching the destination: unsupported input
    # or a predictor failure must not leave a successful-looking partial report.
    report = json.dumps({"schema": "varlens.fusions.v1", "events": events},
                        indent=2, allow_nan=False) + "\n"
    if args.out:
        if args.out.endswith(".gz"):
            with gzip.open(args.out, "wt") as stream:
                stream.write(report)
        else:
            Path(args.out).write_text(report)
    else:
        sys.stdout.write(report)


@input_errors(parser)
def run(raw_args: list[str] | None = None) -> None:
    enable_argcomplete(parser)
    args = parser.parse_args(raw_args)
    run_parsed(args, parser)


if __name__ == "__main__":
    run()
