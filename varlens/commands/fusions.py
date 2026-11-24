"""Reconstruct supplied fusion RNA and optionally predict junction peptides."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import sys

from .arguments import add_version_arg
from ..errors import InputError, input_errors
from ..fusions import analyze_fusion
from ..mhc_binding import (
    DEFAULT_EPITOPE_LENGTHS,
    add_mhc_cli_args,
    has_mhc_cli_args,
    mhc_predictor_from_args,
    peptide_lengths_from_args,
)


parser = argparse.ArgumentParser(description=__doc__)
add_version_arg(parser)
parser.add_argument("inputs", nargs="+", help="Isovar supplied-fusion JSON or JSON.gz files")
parser.add_argument("--out", help="JSON/JSON.gz report; defaults to stdout")
parser.add_argument("--min-fragments", type=int, default=2,
                    help="Minimum direct junction fragments required by Isovar (default: 2)")
add_mhc_cli_args(parser, include_peptide_lengths_alias=True)


@input_errors(parser)
def run(raw_args: list[str] | None = None) -> None:
    args = parser.parse_args(raw_args)
    if args.min_fragments < 1:
        parser.error("--min-fragments must be positive")
    if args.out and not args.out.endswith((".json", ".json.gz")):
        parser.error("--out must end in .json or .json.gz")
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


if __name__ == "__main__":
    run()
