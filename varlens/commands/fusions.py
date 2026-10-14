"""Reconstruct supplied fusion RNA and optionally predict junction peptides."""

import argparse
import gzip
import json
from pathlib import Path
import sys

from ..errors import InputError, input_errors
from ..fusions import analyze_fusion


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="+", help="Isovar supplied-fusion JSON or JSON.gz files")
parser.add_argument("--out", help="JSON/JSON.gz report; defaults to stdout")
parser.add_argument("--hla", nargs="+", default=[],
                    help="MHC alleles; enables NetMHCpan via mhctools")
parser.add_argument("--peptide-lengths", nargs="+", type=int, default=[8, 9, 10, 11])
parser.add_argument("--min-fragments", type=int, default=2,
                    help="Minimum direct junction fragments required by Isovar (default: 2)")


@input_errors(parser)
def run(raw_args=None):
    args = parser.parse_args(raw_args)
    if args.min_fragments < 1 or any(n < 1 for n in args.peptide_lengths):
        parser.error("--min-fragments and --peptide-lengths must be positive")
    if args.out and not args.out.endswith((".json", ".json.gz")):
        parser.error("--out must end in .json or .json.gz")
    events = []
    for filename in args.inputs:
        try:
            if not filename.endswith((".json", ".json.gz")):
                raise InputError("Expected a supplied-fusion .json or .json.gz file")
            opener = gzip.open if filename.endswith(".gz") else open
            with opener(filename, "rt") as stream:
                payload = json.load(stream)
            events.append(analyze_fusion(
                payload, alleles=args.hla, epitope_lengths=args.peptide_lengths,
                min_fragments=args.min_fragments))
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
