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

"""Import-light command-line option builders shared across varlens commands.

This module uses only the Python standard library and ``varlens.version`` /
``varlens.util`` so building any CLI parser or running ``argcomplete`` tab
completion never imports heavy scientific or genomics packages (``pandas``,
``numpy``, ``pysam``, ``varcode``, ``pyensembl``, ``mhctools``, ``topiary``,
or ``pyfaidx``).
"""

from __future__ import annotations

import argparse
from collections.abc import Iterator
from typing import Any

from ..util import PrefixedArgumentParser
from ..version import __version__

READ_BOOLEAN_PROPERTIES: tuple[str, ...] = (
    "is_paired",
    "is_proper_pair",
    "is_qcfail",
    "is_read1",
    "is_read2",
    "is_reverse",
    "is_secondary",
    "is_unmapped",
    "mate_is_reverse",
    "mate_is_unmapped",
    "is_duplicate",
)

READ_STRING_PROPERTIES: tuple[str, ...] = (
    "cigarstring",
    "query_alignment_sequence",
    "query_name",
)

READ_INT_PROPERTIES: tuple[str, ...] = (
    "inferred_length",
    "mapping_quality",
    "query_alignment_length",
    "query_alignment_start",
    "query_length",
    "reference_length",
    "reference_start",
    "template_length",
)

READ_ATTRIBUTE_NAMES: tuple[str, ...] = (
    "cigarstring",
    "flag",
    "inferred_length",
    "is_duplicate",
    "is_paired",
    "is_proper_pair",
    "is_qcfail",
    "is_read1",
    "is_read2",
    "is_reverse",
    "is_secondary",
    "is_unmapped",
    "mapping_quality",
    "mate_is_reverse",
    "mate_is_unmapped",
    "next_reference_id",
    "next_reference_start",
    "query_alignment_end",
    "query_alignment_length",
    "query_alignment_qualities",
    "query_alignment_sequence",
    "query_alignment_start",
    "query_length",
    "query_name",
    "reference_end",
    "reference_id",
    "reference_length",
    "reference_start",
    "template_length",
)


def add_version_arg(parser: Any) -> argparse.Action:
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
        version="%(prog)s " + __version__,
    )


def add_verbose_arg(parser: Any) -> argparse.Action:
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
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Log debug messages to stderr.",
    )


def add_loci_args(parser: Any) -> None:
    """Add ``--locus`` for selecting one or more genomic loci."""
    parser.add_argument(
        "--locus",
        nargs="+",
        default=[],
        help="Genomic locus, like chr1:2342332 or chr1:2342-23423. "
        "Any number of loci may be specified.",
    )


def iter_read_filter_specs() -> Iterator[tuple[str, type[Any], str]]:
    """Yield ``(name, value_type, help_text)`` for every CLI read filter."""
    for prop in READ_BOOLEAN_PROPERTIES:
        yield (prop, bool, "Only reads where %s is True" % prop)
        yield ("not_" + prop, bool, "Only reads where %s is False" % prop)
    for prop in READ_STRING_PROPERTIES:
        yield (prop, str, "Only reads with the specified %s" % prop)
        yield (
            "%s_contains" % prop,
            str,
            "Only reads where %s contains the given string" % prop,
        )
    for prop in READ_INT_PROPERTIES:
        yield (prop, int, "Only reads with the specified %s" % prop)
        yield ("min_%s" % prop, int, "Only reads where %s >=N" % prop)
        yield ("max_%s" % prop, int, "Only reads where %s <=N" % prop)


def add_read_args(parser: argparse.ArgumentParser, positional: bool = False) -> None:
    """Extend a parser with BAM/SAM read-loading and read-filtering options.

    When ``positional`` is True, both optional positional BAM paths and
    ``--reads`` are accepted and combined in order.
    """
    group = parser.add_argument_group("read loading")
    if positional:
        group.add_argument(
            "reads",
            nargs="*",
            default=[],
            metavar="BAM",
            help="Paths to BAM/SAM files. Any number of paths may be specified.",
        )
        group.add_argument(
            "--reads",
            nargs="+",
            action="extend",
            default=[],
            dest="optional_reads",
            metavar="BAM",
            help="Paths to BAM/SAM files (alternative to positional paths).",
        )
    else:
        group.add_argument(
            "--reads",
            nargs="+",
            action="extend",
            default=[],
            metavar="BAM",
            help="Paths to BAM/SAM files. Any number of paths may be specified.",
        )

    group.add_argument(
        "--read-source-name",
        nargs="+",
        help="Names for each read source. The number of names specified "
        "must match the number of bam files. If not specified, filenames are "
        "used for names.",
    )

    filter_group = parser.add_argument_group(
        "read filtering",
        "A number of read filters are available. See the pysam "
        "documentation (http://pysam.readthedocs.org/en/latest/api.html) "
        "for details on what these fields mean. When multiple filter "
        "options are specified, reads must match *all* filters.",
    )

    for name, kind, message in iter_read_filter_specs():
        extra: dict[str, Any] = {}
        if kind is bool:
            extra["action"] = "store_true"
            extra["default"] = None
        elif kind is int:
            extra["type"] = int
            extra["metavar"] = "N"
        elif kind is str:
            extra["metavar"] = "STRING"
        filter_group.add_argument(
            "--" + name.replace("_", "-"),
            help=message,
            **extra,
        )


def add_variant_args(parser: argparse.ArgumentParser, positional: bool = False) -> None:
    """Extend a parser with VCF/CSV variant-loading and variant-filtering options."""
    group = parser.add_argument_group("variant loading")
    group.add_argument(
        "variants" if positional else "--variants",
        nargs=("*" if positional else "+"),
        default=[],
        help="Path to VCF file. Any number of VCF files may be specified. "
        "CSV files in the format outputted by the varlens-variants tool are "
        "also supported.",
    )
    group.add_argument(
        "--genome",
        help="Genome for the variants (e.g. b37). Required when the genome "
        "cannot be guessed from the metadata in the VCF.",
    )
    group.add_argument(
        "--include-failing-variants",
        action="store_true",
        default=False,
        help="Include variants with a non-PASS filter field.",
    )
    group.add_argument(
        "--variant-source-name",
        nargs="+",
        help="Names for variant sources. Must specify one name per variant "
        "source. If not specified, the filenames are used.",
    )
    group.add_argument(
        "--max-variants-per-source",
        type=int,
        metavar="N",
        help="Load at most N variants from each source.",
    )
    group.add_argument(
        "--single-variant",
        nargs=3,
        action="append",
        default=[],
        metavar="X",
        help="Literal variant specified as three arguments: LOCUS REF ALT. "
        "Can be specified any number of times by repeating the "
        "--single-variant option.",
    )

    filter_group = parser.add_argument_group(
        "variant filtering",
        "If multiple filters are specified, the variants must pass *all* "
        "filters. For filtering by loci, any variants that overlap "
        "the specified loci are included. Structural variants match known "
        "breakpoint bases (including remote mates), not the full SV span.",
    )
    filter_group.add_argument(
        "--ref",
        nargs="+",
        help="Include only variants where ref is one of the given values.",
    )
    filter_group.add_argument(
        "--alt",
        nargs="+",
        help="Include only variants where alt is one of the given values.",
    )
    add_loci_args(PrefixedArgumentParser(filter_group, "variant"))


def parse_int_list(string: str) -> list[int]:
    """Parse a comma/space-separated string of integers and inclusive ranges.

    Example: ``parse_int_list("8,9,10-11") == [8, 9, 10, 11]``.
    """
    integers: list[int] = []
    for comma_part in string.split(","):
        for substring in comma_part.split(" "):
            if not substring:
                continue
            if "-" in substring:
                left, right = substring.split("-", 1)
                integers.extend(range(int(left.strip()), int(right.strip()) + 1))
            else:
                integers.append(int(substring.strip()))
    return integers


def parse_mhc_predictor_names(value: str) -> list[str]:
    """Normalize and validate a ``--mhc-predictor`` CLI token against MHCtools."""
    names = [part.strip().lower() for part in value.split(",") if part.strip()]
    from mhctools.cli import args as mhctools_cli_args

    available_fn = getattr(mhctools_cli_args, "_predictor_names", None)
    if available_fn is not None:
        available = set(available_fn())
    else:
        available = set(mhctools_cli_args.mhc_predictors.keys())
    for name in names:
        if name not in available:
            raise argparse.ArgumentTypeError(
                "Unknown predictor %r. Available: %s"
                % (name, ", ".join(sorted(available)))
            )
    return names


class StoreAlleleArguments(argparse.Action):
    """Join multiple ``--mhc-alleles`` tokens into a space-separated string."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Any,
        option_string: str | None = None,
    ) -> None:
        setattr(namespace, self.dest, " ".join(values))


def add_mhc_args(
    parser: argparse.ArgumentParser,
    *,
    include_peptide_lengths_alias: bool = False,
) -> argparse._ArgumentGroup:
    """Add shared MHCtools CLI options and deprecated Varlens aliases.

    Registers the same ``--mhc-predictor``, ``--mhc-alleles``,
    ``--mhc-alleles-file``, and ``--mhc-peptide-lengths`` options as
    ``mhctools.cli.add_mhc_args(parser, predictor_required=False)`` without
    importing ``mhctools`` during parser construction or tab completion.
    """
    group = parser.add_argument_group(
        title="Prediction Options",
        description="Model selection and prediction options",
    )
    group.add_argument(
        "--mhc-predictor",
        nargs="+",
        type=parse_mhc_predictor_names,
        required=False,
        help=(
            "One or more predictor names, space or comma separated. "
            "E.g. '--mhc-predictor mhcflurry' or '--mhc-predictor netmhcpan' "
            "or '--mhc-predictor random'."
        ),
    )
    group.add_argument(
        "--mhc-predictor-path",
        help="Path to executable program or model directory",
        default=None,
        required=False,
    )
    group.add_argument(
        "--mhc-predictor-models-path",
        help="Path to models to use with predictor (e.g. mhcflurry).",
    )
    group.add_argument(
        "--mixmhcpred-python",
        help="Isolated pandas<3 Python executable for MixMHCpred or PRIME; "
        "overrides MIXMHCPRED_PYTHON without changing host dependencies.",
    )
    group.add_argument(
        "--mhc-peptide-lengths",
        type=parse_int_list,
        help="Peptide lengths to predict, e.g. '9' or '8-11' or '8,9,10'. "
        "Defaults to each predictor's own default lengths, which are "
        "narrower than the lengths it supports (the NetMHC family "
        "defaults to 9 alone).",
    )
    group.add_argument(
        "--mhc-epitope-lengths",
        type=parse_int_list,
        help="Deprecated name for --mhc-peptide-lengths",
    )
    group.add_argument(
        "--mhc-alleles-file",
        help="File with one HLA allele per line",
    )
    group.add_argument(
        "--mhc-alleles",
        nargs="+",
        action=StoreAlleleArguments,
        metavar="ALLELE",
        default="",
        help="One or more alleles, as separate quoted values or a quoted "
        "comma/space-separated list. Quote '*' to prevent shell expansion. "
        "E.g. --mhc-alleles 'HLA-A*02:01' 'HLA-B*07:02'. Required (here or via "
        "--mhc-alleles-file) for every allele-specific predictor; "
        "there is no default allele.",
    )
    group.add_argument(
        "--hla",
        nargs="+",
        default=[],
        metavar="ALLELE",
        help=(
            "Deprecated alias for --mhc-alleles. Accepts one or more MHC "
            "alleles as separate arguments or space/comma-separated lists."
        ),
    )
    if include_peptide_lengths_alias:
        group.add_argument(
            "--peptide-lengths",
            nargs="+",
            type=int,
            default=None,
            metavar="N",
            help="Deprecated alias for --mhc-peptide-lengths.",
        )
    return group


def add_effect_include_args(parser: argparse.ArgumentParser) -> None:
    """Add CLI options for Varcode variant effect annotations."""
    group = parser.add_argument_group("variant effect annotations")
    group.add_argument(
        "--include-effect",
        action="store_true",
        default=False,
        help="Include the primary Varcode annotation and its candidate consequences",
    )


def add_gene_include_args(parser: argparse.ArgumentParser) -> None:
    """Add CLI options for overlapping gene annotations."""
    group = parser.add_argument_group("gene annotations")
    group.add_argument(
        "--include-gene",
        action="store_true",
        default=False,
        help="Include gene names",
    )


def add_context_include_args(parser: argparse.ArgumentParser) -> None:
    """Add CLI options for surrounding reference sequence context."""
    group = parser.add_argument_group("variant sequence context")
    group.add_argument(
        "--include-context",
        action="store_true",
        default=False,
        help="Include variant sequence context",
    )
    group.add_argument(
        "--reference",
        help="Path to reference fasta (required for sequence context)",
    )
    group.add_argument(
        "--context-num-bases",
        type=int,
        default=15,
        metavar="N",
        help="Num bases of context to include on each side of the variant",
    )


def add_mhc_binding_include_args(parser: argparse.ArgumentParser) -> None:
    """Add CLI options for MHC binding affinity annotations."""
    group = parser.add_argument_group("MHC binding affinity")
    group.add_argument(
        "--include-mhc-binding",
        action="store_true",
        default=False,
        help="Include MHC binding (tightest affinity and allele)",
    )
    group.add_argument(
        "--hla-file",
        help=(
            "Load per-donor HLA types from the specified CSV file. It must "
            "have columns: 'donor' and 'hla'"
        ),
    )
    add_mhc_args(parser)


def add_read_evidence_include_args(parser: argparse.ArgumentParser) -> None:
    """Add CLI options for BAM read-evidence counts on variants."""
    group = parser.add_argument_group("read evidence")
    group.add_argument(
        "--include-read-evidence",
        action="store_true",
        default=False,
        help="Include counts of supporting / contradicting reads",
    )
    group.add_argument(
        "--read-sources-file",
        help="Load paths to BAMs from the given csv file.",
    )
    group.add_argument(
        "--read-sources-id-column",
        default="source_id",
        help="Column to use to join read sources with the variants dataframe.",
    )
    group.add_argument(
        "--read-sources-column",
        action="append",
        default=[],
        help="Column containing path to reads (e.g. path to a BAM). Can "
        "be specified any number of times. If not specified, all "
        "columns are used.",
    )
    group.add_argument(
        "--always-prefix-column",
        action="store_true",
        default=False,
        help="Always prefix the column names with the source name and "
        "count group, even when there is only one of each.",
    )
    group.add_argument(
        "--survive-errors",
        action="store_true",
        default=False,
        help="If an error is encountered log it and try to continue.",
    )
    add_read_args(parser)


def add_variant_include_args(parser: argparse.ArgumentParser) -> None:
    """Add all ``--include-*`` annotation option groups for ``varlens-variants``."""
    add_effect_include_args(parser)
    add_gene_include_args(parser)
    add_context_include_args(parser)
    add_mhc_binding_include_args(parser)
    add_read_evidence_include_args(parser)
