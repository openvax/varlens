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

from __future__ import annotations

import argparse
from collections.abc import Callable, Iterable
import copy
import re
from typing import Any
import warnings

import pandas
import varcode

DEFAULT_EPITOPE_LENGTHS: tuple[int, ...] = (8, 9, 10, 11)


def parse_hla_alleles(alleles: str | Iterable[str]) -> list[str]:
    """Parse one or more MHC allele strings into a deduplicated list.

    Accepts a single space- or comma-separated string (with optional single
    quotes around legacy ``--hla`` strings) or any iterable of allele strings.
    """
    if isinstance(alleles, str):
        raw_tokens: Iterable[str] = [alleles]
    else:
        raw_tokens = alleles
    parsed: list[str] = []
    seen: set[str] = set()
    for item in raw_tokens:
        if not isinstance(item, str):
            raise TypeError(
                f"MHC alleles must be strings, got {type(item).__name__}: {item!r}"
            )
        cleaned = item.replace("'", "")
        for token in re.split(r"[\s,]+", cleaned):
            if token and token not in seen:
                seen.add(token)
                parsed.append(token)
    return parsed


def validate_epitope_lengths(epitope_lengths: Iterable[int]) -> tuple[int, ...]:
    """Validate and return a tuple of positive integer peptide lengths."""
    lengths = tuple(epitope_lengths)
    if not lengths or any(
        not isinstance(n, int) or isinstance(n, bool) or n <= 0 for n in lengths
    ):
        raise ValueError("Epitope lengths must be positive integers")
    return lengths


def resolve_mhc_predictor_factory(
    mhc_predictor: str | Callable[..., Any] | None = None,
) -> Callable[..., Any]:
    """Resolve an MHC predictor name or callable to an mhctools predictor factory.

    Parameters
    ----------
    mhc_predictor : str, callable, or None, optional
        Predictor name registered in ``mhctools.cli.mhc_predictors`` (such as
        ``"netmhcpan"``, ``"mhcflurry"``, or ``"random"``), a predictor
        class/factory callable, or ``None`` (defaults to ``mhctools.NetMHCpan``
        for backwards compatibility with Python callers).
    """
    import mhctools

    if mhc_predictor is None:
        return mhctools.NetMHCpan
    if isinstance(mhc_predictor, str):
        from mhctools.cli import mhc_predictors

        key = mhc_predictor.strip().lower()
        if key not in mhc_predictors:
            available = ", ".join(sorted(mhc_predictors))
            raise ValueError(
                f"Unknown MHC predictor {mhc_predictor!r}. Available: {available}"
            )
        return mhc_predictors[key]
    if callable(mhc_predictor):
        return mhc_predictor
    raise TypeError(
        "mhc_predictor must be a predictor name, callable, or None, "
        f"got {type(mhc_predictor).__name__}"
    )


def build_mhc_predictor(
    alleles: str | Iterable[str],
    epitope_lengths: Iterable[int] = DEFAULT_EPITOPE_LENGTHS,
    *,
    mhc_predictor: str | Callable[..., Any] | None = None,
) -> Any:
    """Instantiate an mhctools predictor for the given alleles and lengths."""
    parsed_alleles = parse_hla_alleles(alleles)
    lengths = validate_epitope_lengths(epitope_lengths)
    if not parsed_alleles:
        raise ValueError("Specify at least one MHC allele")
    factory = resolve_mhc_predictor_factory(mhc_predictor)
    return factory(alleles=parsed_alleles, default_peptide_lengths=lengths)


def as_topiary_predictor(
    predictor: Any, *, only_novel_epitopes: bool = False
) -> Any:
    """Wrap an mhctools predictor in ``TopiaryPredictor`` if needed."""
    import topiary

    if hasattr(predictor, "predict_from_variants") or hasattr(
        predictor, "predict_from_fragments"
    ):
        return predictor
    return topiary.TopiaryPredictor(
        models=predictor, only_novel_epitopes=only_novel_epitopes
    )


def add_mhc_cli_args(
    parser: argparse.ArgumentParser,
    *,
    include_peptide_lengths_alias: bool = False,
) -> argparse._ArgumentGroup:
    """Add shared MHCtools CLI options and deprecated Varlens aliases.

    Uses ``mhctools.cli.add_mhc_args(parser, predictor_required=False)`` so
    ``varlens-variants`` and ``varlens-fusions`` accept the same
    ``--mhc-predictor``, ``--mhc-alleles``, ``--mhc-alleles-file``, and
    ``--mhc-peptide-lengths`` options as Topiary and Vaxrank. Also registers
    a unified deprecated ``--hla`` alias across both commands.
    """
    from mhctools.cli import add_mhc_args

    group = add_mhc_args(parser, predictor_required=False)
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


def mhc_predictor_names_from_args(args: argparse.Namespace) -> list[str]:
    """Return the flattened list of ``--mhc-predictor`` names from parsed args."""
    raw = getattr(args, "mhc_predictor", None) or []
    return [name for group in raw for name in group]


def has_mhc_cli_args(args: argparse.Namespace) -> bool:
    """Return True if any MHC predictor or allele option was provided."""
    return bool(
        mhc_predictor_names_from_args(args)
        or getattr(args, "mhc_alleles", "")
        or getattr(args, "mhc_alleles_file", None)
        or getattr(args, "hla", None)
    )


def peptide_lengths_from_args(args: argparse.Namespace) -> tuple[int, ...] | None:
    """Resolve and validate peptide lengths from parsed MHC CLI arguments."""
    mhc_lengths = getattr(args, "mhc_peptide_lengths", None)
    if not mhc_lengths:
        mhc_lengths = getattr(args, "mhc_epitope_lengths", None)
    legacy_lengths = getattr(args, "peptide_lengths", None)

    if mhc_lengths and legacy_lengths is not None:
        raise ValueError(
            "Cannot specify both --mhc-peptide-lengths and deprecated --peptide-lengths"
        )
    if legacy_lengths is not None:
        warnings.warn(
            "--peptide-lengths is deprecated; use --mhc-peptide-lengths instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return validate_epitope_lengths(legacy_lengths)
    if mhc_lengths:
        return validate_epitope_lengths(mhc_lengths)
    return None


def mhc_predictor_from_args(
    args: argparse.Namespace,
    *,
    alleles_override: Iterable[str] | None = None,
) -> Any:
    """Build a single mhctools predictor from parsed CLI arguments.

    Delegates predictor construction to
    ``mhctools.cli.mhc_binding_predictor_from_args`` after normalizing
    deprecated Varlens aliases (``--hla`` and ``--peptide-lengths``) and
    optional per-donor allele overrides.
    """
    from mhctools.cli import mhc_binding_predictor_from_args

    effective_args = copy.copy(args)
    legacy_hla = getattr(args, "hla", None) or []
    if legacy_hla:
        if alleles_override is None and (
            getattr(args, "mhc_alleles", None)
            or getattr(args, "mhc_alleles_file", None)
        ):
            raise ValueError(
                "Cannot combine deprecated --hla with --mhc-alleles or --mhc-alleles-file"
            )
        warnings.warn(
            "--hla is deprecated; use --mhc-alleles (and --mhc-predictor) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if alleles_override is None:
            effective_args.mhc_alleles = " ".join(parse_hla_alleles(legacy_hla))
        if not mhc_predictor_names_from_args(effective_args):
            effective_args.mhc_predictor = [["netmhcpan"]]

    if alleles_override is not None:
        effective_args.mhc_alleles = " ".join(parse_hla_alleles(alleles_override))
        effective_args.mhc_alleles_file = None

    lengths = peptide_lengths_from_args(args)
    if lengths is not None:
        effective_args.mhc_peptide_lengths = list(lengths)

    return mhc_binding_predictor_from_args(effective_args)


def fragment_predictions(
    fragments: Iterable[Any],
    alleles: str | Iterable[str] = (),
    epitope_lengths: Iterable[int] = DEFAULT_EPITOPE_LENGTHS,
    *,
    predictor: Any = None,
    mhc_predictor: str | Callable[..., Any] | None = None,
) -> pandas.DataFrame:
    """Predict target-overlapping peptides, retaining all fragment metadata.

    Unlike the small-variant best-affinity summary, this returns Topiary's
    full prediction rows, with score kinds, fragment/sample identity and
    evidence. Target overlap is not a claim of proteome-wide novelty.

    Parameters
    ----------
    fragments : iterable of topiary.ProteinFragment
        Protein fragments to predict over.
    alleles : str or iterable of str, optional
        MHC alleles to predict when ``predictor`` is not supplied.
    epitope_lengths : iterable of int, optional
        Peptide lengths to predict when constructing a predictor from alleles.
    predictor : topiary.TopiaryPredictor or mhctools predictor, optional
        Pre-configured Topiary or mhctools predictor instance.
    mhc_predictor : str or callable, optional
        Predictor name (from ``mhctools.cli.mhc_predictors``) or predictor
        class/factory used when constructing a predictor from ``alleles``.
    """
    fragments = list(fragments)
    if predictor is not None:
        if alleles:
            raise ValueError("Configure alleles on the predictor or pass alleles, not both")
        if mhc_predictor is not None:
            raise ValueError("Pass either predictor or mhc_predictor, not both")
        topiary_predictor = as_topiary_predictor(predictor, only_novel_epitopes=False)
    else:
        model = build_mhc_predictor(
            alleles,
            epitope_lengths=epitope_lengths,
            mhc_predictor=mhc_predictor,
        )
        topiary_predictor = as_topiary_predictor(model, only_novel_epitopes=False)

    predictions = topiary_predictor.predict_from_fragments(fragments)
    if predictions.empty:
        return predictions
    return predictions.loc[predictions["overlaps_target"].eq(True)].reset_index(drop=True)


def binding_affinities(
    variants: Iterable[Any],
    alleles: str | Iterable[str] = (),
    epitope_lengths: Iterable[int] = DEFAULT_EPITOPE_LENGTHS,
    *,
    predictor: Any = None,
    mhc_predictor: str | Callable[..., Any] | None = None,
) -> pandas.DataFrame:
    """Return each variant's lowest mutant-peptide affinity (nM) and allele.

    Variants without an affinity prediction have a missing affinity and allele.
    Topiary's presentation scores are deliberately excluded from this summary.

    Parameters
    ----------
    variants : iterable of varcode.Variant
        Small variants to score.
    alleles : str or iterable of str, optional
        MHC alleles to predict when ``predictor`` is not supplied.
    epitope_lengths : iterable of int, optional
        Peptide lengths to predict when constructing a predictor from alleles.
    predictor : topiary.TopiaryPredictor or mhctools predictor, optional
        Pre-configured Topiary or mhctools predictor instance.
    mhc_predictor : str or callable, optional
        Predictor name (from ``mhctools.cli.mhc_predictors``) or predictor
        class/factory used when constructing a predictor from ``alleles``.
    """
    variants = list(variants)
    from .variants_util import require_small_variants

    require_small_variants(variants, "variant MHC binding summary")

    if predictor is not None:
        if alleles:
            raise ValueError("Configure alleles on the predictor or pass alleles, not both")
        if mhc_predictor is not None:
            raise ValueError("Pass either predictor or mhc_predictor, not both")
        parsed_alleles: list[str] = []
        lengths: tuple[int, ...] = ()
    else:
        parsed_alleles = parse_hla_alleles(alleles)
        lengths = validate_epitope_lengths(epitope_lengths)
        if not parsed_alleles:
            raise ValueError("Specify at least one MHC allele")
        resolve_mhc_predictor_factory(mhc_predictor)

    result = pandas.DataFrame(
        {
            "variant": variants,
            "binding_affinity": float("nan"),
            "binding_allele": pandas.Series([None] * len(variants), dtype=object),
        }
    )
    if not variants:
        return result
    if len({variant.reference_name for variant in variants}) > 1:
        raise ValueError("Cannot combine reference genomes in MHC predictions")

    if predictor is not None:
        topiary_predictor = as_topiary_predictor(predictor, only_novel_epitopes=True)
    else:
        model = build_mhc_predictor(
            parsed_alleles,
            epitope_lengths=lengths,
            mhc_predictor=mhc_predictor,
        )
        topiary_predictor = as_topiary_predictor(model, only_novel_epitopes=True)

    predictions = topiary_predictor.predict_from_variants(
        varcode.VariantCollection(variants)
    )
    if predictions.empty:
        return result

    affinities = predictions.loc[
        predictions["kind"].eq("pMHC_affinity") & predictions["affinity"].notna()
    ]
    # Stable sorting makes ties deterministic without comparing Variant objects.
    best = affinities.sort_values(
        ["affinity", "allele"], kind="stable"
    ).drop_duplicates("variant")
    best = best.set_index("variant")
    # Current Topiary serializes the variant as its short description, while
    # Varlens' public result retains the caller's Variant objects and order.
    keys = result["variant"].map(lambda variant: variant.short_description)
    result["binding_affinity"] = keys.map(best["affinity"])
    result["binding_allele"] = keys.map(best["allele"])

    return result
