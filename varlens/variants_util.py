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

import collections
import json
import warnings

import pandas
import varcode
import logging

from . import util, loci_util
from .locus import Locus
from .errors import InputError
from .read_evidence import pileup_collection

STANDARD_DATAFRAME_COLUMNS = [
    "genome",
    "contig",
    "interbase_start",
    "interbase_end",
    "ref",
    "alt",
]

STRUCTURAL_COLUMN = "structural_variant"
STRUCTURAL_SCHEMA = "varlens.structural_variant.v1"


def require_small_variants(variants, operation):
    for variant in variants:
        if isinstance(variant, varcode.StructuralVariant):
            raise InputError(
                "%s does not support structural variants (%s). Export/filter "
                "the calls with varlens-variants, or use varlens-fusions for "
                "supplied fusion RNA; interval allele counts are not junction "
                "support." % (operation, variant.short_description))


def variant_loci(variant):
    """Small-variant interval, or each known SV breakpoint as a single base."""
    if isinstance(variant, varcode.StructuralVariant):
        return [Locus.from_inclusive_coordinates(contig, position, position)
                for contig, position in variant.breakpoints]
    return [pileup_collection.to_locus(variant)]

def add_args(parser, positional=False):
    group = parser.add_argument_group("variant loading")
    group.add_argument("variants" if positional else "--variants",
        nargs=("*" if positional else "+"), default=[],
        help="Path to VCF file. Any number of VCF files may be specified. "
        "CSV files in the format outputted by the varlens-variants tool are "
        "also supported.")
    group.add_argument("--genome",
        help="Genome for the variants (e.g. b37). Required when the genome "
        "cannot be guessed from the metadata in the VCF.")
    group.add_argument("--include-failing-variants",
        action="store_true",
        default=False,
        help="Include variants with a non-PASS filter field.")
    group.add_argument("--variant-source-name", nargs="+",
        help="Names for variant sources. Must specify one name per variant "
        "source. If not specified, the filenames are used.")
    group.add_argument("--max-variants-per-source", type=int,
        metavar="N",
        help="Load at most N variants from each source.")
    group.add_argument("--single-variant", nargs=3, action="append",
        default=[], metavar="X",
        help="Literal variant specified as three arguments: LOCUS REF ALT. "
        "Can be specified any number of times by repeating the "
        "--single-variant option.")

    # Filters
    group = parser.add_argument_group("variant filtering",
        "If multiple filters are specified, the variants must pass *all* "
        "filters. For filtering by loci, any variants that overlap "
        "the specified loci are included. Structural variants match known "
        "breakpoint bases (including remote mates), not the full SV span.")
    group.add_argument("--ref", nargs="+",
        help="Include only variants where ref is one of the given values.")
    group.add_argument("--alt", nargs="+",
        help="Include only variants where alt is one of the given values.")
    loci_util.add_args(util.PrefixedArgumentParser(group, "variant"))

def load_from_args_as_dataframe(args):
    '''
    Given parsed variant-loading arguments, return a pandas DataFrame.

    If no variant loading arguments are specified, return None.
    '''
    if not args.variants and not args.single_variant:
        return None

    if args.variant_source_name:
        variant_source_names = util.expand(
            args.variant_source_name,
            'variant_source_name',
            'variant source',
            len(args.variants))
    else:
        variant_source_names = util.drop_prefix(args.variants)

    variant_to_sources = collections.defaultdict(list)

    dfs = []
    for i in range(len(args.variants)):
        name = variant_source_names[i]
        prefix = (
            'metadata:' if len(args.variants) == 1 else "metadata:%s:" % name)
        df = load_as_dataframe(
            args.variants[i],
            name=name,
            genome=args.genome,
            max_variants=args.max_variants_per_source,
            only_passing=not args.include_failing_variants,
            metadata_column_prefix=prefix)

        if df.shape[0] == 0:
            logging.warning("No variants loaded from: %s" % args.variants[i])
        else:
            for variant in df.variant:
                variant_to_sources[variant].append(name)
            dfs.append(df)

    if args.single_variant:
        variants = []
        extra_args = {}
        if args.genome:
            extra_args = {
                'genome': args.genome
            }
        for (locus_str, ref, alt) in args.single_variant:
            locus = Locus.parse(locus_str)
            variant = varcode.Variant(
                    locus.contig,
                    locus.inclusive_start,
                    ref,
                    alt,
                    **extra_args)
            variants.append(variant)
            variant_to_sources[variant].append("commandline")
        dfs.append(variants_to_dataframe(variants))

    if not dfs:
        dfs.append(variants_to_dataframe([]))

    # These columns describe the typed object and its source records, not a
    # source-specific annotation. Preserve one object and all source records
    # when the same event occurs in several files.
    structural = {}
    records = collections.defaultdict(list)
    for frame in dfs:
        if STRUCTURAL_COLUMN in frame:
            for variant, value in zip(frame.variant, frame[STRUCTURAL_COLUMN]):
                if pandas.notna(value):
                    structural.setdefault(variant, value)
            frame.drop(columns=STRUCTURAL_COLUMN, inplace=True)
        if "vcf_records" in frame:
            for variant, value in zip(frame.variant, frame.vcf_records):
                if pandas.notna(value):
                    records[variant].extend(json.loads(value))
            frame.drop(columns="vcf_records", inplace=True)

    df = dfs.pop(0)
    for other_df in dfs:
        df = pandas.merge(
            df,
            other_df,
            how='outer',
            on=["variant"] + STANDARD_DATAFRAME_COLUMNS)

    if structural:
        df[STRUCTURAL_COLUMN] = [structural.get(v) for v in df.variant]
    if records:
        df["vcf_records"] = [json.dumps(records[v], sort_keys=True)
                             if v in records else None for v in df.variant]

    genomes = df["genome"].unique()
    if len(genomes) > 1:
        raise InputError(
                "Mixing references is not supported. "
                "Reference genomes: %s" % (", ".join(genomes)))

    df["sources"] = [" ".join(variant_to_sources[v]) for v in df.variant]

    # Apply filters:
    if args.ref:
        df = df.loc[df.ref.isin(args.ref)]
    if args.alt:
        df = df.loc[df.alt.isin(args.alt)]
    loci = loci_util.load_from_args(
        util.remove_prefix_from_parsed_args(args, "variant"))
    if loci is not None:
        df = df.loc[[
            any(loci.intersects(locus) for locus in variant_loci(v))
            for v in df.variant
        ]]
    return df

def load_as_dataframe(
        filename,
        loader=None,
        name=None,
        genome=None,
        max_variants=None,
        only_passing=True,
        metadata_column_prefix=''):

    if name is None:
        name = filename
    filename = str(filename)
    if max_variants is not None and max_variants < 0:
        raise InputError("max_variants must be nonnegative")

    if loader is None:
        if (filename.endswith(".vcf") or filename.endswith(".vcf.gz")):
            # Load from VCF
            def loader(filename):
                # A successful export must not omit unsupported structural
                # alleles. Filtering PASS records and VCF placeholders remains
                # Varcode's normal, documented loading behavior.
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "error", message=r"(Could not parse|Skipped).*symbolic/breakend")
                    try:
                        collection = varcode.load_vcf(
                            filename, genome=genome, max_variants=max_variants,
                            only_passing=only_passing,
                            allow_extended_nucleotides=True,
                            parse_structural_variants=True)
                    except UserWarning as error:
                        raise InputError(
                            "%s No partial output was written." % error) from error
                return variants_to_dataframe(
                    collection,
                    collection.metadata,
                    metadata_column_prefix=metadata_column_prefix,
                    source=filename)

        elif (filename.endswith(".csv") or filename.endswith(".csv.gz")):
            # Load from csv
            def loader(filename):
                # Ignores only_passing
                df = pandas.read_csv(filename, nrows=max_variants,
                                     dtype={"contig": str, "ref": str, "alt": str})
                if "genome" not in df and genome is not None:
                    df["genome"] = genome
                for column in STANDARD_DATAFRAME_COLUMNS:
                    if column not in df:
                        raise InputError("Missing CSV column: %s" % column)
                for column in ['ref', 'alt']:
                    df[column] = df[column].fillna('')
                # CSV cannot contain live Variant objects, even when a caller
                # saved a Python dataframe including its internal column.
                df = df.drop(columns="variant", errors="ignore")
                if "vcf_records" in df:
                    for value in df.vcf_records.dropna():
                        records = json.loads(value)
                        if not isinstance(records, list) or not all(
                                isinstance(record, dict) for record in records):
                            raise InputError("vcf_records must contain a JSON list of records")
                return df
        else:
            raise InputError(
                "Unsupported variant input: %s. Use VCF or Varlens CSV; "
                "for supplied fusion RNA JSON use varlens-fusions." % filename)

    try:
        df = loader(filename)

        if 'genome' not in df:
            df["genome"] = genome

        if "variant" not in df:
            df["variant"] = [
                dataframe_row_to_variant(row) for (i, row) in df.iterrows()
            ]
    except (ValueError, OSError) as error:
        raise InputError("%s: %s" % (filename, error)) from error
    return df

def variants_to_dataframe(
        variants, metadata=None, metadata_column_prefix="", source=None):
    def record(variant):
        d = {
            'variant': variant,
            'genome': str(variant.reference_name),
            'contig': variant.contig,
            'interbase_start': variant.start - 1,
            'interbase_end': variant.end,
            'ref': variant.ref,
            'alt': variant.alt,
        }
        if isinstance(variant, varcode.StructuralVariant):
            values = variant.to_dict()
            # Genome instances are not JSON values; the canonical reference
            # name recreates the same assembly without serializing its cache.
            values["genome"] = variant.reference_name
            d[STRUCTURAL_COLUMN] = json.dumps(
                {"schema": STRUCTURAL_SCHEMA, "variant": values}, sort_keys=True)
            if source is not None:
                d["vcf_records"] = json.dumps([{
                    "source": str(source),
                    "metadata": metadata.get(variant, {}) if metadata else {},
                }], sort_keys=True)
        if metadata:
            for (name, value) in metadata.get(variant, {}).items():
                if name == 'info':
                    for (info_col, value) in (value or {}).items():
                        column = '%sinfo:%s' % (
                            metadata_column_prefix, info_col)
                        d[column] = value
                else:
                    d["%s%s" % (metadata_column_prefix, name.lower())] = value
        return d

    rows = [record(v) for v in variants]
    df = (pandas.DataFrame.from_records(rows) if rows else
          pandas.DataFrame(columns=STANDARD_DATAFRAME_COLUMNS + ["variant"]))
    column_indices = dict(
        (column, i) for (i, column) in enumerate(STANDARD_DATAFRAME_COLUMNS))
    columns = sorted(df.columns, key=lambda col: column_indices.get(col, 100))
    return df[columns]

def dataframe_row_to_variant(row):
    serialized = row.get(STRUCTURAL_COLUMN)
    if serialized is not None and pandas.notna(serialized):
        payload = json.loads(serialized)
        if not isinstance(payload, dict) or payload.get("schema") != STRUCTURAL_SCHEMA:
            raise InputError("Unsupported structural_variant schema")
        try:
            variant = varcode.StructuralVariant(**payload["variant"])
        except (KeyError, TypeError, ValueError) as error:
            raise InputError("Invalid structural_variant payload: %s" % error) from error
        expected = (variant.reference_name, variant.contig, variant.start - 1,
                    variant.end, variant.ref, variant.alt)
        observed = (row.genome, str(row.contig), row.interbase_start,
                    row.interbase_end, row.ref, row.alt)
        if expected != observed:
            raise InputError("CSV columns disagree with structural_variant payload")
        return variant
    if any(character in row.alt for character in "<>[]."):
        raise InputError("Structural CSV alleles require the structural_variant "
                         "column; load the original VCF to preserve event fields")
    return varcode.Variant(
            genome=row.genome,
            contig=row.contig,
            start=row.interbase_start + 1,
            ref=row.ref,
            alt=row.alt,
            allow_extended_nucleotides=True)

def dataframe_to_variants(df):
    for column in STANDARD_DATAFRAME_COLUMNS:
        if column not in df:
            raise ValueError("Missing column: %s" % column)

    extra_columns = [
        c for c in df.columns if c not in STANDARD_DATAFRAME_COLUMNS + ["variant"]
    ]
    metadata = collections.OrderedDict()
    for (i, row) in df.iterrows():
        variant = dataframe_row_to_variant(row)
        metadata[variant] = dict((c, row[c]) for c in extra_columns)

    return varcode.VariantCollection(
        metadata.keys(), source_to_metadata_dict={"dataframe": metadata})

def load_csv(filename, genome=None):
    df = load_as_dataframe(str(filename), genome=genome)
    return dataframe_to_variants(df)
