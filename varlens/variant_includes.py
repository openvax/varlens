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
import collections
from collections.abc import Callable, Iterable, Mapping, Sequence
import logging
import time
from typing import Any

import pandas
import numpy
from varcode import iter_candidate_effects

import pyfaidx

from . import sequence_context
from . import mhc_binding
from . import reads_util
from . import support
from . import read_evidence
from .commands.arguments import (
    add_context_include_args,
    add_effect_include_args,
    add_gene_include_args,
    add_mhc_binding_include_args,
    add_read_evidence_include_args,
)
from .errors import InputError
from .variants_util import require_small_variants

class Includeable(object):
    columns: list[str] | None = None
    string_columns: Sequence[str] = ()
    float_columns: Sequence[str] = ()

    @classmethod
    def from_args(cls, args):
        return cls()

    def process_chunk(self, df):
        raise NotImplementedError()

    def compute(self, df, chunk_rows=None):
        require_small_variants(df.variant, self.name)
        assert self.columns
        for column in self.columns:
            if column not in df.columns:
                dtype = "object" if column in self.string_columns else "float64"
                # Explicit object dtype supports partially populated text
                # columns loaded from CSV.
                df[column] = pandas.Series(index=df.index, dtype=dtype)
            elif column in self.string_columns:
                df[column] = df[column].astype(object)
            elif column in self.float_columns:
                # CSV inference may give a partially annotated affinity column
                # integer dtype; pandas 3 requires room for fractional results.
                df[column] = df[column].astype("float64")
        rows_to_annotate = pandas.isnull(df[self.columns[0]])
        for column in self.columns[1:]:
            rows_to_annotate = rows_to_annotate | pandas.isnull(df[column])

        while rows_to_annotate.sum() > 0:
            if chunk_rows:
                this_chunk_rows = rows_to_annotate & (
                    rows_to_annotate.cumsum() <= chunk_rows)
            else:
                this_chunk_rows = rows_to_annotate

            num_remaining = rows_to_annotate.sum()
            logging.info("%s: %d / %d (%0.1f%%) remaining. Processing %d rows."
                % (
                    self.name,
                    num_remaining,
                    len(rows_to_annotate),
                    num_remaining * 100.0 / len(rows_to_annotate),
                    this_chunk_rows.sum()))

            rows_to_annotate = rows_to_annotate & (~ this_chunk_rows)
            
            if this_chunk_rows.sum() > 0:
                start = time.time()
                df.loc[this_chunk_rows, self.columns] = self.process_chunk(
                    df.loc[this_chunk_rows].copy())[self.columns]
                logging.info("Processed in %0.2f sec" % (time.time() - start))
            yield this_chunk_rows.sum()

class Effect(Includeable):
    name = "variant effect annotations"
    columns = ["effect", "effect_candidates"]
    string_columns = columns

    @staticmethod
    def add_args(parser):
        add_effect_include_args(parser)

    @staticmethod
    def requested(args):
        return args.include_effect

    def process_chunk(self, df):
        effects = [v.effects().top_priority_effect() for v in df["variant"]]
        df["effect"] = [effect.short_description if effect is not None else "unavailable"
                        for effect in effects]
        df["effect_candidates"] = ["; ".join(
            "%s: %s" % (type(candidate).__name__, candidate.short_description)
            if candidate is not None else "Unresolved: unavailable candidate"
            for candidate in iter_candidate_effects(effect)) for effect in effects]
        return df

class Gene(Includeable):
    name = "gene annotations"
    columns = ["gene"]
    string_columns = columns
    
    @staticmethod
    def add_args(parser):
        add_gene_include_args(parser)

    @staticmethod
    def requested(args):
        return args.include_gene

    def process_chunk(self, df):
        df["gene"] = [
            ' '.join(v.gene_names) if v.gene_names else 'None'
            for v in df.variant
        ]
        return df

class Context(Includeable):
    name = "variant sequence context"
    columns = ["context_5_prime", "context_3_prime", "context_mutation"]
    string_columns = columns
    
    @staticmethod
    def add_args(parser):
        add_context_include_args(parser)

    @classmethod
    def from_args(cls, args):
        if not args.reference:
            raise ValueError(
                "The --reference argument is required when including context")
        return cls(
            reference=pyfaidx.Fasta(args.reference),
            context_num_bases=args.context_num_bases)

    def __init__(self, reference, context_num_bases):
        self.reference = reference
        self.context_num_bases = context_num_bases

    @staticmethod
    def requested(args):
        return args.include_context

    def process_chunk(self, df):
        context_5_prime = []
        context_3_prime = []
        context_mutation = []
        for variant in df.variant:
            tpl = sequence_context.variant_context(
                self.reference,
                variant.contig,
                variant.start,
                variant.end,
                variant.alt,
                self.context_num_bases)
            context_5_prime.append(tpl[0])
            context_mutation.append(tpl[1])
            context_3_prime.append(tpl[2])

        df["context_5_prime"] = context_5_prime
        df["context_3_prime"] = context_3_prime
        df["context_mutation"] = context_mutation
        return df
    
class MHCBindingAffinity(Includeable):
    name = "MHC binding affinity"
    columns = ["binding_affinity", "binding_allele"]
    string_columns = ["binding_allele"]
    float_columns = ["binding_affinity"]

    noncoding_effects = set([
        "intergenic",
        "intronic",
        "non-coding-transcript",
        "3' UTR",
        "5' UTR",
        "silent",
    ])
    
    @staticmethod
    def add_args(parser: argparse.ArgumentParser) -> None:
        add_mhc_binding_include_args(parser)

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> MHCBindingAffinity:
        has_global_alleles = bool(
            getattr(args, "mhc_alleles", "")
            or getattr(args, "mhc_alleles_file", None)
            or getattr(args, "hla", None)
        )
        if args.hla_file and has_global_alleles:
            raise InputError(
                "Specify either --hla-file or global MHC alleles "
                "(--mhc-alleles / --mhc-alleles-file / --hla), not both"
            )
        if args.hla_file:
            predictor_names = mhc_binding.mhc_predictor_names_from_args(args)
            if not predictor_names:
                raise InputError("MHC predictor required (use --mhc-predictor)")
            if len(predictor_names) != 1:
                raise InputError(
                    "Expected exactly one MHC predictor, got %d: %s"
                    % (len(predictor_names), ", ".join(predictor_names))
                )
            try:
                mhc_binding.peptide_lengths_from_args(args)
            except ValueError as error:
                raise InputError(str(error)) from error
            return cls(
                hla_dataframe=pandas.read_csv(args.hla_file),
                predictor_factory=lambda alleles: mhc_binding.mhc_predictor_from_args(
                    args, alleles_override=alleles
                ),
            )
        try:
            predictor = mhc_binding.mhc_predictor_from_args(args)
        except ValueError as error:
            raise InputError(str(error)) from error
        return cls(predictor=predictor)

    @staticmethod
    def string_to_hla_alleles(s: str) -> list[str]:
        return mhc_binding.parse_hla_alleles(s)

    def __init__(
        self,
        hla: str | Iterable[str] | None = None,
        hla_dataframe: pandas.DataFrame | None = None,
        donor_to_hla: Mapping[Any, Iterable[str] | None] | None = None,
        *,
        predictor: Any = None,
        mhc_predictor: str | Callable[..., Any] | None = None,
        epitope_lengths: Iterable[int] = mhc_binding.DEFAULT_EPITOPE_LENGTHS,
        predictor_factory: Callable[[Iterable[str]], Any] | None = None,
    ) -> None:
        """
        Specify exactly one of hla, hla_dataframe, donor_to_hla, or predictor.

        Parameters
        -----------
        hla : str or iterable of str
            HLA alleles to use for all donors

        hla_dataframe : pandas.DataFrame with columns 'donor' and 'hla'
            DataFrame giving HLA alleles for each donor. The 'hla' column
            should be a space- or comma-separated list of alleles for that donor.

        donor_to_hla : mapping of donor -> allele list
            Map from donor to HLA alleles for that donor.

        predictor : topiary.TopiaryPredictor or mhctools predictor, optional
            Pre-configured predictor instance to use for all donors.

        mhc_predictor : str or callable, optional
            Predictor name (from ``mhctools.cli.mhc_predictors``) or factory
            used when constructing predictors from alleles.

        epitope_lengths : iterable of int, optional
            Peptide lengths to predict when constructing predictors from alleles.

        predictor_factory : callable, optional
            Factory receiving a donor's allele list and returning a predictor
            instance (used for per-donor CLI ``--hla-file`` runs).
        """
        if (
            bool(hla)
            + (hla_dataframe is not None)
            + bool(donor_to_hla)
            + (predictor is not None)
            != 1
        ):
            raise TypeError(
                "Must specify exactly one of hla, hla_dataframe, donor_to_hla, or predictor"
            )
        if predictor is not None and (
            mhc_predictor is not None or predictor_factory is not None
        ):
            raise ValueError(
                "Pass either predictor or mhc_predictor/predictor_factory, not both"
            )

        self.hla: list[str] | None = (
            self.string_to_hla_alleles(hla)
            if isinstance(hla, str)
            else (list(hla) if hla is not None else None)
        )
        self.donor_to_hla: dict[Any, list[str] | None] | None = (
            {
                donor: (list(alleles) if alleles is not None else None)
                for donor, alleles in donor_to_hla.items()
            }
            if donor_to_hla is not None
            else None
        )
        if hla_dataframe is not None:
            self.donor_to_hla = {}
            for (_i, row) in hla_dataframe.iterrows():
                if row.donor in self.donor_to_hla:
                    raise ValueError("Multiple rows for donor: %s" % row.donor)
                if pandas.isnull(row.hla):
                    self.donor_to_hla[row.donor] = None
                else:
                    self.donor_to_hla[row.donor] = self.string_to_hla_alleles(
                        row.hla)
        self.predictor = predictor
        self.mhc_predictor = mhc_predictor
        self.epitope_lengths = tuple(epitope_lengths)
        self.predictor_factory = predictor_factory
        assert (
            self.hla is not None
            or self.donor_to_hla is not None
            or self.predictor is not None
        )

    @staticmethod
    def requested(args: argparse.Namespace) -> bool:
        return args.include_mhc_binding

    def process_chunk(self, df: pandas.DataFrame) -> pandas.DataFrame:
        if self.predictor is not None:
            rows = (
                ~df.effect.isin(self.noncoding_effects)
                if "effect" in df
                else pandas.Series(True, index=df.index)
            )
            sub_df = df.loc[rows]
            if sub_df.shape[0] > 0:
                result = mhc_binding.binding_affinities(
                    sub_df.variant, predictor=self.predictor
                )
                df.loc[rows, "binding_affinity"] = (
                    result["binding_affinity"].values
                )
                df.loc[rows, "binding_allele"] = (
                    result["binding_allele"].values
                )
            return df

        drop_donor = False
        if 'donor' not in df:
            df["donor"] = "DONOR1"
            drop_donor = True
        for donor in df.donor.unique():
            rows = (df.donor == donor)
            if 'effect' in df:
                rows = rows & (~df.effect.isin(self.noncoding_effects))
            sub_df = df.loc[rows]
            alleles = (
                self.hla
                if self.hla
                else (self.donor_to_hla.get(donor) if self.donor_to_hla else None)
            )
            if alleles and sub_df.shape[0] > 0:
                if self.predictor_factory is not None:
                    donor_predictor = self.predictor_factory(alleles)
                    result = mhc_binding.binding_affinities(
                        sub_df.variant, predictor=donor_predictor
                    )
                else:
                    result = mhc_binding.binding_affinities(
                        sub_df.variant,
                        alleles,
                        epitope_lengths=self.epitope_lengths,
                        mhc_predictor=self.mhc_predictor,
                    )
                df.loc[rows, "binding_affinity"] = (
                    result["binding_affinity"].values)
                df.loc[rows, "binding_allele"] = (
                    result["binding_allele"].values)
        if drop_donor:
            del df["donor"]
        return df

class ReadEvidence(Includeable):
    name = "read evidence"
    default_column_format = "{source}_count_{allele_group}"
    
    @classmethod
    def add_args(cls, parser):
        add_read_evidence_include_args(parser)

    @classmethod
    def from_args(cls, args):
        read_sources = reads_util.load_from_args(args)
        read_sources_df = None
        if args.read_sources_file is not None:
            read_sources_df = pandas.read_csv(
                args.read_sources_file,
                index_col=args.read_sources_id_column)
            if args.read_sources_column:
                read_sources_df = read_sources_df[args.read_sources_column]

        source_names = cls.read_source_names(read_sources, read_sources_df)
        if (args.always_prefix_column or len(source_names) > 1):
            column_format = cls.default_column_format
        else:
            column_format = "{allele_group}"
        return cls(
            read_sources=read_sources,
            read_sources_df=read_sources_df,
            column_format=column_format,
            survive_errors=args.survive_errors)

    def __init__(self,
            read_sources=None,
            read_sources_df=None,
            column_format=default_column_format,
            survive_errors=False):
        """
        
        """
        if sum(x is not None for x in [read_sources, read_sources_df]) != 1:
            raise TypeError(
                "Specify exactly one of read_sources, read_sources_df")
 
        self.read_sources = read_sources
        self.read_sources_df = read_sources_df
        self.column_format = column_format
        self.survive_errors = survive_errors
        self.set_columns()

    @staticmethod
    def read_source_names(read_sources=None, read_sources_df=None):
        if read_sources is not None:
            return [x.name for x in read_sources]
        return read_sources_df.columns.tolist()
   
    def set_columns(self):
        source_names = self.read_source_names(
            read_sources=self.read_sources,
            read_sources_df=self.read_sources_df)
        assert source_names
        self.columns_dict = collections.OrderedDict()
        for source_name in source_names:
            for allele_group in ["num_alt", "num_ref", "total_depth"]:
                column_name = self.column_name(
                    source_name, allele_group)
                self.columns_dict[column_name] = (
                    source_name, allele_group)
        self.columns = list(self.columns_dict)                
        assert self.columns

    def column_name(self, source, allele_group):
        """
        Parameters
        ----------
        source : string
            name of the ReadSource

        allele_group : string
            one of: num_ref, num_alt, total_depth

        Returns
        -------
        column name : string
        """
        return self.column_format.format(
            source=source,
            allele_group=allele_group)

    @staticmethod
    def requested(args):
        return args.include_read_evidence

    def process_chunk(self, df):
        if self.read_sources_df is None:
            def rows_and_read_sources():
                all_rows = numpy.ones(df.shape[0], dtype=bool)
                yield (all_rows, self.read_sources)
        else:
            def rows_and_read_sources():
                join_col = self.read_sources_df.index.name
                for join_value in df[join_col].unique():
                    rows = df[join_col] == join_value
                    read_paths = self.read_sources_df.loc[join_value]
                    read_sources = []
                    for (name, filename) in read_paths.items():
                        if pandas.isnull(filename):
                            continue
                        relevant_columns = [
                            col for (col, (source_name, allele_group))
                            in self.columns_dict.items()
                            if source_name == name
                        ]
                        if (~pandas.isnull(df[relevant_columns].values)).all():
                            logging.info(
                                "Skipping source %s (%s) for %s: data exists" %
                                (name, filename, join_value))
                            continue
                        try:
                            read_sources.append(reads_util.load_bam(
                                filename,
                                name=name))
                        except Exception as e:
                            logging.error("Error loading bam: %s in %s" %
                                (str(e), filename))
                            if not self.survive_errors:
                                raise
                            continue

                    if rows.sum() > 0 and read_sources:
                        logging.info(
                            "Processing %s=%s (%d rows, %d read sources)" % (
                                join_col,
                                join_value,
                                rows.sum(),
                                len(read_sources)))
                        yield (rows, read_sources)
                    else:
                        logging.info(
                            "Skipping %s=%s (%d rows, %d read sources)" % (
                                join_col,
                                join_value,
                                rows.sum(),
                                len(read_sources)))

        for (rows, sources) in rows_and_read_sources():
            variants = df.variant[rows]
            counter = collections.Counter(variants)
            duplicate_variants = dict(
                (v, c) for (v, c) in counter.items() if c > 1)
            if duplicate_variants:
                raise ValueError("Duplicate variant(s) for this source: %s" %
                        duplicate_variants)
            variant_loci = sorted(set(
                read_evidence.pileup_collection.to_locus(variant)
                for variant in variants))

            allele_support_df = support.allele_support_df(
                variant_loci, sources)
            assert set(s.name for s in sources) == set(
                allele_support_df.source.unique())
            variant_support_df = support.variant_support(
                variants, allele_support_df)
            assert set(s.name for s in sources) == set(
                variant_support_df.columns.get_level_values("source"))

            for allele_group in ["num_alt", "num_ref", "total_depth"]:
                sub_panel = variant_support_df[allele_group].loc[list(variants)]
                for source_column in sub_panel.columns:
                    dest_column = self.column_name(
                        source_column, allele_group)
                    assert dest_column in self.columns, (
                            "Bad column: %s not in %s" % (
                                dest_column, " ".join(self.columns)))
                    values = sub_panel[source_column].values
                    assert len(values) == rows.sum(), "%d != %d" % (
                        len(values), rows.sum())
                    df.loc[rows, dest_column] = values
        return df

INCLUDEABLES = Includeable.__subclasses__()
