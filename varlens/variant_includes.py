import logging
import time
import collections

import pandas
import numpy
import typechecks

import pyfaidx

from . import sequence_context
from . import mhc_binding
from . import reads_util
from . import support
from . import read_evidence

class Includeable(object):
    columns = None

    @classmethod
    def from_args(cls, args):
        return cls()

    def process_chunk(self, df):
        raise NotImplementedError()

    def compute(self, df, chunk_rows=None):
        assert self.columns
        for column in self.columns:
            if column not in df.columns:
                df[column] = numpy.nan
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
                df.ix[this_chunk_rows, self.columns] = self.process_chunk(
                    df.ix[this_chunk_rows].copy())[self.columns]
                logging.info("Processed in %f0.2 sec" % (time.time() - start))
            yield this_chunk_rows.sum()

class Effect(Includeable):
    name = "variant effect annotations"
    columns = ["effect"]

    @staticmethod
    def add_args(parser):
        parser = parser.add_argument_group(Effect.name)
        parser.add_argument("--include-effect",
            action="store_true", default=False,
            help="Include varcode effect annotations")

    @staticmethod
    def requested(args):
        return args.include_effect

    def process_chunk(self, df):
        df["effect"] = [
            v.effects().top_priority_effect().short_description
            for v in df["variant"]
        ]
        return df

class Gene(Includeable):
    name = "gene annotations"
    columns = ["gene"]
    
    @staticmethod
    def add_args(parser):
        parser = parser.add_argument_group(Gene.name)
        parser.add_argument("--include-gene",
            action="store_true", default=False,
            help="Include gene names")

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
    
    @staticmethod
    def add_args(parser):
        parser = parser.add_argument_group(Context.name)
        parser.add_argument("--include-context",
            action="store_true", default=False,
            help="Include variant sequence context")
        parser.add_argument("--reference",
            help="Path to reference fasta (required for sequence context)")
        parser.add_argument("--context-num-bases", type=int, default=15,
            metavar="N",
            help="Num bases of context to include on each side of the variant")

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

    noncoding_effects = set([
        "intergenic",
        "intronic",
        "non-coding-transcript",
        "3' UTR",
        "5' UTR",
        "silent",
    ])
    
    @staticmethod
    def add_args(parser):
        parser = parser.add_argument_group(MHCBindingAffinity.name)
        parser.add_argument("--include-mhc-binding",
            action="store_true", default=False,
            help="Include MHC binding (tightest affinity and allele)")
        parser.add_argument("--hla",
            help="Space separated list of MHC alleles, e.g. 'A:02:01 A:02:02'")
        parser.add_argument('--hla-file',
            help="Load HLA types from the specified CSV file. It must have "
            "columns: 'donor' and 'hla'")

    @classmethod
    def from_args(cls, args):
        if bool(args.hla) + bool(args.hla_file) != 1:
            raise ValueError("Must specify exactly one of --hla or --hla-file")
        return cls(
            hla=args.hla,
            hla_dataframe=(
                pandas.read_csv(args.hla_file) if args.hla_file else None))

    @staticmethod
    def string_to_hla_alleles(s):
        return s.replace("'", "").split()

    def __init__(self, hla=None, hla_dataframe=None, donor_to_hla=None):
        """
        Specify exactly one of hla, hla_dataframe, or donor_to_hla.

        Parameters
        -----------
        hla : list of string
            HLA alleles to use for all donors

        hla_dataframe : pandas.DataFrame with columns 'donor' and 'hla'
            DataFrame giving HLA alleles for each donor. The 'hla' column
            should be a space separated list of alleles for that donor.

        donor_to_hla : dict of string -> string list
            Map from donor to HLA alleles for that donor.
        """
        if bool(hla) + (hla_dataframe is not None) + bool(donor_to_hla) != 1:
            raise TypeError(
                "Must specify exactly one of hla, hla_dataframe, donor_to_hla")
        
        self.hla = (
            self.string_to_hla_alleles(hla) if typechecks.is_string(hla)
            else hla)
        self.donor_to_hla = donor_to_hla
        if hla_dataframe is not None:
            self.donor_to_hla = {}
            for (i, row) in hla_dataframe.iterrows():
                if row.donor in self.donor_to_hla:
                    raise ValueError("Multiple rows for donor: %s" % row.donor)
                if pandas.isnull(row.hla):
                    self.donor_to_hla[row.donor] = None
                else:
                    self.donor_to_hla[row.donor] = self.string_to_hla_alleles(
                        row.hla)
        assert self.hla is not None or self.donor_to_hla is not None

    @staticmethod
    def requested(args):
        return args.include_mhc_binding

    def process_chunk(self, df):
        drop_donor = False
        if 'donor' not in df:
            df["donor"] = "DONOR1"
            drop_donor = True
        for donor in df.donor.unique():
            rows = (df.donor == donor)
            if 'effect' in df:
                rows = rows & (~df.effect.isin(self.noncoding_effects))
            sub_df = df.loc[rows]
            alleles = self.hla if self.hla else self.donor_to_hla.get(donor)
            if alleles and sub_df.shape[0] > 0:
                result = mhc_binding.binding_affinities(
                    sub_df.variant, alleles)
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
        group = parser.add_argument_group(cls.name)
        group.add_argument("--include-read-evidence",
            action="store_true", default=False,
            help="Include counts of supporting / contradicting reads")
        group.add_argument("--read-sources-file",
            help="Load paths to BAMs from the given csv file.")
        group.add_argument("--read-sources-id-column",
            default="source_id",
            help="Column to use to join read sources with the variants "
            "dataframe.")
        group.add_argument("--read-sources-column", action="append",
            default=[],
            help="Column containing path to reads (e.g. path to a BAM). Can "
            "be specified any number of times. If not specified, all "
            "columns are used.")
        group.add_argument("--always-prefix-column", action="store_true",
            default=False,
            help="Always prefix the column names with the source name and "
            "count group, even when there is only one of each.")
        group.add_argument("--survive-errors", action="store_true",
            default=False,
            help="If an error is encountered log it and try to continue.")
        
        reads_util.add_args(parser)

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
            read_filters=[],
            column_format=default_column_format,
            survive_errors=False):
        """
        
        """
        if sum(x is not None for x in [read_sources, read_sources_df]) != 1:
            raise TypeError(
                "Specify exactly one of read_sources, read_sources_df")
 
        self.read_sources = read_sources
        self.read_sources_df = read_sources_df
        self.read_filters = read_filters
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
        ----------
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
                    read_paths = self.read_sources_df.ix[join_value]
                    read_sources = []
                    for (name, filename) in read_paths.iteritems():
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
                                name=name,
                                read_filters=self.read_filters))
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
                variant_support_df.minor_axis)

            for allele_group in ["num_alt", "num_ref", "total_depth"]:
                sub_panel = variant_support_df[allele_group, variants]
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
