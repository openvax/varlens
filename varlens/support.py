import collections
import logging

import pandas

from .evaluation import parse_labeled_expression
from . import read_evidence, reads_util

EXPECTED_COLUMNS = [
    "source",
    "contig",
    "interbase_start",
    "interbase_end",
    "allele"
]

def allele_support_df(loci, sources, expressions={}):
    return pandas.DataFrame(
        allele_support_rows(loci, sources, expressions))

def allele_support_rows(loci, sources, expressions=[]):
    extra_columns = collections.OrderedDict()
    for labeled_expression in expressions:
        (label, expression) = parse_labeled_expression(labeled_expression)
        extra_columns[label] = expression
        
    for source in sources:
        logging.info("Reading from: %s" % source.name)
        for locus in loci:
            grouped = dict(source.pileups([locus]).group_by_allele(locus))
            for (allele, group) in grouped.items():
                d = collections.OrderedDict([
                    ("source", source.name),
                    ("contig", locus.contig),
                    ("interbase_start", str(locus.start)),
                    ("interbase_end", str(locus.end)),
                    ("allele", allele),
                    ("count", str(group.num_reads())),
                ])
                for (key, expression) in extra_columns.items():
                    num_reads = len(set(
                        read_evidence.read_key(element.alignment)
                        for pileup in group.pileups.values()
                        for element in pileup
                        if reads_util.evaluate_pileup_element_expression(
                            expression,
                            group,
                            pileup,
                            element)))
                    d[key] = num_reads
                yield pandas.Series(d)

def variant_support(variants, allele_support_df):
    '''
    Collect the read evidence support for the given variants.

    Parameters
    ----------

    variants : iterable of varcode.Variant

    allele_support_df : dataframe
        Allele support dataframe, as output by the varlens-allele-support tool.
        It should have columns: source, contig, interbase_start, interbase_end,
        allele. The remaining columns are interpreted as read counts of various
        subsets of reads (e.g. all reads, non-duplicate reads, etc.)

    Returns
    ----------

    A pandas.Panel4D frame with these axes:

    labels (axis=0) : the type of read being counted, i.e. the read count
        fields in allele_support_df.

    items (axis=1)  : the type of measurement (num_alt, num_ref, num_other,
        total_depth, alt_fraction, any_alt_fraction)

    major axis (axis=2) : the variants

    minor axis (axis=3) : the sources
    '''
    columns = list(allele_support_df.columns)
    if columns[:len(EXPECTED_COLUMNS)] != EXPECTED_COLUMNS:
        raise ValueError("Expected columns: %s" % EXPECTED_COLUMNS)

    sources = sorted(allele_support_df["source"].unique())
    
    count_fields = allele_support_df.columns[len(EXPECTED_COLUMNS):]

    panels = {}
    for field in count_fields:
        allele_support_dict = collections.defaultdict(dict)
        for (i, row) in allele_support_df.iterrows():
            key = (
                row['source'],
                row.contig,
                row.interbase_start,
                row.interbase_end)
            allele_support_dict[key][row.allele] = row[field]

        dataframe_dicts = collections.defaultdict(
            lambda: collections.defaultdict(list))
        
        for variant in variants:
            for source in sources:
                key = (source, variant.contig, variant.start - 1, variant.end)
                alleles = allele_support_dict[key]

                alt = alleles.get(variant.alt, 0)
                ref = alleles.get(variant.ref, 0)
                total = sum(alleles.values())
                other = total - alt - ref

                dataframe_dicts["num_alt"][source].append(alt)
                dataframe_dicts["num_ref"][source].append(ref)
                dataframe_dicts["num_other"][source].append(other)
                dataframe_dicts["total_depth"][source].append(total)
                dataframe_dicts["alt_fraction"][source].append(
                    float(alt) / max(1, total))
                dataframe_dicts["any_alt_fraction"][source].append(
                    float(alt + other) / max(1, total))

        dataframes = dict(
            (label, pandas.DataFrame(value, index=variants))
            for (label, value) in dataframe_dicts.items())

        panels[field] = pandas.Panel(dataframes)

    return pandas.Panel4D(panels)
    