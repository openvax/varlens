# Copyright (c) 2016. Mount Sinai School of Medicine
#
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
import logging

import pandas

EXPECTED_COLUMNS = [
    "source",
    "contig",
    "interbase_start",
    "interbase_end",
    "allele",
    "count",
]

def allele_support_df(loci, sources):
    """
    Returns a DataFrame of allele counts for all given loci in the read sources
    """
    return pandas.DataFrame(
        allele_support_rows(loci, sources),
        columns=EXPECTED_COLUMNS)

def allele_support_rows(loci, sources):
    for source in sources:
        logging.info("Reading from: %s (%s)" % (source.name, source.filename))
        for locus in loci:
            grouped = dict(source.pileups([locus]).group_by_allele(locus))
            if grouped:
                items = grouped.items()
            else:
                items = [("N" * (locus.end - locus.start), None)]
            for (allele, group) in items:
                d = collections.OrderedDict([
                    ("source", source.name),
                    ("contig", locus.contig),
                    ("interbase_start", str(locus.start)),
                    ("interbase_end", str(locus.end)),
                    ("allele", allele),
                    ("count", group.num_reads() if group is not None else 0),
                ])
                yield pandas.Series(d)

def variant_support(variants, allele_support_df, ignore_missing=False):
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

    ignore_missing : boolean
        If True, then varaints with no allele counts will be interpreted as
        having 0 depth. If False, then an exception will be raised if any
        variants have no allele counts.

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
    missing = [
        c for c in EXPECTED_COLUMNS if c not in allele_support_df.columns
    ]
    if missing:
        raise ValueError("Missing columns: %s" % " ".join(missing))

    # Ensure our start and end fields are ints.
    allele_support_df[["interbase_start", "interbase_end"]] = (
        allele_support_df[["interbase_start", "interbase_end"]].astype(int))

    sources = sorted(allele_support_df["source"].unique())

    allele_support_dict = collections.defaultdict(dict)
    for (i, row) in allele_support_df.iterrows():
        key = (
            row['source'],
            row.contig,
            row.interbase_start,
            row.interbase_end)
        allele_support_dict[key][row.allele] = row["count"]

    # We want an exception on bad lookups, so convert to a regular dict.
    allele_support_dict = dict(allele_support_dict)

    dataframe_dicts = collections.defaultdict(
        lambda: collections.defaultdict(list))

    for variant in variants:
        for source in sources:
            key = (source, variant.contig, variant.start - 1, variant.end)
            try:
                alleles = allele_support_dict[key]
            except KeyError:
                message = (
                    "No allele counts in source %s for variant %s" % (
                        source, str(variant)))
                if ignore_missing:
                    logging.warning(message)
                    alleles = {}
                else:
                    raise ValueError(message)

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

    return pandas.Panel(dataframes)
