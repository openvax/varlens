from __future__ import with_statement

import logging
import numpy
import collections
import math
try:
    import cPickle as pickle
except ImportError:
    import pickle

from matplotlib import pyplot

import varcode
from varcode import read_evidence

from .. import load_variants

PIE_CHART_COLORS = [
    "#37AE37",
    "#EEDB00",
    "#DA4545",
    "#0D5555",
    "#5EC1C1",
    'brown',
    'black',
    'beige',
]


def plot(
        variants=[],
        loci=[],
        read_inputs=None,
        evidence=None,
        evidence_out=None,
        neighboring_loci_offsets=[],
        allele_min_percent=1.0,
        max_alleles=5,
        radius_function=lambda values: math.pow(sum(values) / 1300., 1 / 3.0),
        cell_label_function=lambda read_label, locus: "",
        draw_kwargs={}):
    '''
    High level plotting function.

    Parameters
    -----------

    variants : iterable of Variant

    evidence : dict read label -> locus -> allele -> count [optional]
        If not specified, collect_evidence will be used to get this.

    evidence_out : dict [optional]
        If specified, the given dict is mutated to contain the collected
        evidence. Next time you call plot(), you can pass this in as evidence
        for a performance improvement.

    cell_labels : [optional] dict of (read label, locus label) -> string
        Extra info to include below a pie chart.

    '''

    # rows is a list of (locus, description, first_alleles)
    rows = [(locus, locus_label(locus), ()) for locus in loci]

    for (variant_num, variant) in enumerate(variants):
        locus = varcode.read_evidence.pileup_collection.to_locus(variant)
        rows.append(
            (locus,
                variant_label(variant_num + 1, variant),
                (variant.ref, variant.alt)))
        for offset in neighboring_loci_offsets:
            nearby_locus = varcode.Locus(
                contig=locus.contig,
                start=locus.start + offset,
                end=locus.end + offset)
            rows.append((locus, locus_label(nearby_locus), ()))

    if not rows:
        if not evidence:
            raise ValueError(
                "Must specify at least one of loci, variants, and evidence")

        # Use the loci directly from evidence.
        all_loci = set()
        for loci in evidence.values():
            all_loci.update(loci)
        all_loci = sorted(all_loci)
        rows = [(x, locus_label(x), ()) for x in all_loci]

    if evidence is None:
        if not read_inputs:
            raise ValueError("Must specify one of evidence, read_inputs")
        evidence = collect_evidence([x for (x, _, _) in rows], read_inputs)

    if evidence_out is not None:
        evidence_out.clear()
        evidence_out.update(evidence)

    # For each row, a list of (allele, [counts for each col]) pairs.
    allele_vectors_matrix = []
    radii = []
    cell_labels = []
    for (locus, description, first_alleles) in rows:
        # Append to allele_vectors_matrix
        complete_allele_dicts = [
            value.get(locus, {})
            for value in evidence.values()
        ]
        allele_vectors = allele_dicts_to_allele_vector_pairs(
                complete_allele_dicts,
                first_alleles=first_alleles,
                min_percent=allele_min_percent,
                max_alleles=max_alleles)
        allele_vectors_matrix.append(allele_vectors)

        # Append to radii
        radii.append([])
        for (col_num, value) in enumerate(evidence.values()):
            values = [vector[col_num] for (_, vector) in allele_vectors]
            radius = radius_function(values)
            radii[-1].append(radius)

        # Append to cell labels.
        cell_labels.append([
            cell_label_function(reads_label, locus)
            for reads_label in evidence
        ])

    print(radii)
    return draw_plot(
        allele_vectors_matrix,
        radii,
        [description for (_, description, _) in rows],
        evidence.keys(),
        cell_labels,
        **draw_kwargs)

def abbreviate_allele(allele):
    if allele == "":
        return "del"
    if len(allele) > 5:
        return "%s+%d" % (allele[0], len(allele))
    return allele

def allele_dicts_to_allele_vector_pairs(
        full_dicts,
        first_alleles=(),
        min_percent=None,
        max_alleles=None,
        abbreviate_allele_function=abbreviate_allele):

    total_counts = collections.Counter()
    alleles_with_sufficient_evidence = set()
    for full_dict in full_dicts:
        counts = collections.Counter(full_dict)
        total_counts += counts
        total = sum(counts.values())
        if min_percent is not None:
            alleles_with_sufficient_evidence.update(
                allele for (allele, count)
                in counts.most_common()
                if count * 100.0 / (total if total > 0 else 1) >= min_percent)

    if min_percent is None:
        # All alleles have sufficient evidence.
        alleles_with_sufficient_evidence = set(total_counts)

    alleles = list(first_alleles) + [
                allele for (allele, _) in total_counts.most_common()
                if (allele not in first_alleles and
                    allele in alleles_with_sufficient_evidence)]

    if max_alleles is not None:
        alleles = alleles[:max_alleles]

    skipped_alleles = set(total_counts) - set(alleles)

    allele_vector_pairs = [
        (abbreviate_allele_function(allele),
            [full_dict.get(allele, 0) for full_dict in full_dicts])
        for allele in alleles
    ]
    if skipped_alleles:
        # Add "other category".
        vector = [
            sum(full_dict.get(allele, 0) for allele in skipped_alleles)
            for full_dict in full_dicts
        ]
        allele_vector_pairs.append(
            ("%d others" % len(skipped_alleles), vector))

    return allele_vector_pairs

def locus_label(locus):
    return "{locus.contig}:{locus.position}".format(locus=locus)

def variant_label(num, variant):
    effects = variant.effects()
    return "\n".join([
        "#{num}",
        "{variant.contig}:{variant.start}{optional_end}",
        "{gene_names}",
        "{variant.ref}->{variant.alt}",
        "{effect}",
    ]).format(
        num=num,
        variant=variant,
        optional_end=(
            "" if variant.start == variant.end else "-%s" % variant.end),
        gene_names=" ".join(variant.gene_names),
        effect=effects.top_priority_effect().short_description)




def collect_evidence(loci, read_inputs):
    '''
    Parameters
    ----------
    loci : list of varcode.Locus or varcode.Variant instances

    read_inputs : list of ReadInput

    Returns
    ----------
    dict : read label -> locus -> allele -> count
    '''
    print("Loci", loci)
    result = collections.OrderedDict()
    for (i, read_input) in enumerate(read_inputs):
        logging.info("Reading pileups for %d / %d: %s" % (
            (i + 1), len(read_inputs), read_input.path))
        evidence = read_evidence.PileupCollection.from_bam(
            read_input.path, loci)
        logging.info("Done loading evidence")
        result[read_input.name] = collections.OrderedDict()
        for (j, locus) in enumerate(loci):
            if j % 100 == 0:
                logging.info("Locus %d / %d" % ((j + 1), len(loci)))
            result[read_input.name][locus] = dict(
                evidence.allele_summary(locus))
    return result

def text_args(user_args, **default_args):
    if not isinstance(user_args, dict):
        user_args = {"s": user_args}
    default_args.update(user_args)
    return default_args

def draw_plot(
        allele_vectors_matrix,
        radii,
        row_labels,
        col_labels,
        cell_labels,
        rows_per_figure=10,
        col_label_height=8.0,
        row_label_width=3.5,
        cell_width=2.5,
        cell_height=3.5,
        pie_chart_colors=PIE_CHART_COLORS,
        legend_extra_properties={},
        subplots_adjust_properties={'hspace': 0.1}):
    """

    Low level plotting function.

    Parameters
    ---


    radii
        Num rows x num cols array

    """

    num_rows = len(radii)
    num_cols = len(radii[0])

    assert len(allele_vectors_matrix) == num_rows
    assert len(row_labels) == num_rows
    assert len(col_labels) == num_cols

    legend_properties = {
        'bbox_to_anchor': (-.5, .8),
        'fontsize': 'large',
        'loc': 'upper left',
        'frameon': False,
    }
    legend_properties.update(legend_extra_properties)

    figure_row_plans = []
    current_figure_row_end = 0  # exclusive
    while current_figure_row_end < num_rows:
        current_figure_rows = min(
            rows_per_figure, (num_rows - current_figure_row_end))
        current_figure_row_start = current_figure_row_end
        current_figure_row_end += current_figure_rows
        figure_row_plans.append((
            current_figure_rows,
            current_figure_row_start,
            current_figure_row_end))

    for (figure_num, plan) in enumerate(figure_row_plans):
        logging.info("Creating figure %d / %d" % (
            figure_num + 1, len(figure_row_plans)))
        (current_figure_rows,
            current_figure_row_start,
            current_figure_row_end) = plan
   
        figure = pyplot.figure(
            figsize=(
                row_label_width + cell_width * num_cols,
                col_label_height + cell_height * current_figure_rows))

        label_height_ratio = col_label_height / float(cell_height)
        label_width_ratio = row_label_width / float(cell_width)
        gridspec = pyplot.GridSpec(
            current_figure_rows + 1,
            num_cols + 1,
            height_ratios=([label_height_ratio] +
                [1] * current_figure_rows),
            width_ratios=[label_width_ratio] + [1] * num_cols)

        # The first row of each figure gives the column labels.
        for col_num in range(num_cols):
            ax = pyplot.subplot(gridspec.new_subplotspec((0, col_num + 1)))
            ax.axis("off")
            ax.text(**text_args(
                col_labels[col_num],
                x=0.5,
                y=1.0,
                rotation=90,
                fontsize='x-large'))
            
        # Subsequent rows are the variants.
        for row_num in range(current_figure_row_start, current_figure_row_end):
            local_row_num = row_num - current_figure_row_start

            # First column gives row label.
            ax_row_label = pyplot.subplot(
                gridspec.new_subplotspec((local_row_num + 1, 0)))
            ax_row_label.text(**text_args(
                row_labels[row_num],
                x=0.5,
                y=0.5,
                transform=ax_row_label.transAxes,
                verticalalignment="center"))
            ax_row_label.axis("off")

            # List of (display name, [col1 count, col2 count, ...]) pairs.
            allele_vectors = allele_vectors_matrix[row_num]
            display_alleles = [allele for (allele, _) in allele_vectors]
            col_to_counts_vector = numpy.array(
                [vector for (_, vector) in allele_vectors]).T
            axis_with_a_pie_chart = None
            for col_num in range(num_cols):
                ax = pyplot.subplot(
                    gridspec.new_subplotspec((local_row_num + 1, col_num + 1)))
                ax.axis("off")

                vector = col_to_counts_vector[col_num]

                radius = radii[row_num][col_num]
                if radius > 0:
                    ax.pie(vector, radius=radius, colors=pie_chart_colors)
                    axis_with_a_pie_chart = ax

                ax.text(**text_args(
                    cell_labels[row_num][col_num],
                    x=0.5,
                    y=-0.2,
                    transform=ax.transAxes,
                    fontsize="large",
                    horizontalalignment="center"))

            if axis_with_a_pie_chart:
                axis_with_a_pie_chart.legend(
                    display_alleles,
                    bbox_transform=ax_row_label.transAxes,
                    **legend_properties)

        pyplot.subplots_adjust(**subplots_adjust_properties)
        yield figure

def load_variants_dict(variant_inputs, filter=None, ensembl_version=None):
    '''

    Returns
    --------
    dict : Variant -> set(VariantInput) giving the inputs that have the variant
    '''
    result = collections.defaultdict(set)
    for (i, variant_input) in enumerate(variant_inputs):
        logging.info("Loading vcf %d / %d" % ((i + 1), len(variant_inputs)))
        vc = load_variants.load_vcf(
            variant_input.path,
            filter=filter,
            ensembl_version=ensembl_version)
        for variant in vc:
            result[variant].add(variant_input.name)
    return result