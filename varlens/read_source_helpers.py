import typechecks

from . import evaluation


def evaluate_read_expression(
        expression,
        alignment,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        bindings = evaluation.EvaluationEnvironment(
            [alignment],
            extra=extra_bindings)
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(alignment)


def evaluate_pileup_element_expression(
        expression,
        collection,
        pileup,
        element,
        error_value=evaluation.RAISE,
        extra_bindings={}):

    if typechecks.is_string(expression):
        bindings = evaluation.EvaluationEnvironment(
            [element, element.alignment, pileup],
            extra={
                'element': element,
                'pileup': pileup,
                'collection': collection,
            })
        return evaluation.evaluate_expression(
            expression,
            bindings,
            error_value=error_value)
    else:
        return expression(pileup)


def alignment_key(pysam_alignment_record):
    '''
    Return the identifying attributes of a `pysam.AlignedSegment` instance.
    This is necessary since these objects do not support a useful notion of
    equality (they compare on identify by default).
    '''
    return (
        read_key(pysam_alignment_record),
        pysam_alignment_record.query_alignment_start,
        pysam_alignment_record.query_alignment_end,
    )


def read_key(pysam_alignment_record):
    '''
    Given a `pysam.AlignedSegment` instance, return the attributes identifying
    the *read* it comes from (not the alignment). There may be more than one
    alignment for a read, e.g. chimeric and secondary alignments.
    '''
    return (
        pysam_alignment_record.query_name,
        pysam_alignment_record.is_duplicate,
        pysam_alignment_record.is_read1,
        pysam_alignment_record.is_read2,
    )


def flatten_header(header):
    for (group, rows) in header.items():
        for (index, row) in enumerate(rows):
            if not isinstance(row, dict):
                key_values = [(row, "")]
            else:
                key_values = row.items()
            for (key, value) in key_values:
                yield (str(group), index, str(key), str(value))
