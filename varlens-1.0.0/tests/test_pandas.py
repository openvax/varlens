from types import SimpleNamespace

import pandas
from varcode import Variant

from varlens.variant_includes import Gene, ReadEvidence
from . import data_path


def test_chunked_text_annotations_with_missing_csv_values():
    frame = pandas.DataFrame(
        {
            "variant": [
                SimpleNamespace(gene_names=["TP53"]),
                SimpleNamespace(gene_names=["BRAF"]),
            ],
            "gene": [float("nan"), float("nan")],
        },
        index=[4, 9],
    )
    assert list(Gene().compute(frame, chunk_rows=1)) == [1, 1]
    assert frame["gene"].tolist() == ["TP53", "BRAF"]
    assert list(Gene().compute(frame, chunk_rows=1)) == []


def test_read_source_table_with_nonconsecutive_rows():
    frame = pandas.DataFrame(
        {
            "variant": [Variant("22", 46931060, "A", "C", genome="GRCh37")],
            "sample": ["patient"],
        },
        index=[7],
    )
    sources = pandas.DataFrame(
        {
            "tumor": [data_path("CELSR1/bams/bam_0.bam")],
            "missing": [float("nan")],
        },
        index=pandas.Index(["patient"], name="sample"),
    )
    assert list(ReadEvidence(read_sources_df=sources).compute(frame, chunk_rows=1)) == [
        1
    ]
    assert frame.loc[7, "tumor_count_num_ref"] == 50
    assert frame.loc[7, "tumor_count_num_alt"] == 0
    assert pandas.isna(frame.loc[7, "missing_count_num_ref"])
