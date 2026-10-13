import pandas
import pytest
from varcode import Variant

from varlens.support import EXPECTED_COLUMNS, variant_support


def test_variant_support_dataframe():
    first = Variant("22", 101, "A", "C", genome="GRCh37")
    second = Variant("22", 102, "G", "T", genome="GRCh37")
    counts = pandas.DataFrame(
        [
            ["tumor", "chr22", "100", "101", "A", 7],
            ["tumor", "chr22", "100", "101", "C", 3],
            ["tumor", "chr22", "100", "101", "G", 2],
            ["normal", "22", "100", "101", "A", 5],
            ["tumor", "22", "101", "102", "G", 4],
            ["normal", "22", "101", "102", "N", 0],
        ],
        columns=EXPECTED_COLUMNS,
    )
    original = counts.copy(deep=True)
    result = variant_support(iter([second, first]), counts)
    assert result.index.tolist() == [second, first]
    assert result.columns.names == ["measurement", "source"]
    assert result["num_alt"]["tumor"].tolist() == [0, 3]
    assert result["num_ref"]["normal"].tolist() == [0, 5]
    assert result.loc[first, ("num_other", "tumor")] == 2
    assert result.loc[first, ("total_depth", "tumor")] == 12
    assert result.loc[first, ("alt_fraction", "tumor")] == pytest.approx(3 / 12)
    assert result.loc[first, ("any_alt_fraction", "tumor")] == pytest.approx(5 / 12)
    assert result.loc[second, ("alt_fraction", "normal")] == 0
    pandas.testing.assert_frame_equal(counts, original)


def test_missing_and_empty_variant_support():
    variant = Variant("22", 101, "A", "C", genome="GRCh37")
    counts = pandas.DataFrame(
        [
            ["tumor", "22", 200, 201, "A", 7],
        ],
        columns=EXPECTED_COLUMNS,
    )
    with pytest.raises(ValueError, match="No allele counts"):
        variant_support([variant], counts)
    result = variant_support([variant], counts, ignore_missing=True)
    assert (result == 0).all().all()
    assert variant_support([], counts).empty
    assert variant_support([variant], counts.iloc[:0]).empty
