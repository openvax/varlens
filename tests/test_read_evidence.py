# Copyright (c) 2015. Mount Sinai School of Medicine
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
import pytest
import pysam

from pysam import Samfile

from varcode import Variant as VarcodeVariant

from . import data_path

from varlens.locus import Locus
from varlens.read_evidence import PileupCollection

Variant = collections.namedtuple("Variant", "locus ref alt")


def filtered_read_names(filtered_evidence):
    assert filtered_evidence.parent is not None
    full = set(filtered_evidence.parent.read_attribute("query_name"))
    filtered = set(filtered_evidence.read_attribute("query_name"))
    assert filtered - full == set()
    return full - filtered


def test_raw_contig_names_remain_distinct(tmp_path):
    path = tmp_path / "raw-contigs.bam"
    sequences = {"chr1": "AC", "1": "TG"}
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": contig, "LN": 10} for contig in sequences],
    }
    with Samfile(str(path), "wb", header=header) as bam:
        for reference_id, (contig, sequence) in enumerate(sequences.items()):
            read = pysam.AlignedSegment()
            read.query_name = contig
            read.query_sequence = sequence
            read.reference_id = reference_id
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigar = [(0, len(sequence))]
            read.query_qualities = [30] * len(sequence)
            bam.write(read)
    pysam.index(str(path))

    loci = [Locus(contig, 0, 2) for contig in sequences]
    evidence = PileupCollection.from_bam(str(path), loci, normalized_contig_names=False)
    assert set(evidence.loci()) == {
        Locus(contig, position, position + 1)
        for contig in sequences
        for position in range(2)
    }
    for locus in loci:
        sequence = sequences[locus.contig]
        first_base = Locus(locus.contig, 0, 1)
        assert evidence.allele_summary(locus) == [(sequence, 1)]
        assert evidence.at(first_base).allele_summary(first_base) == [(sequence[0], 1)]
        assert evidence.pileup(first_base).locus == first_base
        assert evidence.allele_summary(Locus(locus.contig, 1, 1)) == [("", 1)]
        assert evidence.match_summary(Variant(locus, sequence, "GG")) == [
            (sequence, 1),
            ("GG", 0),
        ]

    # Subsets must retain raw-name lookup even when only the alias is present.
    subset = evidence.at(Locus("1", 0, 1))
    with pytest.raises(KeyError):
        subset.allele_summary(Locus("chr1", 0, 1))
    with pytest.raises(KeyError):
        subset.pileup(Locus("chr1", 0, 1))


def test_normalized_contig_alias_queries():
    raw_locus = Locus("chr17", 41244935, 41244936)
    normalized_locus = Locus("17", 41244935, 41244936)
    evidence = PileupCollection.from_bam(
        data_path("rna_chr17_41244936.bam"), [raw_locus]
    )
    assert evidence.loci() == [normalized_locus]
    assert evidence.pileup(raw_locus) is evidence.pileup(normalized_locus)
    assert evidence.allele_summary(raw_locus) == evidence.allele_summary(
        normalized_locus
    )
    assert evidence.at(raw_locus).loci() == [normalized_locus]
    assert evidence.at(raw_locus).pileup(raw_locus) is evidence.pileup(normalized_locus)


def test_read_evidence_rna1_single_base_loci():
    loci = [
        Locus.from_inclusive_coordinates("17", 41244936, 41244936),  # 0
        Locus.from_inclusive_coordinates("17", 41244937, 41244937),  # 1
        Locus.from_inclusive_coordinates("17", 41244935, 41244935),  # 2
        Locus.from_inclusive_coordinates("17", 41244933, 41244933),  # 3
        Locus.from_inclusive_coordinates("17", 41244853, 41244853),  # 4
        Locus.from_inclusive_coordinates("17", 41244857, 41244857),  # 5
        Locus.from_inclusive_coordinates("17", 41244864, 41244864),  # 6
        Locus.from_inclusive_coordinates("17", 41244879, 41244879),  # 7
        Locus.from_inclusive_coordinates("17", 41244901, 41244901),  # 8
        Locus.from_inclusive_coordinates("17", 41244910, 41244910),  # 9
        Locus.from_inclusive_coordinates("17", 41244917, 41244917),  # 10
        Locus.from_inclusive_coordinates("17", 41244972, 41244972),  # 11
        Locus.from_inclusive_coordinates("17", 41244973, 41244973),  # 12
        Locus.from_inclusive_coordinates("17", 41245026, 41245026),  # 13
        Locus.from_inclusive_coordinates("17", 41245027, 41245027),  # 14
        Locus.from_inclusive_coordinates("17", 41245019, 41245019),  # 15
        Locus.from_inclusive_coordinates("17", 41245018, 41245018),  # 16
    ]
    evidence = PileupCollection.from_bam(data_path("rna_chr17_41244936.bam"), loci)

    assert (evidence.allele_summary(loci[0])) == ([("A", 11), ("G", 7)])
    assert (evidence.allele_summary(loci[1])) == ([("G", 17), ("A", 1)])
    assert (evidence.allele_summary(loci[2])) == ([("C", 18)])
    assert (evidence.allele_summary(loci[3])) == ([("A", 17), ("G", 1)])
    assert (evidence.allele_summary(loci[4])) == ([("C", 1)])
    assert (evidence.allele_summary(loci[5])) == ([("T", 2)])
    assert (evidence.allele_summary(loci[6])) == ([("T", 4)])
    assert (evidence.allele_summary(loci[7])) == ([("C", 8)])
    assert (evidence.allele_summary(loci[8])) == ([("C", 8)])
    assert (evidence.allele_summary(loci[9])) == ([("C", 9)])
    assert (evidence.allele_summary(loci[10])) == ([("A", 10)])
    assert (evidence.allele_summary(loci[11])) == ([("T", 11)])
    assert (evidence.allele_summary(loci[12])) == ([("T", 11)])
    assert (evidence.allele_summary(loci[13])) == ([("C", 1)])
    assert (evidence.allele_summary(loci[14])) == ([("G", 1)])
    assert (evidence.allele_summary(loci[15])) == ([("T", 8)])
    assert (evidence.allele_summary(loci[16])) == ([("T", 8)])


def test_read_evidence_rna1_multi_base_loci():
    loci = [
        Locus.from_inclusive_coordinates("17", 41244853, 41244854),  # 0
        Locus.from_inclusive_coordinates("17", 41244853, 41244857),  # 1
        Locus.from_inclusive_coordinates("17", 41244854, 41244857),  # 2
        Locus.from_inclusive_coordinates("17", 41244852, 41244857),  # 3
        Locus.from_inclusive_coordinates("17", 41244933, 41244936),  # 4
        Locus.from_inclusive_coordinates("17", 41244933, 41244937),  # 5
        Locus.from_inclusive_coordinates("17", 41244971, 41244973),  # 6
        Locus.from_inclusive_coordinates("17", 41265063, 41265067),  # 7
    ]
    evidence = PileupCollection.from_bam(data_path("rna_chr17_41244936.bam"), loci)
    assert (evidence.allele_summary(loci[0])) == ([("CT", 1)])
    assert (evidence.allele_summary(loci[1])) == ([("CTTTT", 1)])
    assert (evidence.allele_summary(loci[2])) == ([("TTTT", 1)])
    assert (evidence.allele_summary(loci[3])) == ([])
    assert (evidence.allele_summary(loci[4])) == (
        [("AACA", 11), ("AACG", 6), ("GACG", 1)]
    )
    assert (evidence.allele_summary(loci[5])) == (
        [("AACAG", 10), ("AACGG", 6), ("AACAA", 1), ("GACGG", 1)]
    )
    assert (evidence.allele_summary(loci[6])) == ([("ATT", 11)])
    assert (evidence.allele_summary(loci[7])) == ([("ACCCG", 1)])


def test_read_evidence_gatk_mini_bundle_extract():
    loci = [
        Locus.from_inclusive_coordinates("20", 9999996, 9999996),  # 0
        Locus.from_inclusive_coordinates("20", 10260442),  # 1
        Locus.from_inclusive_coordinates("20", 10006823),  # 2
        Locus.from_inclusive_coordinates("20", 10006819, 10006823),  # 3
        Locus.from_inclusive_coordinates("20", 10006819, 10006825),  # 4
        Locus.from_inclusive_coordinates("20", 10006822, 10006827),  # 5
        Locus.from_inclusive_coordinates("20", 10007175),  # 6
        Locus.from_inclusive_coordinates("20", 10007174, 10007176),  # 7
        Locus.from_inclusive_coordinates("20", 1, 3),  # 8
        Locus.from_inclusive_coordinates("20", 10008796),  # 9
        Locus.from_inclusive_coordinates("20", 10008921),  # 10
    ]
    handle = Samfile(data_path("gatk_mini_bundle_extract.bam"))
    evidence = PileupCollection.from_bam(handle, loci)

    assert (evidence.allele_summary(loci[0])) == ([("ACT", 9)])
    assert (evidence.filter(drop_duplicates=True).allele_summary(loci[0])) == (
        [("ACT", 8)]
    )
    assert (evidence.allele_summary(loci[1])) == ([("T", 7)])
    assert (evidence.filter().allele_summary(loci[2])) == ([("", 6), ("C", 2)])
    assert (
        evidence.filter(drop_duplicates=True, min_base_quality=50).allele_summary(
            loci[2]
        )
    ) == ([])
    assert (evidence.filter(drop_duplicates=True).allele_summary(loci[2])) == (
        [("", 5), ("C", 1)]
    )
    assert (
        evidence.filter(drop_duplicates=True, min_mapping_quality=60).allele_summary(
            loci[2]
        )
    ) == ([("", 5), ("C", 1)])
    assert (
        evidence.filter(drop_duplicates=True, min_mapping_quality=61).allele_summary(
            loci[2]
        )
    ) == ([("", 2)])
    assert (
        evidence.filter(drop_duplicates=True, min_mapping_quality=61).allele_summary(
            loci[3]
        )
    ) == ([("A", 2)])
    assert (
        evidence.filter(drop_duplicates=True, min_mapping_quality=61).allele_summary(
            loci[4]
        )
    ) == ([("AAA", 2)])
    assert (
        evidence.filter(drop_duplicates=True, min_mapping_quality=61).allele_summary(
            loci[5]
        )
    ) == ([("AAAC", 2)])
    assert (evidence.filter().allele_summary(loci[6])) == ([("T", 5), ("C", 3)])
    assert (evidence.filter(min_base_quality=30).allele_summary(loci[6])) == (
        [("T", 4), ("C", 3)]
    )
    assert (evidence.filter().allele_summary(loci[7])) == ([("CTT", 5), ("CCT", 3)])
    assert (evidence.filter(min_base_quality=30).allele_summary(loci[7])) == (
        [("CTT", 3), ("CCT", 2)]
    )
    assert (evidence.filter(min_base_quality=32).allele_summary(loci[2])) == (
        [("", 6), ("C", 1)]
    )
    assert (filtered_read_names(evidence.at(loci[2]).filter(min_base_quality=32))) == (
        {"20GAVAAXX100126:4:3:18352:43857"}
    )
    assert (evidence.allele_summary(loci[8])) == ([])
    assert (evidence.filter(drop_duplicates=True).allele_summary(loci[8])) == ([])
    with pytest.raises(KeyError):
        evidence.allele_summary(
            Locus.from_inclusive_coordinates("20", 10009174, 10009176)
        )
    assert (
        filtered_read_names(evidence.at(loci[9]).filter(drop_improper_mate_pairs=True))
    ) == ({"20FUKAAXX100202:8:68:1530:49310"})
    assert (len(evidence.at(loci[8]).read_attribute("mapping_quality"))) == (0)
    assert (list(evidence.at(loci[9]).read_attribute("mapping_quality"))) == (
        list(evidence.at(loci[9]).read_attributes().mapping_quality)
    )
    assert (evidence.filter(drop_duplicates=True).allele_summary(loci[10])) == (
        [("C", 2), ("CA", 1), ("CAA", 1)]
    )
    assert (
        evidence.filter(drop_duplicates=True).allele_summary(
            Locus.from_interbase_coordinates(
                loci[10].contig, loci[10].end, loci[10].end
            )
        )
    ) == ([("", 2), ("A", 1), ("AA", 1)])


def test_read_evidence_variant_matching_gatk_mini_bundle_extract():
    handle = Samfile(data_path("gatk_mini_bundle_extract.bam"))

    loci = [
        Locus.from_inclusive_coordinates("20", 10008951),  # 0
        Locus.from_inclusive_coordinates("20", 10009053),  # 1
        Locus.from_inclusive_coordinates("20", 10009053, 10009054),  # 2
        Locus.from_inclusive_coordinates("20", 10006822),  # 3
        Locus.from_inclusive_coordinates("20", 10006822, 10006823),  # 4
    ]
    evidence = PileupCollection.from_bam(handle, loci)

    assert (evidence.match_summary(Variant(loci[0], "A", "C"))) == (
        [("A", 1), ("C", 4)]
    )
    assert (
        evidence.filter(drop_duplicates=True).match_summary(Variant(loci[0], "A", "C"))
    ) == ([("A", 0), ("C", 3)])
    assert (evidence.match_summary(Variant(loci[1], "A", "C"))) == (
        [("A", 3), ("C", 0)]
    )
    assert (evidence.match_summary(Variant(loci[1], "A", "CC"))) == (
        [("A", 3), ("CC", 0)]
    )
    assert (evidence.match_summary(Variant(loci[1], "A", ""))) == ([("A", 3), ("", 0)])
    assert (evidence.match_summary(Variant(loci[1], "A", ""))) == ([("A", 3), ("", 0)])
    assert (evidence.match_summary(Variant(loci[2], "AT", ""))) == (
        [("AT", 3), ("", 0)]
    )
    assert (evidence.match_summary(Variant(loci[3], "A", ""))) == ([("A", 2), ("", 6)])
    assert (evidence.match_summary(Variant(loci[4], "AC", ""))) == (
        [("AC", 2), ("", 6)]
    )
    assert (
        evidence.match_summary(
            Variant(loci[4], "AC", ""),
            lambda e: e.read_attributes().mapping_quality.mean(),
        )
    ) == ([("AC", 60.0), ("", 65.0)])


def test_read_evidence_variant_matching_gatk_bundle_native_varcode_variant():
    # Try native varcode Variant.
    handle = Samfile(data_path("gatk_mini_bundle_extract.bam"))
    locus = Locus.from_inclusive_coordinates("20", 10008951)
    variant = VarcodeVariant(
        locus.contig,
        locus.position + 1,  # inclusive not interbase
        "A",
        "C",
    )
    evidence = PileupCollection.from_bam(handle, [variant])
    assert (evidence.match_summary(variant)) == ([("A", 1), ("C", 4)])


def test_read_evidence_variant_matching_gatk_mini_bundle_extract_warning():
    filename = data_path("gatk_mini_bundle_extract.bam")

    # Should log a warning but pass.
    loci = [
        Locus.from_inclusive_coordinates("20", 10009053, 10009054),  # 0
    ]
    evidence = PileupCollection.from_bam(filename, loci)
    assert (evidence.match_summary(Variant(loci[0], "A", ""))) == (
        [("A", 0), ("", 0), ("AT", 3)]
    )
