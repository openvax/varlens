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

import pandas
import varcode


def binding_affinities(variants, alleles, epitope_lengths=(8, 9, 10, 11)):
    """Return each variant's lowest mutant-peptide affinity (nM) and allele.

    Variants without an affinity prediction have a missing affinity and allele.
    Topiary's presentation scores are deliberately excluded from this summary.
    """
    import mhctools
    import topiary

    variants = list(variants)
    alleles = list(dict.fromkeys(alleles))
    lengths = tuple(epitope_lengths)
    if not alleles:
        raise ValueError("Specify at least one MHC allele")
    if not lengths or any(length <= 0 for length in lengths):
        raise ValueError("Epitope lengths must be positive")

    result = pandas.DataFrame(
        {
            "variant": variants,
            "binding_affinity": float("nan"),
            "binding_allele": pandas.Series([None] * len(variants), dtype=object),
        }
    )
    if not variants:
        return result
    if len({variant.reference_name for variant in variants}) > 1:
        raise ValueError("Cannot combine reference genomes in MHC predictions")

    model = mhctools.NetMHCpan(alleles=alleles, default_peptide_lengths=lengths)
    predictor = topiary.TopiaryPredictor(models=model, only_novel_epitopes=True)
    predictions = predictor.predict_from_variants(varcode.VariantCollection(variants))
    if predictions.empty:
        return result

    affinities = predictions.loc[
        predictions["kind"].eq("pMHC_affinity") & predictions["affinity"].notna()
    ]
    # Stable sorting makes ties deterministic without comparing Variant objects.
    best = affinities.sort_values(
        ["affinity", "allele"], kind="stable"
    ).drop_duplicates("variant")
    best = best.set_index("variant")
    # Current Topiary serializes the variant as its short description, while
    # Varlens' public result retains the caller's Variant objects and order.
    keys = result["variant"].map(lambda variant: variant.short_description)
    result["binding_affinity"] = keys.map(best["affinity"])
    result["binding_allele"] = keys.map(best["allele"])

    return result
