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

import math

import pandas
import varcode

# Peptide lengths affect predictions and must be part of both cache keys.
CACHED_BINDING_AFFINITIES = {}
BINDING_PREDICTORS = {}


def binding_affinities(variants, alleles, epitope_lengths=(8, 9, 10, 11)):
    """Return the best mutant-peptide affinity (nM) and input allele per variant."""
    import mhctools
    import topiary

    variants = list(variants)
    alleles = list(alleles)
    lengths = tuple(sorted(set(epitope_lengths)))
    if not alleles or not lengths:
        raise ValueError("At least one allele and peptide length are required")

    for allele in alleles:
        pending = [v for v in variants
                   if (v, allele, lengths) not in CACHED_BINDING_AFFINITIES]
        if not pending:
            continue
        predictor_key = (allele, lengths)
        if predictor_key not in BINDING_PREDICTORS:
            BINDING_PREDICTORS[predictor_key] = mhctools.NetMHCpan(
                [allele], default_peptide_lengths=list(lengths))
        predictions = topiary.TopiaryPredictor(
            models=BINDING_PREDICTORS[predictor_key],
            only_novel_epitopes=True,
        ).predict_from_variants(varcode.VariantCollection(pending))
        # Topiary returns long-form prediction rows, or an empty DataFrame
        # when no effects yield peptides.
        values = {}
        if len(predictions):
            frame = predictions
            affinity_rows = frame.loc[frame["kind"] == "pMHC_affinity"]
            values = affinity_rows.groupby("variant")["value"].min().to_dict()
        for variant in pending:
            CACHED_BINDING_AFFINITIES[(variant, allele, lengths)] = values.get(
                variant.short_description, float("nan"))

    rows = []
    for variant in variants:
        available = [
            (CACHED_BINDING_AFFINITIES[(variant, allele, lengths)], allele)
            for allele in alleles
            if math.isfinite(CACHED_BINDING_AFFINITIES[(variant, allele, lengths)])
        ]
        affinity, allele = min(available) if available else (float("nan"), None)
        rows.append((variant, affinity, allele))
    return pandas.DataFrame(
        rows, columns=["variant", "binding_affinity", "binding_allele"])
