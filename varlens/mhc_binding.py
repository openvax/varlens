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

import pandas
import varcode

CACHED_BINDING_AFFINITIES = {}  # (variant, allele -> nm affinity)
BINDING_PREDICTORS = {}
def binding_affinities(variants, alleles, epitope_lengths=[8, 9, 10, 11]):
    # We import these here so we don't depend on these libraries unless this
    # function is called.
    import mhctools
    import topiary

    for allele in alleles:
        if allele not in BINDING_PREDICTORS:
            BINDING_PREDICTORS[allele] = mhctools.NetMHCpan(
                [allele], epitope_lengths=epitope_lengths)
        predictor = BINDING_PREDICTORS[allele]
        predictions = topiary.predict_epitopes_from_variants(
            varcode.VariantCollection([
                v for v in variants
                if (v, allele) not in CACHED_BINDING_AFFINITIES
            ]),
            predictor,
            ic50_cutoff=float('inf'),
            percentile_cutoff=100)
        if len(predictions) > 0:
            predictions_df = pandas.DataFrame(
                predictions.elements, columns=predictions[0]._fields)
            values = predictions_df.groupby("variant")["value"].min()
            for (variant, value) in zip(values.index, values):
                CACHED_BINDING_AFFINITIES[(variant, allele)] = value
        
    result_df = collections.defaultdict(list)
    for variant in variants:
        (binding_affinity, binding_allele) = min(
            (CACHED_BINDING_AFFINITIES.get((variant, allele), float('nan')),
                allele)
            for allele in alleles)
        if pandas.isnull(binding_affinity):
            binding_allele = None
        result_df["variant"].append(variant)
        result_df["binding_affinity"].append(binding_affinity)
        result_df["binding_allele"].append(binding_allele)

    return pandas.DataFrame(result_df)
    