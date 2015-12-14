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

import pyfaidx

def variant_context(
        reference_fasta,
        contig,
        inclusive_start,
        inclusive_end,
        alt,
        context_length):  
    """
    Retrieve the surronding reference region from a variant.

    SNVs are canonicalized so the reference base is a pyrmidine (C/T). For
    indels the reverse complement will still be taken if the first base of
    the reference is not a pyrmidine, but since the reference will also be
    reversed, that doesn't guarantee it will start with a pyrmidine.

    Parameters
    ----------
    reference_fasta : FastaReference
        reference sequence from pyfaidx package

    contig : str
        Chromosome of the variant

    inclusive_start : int
        start of the variant in 1-based inclusive coordinates

    inclusive_end : int
        end of the variant in 1-based inclusive coordinates

    alt : string
        alt sequence

    context_length : int
        number of bases on either side of the variant to return

    Returns
    ---------
    A tuple of (5', mutation, 3') where
        5' - bases immediately 5 prime to the mutation
        
        3' - bases immediately 3 prime to the mutation
        
        mutation - the ref sequence followed by a > character followed by the
            the alt sequence
    """

    # Move from 1-base coorindates to 0-base coordinates
    start = int(inclusive_start) - 1
    end = int(inclusive_end)

    full_sequence = reference_fasta[contig]

    left = str(full_sequence[start - context_length:start].seq).upper()
    middle = str(full_sequence[start: end].seq).upper()
    right = str(full_sequence[end: end + context_length].seq).upper()

    # Complement and reverse the context if necessary so the ref base is a
    # pyrmidine (C/T)
    if middle[0] in ('A', 'G'):
        context_5prime = pyfaidx.complement(right)[::-1]
        context_3prime = pyfaidx.complement(left)[::-1]
        context_mutation = "%s>%s" % (
            pyfaidx.complement(middle)[::-1], pyfaidx.complement(alt)[::-1])
    else:
        context_5prime = left
        context_3prime = right
        context_mutation = "%s>%s" % (middle, alt)

    return (context_5prime, context_mutation, context_3prime)

