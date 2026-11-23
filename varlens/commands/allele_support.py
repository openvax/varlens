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

'''
Given one or more BAMs and some genomic sites to consider, write a csv file
giving counts of reads supporting each allele at each site for each BAM.

The genomic sites to consider may be specified by locus (--locus option), or via
one or more VCF files.

The positions outputted by this command are in *interbase coordinates*, i.e.
starting at 0, inclusive on first index, exclusive on second (as opposed to
the one-based inclusive coordinates used in VCF files).

Example:

%(prog)s \\
    --reads tests/data/CELSR1/bams/bam_1.bam \\
    --locus 22:46931061 22:46931063

'''

import argparse
import csv
import sys
import logging

from .. import loci_util
from .. import reads_util
from .. import variants_util

from . import configure_logging
from .arguments import add_verbose_arg, add_version_arg
from .. import support
from ..errors import input_errors

parser = argparse.ArgumentParser(usage=__doc__)
add_version_arg(parser)
group = parser.add_argument_group("output arguments")
group.add_argument("--out",
    help="Output CSV file. If not specified, the CSV is written to stdout.")
add_verbose_arg(group)
loci_util.add_args(parser.add_argument_group("loci specification"))
variants_util.add_args(parser)
reads_util.add_args(parser)

@input_errors(parser)
def run(raw_args=None):
    args = parser.parse_args(raw_args)
    configure_logging(args)

    loci = loci_util.load_from_args(args)  # may be None
    variants_df = variants_util.load_from_args_as_dataframe(args)
    if variants_df is not None:
        loci = variants_util.union_with_variant_loci(
            loci, variants_df.variant, "allele support")

    if not loci:
        if variants_df is not None:
            parser.error("No loci: variants specified but none remained "
                "after filtering")
        else:
            parser.error("No genomic loci or variants specified.")

    logging.info("Loaded %d genomic loci." % len(loci))

    read_sources = reads_util.load_from_args(args)

    if read_sources is None:
        parser.error("No read sources (--reads argument) specified.")

    out_fd = open(args.out, "w") if args.out else sys.stdout
    writer = csv.writer(out_fd)

    rows_generator = support.allele_support_rows(loci, read_sources)
    for (i, row) in enumerate(rows_generator):
        if i == 0:
            writer.writerow(row.index.tolist())
        writer.writerow([str(x) for x in row])

    if out_fd is not sys.stdout:
        out_fd.close()
        print("Wrote: %s" % args.out)
