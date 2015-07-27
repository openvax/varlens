'''
Convert cufflinks fpkm_tracking files to txt files usable in the Broad's GSEA
tool.

%(prog)s \
    --input file1.fpkm_tracking [file2.fpkm_tracking ...]
    --name name1 [name2 ...]
    --out out.txt
'''

import argparse
import os
import collections

import pandas

import varcode
import varcode.cufflinks

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("--input", required=True, nargs="+",
    help="Input fpkm_tracking files from cufflinks, one per sample")
parser.add_argument("--name", nargs="+",
    help="Sample names. Default: filenames of input files.")
parser.add_argument("--out", required=True, help="Output path.")

def named_dataframes(inputs, names=None):
    if not names:
        prefix = os.path.commonprefix(inputs)
        suffix = os.path.commonprefix([x[::-1] for x in inputs])[::-1]
        names = [x[len(prefix):][:-len(suffix)] for x in inputs]

    seen = set()
    for i in range(len(names)):
        if names[i] in seen:
            num = 2
            while names[i] in seen:
                names[i] = "%s_%d" % (names[i], num)
                num += 1
        seen.add(names[i])

    if len(names) != len(inputs):
        raise ValueError("%d inputs but %d names" % (len(names, len(inputs))))

    result = collections.OrderedDict()
    for (i, (name, filename)) in enumerate(zip(names, inputs)):
        print("Opening file %d / %d: %s = %s" %
            (i + 1, len(names), name, filename))
        result[name] = varcode.cufflinks.load_cufflinks_tracking_file(filename)
    return result

def run():
    args = parser.parse_args()

    dataframes = named_dataframes(args.input, args.name)
    for (name, df) in dataframes.items():
        df["SAMPLE"] = name

        # We use just the first gene from the list. May want to change this!
        df["NAME"] = [lst[0] for lst in df.gene_names]

    print("Concatenating dataframes.")
    concat = pandas.concat(dataframes.values())

    print("Pivoting.")
    pivot = pandas.pivot_table(
        concat, values="fpkm", index="NAME", columns="SAMPLE")

    print("Writing.")
    pivot.to_csv(args.out, sep='\t')
    print("Wrote: %s" % args.out)
