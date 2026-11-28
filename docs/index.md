# varlens

Varlens provides command-line and Python utilities to combine and annotate
genomic variants from VCF/CSV files, filter and inspect aligned reads from
BAM/SAM files, count allele-supporting reads at genomic sites, and reconstruct
supplied fusion RNA transcripts.

It builds on [Varcode](https://github.com/openvax/varcode) (variants and
effects) and [pysam](https://github.com/pysam-developers/pysam) (reads), with
[Topiary](https://github.com/openvax/topiary) and
[MHCtools](https://github.com/openvax/mhctools) for MHC binding predictions and
[Isovar](https://github.com/openvax/isovar) for fusion RNA reconstruction.

## Installation

Varlens requires Python 3.11 or newer, pandas 3, Varcode 11.2+ (<14),
DataCache 1.38+, MHCtools 4.1.1+, and Topiary 5.109+:

```bash
python -m pip install varlens
```

Install the optional `fusion` extra to enable `varlens fusions` with Isovar:

```bash
python -m pip install 'varlens[fusion]'
```

Effect and gene annotations require a PyEnsembl release matching your genome
build (for example `pyensembl install --release 75 --species human` for GRCh37).

## Command-line interface

Varlens ships a unified `varlens` command with four subcommands:

| Subcommand | Summary |
| --- | --- |
| `varlens variants` | Combine, filter, and annotate variants from VCF/CSV files |
| `varlens reads` | Display, filter, and copy reads from BAM/SAM files |
| `varlens allele-support` | Count reads supporting each allele at genomic sites |
| `varlens fusions` | Reconstruct supplied fusion RNA and predict junction peptides |

See the [Command guide](commands.md) for worked examples,
[Structural variants](structural_variants.md) for symbolic SV and breakend
handling, and [Fusion RNA](fusions.md) for `varlens fusions`.

### Quick start

Combine variants from two VCF files with effect and gene annotations:

```bash
varlens variants tests/data/CELSR1/vcfs/vcf_1.vcf tests/data/CELSR1/vcfs/vcf_2.vcf \
    --include-effect \
    --include-gene
```

Count reads supporting each variant's reference and alternate alleles in a BAM:

```bash
varlens variants tests/data/CELSR1/vcfs/vcf_1.vcf \
    --include-read-evidence \
    --reads tests/data/CELSR1/bams/bam_1.bam \
    --min-mapping-quality 10
```

Filter reads overlapping a locus (1-based inclusive) and write a sorted BAM:

```bash
varlens reads tests/data/CELSR1/bams/bam_0.bam \
    --min-mapping-quality 30 \
    --locus 22:46932040-46932050 \
    --out /tmp/result.bam
```

Count reads supporting each allele at explicit genomic sites:

```bash
varlens allele-support tests/data/CELSR1/bams/bam_1.bam \
    --locus 22:46931061 22:46931063
```

Reconstruct a supplied fusion transcript (requires `varlens[fusion]`):

```bash
varlens fusions fusion.json.gz --out reconstruction.json
```

### Shell tab completion

`varlens` supports tab completion via
[argcomplete](https://kislyuk.github.io/argcomplete/). Enable it globally once
or register `varlens` in your shell startup file:

```bash
activate-global-python-argcomplete
# or:
eval "$(register-python-argcomplete varlens)"
```

## Coordinate systems

Varlens outputs 0-based half-open (*interbase*) coordinates (`interbase_start`,
`interbase_end`). On the command line, `--locus chr22:10-20` uses 1-based
inclusive coordinates (matching samtools and VCF), while `--locus chr22/9-20`
uses 0-based half-open interbase coordinates for the same 11-base span.

Insertions use an empty interbase interval at the insertion boundary (for
example VCF `1:100 A>AT` becomes `interbase_start=100, interbase_end=100` with
`ref=""` and `alt="T"`).

## Python API

`varlens.support.variant_support(variants, counts)` returns a pandas DataFrame
indexed by variant with two column levels (`measurement` and `source`) across
`num_alt`, `num_ref`, `num_other`, `total_depth`, `alt_fraction`, and
`any_alt_fraction`:

```python
from varlens.support import variant_support

evidence = variant_support(variants, counts)
alt_counts = evidence["num_alt"]
tumor_depth = evidence[("total_depth", "tumor")]
```

MHC binding helpers (`varlens.mhc_binding.binding_affinities`,
`varlens.mhc_binding.fragment_predictions`, and `varlens.fusions.analyze_fusion`)
use the Topiary API with any configured `mhctools` predictor. See the
[Reference](reference.md) for the full CLI and Python module APIs.
