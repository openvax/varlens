# Fusion RNA with `varlens fusions`

`varlens fusions` reconstructs a fusion transcript that another tool has
already assembled, using Isovar, and can predict MHC binding of the
junction-overlapping peptides.

Install the optional Isovar adapter and reconstruct supplied fusion RNA:

```bash
pip install 'varlens[fusion]'
varlens fusions fusion.json.gz --out reconstruction.json
```

Optionally predict junction-overlapping peptides with Topiary and any mhctools
predictor:

```bash
varlens fusions fusion.json.gz --mhc-predictor mhcflurry \
    --mhc-alleles 'HLA-A*02:01' 'HLA-B*07:02' \
    --mhc-peptide-lengths 8-11 --out predictions.json
```

Any predictor supported by `mhctools` (such as `mhcflurry` or `netmhcpan`) can
be selected with `--mhc-predictor`; its runtime/model dependencies must be
installed for prediction runs. The reconstruction-only command does not invoke
any predictor. Multiple JSON/JSON.gz input files can be supplied, and `--out`
also accepts JSON.gz; omitted output goes to stdout. A failed input or
prediction prevents writing a partial report.

## Input and report schema

Each input is an [Isovar supplied-fusion input](https://github.com/openvax/isovar/blob/master/docs/fusion.md#input):
a JSON object with `fusion`, `references`, and `reads` (plus optional
`reference_names`). `fusion` supplies event identity, assembly, oriented
donor/acceptor breakpoints, RNA sequence, mapped blocks, half-open junction
coordinates, and provenance with a sample ID. References supply annotated
transcript sequence/mapping and CDS coordinates; reads supply observed RNA
sequence/mapping and sample/library/read/fragment identity. Isovar validates
these fields and transfers frames only from compatible annotated donors. The
`note` and `original_records` fields in audit envelopes are retained in the
report, but are not counted as additional RNA evidence. Unknown fields are
rejected.

The output `varlens.fusions.v1` report contains an `events` list. Each
`varlens.fusion.v1` event retains the entire input, the unmodified
`isovar.fusion_rna.v3` reconstruction, software versions, Topiary fragments,
and full peptide prediction rows. `prediction_status` distinguishes
`not_requested`, `no_coding_products`, and `completed`; a completed run may
have no target-overlapping predictions.

Unresolved and insufficient-support outcomes keep their reasons and measured
evidence. Ambiguous or partial-CDS coding hypotheses remain separate fragments
with explicit status on every prediction; Varlens does not choose an isoform.
Junction read/fragment counts retain their supplied-observation scope. Unknown
reference counts, protein-observation evidence, UMI/cell counts, and comparator
sequences are not filled with zero or inferred from alignment counts.
Junction overlap is not a claim of proteome-wide novelty or tumor specificity.
This command does not discover fusions, assemble arbitrary BAMs, translate
structural VCFs, or establish whole-sample expression from selected reads.

## Python API

Python callers can use the same workflow via [`analyze_fusion`](reference.md):

```python
from varlens.fusions import analyze_fusion

result = analyze_fusion(
    payload,
    alleles=["HLA-A*02:01"],
    mhc_predictor="mhcflurry",
)
```

`analyze_fusion(payload)` performs reconstruction only. A configured
`topiary.TopiaryPredictor` or `mhctools` predictor can be passed as `predictor=`
instead of `alleles` / `mhc_predictor`. `varlens.mhc_binding.fragment_predictions`
also accepts existing `topiary.ProteinFragment` objects and returns
target-overlapping prediction rows with source, sample, and evidence metadata
intact.
