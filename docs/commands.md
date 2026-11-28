# Command guide

Worked examples for `varlens variants`, `varlens reads`, and
`varlens allele-support`, with the output each prints. Run them from a git
checkout: they use the test data in `tests/data`, and the effect and gene
annotations need `pyensembl install --release 75 --species human` (GRCh37).
Run `varlens <subcommand> --help` for the full option list.

## `varlens variants`

Given variants from one or more VCF or CSV files, apply filters, add additional
columns, and output to CSV. (Currently only CSV output is supported, not VCF.)

Additional annotations can be added with `--include-*` options such as
`--include-gene` or `--include-effect`.

Print basic info for the variants found in two VCF files. Variants found in
both files are listed in one row, and the `sources` column lists the files each
variant was found in:

```bash
$ varlens variants tests/data/CELSR1/vcfs/vcf_1.vcf tests/data/CELSR1/vcfs/vcf_2.vcf

genome,contig,interbase_start,interbase_end,ref,alt,sources
GRCh37,22,21829554,21829555,T,G,1
GRCh37,22,45309892,45309893,T,G,2
GRCh37,22,46931059,46931060,A,C,1
GRCh37,22,46931061,46931062,G,A,1 2
GRCh37,22,50636217,50636218,A,C,1
GRCh37,22,50875932,50875933,A,C,1
```

Include Varcode variant effect annotations and overlapping gene names, and
write to a file:

```bash
$ varlens variants tests/data/CELSR1/vcfs/vcf_1.vcf tests/data/CELSR1/vcfs/vcf_2.vcf \
    --include-effect \
    --include-gene \
    --out /tmp/result.csv

Wrote: /tmp/result.csv

$ cat /tmp/result.csv
genome,contig,interbase_start,interbase_end,ref,alt,sources,effect,effect_candidates,gene
GRCh37,22,21829554,21829555,T,G,1,non-coding-transcript,NoncodingTranscript: non-coding-transcript,PI4KAP2
GRCh37,22,45309892,45309893,T,G,2,p.T214P,Substitution: p.T214P,PHF21B
GRCh37,22,46931059,46931060,A,C,1,p.S670A,Substitution: p.S670A,CELSR1
GRCh37,22,46931061,46931062,G,A,1 2,p.S669F,Substitution: p.S669F,CELSR1
GRCh37,22,50636217,50636218,A,C,1,intronic,Intronic: intronic,TRABD
GRCh37,22,50875932,50875933,A,C,1,p.A562fs,FrameShift: p.A562fs; Unresolved: unresolved-intron-retention; Unresolved: unresolved-cryptic; Intronic: intronic,PPP6R2
```

`effect` is the producer-preferred consequence on the selected transcript.
`effect_candidates` lists that consequence and its alternatives in Varcode's
order, including unresolved outcomes. For example, an intronic primary can
retain a protein-changing exon-skip deletion. These are predicted alternatives,
not observed products or calibrated probabilities. MHC predictions use Topiary's
shared candidate-aware selection and can therefore use a coding alternative to
that primary.

Print counts of reads supporting reference/variant/other alleles from a BAM,
counting only reads with mapping quality >= 10:

```bash
$ varlens variants tests/data/CELSR1/vcfs/vcf_1.vcf \
    --include-read-evidence \
    --reads tests/data/CELSR1/bams/bam_1.bam \
    --min-mapping-quality 10

genome,contig,interbase_start,interbase_end,ref,alt,sources,num_alt,num_ref,total_depth
GRCh37,22,21829554,21829555,T,G,vcf_1,0.0,0.0,0.0
GRCh37,22,46931059,46931060,A,C,vcf_1,0.0,216.0,320.0
GRCh37,22,46931061,46931062,G,A,vcf_1,0.0,321.0,321.0
GRCh37,22,50636217,50636218,A,C,vcf_1,0.0,0.0,0.0
GRCh37,22,50875932,50875933,A,C,vcf_1,0.0,0.0,0.0
```

## `varlens reads`

Filter reads from one or more BAMs and output a CSV or a new BAM. Loci and VCF
files may be specified to restrict reads to those overlapping specific loci or
variants.

Print basic fields for the reads in a BAM (`source` names the BAM each read
came from):

```bash
$ varlens reads tests/data/CELSR1/bams/bam_0.bam

source,query_name,reference_start,reference_end,cigarstring
bam_0,HISEQ:142:C5822ANXX:3:2116:16538:101199,46929962,46930062,100M
bam_0,HISEQ:142:C5822ANXX:3:1106:18985:32932,46929964,46930064,100M
...
```

Filter to reads aligned on the (-) strand, write to a CSV file, and include
mapping quality and aligned sequence:

```bash
$ varlens reads tests/data/CELSR1/bams/bam_0.bam \
    --is-reverse \
    --field mapping_quality query_alignment_sequence \
    --out /tmp/result.csv

Wrote: /tmp/result.csv
```

Write a BAM consisting of reads with mapping quality >= 30 overlapping a locus:

```bash
$ varlens reads tests/data/CELSR1/bams/bam_0.bam \
    --min-mapping-quality 30 \
    --locus 22:46932040-46932050 \
    --out /tmp/result.bam

Sorting read file /tmp/result.bam
Wrote 65 reads: /tmp/result.bam
```

Write a BAM consisting of reads overlapping variants from a VCF:

```bash
$ varlens reads tests/data/CELSR1/bams/bam_0.bam \
    --variants tests/data/CELSR1/vcfs/vcf_1.vcf \
    --out /tmp/result.bam

Sorting read file /tmp/result.bam
Wrote 51 reads: /tmp/result.bam
```

Print the header for a BAM in CSV format:

```bash
$ varlens reads tests/data/CELSR1/bams/bam_0.bam --header

read_source,group,index,key,value
bam_0,HD,0,VN,
bam_0,HD,1,GO,
bam_0,HD,2,SO,
bam_0,SQ,0,SN,chrM
...
```

## `varlens allele-support`

Given one or more BAMs and genomic sites (`--locus` or `--variants`), write a
CSV giving counts of reads supporting each allele at each site for each BAM.

Output positions are in *interbase coordinates* (0-based half-open; see
[Coordinate systems](index.md#coordinate-systems)). Insertions use an empty
interval at the insertion boundary: VCF `1:100 A>AT` becomes
`interbase_start=100, interbase_end=100` with `ref=""` and `alt="T"`. The
pileup reads the base immediately to the left and counts only the inserted
sequence.

```bash
$ varlens allele-support \
    tests/data/CELSR1/bams/bam_1.bam \
    --locus 22:46931061 22:46931063

source,contig,interbase_start,interbase_end,allele,count
bam_1,22,46931060,46931061,G,329
bam_1,22,46931060,46931061,,1
bam_1,22,46931062,46931063,A,327
bam_1,22,46931062,46931063,AG,2
bam_1,22,46931062,46931063,AC,1
```

BAM paths may also be passed with `--reads`.
