# varlens
tools for looking at somatic variation across multiple samples

## variant support

Under development, not working yet. The `varlens-support` tool gives read counts supporting variants across a number of BAM files.

## convert cufflinks fpkm_tracking files to a format accepted by GSEA

The `varlens-fpkm2gsea` tool combines cufflinks fpkm_tracking files that look like this:

```
tracking_id class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	FPKM	FPKM_conf_lo	FPKM_conf_hi	FPKM_status
ENSG00000186092	-	-	ENSG00000186092	OR4F5	-	chr1:69090-70008	-	-	0       	0      	0       	OK
ENSG00000237613	-	-	ENSG00000237613	FAM138A	-	chr1:34553-36081	-	-	0.00847125	0	    0.0807049	OK
...
```

into a single file that gives the FPKM expression for each sample at each gene:

```
NAME	189_11_13_RNA	189_1_13_RNA_Replacement	189_6_14_RNA_Replacement	189_9_12_RNA_Replacement	189_Omentum_RNA	189_Right_Ovary_RNA	189_Tumor_Colon_RN	189_Uterus_RNA	_189_Left_Ovary_RNA
5S_rRNA	0.0         	0.0	                        0.0	                        0.0	                        0.0	            0.0	                0.0	                0.0             0.0
7SK	    0.344282341667	0.28889425	                0.313886916667          	0.290836583333	            0.238506258333	0.163143416667  	0.23626425      	0.443587633333	0.2186565
...
```

This can be given to the Broad's [GSEA](http://www.broadinstitute.org/gsea/index.jsp) tool as a ".txt" format file.

Example:

```
varlens-fpkm2gsea \
    --input cufflinks-all-samples/*/genes.fpkm_tracking \
    --out cufflinks_fpkms_all_samples.txt
```
