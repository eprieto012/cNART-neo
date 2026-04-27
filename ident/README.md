# Dominant Clonotype Discovery

## Preprocessing

Raw 10x Genomics 5′ gene-expression and V(D)J FASTQ files were processed using `cellranger multi` (Cell Ranger v9.0.0) within a Singularity container. Sample demultiplexing was performed using On-Chip Multiplexing (OCM) barcodes, allowing pooled cells to be reassigned to their original biological samples.

The following reference datasets were used:

- `refdata-gex-GRCh38-2024-A` (gene expression)
- `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` (V(D)J annotation)

See `crsing_multi_9.0.0.sh` and `multiconfig.csv` for additional details on the Cell Ranger workflow.

## Analysis Workflow

Run `ITAG_TIL15_S02_m1.Rmd` to reproduce the analysis from Cell Ranger outputs to final clonotype frequency tables and plots.

This workflow includes:

- ambient RNA assessment with [SoupX](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html)
- import of V(D)J annotations into Seurat using [djvdj](https://rnabioco.github.io/djvdj/)
- sample merging
- quality-control assessment
- downstream filtering and cell annotation

## TCR Rescue with MiXCR

Cell Ranger did not recover TCR information for a subset of T cells despite detectable TCR gene expression. Missing assignments were rescued using the [MiXCR](https://mixcr.com/) `analyze 10x-sc-xcr-vdj-gemx-v3` workflow, executed with `script_mixcr_gex.sh`.

This workflow requires a `10x-sample-sheet.tsv` file, generated automatically within the MiXCR section of the main analysis script, to restrict analysis to Cell Ranger-defined cell-containing GEMs.

## Filtering

The following cells or features were excluded:

- doublets
- ribosomal genes
- TCR genes
- genes with very low expression
- cells with >20% mitochondrial reads (`percent.mt`)
- cells without recovered TCR information

## CD4+/CD8+ T-Cell Annotation

[MAGIC](https://magic.readthedocs.io/en/stable/) imputation was first applied to gene-expression data. CD4 and CD8A/CD8B expression thresholds were then estimated from imputed values using nonparametric density estimation and clustering with [`pdfCluster`](https://github.com/cran/pdfCluster/tree/master).

Final CD4+ and CD8+ annotations were assigned using [scGate](https://github.com/carmonalab/scGate) as the primary lineage classifier, with MAGIC-derived expression used as complementary evidence. When classifications disagreed, `scGate` was prioritised.

## TCR Analysis

Dominant clonotypes were identified by comparing TCR frequencies across conditions and cell subsets. Frequency plots and summary tables were generated for each sample.

Outputs include candidate clonotype sets covering 70% and 100% of cumulative clonotype frequency:

- `TCR_interrogated_top100.tsv`
- `TCR_interrogated_top70.tsv`