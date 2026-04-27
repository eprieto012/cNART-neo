# Antigen-Associated Clonotype Discovery

## Overview

This module identifies expanded antigen-associated TCR clonotypes from antigen-specific co-culture assays using single-cell 10x Genomics 5′ gene-expression and V(D)J data.

## Preprocessing

- Processes raw 10x Genomics 5′ gene-expression and V(D)J FASTQ files using `cellranger multi` (Cell Ranger v9.0.0) within a Singularity container.
- Demultiplexes pooled samples using On-Chip Multiplexing (OCM).

### Reference Datasets

- `refdata-gex-GRCh38-2024-A` — gene-expression reference.
- `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` — V(D)J reference.

Additional Cell Ranger details are provided in `crsing_multi_9.0.0.sh` and `multiconfig.csv`.

## Single-Cell Analysis: `ITAG_TIL15_S02_m1.Rmd`

- Assesses ambient RNA contamination with [SoupX](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html).
- Imports V(D)J annotations into Seurat using [djvdj](https://rnabioco.github.io/djvdj/).
- Merges samples and performs quality-control assessment.
- Filters low-quality cells, genes, doublets ([scDblFinder](https://github.com/plger/scDblFinder)), and cells without recovered TCR information.
- Identifies expanded antigen-associated clonotypes across conditions and cell subsets.

### TCR Rescue with MiXCR

Cell Ranger did not recover TCR information for a subset of T cells despite detectable TCR gene expression. Missing assignments were rescued using [MiXCR](https://mixcr.com/) through the `analyze 10x-sc-xcr-vdj-gemx-v3` workflow.

A `10x-sample-sheet.tsv` file, generated automatically within the pipeline, was used to restrict rescue analysis to Cell Ranger-defined cell-containing GEMs.

### CD4+/CD8+ T-Cell Annotation

CD4+ and CD8+ lineage labels were assigned using a combined strategy integrating expression imputation, unsupervised thresholding, and reference-based gating.

- [MAGIC](https://magic.readthedocs.io/en/stable/) was used for gene-expression imputation.
- CD4 and CD8A/CD8B thresholds were estimated from imputed values using [pdfCluster](https://github.com/cran/pdfCluster/tree/master).
- Final lineage labels were assigned using [scGate](https://github.com/carmonalab/scGate), prioritised when classifications disagreed.

## Output

Candidate clonotype sets covering 70% and 100% of cumulative clonotype frequency, subsequently used for cross-compartment clonotype tracking in Phase 2 (`track/`):

- `TCR_interrogated_top100.tsv`
- `TCR_interrogated_top70.tsv`