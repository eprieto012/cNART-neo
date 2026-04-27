# Cross-Compartment Clonotype Tracking

## Overview

This module tracks antigen-associated TCR clonotypes across peripheral blood T-cell compartments and maps experimentally validated reactive clonotypes onto CD8+ T-cell phenotypic states.

## Preprocessing

- Processes raw 10x Genomics 5′ gene-expression and V(D)J FASTQ files using `cellranger multi` (Cell Ranger v9.0.0) within a Singularity container.
- Demultiplexes pooled samples using On-Chip Multiplexing (OCM).

### Reference Datasets

- `refdata-gex-GRCh38-2024-A` — gene-expression reference.
- `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` — V(D)J reference.

Additional Cell Ranger details are provided in `crsing_multi_9.0.0.sh` and `multiconfig.csv`.

## Single-Cell Analysis: `ITAG_TIL15_S01_m1.Rmd`

- Assesses ambient RNA contamination with [SoupX](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html).
- Imports V(D)J annotations into Seurat using [djvdj](https://rnabioco.github.io/djvdj/).
- Merges samples and performs quality-control assessment.
- Filters low-quality cells, genes, doublets ([scDblFinder](https://github.com/plger/scDblFinder)), and cells without recovered TCR information.
- Generates filtered and annotated single-cell objects for downstream tracking analyses.

### TCR Rescue with MiXCR

Cell Ranger did not recover TCR information for a subset of T cells despite detectable TCR gene expression. Missing assignments were rescued using [MiXCR](https://mixcr.com/) through the `analyze 10x-sc-xcr-vdj-gemx-v3` workflow.

A `10x-sample-sheet.tsv` file, generated automatically within the pipeline, was used to restrict rescue analysis to Cell Ranger-defined cell-containing GEMs.

### CD4+/CD8+ T-Cell Annotation

CD4+ and CD8+ lineage labels were assigned using a combined strategy integrating expression imputation, unsupervised thresholding, and reference-based gating.

- [MAGIC](https://magic.readthedocs.io/en/stable/) was used for gene-expression imputation.
- CD4 and CD8A/CD8B thresholds were estimated from imputed values using [`pdfCluster`](https://github.com/cran/pdfCluster/tree/master).
- Independent lineage labels were generated using [scGate](https://github.com/carmonalab/scGate).
- Final CD4+ or CD8+ annotations were retained only when both approaches were fully concordant.

## Cross-Compartment Analysis: `ITAG_TIL15_S01_m2.Rmd`

### Compartment Annotation

Cells were assigned to the following compartments:

- `CD4+Bulk`
- `CD4+CD39+`
- `CD4+CD39-`
- `CD8+Bulk`
- `CD8+CD39+`
- `CD8+CD39-`

### PDCD1 Thresholding

PDCD1 expression thresholds were estimated independently within each compartment using MAGIC-imputed values and density-based clustering with [`pdfCluster`](https://github.com/cran/pdfCluster/tree/master).

These thresholds were used to define PDCD1-high and PDCD1-low control clonotypes.

### CD8+ T-Cell State Analysis

CD8+ T cells were subsetted and reclustered to refine cell-state annotation.

Clusters were manually annotated using marker-gene expression and aggregate cluster-level profiles.

### TCR Candidate Prioritisation and Validation

Neoantigen-associated clonotypes were tracked across `CD8+Bulk`, `CD8+CD39+`, and `CD8+CD39-` compartments.

Candidate TCRs were prioritised based on:

- Presence in NeoAg stimulation pools from Phase 1 (`disc/`).
- Expansion across CD8+ compartments.
- PDCD1 expression levels.

A total of 26 candidate TCRs were reconstructed, transduced into healthy-donor PBMCs, and tested experimentally against patient-specific candidate (neo)antigens.

Validated TCR annotations are provided in `validated-tcrs.txt`.

## Output

Cross-compartment clonotype rankings, CD8+ state projections, and visualisation of validated reactive clonotypes.

<p align="center">
  <img src="plot.png" width="1000">
</p>