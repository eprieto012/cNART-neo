# Cross-Compartment Tracking and Validation

## Preprocessing

Raw 10x Genomics 5′ gene-expression and V(D)J FASTQ files were processed with `cellranger multi` (Cell Ranger v9.0.0) using a Singularity container. Sample demultiplexing was performed with On-Chip Multiplexing (OCM) barcodes, allowing pooled cells to be reassigned to their original biological samples.

The following references were used:

- `refdata-gex-GRCh38-2024-A` (gene expression)
- `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` (V(D)J annotation)

See `crsing_multi_9.0.0.sh` and `multiconfig.csv` for additional details on the Cell Ranger workflow.

## Analysis Workflow

Run `ITAG_TIL15_S02_m1.Rmd` to reproduce the preprocessing workflow from Cell Ranger outputs to filtered annotated single-cell objects.

Ambient RNA contamination was assessed with [SoupX](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html). V(D)J annotations were imported into Seurat using [djvdj](https://rnabioco.github.io/djvdj/), samples were merged, and quality-control metrics were evaluated before downstream filtering and annotation.

Run `ITAG_TIL15_S02_m2.Rmd` to reproduce the downstream cross-compartment clonotype tracking and prioritisation workflow.

## TCR Rescue with MiXCR

Cell Ranger failed to recover TCR information for some T cells despite detectable TCR gene expression. Missing TCR assignments were rescued using the [MiXCR](https://mixcr.com/) `analyze 10x-sc-xcr-vdj-gemx-v3` workflow, executed with `script_mixcr_gex.sh`.

This script requires a `10x-sample-sheet.tsv` file, generated automatically in the MiXCR section of the main analysis script, to restrict analysis to cell-containing GEMs identified by Cell Ranger.

## Filtering

The following cells or features were excluded:

- Doublets
- Ribosomal genes
- TCR genes
- Genes with very low expression
- Cells with >20% mitochondrial reads (`percent.mt`)
- Cells without recovered TCR information

## CD4+/CD8+ T-Cell Annotation

[MAGIC](https://magic.readthedocs.io/en/stable/) imputation was first applied to gene-expression data. CD4 and CD8A/CD8B expression thresholds were then estimated from the imputed values using nonparametric density estimation and clustering with [`pdfCluster`](https://github.com/cran/pdfCluster/tree/master).

Final CD4+/CD8+ labels were assigned using a consensus-based rule combining MAGIC-derived thresholding and [scGate](https://github.com/carmonalab/scGate) annotation.

## Compartment Annotation

Cells were assigned to compartments based on sample identity and CD4+/CD8+ annotation. When required to meet minimum cell numbers for downstream analyses, sorted populations were supplemented with cells from the opposite lineage and/or CD39 status.

The following compartments were considered:

- `CD4+Bulk`
- `CD4+CD39+`
- `CD4+CD39-`
- `CD8+Bulk`
- `CD8+CD39+`
- `CD8+CD39-`

## PDCD1 Thresholding

PDCD1 expression thresholds were estimated independently within each compartment using MAGIC-imputed expression values and density-based clustering with `pdfCluster`.

These thresholds were used to define PDCD1-high and PDCD1-low control clonotypes.

## CD8+ T-Cell State Analysis

CD8+ T cells were subsetted and reclustered to refine cell-state annotation. Clusters were manually annotated using marker-gene expression and aggregate cluster-level profiles.

## TCR Candidate Prioritisation and Validation

NeoAg-associated clonotypes were tracked across CD8+ compartments, including bulk, CD39+, and CD39− populations.

Candidate TCRs were prioritised based on:

- Presence in NeoAg stimulation pools
- Detection in the top 70% or top 100% cumulative clonotype-frequency sets
- Expansion across CD8+ compartments
- PDCD1 expression
- Availability of matched TCRβ repertoire information

At this stage, 26 TCRs were reconstructed, transduced into healthy-donor PBMCs, and tested for reactivity against the patient’s candidate (neo)antigens.

After experimental validation, TCR reactivity status was projected back onto the single-cell dataset for downstream visualization and phenotypic interpretation.

## Output

<p align="center">
  <img src="plot.png" width="1000">
</p>