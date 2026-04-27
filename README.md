# Single-Cell Profiling of (Neo)Antigen-Reactive TCRs Across Blood Compartments

> Computational analyses supporting the manuscript  
> ***Peripheral Blood as the Sole Source for Cancer Neoantigen and Reactive T Cell Discovery***  
> Garcia-Garijo A et al., 2026

## Overview

This repository contains the computational workflows used to identify dominant neoantigen-associated TCR clonotypes from antigen-specific co-culture experiments and to track their phenotypic distribution across peripheral blood T-cell compartments.

The project is organised into three sequential analysis phases.

---

## Workflow Structure

### Phase 1 — Dominant Clonotype Discovery

Run:

`ITAG_TIL15_S02_m1.Rmd`  
[Module documentation](https://github.com/eprieto012/cNART-neo/blob/main/ident/README.md)

This workflow:

- preprocesses 10x Genomics 5′ gene-expression and V(D)J data
- rescues missing TCR assignments using MiXCR
- filters low-quality cells and features
- annotates CD4+ and CD8+ T cells
- identifies dominant clonotypes

Outputs include cumulative clonotype-frequency candidate sets:

- `TCR_interrogated_top100.tsv`
- `TCR_interrogated_top70.tsv`

---

### Phase 2 — Cross-Compartment Tracking and Validation

Run:

`ITAG_TIL15_S02_m2.Rmd`  
[Module documentation](https://github.com/eprieto012/cNART-neo/blob/main/track/README.md)

This workflow uses Phase 1 outputs to:

- track clonotypes across CD8+ bulk / CD39+ / CD39− compartments
- prioritise candidates based on expansion and PDCD1 status
- map clonotypes onto CD8+ T-cell states
- project experimentally validated reactive TCRs onto the single-cell atlas

---

### Phase 3 — Compartment Frequency Analysis

Run:

`ITAG_TIL15_Bulk_m1.Rmd`  
[Module documentation](https://github.com/eprieto012/cNART-neo/blob/main/bulk/README.md)

This workflow compares candidate clonotype frequencies across tumor biopsy, peripheral blood, and TIL infusion product samples.

---

## Repository Modules

- **ident/** — Dominant clonotype discovery following antigen-specific co-culture with TMG-loaded autologous B cells.
- **track/** — Cross-compartment tracking and validation of dominant clonotypes by PD-1/CD39 phenotype.
- **bulk/** — Frequency analysis across tumor, blood, and infusion-product compartments.