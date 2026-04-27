# Single-Cell Profiling of (Neo)Antigen-Reactive TCRs Across Blood Compartments

> Computational analyses supporting the manuscript ***Peripheral Blood as the Sole Source for Cancer Neoantigen and Reactive T Cell Discovery*** by Garcia-Garijo A et al., 2026.

## Overview

This repository contains computational workflows to identify dominant TCR clonotypes from antigen-specific co-culture assays of TILs and PBMCs with autologous B cells presenting candidate (neo)antigens in tandem minigene (TMG) format, track their phenotypic distribution across peripheral blood CD8+ T-cell subsets defined by PD-1 and CD39 expression, and compare their frequencies across tumor biopsy, peripheral blood, and TIL infusion product samples.

The project is organised into three sequential analysis phases.

## Workflow

### Phase 1: Antigen-Associated Clonotype Discovery
**Module:** [`disc/`](https://github.com/eprieto012/cNART-neo/tree/main/disc)

Script m1 `ITAG_TIL15_S02_m1.Rmd`:

- Preprocesses 10x Genomics 5′ gene-expression and V(D)J data.
- Rescues missing TCR assignments using MiXCR.
- Filters low-quality cells and features.
- Annotates CD4+ and CD8+ T cells.
- Prioritises expanded antigen-associated clonotypes.

### Phase 2: Cross-Compartment Clonotype Tracking
**Module:** [`track/`](https://github.com/eprieto012/cNART-neo/tree/main/track)

Script m1 `ITAG_TIL15_S01_m1.Rmd` reproduces the shared preprocessing workflow described in Phase 1.

Script m2 `ITAG_TIL15_S01_m2.Rmd` uses the Phase 1 outputs to:

- Track clonotypes across CD8+ bulk, CD39+, and CD39− compartments.
- Prioritise candidates for validation based on expansion, PDCD1 expression, and compartment of origin.
- Map clonotypes and validated reactivity status onto CD8+ T-cell states.

### Phase 3: Cross-Sample Clonotype Frequency Analysis
**Module:** [`bulk/`](https://github.com/eprieto012/cNART-neo/tree/main/bulk)

Script `ITAG_TIL15_S03_m1.Rmd`:

- Compares candidate clonotype frequencies across tumor biopsy, peripheral blood, and TIL infusion product samples.

### Setup

**Module:** [`setup/`](https://github.com/eprieto012/cNART-neo/tree/main/setup)

These scripts provide shared R/Python configuration, plotting themes, random seeds, and helper functions for the main workflows.

- `setup-m1.R` — common setup for the core single-cell workflows, including preprocessing, quality control, clonotype recovery, filtering, annotation, and clonotype-frequency analysis.
- `setup-m2.R` — additional setup for cross-compartment tracking, candidate prioritisation, and validation analyses.