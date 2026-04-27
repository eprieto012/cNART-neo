# Compartment Frequency Analysis

## Analysis Workflow

Run `ITAG_TIL15_Bulk_m1.Rmd` to compare TCRβ clonotype frequencies across bulk repertoire samples from tumor, peripheral blood, and TIL infusion product.

## TCRβ Processing

Only productive in-frame clonotypes were retained. Clonotypes were defined by CDR3β amino-acid sequence together with TRBV and TRBJ genes, and frequencies were calculated within each sample.

## Reactive TCR Annotation

Validated clonotypes were annotated using `reactive-clono.txt`, including TCR identifier and experimental reactivity status.

## Output

<p align="center">
  <img src="plot.png" width="1000">
</p>

Frequency plots were generated to compare reactive and non-reactive clonotypes across compartments.