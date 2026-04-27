# Compartment Frequency Analysis

## Analysis Workflow

Run `ITAG_TIL15_Bulk_m1.Rmd` to compare TCRβ clonotype frequencies across bulk repertoire samples.

## TCRβ Processing

Only productive in-frame clonotypes were retained. TCRβ clonotypes were defined using CDR3β amino-acid sequence together with TRBV and TRBJ genes. Frequencies were calculated within each sample.

## Reactive TCR Annotation

Validated clonotypes were annotated using `Reactive_IDs_TIL15.txt`, including TCR identifier and reactivity status.

## Output

Frequency plots were generated to compare reactive and non-reactive TCRβ clonotypes across tumor, blood, and infusion-product compartments.