# Resolution of reviewer-facing consistency and heterogeneity issues

## 1. M04 FDR discrepancy

The source file `results/WGCNA/GSE29272_WGCNA_module_trait_correlations.csv` was treated as authoritative. For M04 it records:

- correlation: -0.756230991028608;
- P value: 5.39409408329265e-50;
- Benjamini-Hochberg FDR: 3.23645644997559e-49.

The `10^-29` value in the Figure 1 legend was a transcription error. The Results, Figure 1 plot subtitle and source table already carried the correct order of magnitude. The legend was corrected to `P=5.39×10^-50, FDR=3.24×10^-49`. The manuscript audit now checks the source value and both textual locations and explicitly rejects the former erroneous value.

## 2. Two-cohort statistical comparison

The fixed-effect and random-effects summaries, tau-squared and I-squared were removed from the manuscript, Figure 5 and reproducible output. The revised analysis reports only the two cohort-specific Hedges' gz estimates and confidence intervals. They are displayed side by side as a descriptive comparison. No common effect is estimated.

This change avoids implying that two cohorts provide a stable estimate of between-study heterogeneity or a meaningful pooled effect.

## 3. Investigation of discordant single-cell estimates

Published patient-characteristic tables were retrieved from the official supplementary files for both source articles. The analysis was restricted to the 10 and 11 patients who passed the prespecified epithelial-cell threshold.

Clinical metadata coverage was 9/10 for GSE206785 and 11/11 for GSE270680. Age, sex, stage, broad histologic group and sequencing platform were compared descriptively and with explicitly exploratory tests. Gastric anatomic position was reported for GSE206785 but not GSE270680 and therefore could not be compared.

A concrete experimental-design difference was identified. GSE206785 used EPCAM-positive cell depletion before sequencing to enrich the tumor microenvironment; its epithelial cells are the incompletely depleted remainder. GSE270680 reported viable-cell sorting with CD235a exclusion and did not report EPCAM depletion. This difference can alter epithelial lineage representation and eligibility under the 20-cell rule. It is now discussed as a plausible selection mechanism, not as a proven cause.

The new reproducible outputs are:

- `singlecell_cohort_characteristics.csv`;
- `singlecell_cohort_comparison_tests.csv`;
- `singlecell_eligible_patient_clinical_metadata.csv`;
- `singlecell_clinical_stratified_effects.csv`;
- `singlecell_cohort_heterogeneity_report.md`.

Supplementary Tables S6 and S7 embed the cohort characteristics and exploratory comparison results directly in the supplementary document.
