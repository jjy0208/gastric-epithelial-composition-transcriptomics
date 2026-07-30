# Data dictionary

## Figures

`figures/` contains the six frozen main figures in PDF and PNG formats. These
files are publication-facing outputs rather than independent raw measurements.

## Figure source data

`source_data/` contains CSV files mapped to main-figure panels. Filenames begin
with the corresponding figure number. Statistical source tables retain the
units used in the manuscript: paired patients for patient-level comparisons,
sections or donors for spatial summaries, and cohorts for bulk validation.

## Supplementary tables

`supplementary_tables/` contains:

- discovery differential-expression and enrichment tables;
- external bulk validation summaries;
- composition-adjustment summaries;
- patient-level single-cell summaries for GSE206785 and GSE270680;
- descriptive cross-cohort standardized effects, without a pooled estimate;
- spatial section manifests, gene coverage and concordance tests for GSE270678;
- download inventories for public inputs.

`singlecell_eligible_patient_clinical_metadata.csv` uses public study patient
codes, not names or newly collected identifiers. Blank fields represent
information not reported in the public source.

## Audit and validation records

`metadata/audit/` contains the frozen version 2.2.0 claim-evidence map,
manuscript numeric audit, reference audit and terminology records.

`metadata/validation_summaries/` contains cohort-specific validation reports and
R session-information files.

Historical version 1.0.0 audit files are retained separately in
`metadata/audit_v1.0.0/`.

## Missing values and exclusions

Blank fields or `NA` represent unavailable quantities. Exclusion rules and
minimum-cell thresholds are documented in the R scripts and validation reports.
Cells and spatial spots were not treated as independent biological replicates.

## Provenance and licensing

All outputs derive from the public sources recorded in `data_sources.tsv`.
Original GEO files, author-provided objects and publication supplements are not
redistributed. The public release integrity is recorded in
`file_manifest_sha256.csv`.
