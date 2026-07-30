# Gastric epithelial composition-aware transcriptomics

Version 2.2.0

[![Concept DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21677617.svg)](https://doi.org/10.5281/zenodo.21677617)

This repository contains the reproducibility code and frozen, publication-facing
derived outputs supporting the manuscript:

> Independent single-cell and spatial validation reveals mixed compositional
> and epithelial contributions to a gastric cancer transcriptomic module

## Study overview

The analysis used three paired bulk microarray cohorts (GSE29272, GSE79973 and
GSE19826), two independent single-cell cohorts (GSE206785 and GSE270680), and
one spatial transcriptomic cohort (GSE270678).

A 73-gene module identified in GSE29272 was tested in two independent bulk
cohorts, localized and evaluated at patient level in two single-cell cohorts,
examined with composition-aware bulk models, and assessed for spatial
colocalization across 33 sections from 19 donors.

The bounded conclusion is that the bulk module is sensitive to gastric
epithelial composition, while its within-epithelial component differs between
single-cell cohorts. The analysis does not establish a tumor-cell-intrinsic
mechanism, a causal mediator, a clinical biomarker, or protein-level validation.

## Repository contents

- `scripts/`: 23 R scripts covering data retrieval, preprocessing, statistical
  analysis, validation, auditing and release assembly.
- `figures/`: frozen main figures in PDF and PNG formats.
- `source_data/`: CSV source data for the main figure panels.
- `supplementary_tables/`: frozen publication-facing statistical outputs and
  download inventories.
- `metadata/audit/`: version 2.2.0 manuscript and claim-consistency audits.
- `metadata/validation_summaries/`: validation reports and R session records.
- `metadata/audit_v1.0.0/`: retained historical audit records from version
  1.0.0.
- `metadata/data_sources.tsv`: public source accessions, downloaded filenames
  and recorded checksums.
- `metadata/file_manifest_sha256.csv`: SHA-256 manifest for the public release.
- `metadata/ZENODO_METADATA_TEMPLATE.md`: fields to verify in the Zenodo record.
- `LICENSE_CODE`: MIT licence for the R code.
- `LICENSE_DATA.md`: CC BY 4.0 licence for original derived outputs, figures and
  documentation.

The manuscript, cover letter, submission-only DOCX/PDF files, original GEO
files, third-party publication supplements and local analysis caches are not
included in this repository.

## Reproduction outline

1. Install R and the packages reported in the included session-information
   files.
2. Obtain the public inputs recorded in `metadata/data_sources.tsv`. The
   download scripts cover the GEO inputs used by the workflow.
3. Run the bulk preprocessing and discovery scripts `01` through `06`.
4. Run the first paper-level validation scripts `07` through `10`.
5. Run `13_download_GSE270680_GSE270678.R`,
   `14_GSE270680_singlecell_replication.R`, and
   `15_GSE270678_spatial_validation.R` for the second single-cell and spatial
   cohorts.
6. Run `20_singlecell_cohort_heterogeneity_audit.R` for the prespecified
   cross-cohort descriptive comparison. No pooled two-cohort meta-analysis is
   reported.
7. Compare regenerated outputs with the frozen audit and source-data files.

Scripts locate the project directory relative to their own location rather than
requiring the original Windows drive path. Some manuscript-assembly utilities
(`11`, `12`, and `16` through `19`) require submission-source files that are
intentionally excluded from the public repository; they are not required to
rerun the biological analyses.

## Data and code availability

The stable Zenodo concept DOI for all versions is
[10.5281/zenodo.21677617](https://doi.org/10.5281/zenodo.21677617). The archived
version 1.0.0 record is
[10.5281/zenodo.21677618](https://doi.org/10.5281/zenodo.21677618).

The version-specific Zenodo DOI for version 2.2.0 will be added only after the
new Zenodo version is actually published.

Original public datasets remain available from NCBI GEO under accession
numbers GSE29272, GSE79973, GSE19826, GSE206785, GSE270680 and GSE270678.
They are not redistributed here.

## Reuse and attribution

The R scripts are released under the MIT licence. The authors' original derived
tables, statistical outputs, figures and documentation are released under
CC BY 4.0. Third-party GEO and publication-source files remain subject to their
original terms and citation requirements; no additional licence is asserted
over those files.

## Contact

Corresponding author: Changlong Yang  
The Third Affiliated Hospital of Kunming Medical University/Yunnan Cancer
Hospital  
Email: 965084451@qq.com  
ORCID: https://orcid.org/0009-0001-6553-4711
