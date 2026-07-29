# Gastric epithelial composition-aware transcriptomics

Version 1.0.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21677618.svg)](https://doi.org/10.5281/zenodo.21677618)

Archived release: [Zenodo version 1.0.0](https://doi.org/10.5281/zenodo.21677618).

This repository contains the R code, processed expression matrices, statistical
outputs and figure source data supporting the manuscript:

> A reproducible gastric epithelial transcriptomic module is largely shaped by
> tissue composition in gastric cancer

## Study overview

The analysis used three paired bulk microarray cohorts (GSE29272, GSE79973 and
GSE19826) and one paired single-cell cohort (GSE206785). A 73-gene module
identified in GSE29272 was tested in two independent bulk cohorts and then
evaluated using patient-level single-cell pseudobulk and composition-aware
models.

The repository supports the bounded conclusion that the reproducible bulk
module is sensitive to differentiated gastric epithelial representation. It
does not establish a tumor-cell-intrinsic mechanism or a clinical biomarker.

## Repository contents

- `scripts/`: reproducible R scripts.
- `clean_data/`: processed gene-by-sample matrices and sample metadata.
- `results/`: statistical outputs, reports, figures and figure source data.
- `metadata/data_sources.tsv`: public source accessions, downloaded filenames
  and SHA-256 checksums.
- `metadata/file_manifest_sha256.csv`: checksum manifest for this release.
- `metadata/DATA_DICTIONARY.md`: file and variable-level guidance.
- `metadata/ZENODO_METADATA_TEMPLATE.md`: fields used for the Zenodo record.
- `LICENSE_CODE`: MIT licence for the R code.
- `LICENSE_DATA.md`: CC BY 4.0 licence for original processed data, derived
  outputs, figures and documentation.

The original GEO files are not redistributed. They remain available from NCBI
GEO under their original accessions. The download inventory and checksums are
provided so that the exact inputs can be verified.

## Reproduction outline

1. Install R and the packages reported in the included `sessionInfo` files.
2. Obtain the public GEO inputs listed in `metadata/data_sources.tsv`, or use
   `scripts/01_download_GEO_data.R` where applicable.
3. Run the preprocessing scripts in numerical order.
4. Run the paper-level validation scripts:
   - `07_paper_bulk_validation.R`
   - `08_paper_singlecell_composition.R`
   - `09_paper_bulk_deconvolution.R`
   - `10_paper_figures.R`
5. Compare the regenerated statistical summaries with the frozen
   `metadata/audit/manuscript_numeric_audit.csv`. The pre-publication
   consistency audit passed 25 of 25 checks; its report is retained in
   `metadata/audit/`.

Scripts locate the project directory relative to their own location and should
not require the original Windows drive path.

## Main data sources

- GSE29272: discovery paired bulk microarray cohort.
- GSE79973: external paired bulk validation cohort.
- GSE19826: external paired bulk validation cohort.
- GSE206785: paired single-cell localization and composition analysis.

## Data and code availability

The R scripts, processed expression matrices, statistical outputs, figure source
data and reproducibility materials are publicly available in this GitHub
repository and archived as version 1.0.0 in Zenodo at
[https://doi.org/10.5281/zenodo.21677618](https://doi.org/10.5281/zenodo.21677618).

Original public datasets are available from the NCBI Gene Expression Omnibus
(GEO) under accession numbers GSE29272, GSE79973, GSE19826 and GSE206785.

## Reuse and attribution

The R scripts are released under the MIT licence. The authors' original
processed data tables, derived statistical outputs, figures and documentation
are released under CC BY 4.0. The public GEO source datasets remain subject to
the terms and citation requirements of their original records and
publications; no additional licence is asserted over those third-party data.

## Contact

Corresponding author: Changlong Yang  
The Third Affiliated Hospital of Kunming Medical University/Yunnan Cancer
Hospital  
Email: 965084451@qq.com  
ORCID: https://orcid.org/0009-0001-6553-4711
