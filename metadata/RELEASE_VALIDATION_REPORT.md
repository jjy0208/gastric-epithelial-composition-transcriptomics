# Public release validation report

Validation date: 2026-07-30

Release version: v2.2.0

## Intended public scope

- 23 R scripts
- 12 frozen main-figure files
- 23 figure source-data tables
- 21 supplementary result and provenance tables
- 8 current audit files
- 6 validation summaries and R session records
- repository documentation, citation metadata and licence files

## Explicit exclusions

- original GEO archives and matrices
- third-party publication supplements
- author-provided large single-cell objects
- manuscript and supplementary-information DOCX/PDF files
- cover letter and submission-only documents
- local caches, temporary files and manuscript render directories

## Validation status

- Source package SHA-256 confirmed before staging:
  `db2b751e72605228513261333d11251b0601493aba0ae7b96e9b3e27d74d2c07`.
- Source-to-staging comparison covered 93 files across all six public content
  groups: 89 were byte-identical and four had documented public-release safety
  edits.
- The four intentional edits were:
  - conversion of local absolute paths to project-relative paths in three
    provenance tables;
  - omission of non-corresponding-author email addresses from the public copy
    of `12_build_manuscript_docx.R`.
- Missing source files: 0.
- Unexpected source-to-staging hash mismatches: 0.
- R syntax parsing: 23/23 scripts passed.
- `CITATION.cff`: valid YAML, version 2.2.0, five creators.
- CSV structural loading: all CSV files passed with at least one data row.
- Forbidden raw/submission file types or directories: 0.
- Files at or above 25 MB: 0.
- Local absolute project-path hits: 0.
- credential/private-key pattern hits: 0.
- Public email exposure is limited to the already-public corresponding-author
  address in README, citation metadata and the sanitized manuscript builder.
- Public-release SHA-256 manifest generated and verified with zero mismatches.
- GitHub authentication and repository administrator permission confirmed.
- Reviewed content commit:
  `7e257ec4621e35d053dfed9fe1115c22ab2eb581`.
- PR #1 merge commit:
  `fda9b0a8d053399f08a09855f54e55c8a97005bd`.
- Git tag, GitHub Release and Zenodo publication: pending.

This report must be updated only from actual validation results. A pending item
must not be represented as passed.
