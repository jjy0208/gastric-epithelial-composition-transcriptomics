# Data dictionary

## Processed bulk expression matrices

Files:

- `clean_data/GSE29272_clean_expression_matrix.csv`
- `clean_data/GSE79973_clean_expression_matrix.csv`
- `clean_data/GSE19826_clean_expression_matrix.csv`

Rows are HGNC gene symbols and columns are GEO samples. Values are processed,
gene-level microarray expression values. Multiple probes assigned to the same
gene were averaged after annotation. The GSE19826 positive signal matrix was
transformed as `log2(signal + 1)` before probe aggregation.

## Sample metadata

Files ending in `_sample_info.csv` link expression-matrix columns to analysis
groups. Core fields include:

- `sample`: GEO sample identifier.
- `group`: analysis group, standardized as tumor or normal.
- Additional columns retain available pairing, tissue, platform or source
  metadata.

Only complete tumor-normal patient pairs entered the paired primary analyses.

## Paper validation outputs

`results/PaperValidation/` contains:

- external bulk gene-direction and module-score validation;
- patient-level single-cell pseudobulk scores;
- whole-tissue composition sensitivity analyses;
- independent gastric-lineage decomposition scores;
- statistical summaries, reports and software environment records.

Single-cell patient and sample codes are retained from the already public
GSE206785 metadata. No newly collected direct identifiers are included.

## Figure source data

`results/ManuscriptFigures/SourceData/` maps one CSV file to each main-figure
panel. Filenames begin with the corresponding figure and panel number.

## Missing values and exclusions

Blank fields or `NA` represent unavailable quantities. Exclusion rules and
minimum-cell thresholds are documented in the R scripts and analysis reports.
Cells were not treated as independent biological replicates; paired patients
were the statistical unit where patient pairing was available.

## Provenance

All processed outputs derive from the public GEO accessions recorded in
`data_sources.tsv`. Raw third-party GEO files are not included in this release.
Release file integrity can be checked with
`metadata/file_manifest_sha256.csv`.

