# 构建不含原始大文件的投稿与可复现性发布包，并生成 SHA-256 清单。
suppressPackageStartupMessages(library(digest))

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) != 1L) stop("请使用 Rscript 运行本脚本。")
  normalizePath(sub("^--file=", "", f), winslash = "/", mustWork = TRUE)
}

project_dir <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  winslash = "/",
  mustWork = TRUE
)
package_name <- "gastric-epithelial-composition-transcriptomics-v2.2.0-submission"
package_dir <- file.path(project_dir, "release", package_name)
if (dir.exists(package_dir)) {
  stop("目标发布目录已存在，为避免覆盖请先人工核对：", package_dir)
}

dirs <- file.path(
  package_dir,
  c(
    "submission", "figures", "source_data", "supplementary_tables",
    "scripts", "audit", "validation_summaries"
  )
)
dir.create(package_dir, recursive = TRUE, showWarnings = FALSE)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

copy_checked <- function(files, destination) {
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("打包时缺少文件：", paste(missing, collapse = ", "))
  }
  ok <- file.copy(files, destination, overwrite = FALSE, copy.date = TRUE)
  if (!all(ok)) stop("文件复制失败：", paste(files[!ok], collapse = ", "))
}

manuscript_dir <- file.path(project_dir, "report", "manuscript")
submission_files <- file.path(
  manuscript_dir,
  c(
    "gastric_epithelial_composition_manuscript_v2.docx",
    "gastric_epithelial_composition_manuscript_v2.pdf",
    "cover_letter_scientific_reports.docx",
    "cover_letter_scientific_reports.pdf",
    "supplementary_information_v2.docx",
    "supplementary_information_v2.pdf",
    "submission_readiness_checklist.md",
    "submission_package_readme.md",
    "data_code_availability.md"
  )
)
copy_checked(submission_files, file.path(package_dir, "submission"))

figure_dir <- file.path(project_dir, "results", "ManuscriptFigures")
figure_names <- c(
  "Figure1_discovery", "Figure2_bulk_replication",
  "Figure3_single_cell_localization", "Figure4_composition_adjustment",
  "Figure5_independent_single_cell_validation", "Figure6_spatial_validation"
)
figure_files <- unlist(lapply(
  figure_names,
  function(x) file.path(figure_dir, paste0(x, c(".pdf", ".png")))
))
copy_checked(figure_files, file.path(package_dir, "figures"))

source_files <- list.files(
  file.path(figure_dir, "SourceData"),
  pattern = "\\.csv$", full.names = TRUE
)
copy_checked(source_files, file.path(package_dir, "source_data"))

supplementary_tables <- c(
  file.path(project_dir, "results", "DEG", "GSE29272_DEG_all.csv"),
  file.path(project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"),
  file.path(project_dir, "results", "GSEA", "GSE29272_GSEA_result.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCell",
            "epithelial_lineage_paired_tests.csv"),
  file.path(project_dir, "raw_data", "download_log.txt"),
  file.path(project_dir, "raw_data", "GSE270680_GSE270678_download_log.tsv"),
  file.path(project_dir, "results", "PaperValidation",
            "external_bulk_validation_summary.csv"),
  file.path(project_dir, "results", "PaperValidation", "Deconvolution",
            "bulk_deconvolution_adjustment_summary.csv"),
  file.path(project_dir, "results", "PaperValidation", "Deconvolution",
            "bulk_lineage_marker_panel.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCell",
            "singlecell_analysis_summary.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "GSE270680_singlecell_replication_summary.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "singlecell_cross_cohort_standardized_effects.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "GSE270680_epithelial_lineage_paired_tests.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "singlecell_cohort_characteristics.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "singlecell_cohort_comparison_tests.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "singlecell_eligible_patient_clinical_metadata.csv"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "singlecell_clinical_stratified_effects.csv"),
  file.path(project_dir, "raw_data", "publication_supplements",
            "download_log.tsv"),
  file.path(project_dir, "results", "PaperValidation", "Spatial",
            "GSE270678_section_manifest.csv"),
  file.path(project_dir, "results", "PaperValidation", "Spatial",
            "GSE270678_spatial_concordance_tests.csv"),
  file.path(project_dir, "results", "PaperValidation", "Spatial",
            "GSE270678_signature_gene_coverage.csv")
)
copy_checked(supplementary_tables, file.path(package_dir, "supplementary_tables"))

script_files <- list.files(
  file.path(project_dir, "scripts"),
  pattern = "\\.R$", full.names = TRUE
)
copy_checked(script_files, file.path(package_dir, "scripts"))

audit_files <- file.path(
  manuscript_dir,
  c(
    "manuscript_numeric_audit_v2.csv",
    "manuscript_audit_report_v2.md",
    "reference_recency_and_doi_audit.csv",
    "reference_recency_and_doi_audit.md",
    "claim_evidence_map.csv",
    "terminology_ledger.csv",
    "scientific_editing_and_humanization_report.md",
    "reviewer_issues_resolution_report.md"
  )
)
copy_checked(audit_files, file.path(package_dir, "audit"))

validation_files <- c(
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "GSE270680_singlecell_replication_report.md"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "singlecell_cohort_heterogeneity_report.md"),
  file.path(project_dir, "results", "PaperValidation", "Spatial",
            "GSE270678_spatial_validation_report.md"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "14_sessionInfo.txt"),
  file.path(project_dir, "results", "PaperValidation", "SingleCellReplication",
            "20_sessionInfo.txt"),
  file.path(project_dir, "results", "PaperValidation", "Spatial",
            "15_sessionInfo.txt")
)
copy_checked(validation_files, file.path(package_dir, "validation_summaries"))

all_files <- list.files(package_dir, recursive = TRUE, full.names = TRUE)
manifest <- data.frame(
  relative_path = substring(
    normalizePath(all_files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(package_dir, winslash = "/", mustWork = TRUE)) + 2L
  ),
  bytes = file.info(all_files)$size,
  sha256 = vapply(
    all_files,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$relative_path), ]
write.csv(
  manifest,
  file.path(package_dir, "SHA256_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("投稿发布目录已生成：", package_dir, "\n")
cat("清单文件数：", nrow(manifest), "\n")
