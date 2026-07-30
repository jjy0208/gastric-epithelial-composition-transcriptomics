# 稿件终审：自动核对关键数值、文件、措辞边界和投稿材料完整性。
suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

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
report_dir <- file.path(project_dir, "report", "manuscript")
fig_dir <- file.path(project_dir, "results", "ManuscriptFigures")

checks <- list()
add_check <- function(id, description, observed, expected, pass, evidence) {
  checks[[length(checks) + 1L]] <<- data.table(
    id = id,
    description = description,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = isTRUE(pass),
    evidence = evidence
  )
}
metric <- function(dt, name) as.numeric(dt[metric == name, value])
near <- function(x, y, tol = 1e-8) is.finite(x) && abs(x - y) <= tol

deg <- fread(file.path(project_dir, "results", "DEG", "GSE29272_DEG_summary.csv"))
wgcna <- fread(file.path(project_dir, "results", "WGCNA", "GSE29272_WGCNA_summary.csv"))
candidate <- fread(file.path(project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"))
external <- fread(file.path(
  project_dir, "results", "PaperValidation",
  "external_bulk_validation_summary.csv"
))
sc1 <- fread(file.path(
  project_dir, "results", "PaperValidation", "SingleCell",
  "singlecell_analysis_summary.csv"
))
sc2 <- fread(file.path(
  project_dir, "results", "PaperValidation", "SingleCellReplication",
  "GSE270680_singlecell_replication_summary.csv"
))
cross_cohort_effects <- fread(file.path(
  project_dir, "results", "PaperValidation", "SingleCellReplication",
  "singlecell_cross_cohort_standardized_effects.csv"
))
spatial <- fread(file.path(
  project_dir, "results", "PaperValidation", "Spatial",
  "GSE270678_spatial_concordance_tests.csv"
))
spatial_manifest <- fread(file.path(
  project_dir, "results", "PaperValidation", "Spatial",
  "GSE270678_section_manifest.csv"
))
deconv <- fread(file.path(
  project_dir, "results", "PaperValidation", "Deconvolution",
  "bulk_deconvolution_adjustment_summary.csv"
))

add_check("A01", "Discovery significant DEG count", deg$significant_genes, 122,
          deg$significant_genes == 122, "results/DEG/GSE29272_DEG_summary.csv")
add_check("A02", "Discovery up/down counts",
          paste(deg$upregulated, deg$downregulated, sep = "/"), "47/75",
          deg$upregulated == 47 && deg$downregulated == 75,
          "results/DEG/GSE29272_DEG_summary.csv")
add_check("A03", "Candidate intersection count", nrow(candidate), 73,
          nrow(candidate) == 73,
          "results/HubGene/GSE29272_candidate_genes.csv")
add_check("A04", "All candidate genes lower in tumor",
          sum(candidate$DEG_direction == "Down"), 73,
          all(candidate$DEG_direction == "Down"),
          "results/HubGene/GSE29272_candidate_genes.csv")
add_check("A05", "WGCNA key-module gene count",
          metric(wgcna, "key_module_genes"), 192,
          metric(wgcna, "key_module_genes") == 192,
          "results/WGCNA/GSE29272_WGCNA_summary.csv")
add_check("A06", "WGCNA scale-free R-squared criterion not reached",
          wgcna[metric == "R2_criterion_reached", value], "FALSE",
          wgcna[metric == "R2_criterion_reached", value] == "FALSE",
          "results/WGCNA/GSE29272_WGCNA_summary.csv")

for (ds in c("GSE79973", "GSE19826")) {
  row <- external[dataset == ds]
  add_check(
    paste0("B_", ds, "_direction"),
    paste(ds, "assayed/concordant candidate genes"),
    paste(row$available_module_genes, row$down_direction_genes, sep = "/"),
    "72/72",
    row$available_module_genes == 72 && row$down_direction_genes == 72,
    "results/PaperValidation/external_bulk_validation_summary.csv"
  )
  add_check(
    paste0("B_", ds, "_paired"),
    paste(ds, "paired module test"),
    row$paired_t_pvalue, "<0.05",
    row$paired_t_pvalue < 0.05,
    "results/PaperValidation/external_bulk_validation_summary.csv"
  )
}

add_check("C01", "GSE206785 total cells", metric(sc1, "cells"), 111140,
          metric(sc1, "cells") == 111140,
          "results/PaperValidation/SingleCell/singlecell_analysis_summary.csv")
add_check("C02", "GSE206785 eligible epithelial pairs",
          metric(sc1, "eligible_complete_epithelial_pairs"), 10,
          metric(sc1, "eligible_complete_epithelial_pairs") == 10,
          "results/PaperValidation/SingleCell/singlecell_analysis_summary.csv")
add_check("C03", "GSE206785 paired test did not detect a difference",
          metric(sc1, "epithelial_score_paired_t_pvalue"), ">=0.05",
          metric(sc1, "epithelial_score_paired_t_pvalue") >= 0.05,
          "results/PaperValidation/SingleCell/singlecell_analysis_summary.csv")
add_check("C04", "GSE270680 total cells", metric(sc2, "cells_total"), 470609,
          metric(sc2, "cells_total") == 470609,
          "results/PaperValidation/SingleCellReplication/GSE270680_singlecell_replication_summary.csv")
add_check("C05", "GSE270680 patients/libraries",
          paste(metric(sc2, "patients_total"), metric(sc2, "libraries_total"), sep = "/"),
          "27/77",
          metric(sc2, "patients_total") == 27 && metric(sc2, "libraries_total") == 77,
          "results/PaperValidation/SingleCellReplication/GSE270680_singlecell_replication_summary.csv")
add_check("C06", "GSE270680 eligible epithelial pairs",
          metric(sc2, "eligible_complete_epithelial_pairs"), 11,
          metric(sc2, "eligible_complete_epithelial_pairs") == 11,
          "results/PaperValidation/SingleCellReplication/GSE270680_singlecell_replication_summary.csv")
add_check("C07", "GSE270680 paired mean difference",
          metric(sc2, "mean_tumor_minus_normal"), "-0.6997324",
          near(metric(sc2, "mean_tumor_minus_normal"), -0.699732405464963),
          "results/PaperValidation/SingleCellReplication/GSE270680_singlecell_replication_summary.csv")
add_check("C08", "GSE270680 paired P value",
          metric(sc2, "paired_t_pvalue"), "<0.001",
          metric(sc2, "paired_t_pvalue") < 0.001,
          "results/PaperValidation/SingleCellReplication/GSE270680_singlecell_replication_summary.csv")
add_check("C09", "GSE270680 lineage-adjusted P value",
          metric(sc2, "lineage_adjusted_cluster_robust_pvalue"), "<0.05",
          metric(sc2, "lineage_adjusted_cluster_robust_pvalue") < 0.05,
          "results/PaperValidation/SingleCellReplication/GSE270680_singlecell_replication_summary.csv")

effect_206785 <- cross_cohort_effects[cohort == "GSE206785"]
effect_270680 <- cross_cohort_effects[cohort == "GSE270680"]
add_check("C10", "GSE206785 cohort-specific standardized effect",
          effect_206785$hedges_gz, "-0.03483035",
          near(effect_206785$hedges_gz, -0.0348303519967794),
          "results/PaperValidation/SingleCellReplication/singlecell_cross_cohort_standardized_effects.csv")
add_check("C11", "GSE270680 cohort-specific standardized effect",
          effect_270680$hedges_gz, "-1.4112448",
          near(effect_270680$hedges_gz, -1.41124480123282),
          "results/PaperValidation/SingleCellReplication/singlecell_cross_cohort_standardized_effects.csv")

primary_sections <- sum(spatial_manifest$primary_section %in% TRUE)
add_check("D01", "Primary spatial sections", primary_sections, 33,
          primary_sections == 33,
          "results/PaperValidation/Spatial/GSE270678_section_manifest.csv")
add_check("D02", "Spatial donors",
          spatial[comparison == "Epithelial", n_donors], 19,
          spatial[comparison == "Epithelial", n_donors] == 19,
          "results/PaperValidation/Spatial/GSE270678_spatial_concordance_tests.csv")
add_check("D03", "Epithelial spatial localization",
          spatial[comparison == "Epithelial", median_rho], "positive, FDR<0.05",
          spatial[comparison == "Epithelial", median_rho] > 0 &&
            spatial[comparison == "Epithelial", wilcoxon_fdr] < 0.05,
          "results/PaperValidation/Spatial/GSE270678_spatial_concordance_tests.csv")
add_check("D04", "Stromal association not detected",
          spatial[comparison == "Stromal", wilcoxon_fdr], ">=0.05",
          spatial[comparison == "Stromal", wilcoxon_fdr] >= 0.05,
          "results/PaperValidation/Spatial/GSE270678_spatial_concordance_tests.csv")
add_check("D05", "Immune association not detected",
          spatial[comparison == "Immune", wilcoxon_fdr], ">=0.05",
          spatial[comparison == "Immune", wilcoxon_fdr] >= 0.05,
          "results/PaperValidation/Spatial/GSE270678_spatial_concordance_tests.csv")

for (ds in c("GSE79973", "GSE19826")) {
  row <- deconv[dataset == ds]
  add_check(
    paste0("E_", ds, "_attenuation"),
    paste(ds, "composition-adjusted attenuation"),
    row$coefficient_attenuation, ">0.5",
    row$coefficient_attenuation > 0.5,
    "results/PaperValidation/Deconvolution/bulk_deconvolution_adjustment_summary.csv"
  )
}

figure_names <- c(
  "Figure1_discovery", "Figure2_bulk_replication",
  "Figure3_single_cell_localization", "Figure4_composition_adjustment",
  "Figure5_independent_single_cell_validation", "Figure6_spatial_validation"
)
figure_files <- as.vector(outer(
  file.path(fig_dir, figure_names),
  c(".pdf", ".png"),
  paste0
))
add_check("F01", "All 12 main-figure files exist",
          sum(file.exists(figure_files)), 12,
          all(file.exists(figure_files)), "results/ManuscriptFigures")
source_files <- list.files(
  file.path(fig_dir, "SourceData"),
  pattern = "\\.csv$", full.names = TRUE
)
add_check("F02", "Figure source-data files", length(source_files), ">=22",
          length(source_files) >= 22,
          "results/ManuscriptFigures/SourceData")
new_scripts <- file.path(
  project_dir, "scripts",
  c(
    "13_download_GSE270680_GSE270678.R",
    "14_GSE270680_singlecell_replication.R",
    "15_GSE270678_spatial_validation.R"
  )
)
add_check("F03", "New validation scripts exist", sum(file.exists(new_scripts)), 3,
          all(file.exists(new_scripts)), "scripts")

manuscript_path <- file.path(report_dir, "manuscript_en.md")
manuscript <- paste(
  readLines(manuscript_path, warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)
unfinished_patterns <- c(
  "Analysis freeze", "automated 25-item audit",
  "should confirm", "require confirmation and approval",
  "DOI to be assigned", "Draft for author review"
)
unfinished_hits <- unfinished_patterns[vapply(
  unfinished_patterns,
  function(x) grepl(x, manuscript, fixed = TRUE),
  logical(1)
)]
add_check("G01", "No unfinished-template language in manuscript",
          paste(unfinished_hits, collapse = "; "), "none",
          length(unfinished_hits) == 0,
          "report/manuscript/manuscript_en.md")
required_phrases <- c(
  "descriptive rather than a causal mediation analysis",
  "No perturbation or protein-level validation was performed",
  "no source publication was linked from the GSE79973 GEO record",
  "No fixed-effect or random-effects summary, between-study variance or I² was calculated or interpreted."
)
missing_phrases <- required_phrases[!vapply(
  required_phrases,
  function(x) grepl(x, manuscript, fixed = TRUE),
  logical(1)
)]
add_check("G02", "Key limitation statements retained",
          paste(missing_phrases, collapse = "; "), "none missing",
          length(missing_phrases) == 0,
          "report/manuscript/manuscript_en.md")
metadata_complete <-
  grepl("Junyi Jia", manuscript, fixed = TRUE) &&
  grepl("Changlong Yang", manuscript, fixed = TRUE) &&
  grepl(
    "The Third Affiliated Hospital of Kunming Medical University/Yunnan Cancer Hospital",
    manuscript, fixed = TRUE
  ) &&
  grepl("2024XKTDTS06", manuscript, fixed = TRUE) &&
  grepl("The authors declare no competing interests.", manuscript, fixed = TRUE)
add_check("G03", "Author, affiliation, funding and conflict metadata complete",
          metadata_complete, "TRUE", metadata_complete,
          "report/manuscript/manuscript_en.md")
manuscript_lines <- readLines(
  manuscript_path, warn = FALSE, encoding = "UTF-8"
)
count_words <- function(x) {
  x <- trimws(gsub("[[:space:]]+", " ", paste(x, collapse = " ")))
  if (!nzchar(x)) return(0L)
  length(strsplit(x, " ", fixed = TRUE)[[1]])
}
title_words <- count_words(sub("^# ", "", manuscript_lines[1]))
abstract_start <- which(manuscript_lines == "## Abstract")
introduction_start <- which(manuscript_lines == "## Introduction")
abstract_words <- count_words(
  manuscript_lines[(abstract_start + 1L):(introduction_start - 1L)]
)
add_check("G04", "Scientific Reports title length", title_words, "<=20",
          title_words <= 20,
          "report/manuscript/manuscript_en.md")
add_check("G05", "Scientific Reports abstract length", abstract_words, "<=200",
          abstract_words <= 200,
          "report/manuscript/manuscript_en.md")
body_without_affiliations <- paste(manuscript_lines[-seq_len(15)], collapse = "\n")
old_citation_hits <- gregexpr(
  "\\^[0-9]+([,–-][0-9]+)*\\^",
  body_without_affiliations,
  perl = TRUE
)[[1]]
add_check("G06", "No superscript-style citations remain in article body",
          if (old_citation_hits[1] == -1) 0 else length(old_citation_hits), 0,
          old_citation_hits[1] == -1,
          "report/manuscript/manuscript_en.md")
reference_audit <- fread(file.path(
  report_dir, "reference_recency_and_doi_audit.csv"
))
add_check("G07", "Reference list size after removing uncited legacy items",
          nrow(reference_audit), 22,
          nrow(reference_audit) == 22,
          "report/manuscript/reference_recency_and_doi_audit.csv")
add_check("G08", "References from 2021–2026",
          sum(reference_audit$recent_2021_2026), 15,
          sum(reference_audit$recent_2021_2026) == 15,
          "report/manuscript/reference_recency_and_doi_audit.csv")
add_check("G09", "Every retained older reference has an essential-source rationale",
          sum(
            !reference_audit$recent_2021_2026 &
              nzchar(reference_audit$old_reference_rationale)
          ),
          sum(!reference_audit$recent_2021_2026),
          all(
            reference_audit$recent_2021_2026 |
              nzchar(reference_audit$old_reference_rationale)
          ),
          "report/manuscript/reference_recency_and_doi_audit.csv")
add_check("G10", "All reference DOIs resolve through Crossref",
          sum(reference_audit$crossref_status == "resolved"), nrow(reference_audit),
          all(reference_audit$crossref_status == "resolved"),
          "report/manuscript/reference_recency_and_doi_audit.csv")

reference_start <- which(manuscript_lines == "## References")
article_prose <- paste(
  manuscript_lines[seq_len(reference_start - 1L)],
  collapse = "\n"
)
prose_dash_hits <- gregexpr("[—–]", article_prose, perl = TRUE)[[1]]
add_check("G11", "No em or en dash remains in article prose",
          if (prose_dash_hits[1] == -1) 0 else length(prose_dash_hits), 0,
          prose_dash_hits[1] == -1,
          "report/manuscript/manuscript_en.md")
template_phrases <- c(
  "Here, we", "In conclusion", "These findings establish",
  "The most defensible interpretation", "Several additional limitations",
  "We next", "Taken together", "It is important to note"
)
template_hits <- template_phrases[vapply(
  template_phrases,
  function(x) grepl(x, article_prose, fixed = TRUE),
  logical(1)
)]
add_check("G12", "No targeted formulaic transition remains in article prose",
          paste(template_hits, collapse = "; "), "none",
          length(template_hits) == 0,
          "report/manuscript/manuscript_en.md")

figure_legend <- paste(
  readLines(
    file.path(report_dir, "figure_legends.md"),
    warn = FALSE, encoding = "UTF-8"
  ),
  collapse = "\n"
)
m04_result <- fread(file.path(
  project_dir, "results", "WGCNA",
  "GSE29272_WGCNA_module_trait_correlations.csv"
))[module == "M04"]
m04_consistent <-
  nrow(m04_result) == 1L &&
  isTRUE(all.equal(m04_result$correlation, -0.756230991028608, tolerance = 1e-12)) &&
  isTRUE(all.equal(m04_result$p_value, 5.39409408329265e-50, tolerance = 1e-12)) &&
  isTRUE(all.equal(m04_result$FDR, 3.23645644997559e-49, tolerance = 1e-12)) &&
  grepl("P=5.39×10^−50^, FDR=3.24×10^−49^", manuscript, fixed = TRUE) &&
  grepl("P=5.39×10^−50, FDR=3.24×10^−49", figure_legend, fixed = TRUE) &&
  !grepl("3.24×10^−29", paste(manuscript, figure_legend), fixed = TRUE)
add_check("G13", "M04 P value and FDR agree across source, text and legend",
          m04_consistent, "TRUE", m04_consistent,
          "WGCNA source, manuscript and Figure 1 legend")

forbidden_pooling <- c(
  "random-effects estimate", "fixed- and random-effects summaries",
  "I²=86.7%", "random-effects confidence interval"
)
pooling_hits <- forbidden_pooling[vapply(
  forbidden_pooling,
  function(x) grepl(
    x, paste(manuscript, figure_legend), fixed = TRUE
  ),
  logical(1)
)]
add_check("G14", "No formal two-cohort pooled-effect language remains",
          paste(pooling_hits, collapse = "; "), "none",
          length(pooling_hits) == 0,
          "manuscript and Figure 5 legend")

clinical_file <- file.path(
  project_dir, "results", "PaperValidation", "SingleCellReplication",
  "singlecell_eligible_patient_clinical_metadata.csv"
)
clinical_dt <- if (file.exists(clinical_file)) fread(clinical_file) else data.table()
clinical_ok <-
  nrow(clinical_dt) == 21L &&
  clinical_dt[dataset == "GSE206785", .N] == 10L &&
  clinical_dt[dataset == "GSE270680", .N] == 11L &&
  clinical_dt[
    dataset == "GSE206785" & clinical_metadata_available == TRUE,
    .N
  ] == 9L &&
  clinical_dt[
    dataset == "GSE270680" & clinical_metadata_available == TRUE,
    .N
  ] == 11L
add_check("G15", "Eligible-patient clinical audit has expected coverage",
          if (nrow(clinical_dt)) {
            paste0(
              nrow(clinical_dt), " rows; ",
              clinical_dt[clinical_metadata_available == TRUE, .N],
              " with clinical metadata"
            )
          } else "missing",
          "21 rows; 20 with clinical metadata",
          clinical_ok,
          "singlecell_eligible_patient_clinical_metadata.csv")

audit <- rbindlist(checks, use.names = TRUE)
fwrite(audit, file.path(report_dir, "manuscript_numeric_audit_v2.csv"))
status <- if (all(audit$pass)) "PASS" else "FAIL"
failed <- audit[pass == FALSE]
md <- c(
  "# Manuscript numerical and file audit v2",
  "",
  paste0("- Audit time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- Overall status: **", status, "**"),
  paste0("- Checks passed: ", sum(audit$pass), "/", nrow(audit)),
  "",
  "## Failed checks",
  "",
  if (nrow(failed) == 0) {
    "None."
  } else {
    paste0(
      "- ", failed$id, ": ", failed$description,
      " (observed ", failed$observed, "; expected ", failed$expected, ")"
    )
  },
  "",
  "## Remaining author-side submission actions",
  "",
  "- All authors must review and approve the exact submitted version.",
  "- Upload the v2 scripts and derived source data to the existing repository/Zenodo record before submission.",
  "- Confirm journal-specific article type, word limits and portal declarations.",
  "- If requested by the journal or institution, obtain a written secondary-analysis ethics determination."
)
writeLines(
  md,
  file.path(report_dir, "manuscript_audit_report_v2.md"),
  useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(report_dir, "11_sessionInfo_v2.txt"),
  useBytes = TRUE
)

if (!all(audit$pass)) {
  print(failed)
  stop("稿件审计存在失败项目，请查看 manuscript_numeric_audit_v2.csv。")
}
cat("稿件审计通过：", nrow(audit), "项检查全部 PASS。\n")
