# 比较两个单细胞队列的临床构成与实验设计，并评估可解释的异质性来源。
# 仅使用作者公开的元数据、论文补充表和既有患者配对结果。
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Seurat)
})

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) != 1L) stop("请使用 Rscript 运行本脚本。")
  normalizePath(sub("^--file=", "", f), winslash = "/", mustWork = TRUE)
}

project_dir <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  winslash = "/", mustWork = TRUE
)
raw_dir <- file.path(project_dir, "raw_data", "publication_supplements")
out_dir <- file.path(
  project_dir, "results", "PaperValidation", "SingleCellReplication"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A")] <- NA_character_
  x
}

stage_group <- function(x) {
  x <- normalize_stage(x)
  fifelse(x %chin% c("I", "II"), "Early (I-II)",
          fifelse(x %chin% c("III", "IV"), "Advanced (III-IV)", NA_character_))
}

histology_group <- function(x) {
  x <- tolower(trimws(as.character(x)))
  fifelse(x == "intestinal", "Intestinal",
          fifelse(!is.na(x) & nzchar(x), "Non-intestinal", NA_character_))
}

collapse_values <- function(x) {
  x <- sort(unique(na.omit(as.character(x))))
  if (length(x) == 0L) NA_character_ else paste(x, collapse = ";")
}

# GSE206785：论文补充表 Table S1 的第3行为字段名，第4行起为患者。
gse206_file <- file.path(
  raw_dir, "GSE206785_publication_TableS1_patient_characteristics.xlsx"
)
gse206_raw <- as.data.table(
  read_excel(gse206_file, sheet = "Table S1", col_names = FALSE)
)
headers206 <- as.character(gse206_raw[3, 2:27])
clinical206 <- copy(gse206_raw[4:nrow(gse206_raw), 2:27])
setnames(clinical206, make.unique(headers206))
clinical206 <- clinical206[
  !is.na(`Donor ID`) & nzchar(as.character(`Donor ID`)),
  .(
    patient = paste0("P", sub("^D", "", as.character(`Donor ID`))),
    age = suppressWarnings(as.numeric(Age)),
    sex = fifelse(Sex == "M", "Male", fifelse(Sex == "F", "Female", NA_character_)),
    stage = normalize_stage(Stage),
    histologic_subtype = as.character(`Lauren Classification`),
    anatomic_position = as.character(Position),
    differentiation = as.character(Differentiation)
  )
]

meta206 <- fread(file.path(project_dir, "raw_data", "GSE206785_metadata.txt.gz"))
platform206 <- meta206[, .(
  platform = collapse_values(sub("^SC", "", Platform))
), by = Patient]
setnames(platform206, "Patient", "patient")

pairs206_raw <- fread(file.path(
  project_dir, "results", "PaperValidation", "SingleCell",
  "epithelial_complete_pair_scores.csv"
))[eligible == TRUE]
pairs206 <- pairs206_raw[, .(
  n_cells_normal = n_cells[Tissue == "Normal"][1],
  n_cells_tumor = n_cells[Tissue == "Tumor"][1],
  module_score_normal = module_score[Tissue == "Normal"][1],
  module_score_tumor = module_score[Tissue == "Tumor"][1]
), by = Patient]
setnames(pairs206, "Patient", "patient")
pairs206[, module_difference := module_score_tumor - module_score_normal]

# GSE270680：临床字段取自作者公开Seurat对象，并与论文Supplementary Data 2核对。
gse270_file <- file.path(
  raw_dir, "GSE270680_publication_Supplementary_Data_1-10.xlsx"
)
published270 <- as.data.table(
  read_excel(gse270_file, sheet = "Data 2")
)
setnames(
  published270, c("Donor ID", "Subtype", "AJCC"),
  c("patient", "published_subtype", "published_stage")
)

obj270 <- readRDS(file.path(project_dir, "raw_data", "GSE270680", "sc.rds"))
meta270 <- as.data.table(obj270[[]])
clinical270 <- meta270[, .(
  age = suppressWarnings(as.numeric(as.character(Age[1]))),
  sex = as.character(Gender[1]),
  stage = normalize_stage(AJCC[1]),
  histologic_subtype = as.character(Subtype[1]),
  platform = collapse_values(platform)
), by = patient]
clinical270 <- merge(
  clinical270, published270, by = "patient", all.x = TRUE, sort = FALSE
)
if (clinical270[
  !is.na(published_subtype) & histologic_subtype != published_subtype,
  .N
] > 0L || clinical270[
  !is.na(published_stage) & stage != normalize_stage(published_stage),
  .N
] > 0L) {
  stop("GSE270680对象与论文Supplementary Data 2的临床字段不一致。")
}
clinical270[, c("published_subtype", "published_stage") := NULL]
rm(obj270)
gc()

pairs270 <- fread(file.path(
  out_dir, "GSE270680_epithelial_complete_pair_scores.csv"
))
setnames(
  pairs270,
  c("module_score_Normal", "module_score_Tumor",
    "n_cells_Normal", "n_cells_Tumor", "difference"),
  c("module_score_normal", "module_score_tumor",
    "n_cells_normal", "n_cells_tumor", "module_difference")
)

# 形成21名实际进入配对上皮分析的患者级表。
patients206 <- merge(pairs206, clinical206, by = "patient", all.x = TRUE)
patients206 <- merge(patients206, platform206, by = "patient", all.x = TRUE)
patients206[, `:=`(
  dataset = "GSE206785",
  clinical_metadata_available = !is.na(age),
  anatomic_position = fifelse(
    is.na(anatomic_position), "Not reported for this patient", anatomic_position
  ),
  differentiation = fifelse(
    is.na(differentiation), "Not reported for this patient", differentiation
  ),
  epithelial_selection = "EPCAM-positive cell depletion before sequencing"
)]

patients270 <- merge(pairs270, clinical270, by = "patient", all.x = TRUE)
patients270[, `:=`(
  dataset = "GSE270680",
  clinical_metadata_available = !is.na(age),
  anatomic_position = "Not reported",
  differentiation = "Not reported",
  epithelial_selection =
    "No EPCAM depletion reported; viable CD235a-negative cells selected by FACS"
)]

patient_table <- rbindlist(
  list(patients206, patients270),
  use.names = TRUE, fill = TRUE
)
patient_table[, stage_group := stage_group(stage)]
patient_table[, histology_group := histology_group(histologic_subtype)]
setcolorder(
  patient_table,
  c(
    "dataset", "patient", "age", "sex", "stage", "stage_group",
    "histologic_subtype", "histology_group", "anatomic_position",
    "differentiation", "platform", "epithelial_selection",
    "clinical_metadata_available", "n_cells_normal", "n_cells_tumor",
    "module_score_normal", "module_score_tumor", "module_difference"
  )
)
fwrite(
  patient_table,
  file.path(out_dir, "singlecell_eligible_patient_clinical_metadata.csv")
)

fmt_count <- function(dt, variable, level) {
  n_complete <- dt[!is.na(get(variable)), .N]
  n_level <- dt[get(variable) == level, .N]
  paste0(n_level, "/", n_complete)
}

fmt_age <- function(dt) {
  z <- dt[!is.na(age), age]
  if (!length(z)) return("Not reported")
  paste0(
    sprintf("%.1f", median(z)), " [",
    sprintf("%.1f", quantile(z, 0.25)), "-",
    sprintf("%.1f", quantile(z, 0.75)), "]"
  )
}

d206 <- patient_table[dataset == "GSE206785"]
d270 <- patient_table[dataset == "GSE270680"]
summary_table <- rbindlist(list(
  data.table(
    characteristic = c(
      "Full source cohort, donors",
      "Eligible paired epithelial analysis",
      "Eligible patients with published clinical metadata",
      "Age, median [IQR], years",
      "Male sex",
      "Stage I", "Stage II", "Stage III", "Stage IV",
      "Intestinal histology", "Non-intestinal histology",
      "Anatomic gastric position",
      "Single-cell platform",
      "Epithelial selection before sequencing"
    ),
    GSE206785 = c(
      "24", "10", paste0(sum(d206$clinical_metadata_available), "/10"),
      fmt_age(d206),
      fmt_count(d206, "sex", "Male"),
      fmt_count(d206, "stage", "I"),
      fmt_count(d206, "stage", "II"),
      fmt_count(d206, "stage", "III"),
      fmt_count(d206, "stage", "IV"),
      fmt_count(d206, "histology_group", "Intestinal"),
      fmt_count(d206, "histology_group", "Non-intestinal"),
      paste(sort(unique(d206$anatomic_position)), collapse = "; "),
      paste(sort(unique(d206$platform)), collapse = "; "),
      "EPCAM-positive depletion"
    ),
    GSE270680 = c(
      "27", "11", paste0(sum(d270$clinical_metadata_available), "/11"),
      fmt_age(d270),
      fmt_count(d270, "sex", "Male"),
      fmt_count(d270, "stage", "I"),
      fmt_count(d270, "stage", "II"),
      fmt_count(d270, "stage", "III"),
      fmt_count(d270, "stage", "IV"),
      fmt_count(d270, "histology_group", "Intestinal"),
      fmt_count(d270, "histology_group", "Non-intestinal"),
      "Not reported",
      paste(sort(unique(d270$platform)), collapse = "; "),
      "No EPCAM depletion reported"
    )
  )
))
fwrite(
  summary_table,
  file.path(out_dir, "singlecell_cohort_characteristics.csv")
)

# 对两个实际分析人群的共同临床字段做小样本、探索性比较。
safe_wilcox <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  if (sum(keep) < 4L || length(unique(g[keep])) < 2L) return(NA_real_)
  suppressWarnings(wilcox.test(x[keep] ~ g[keep], exact = FALSE)$p.value)
}

safe_fisher <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  tab <- table(x[keep], g[keep])
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  fisher.test(tab)$p.value
}

cohort_tests <- data.table(
  variable = c("Age", "Sex", "Stage group", "Histologic group", "Platform"),
  test = c(
    "Wilcoxon rank-sum", "Fisher exact", "Fisher exact",
    "Fisher exact", "Fisher exact"
  ),
  n_GSE206785 = c(
    sum(!is.na(d206$age)), sum(!is.na(d206$sex)),
    sum(!is.na(d206$stage_group)), sum(!is.na(d206$histology_group)),
    sum(!is.na(d206$platform))
  ),
  n_GSE270680 = c(
    sum(!is.na(d270$age)), sum(!is.na(d270$sex)),
    sum(!is.na(d270$stage_group)), sum(!is.na(d270$histology_group)),
    sum(!is.na(d270$platform))
  ),
  p_value = c(
    safe_wilcox(patient_table$age, patient_table$dataset),
    safe_fisher(patient_table$sex, patient_table$dataset),
    safe_fisher(patient_table$stage_group, patient_table$dataset),
    safe_fisher(patient_table$histology_group, patient_table$dataset),
    safe_fisher(patient_table$platform, patient_table$dataset)
  )
)
cohort_tests[, FDR := p.adjust(p_value, method = "BH")]
cohort_tests[, interpretation :=
  "Exploratory only; small selected analysis subsets and incomplete metadata"]
fwrite(
  cohort_tests,
  file.path(out_dir, "singlecell_cohort_comparison_tests.csv")
)

# 各队列内部按共同临床字段分层展示效应，不作因果归因。
stratify_effect <- function(dt, variable) {
  dt[!is.na(get(variable)), .(
    n = .N,
    mean_difference = mean(module_difference),
    sd_difference = if (.N > 1L) sd(module_difference) else NA_real_,
    median_difference = median(module_difference),
    min_difference = min(module_difference),
    max_difference = max(module_difference)
  ), by = c("dataset", variable)][
    , variable_name := variable
  ][
    , level := as.character(get(variable))
  ][
    , (variable) := NULL
  ]
}

stratified <- rbindlist(lapply(
  c("sex", "stage_group", "histology_group", "platform"),
  function(v) stratify_effect(patient_table, v)
), use.names = TRUE, fill = TRUE)
setcolorder(
  stratified,
  c(
    "dataset", "variable_name", "level", "n", "mean_difference",
    "sd_difference", "median_difference", "min_difference", "max_difference"
  )
)
fwrite(
  stratified,
  file.path(out_dir, "singlecell_clinical_stratified_effects.csv")
)

report <- c(
  "# Cross-cohort single-cell heterogeneity audit",
  "",
  paste0("- Analysis date: ", format(Sys.Date(), "%Y-%m-%d")),
  "- Unit of comparison: patients eligible for the prespecified paired epithelial analysis.",
  "- GSE206785: 10 eligible pairs; clinical characteristics were available for 9.",
  "- GSE270680: 11 eligible pairs; clinical characteristics were available for all 11.",
  "",
  "## Directly observed design difference",
  "",
  paste(
    "GSE206785 used EPCAM-positive cell depletion before single-cell sequencing.",
    "The retained epithelial cells therefore represent an incompletely depleted subset.",
    "GSE270680 reported viable-cell sorting with CD235a exclusion and no EPCAM depletion."
  ),
  "",
  "## Clinical comparability",
  "",
  paste(
    "Age, sex, stage and histologic subtype could be compared descriptively.",
    "The published GSE206785 table also reported gastric position and differentiation;",
    "corresponding fields were not reported for GSE270680, so these factors could not be tested."
  ),
  "",
  "## Interpretation",
  "",
  paste(
    "The tables document real differences in available clinical composition and epithelial",
    "sampling protocols. Because the analyzed subsets contain only 10 and 11 patients,",
    "the exploratory tests cannot identify a cause of the discordant transcriptomic effects.",
    "The protocol difference is a concrete source of selection heterogeneity, whereas",
    "stage, histology and position remain incompletely evaluated explanations."
  )
)
writeLines(
  report,
  file.path(out_dir, "singlecell_cohort_heterogeneity_report.md"),
  useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "20_sessionInfo.txt"),
  useBytes = TRUE
)

cat("单细胞队列异质性审计完成：", out_dir, "\n")
