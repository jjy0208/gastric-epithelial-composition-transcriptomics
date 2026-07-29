# GEO 芯片表达矩阵清洗
# 数据集：GSE79973、GSE29272
# 作用：读取本地 Series Matrix 和 GEO 平台注释，完成探针注释、聚合、
#       缺失/低表达过滤、结果与质量报告输出。
# 注意：不执行差异表达分析。

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("请使用 Rscript 运行本脚本，以便自动定位项目目录。")
  }
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_dir <- normalizePath(
  file.path(dirname(script_path), ".."),
  winslash = "/",
  mustWork = TRUE
)
raw_dir <- file.path(project_dir, "raw_data")
clean_dir <- file.path(project_dir, "clean_data")
report_dir <- file.path(project_dir, "report")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

run_log_path <- file.path(report_dir, "02_clean_GEO_microarray_run.log")
run_log_con <- file(run_log_path, open = "wt", encoding = "UTF-8")
sink(run_log_con, type = "output")
sink(run_log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(run_log_con)
}, add = TRUE)

cat("开始时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
cat("R：", R.version.string, "\n", sep = "")
cat("项目目录：", project_dir, "\n\n", sep = "")

specs <- data.frame(
  accession = c("GSE79973", "GSE29272"),
  platform = c("GPL570", "GPL96"),
  matrix_file = c(
    "GSE79973_series_matrix.txt.gz",
    "GSE29272_series_matrix.txt.gz"
  ),
  annotation_file = c(
    "GPL570.annot.gz",
    "GPL96.soft.gz"
  ),
  source_preprocessing = c(
    "GEO提交者说明：MAS 5.0标准化；上传值呈log2样强度分布",
    "GEO提交者说明：RMA标准化；上传值呈log2样强度分布"
  ),
  stringsAsFactors = FALSE
)

required_files <- unique(c(specs$matrix_file, specs$annotation_file))
missing_files <- required_files[
  !file.exists(file.path(raw_dir, required_files))
]
if (length(missing_files) > 0L) {
  stop("缺少必需文件：", paste(missing_files, collapse = ", "))
}

find_col <- function(tab, candidates) {
  idx <- match(tolower(candidates), tolower(colnames(tab)))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    stop(
      "平台注释中未找到候选列：",
      paste(candidates, collapse = ", "),
      "；实际列：",
      paste(colnames(tab), collapse = ", ")
    )
  }
  colnames(tab)[idx[1L]]
}

classify_symbols <- function(x) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | x == "" | x == "---" |
    toupper(x) %in% c("NA", "N/A", "NULL")
  ambiguous <- !missing & (
    grepl("///|//|;|,", x) |
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", x)
  )
  clean <- x
  clean[missing | ambiguous] <- NA_character_
  list(clean = clean, missing = missing, ambiguous = ambiguous)
}

impute_row_median <- function(mat) {
  if (!anyNA(mat)) {
    return(mat)
  }
  affected <- which(rowSums(is.na(mat)) > 0L)
  for (i in affected) {
    mat[i, is.na(mat[i, ])] <- median(mat[i, ], na.rm = TRUE)
  }
  mat
}

process_one <- function(spec) {
  accession <- spec$accession
  cat("==== ", accession, " ====\n", sep = "")

  matrix_path <- file.path(raw_dir, spec$matrix_file)
  annotation_path <- file.path(raw_dir, spec$annotation_file)

  eset <- getGEO(filename = matrix_path)
  expression <- Biobase::exprs(eset)
  storage.mode(expression) <- "double"
  sample_info <- Biobase::pData(eset)

  if (ncol(expression) != nrow(sample_info)) {
    stop(accession, "：表达矩阵样本数与样本信息行数不一致。")
  }
  if (!identical(colnames(expression), rownames(sample_info))) {
    stop(accession, "：表达矩阵样本顺序与样本信息不一致。")
  }
  if (!identical(Biobase::annotation(eset), spec$platform)) {
    stop(
      accession, "：Series Matrix平台为 ", Biobase::annotation(eset),
      "，但预期为 ", spec$platform, "。"
    )
  }

  raw_probe_n <- nrow(expression)
  raw_sample_n <- ncol(expression)
  raw_na_n <- sum(!is.finite(expression))
  expression[!is.finite(expression)] <- NA_real_

  missing_fraction <- rowMeans(is.na(expression))
  keep_missing <- missing_fraction <= 0.20
  removed_missing_n <- sum(!keep_missing)
  expression_after_missing <- expression[keep_missing, , drop = FALSE]
  expression_after_missing <- impute_row_median(expression_after_missing)

  gpl <- getGEO(filename = annotation_path)
  annotation_table <- GEOquery::Table(gpl)
  id_col <- find_col(annotation_table, c("ID", "ID_REF"))
  symbol_col <- find_col(
    annotation_table,
    c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol")
  )

  matched_index <- match(
    rownames(expression_after_missing),
    as.character(annotation_table[[id_col]])
  )
  raw_symbols <- annotation_table[[symbol_col]][matched_index]
  symbol_status <- classify_symbols(raw_symbols)
  gene_symbols <- symbol_status$clean

  unannotated_n <- sum(symbol_status$missing)
  ambiguous_n <- sum(symbol_status$ambiguous)
  annotated_probe_n <- sum(!is.na(gene_symbols))
  annotation_rate <- annotated_probe_n / raw_probe_n

  annotated_expression <- expression_after_missing[
    !is.na(gene_symbols), ,
    drop = FALSE
  ]
  annotated_symbols <- gene_symbols[!is.na(gene_symbols)]

  sample_columns <- colnames(annotated_expression)
  expression_dt <- as.data.table(annotated_expression)
  expression_dt[, GeneSymbol := annotated_symbols]
  setcolorder(expression_dt, c("GeneSymbol", sample_columns))

  collapsed_dt <- expression_dt[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = GeneSymbol,
    .SDcols = sample_columns
  ]
  collapsed_matrix <- as.matrix(collapsed_dt[, ..sample_columns])
  rownames(collapsed_matrix) <- collapsed_dt$GeneSymbol
  storage.mode(collapsed_matrix) <- "double"

  genes_before_filter_n <- nrow(collapsed_matrix)
  gene_sd <- apply(collapsed_matrix, 1L, sd)
  keep_nonconstant <- is.finite(gene_sd) & gene_sd > 0
  removed_constant_n <- sum(!keep_nonconstant)
  filtered_matrix <- collapsed_matrix[keep_nonconstant, , drop = FALSE]

  probe_row_means <- rowMeans(expression_after_missing)
  centered_scale <- (
    as.numeric(quantile(expression_after_missing, 0.01, na.rm = TRUE)) < 0 &&
      median(abs(probe_row_means), na.rm = TRUE) < 0.20
  )

  if (centered_scale) {
    low_expression_method <- paste0(
      "未执行绝对低表达过滤：矩阵以0为中心，数值不代表绝对表达强度；",
      "仅移除零方差基因"
    )
    low_expression_cutoff <- NA_real_
    removed_low_expression_n <- 0L
  } else {
    gene_medians <- apply(filtered_matrix, 1L, median, na.rm = TRUE)
    low_expression_cutoff <- as.numeric(
      quantile(gene_medians, 0.10, na.rm = TRUE, names = FALSE)
    )
    keep_expression <- gene_medians > low_expression_cutoff
    removed_low_expression_n <- sum(!keep_expression)
    filtered_matrix <- filtered_matrix[keep_expression, , drop = FALSE]
    low_expression_method <- paste0(
      "移除基因中位表达不高于第10百分位的最低表达基因；阈值=",
      format(low_expression_cutoff, digits = 6)
    )
  }

  filtered_matrix <- filtered_matrix[
    order(rownames(filtered_matrix)),
    ,
    drop = FALSE
  ]
  if (anyDuplicated(rownames(filtered_matrix)) > 0L) {
    stop(accession, "：输出Gene Symbol仍有重复。")
  }
  if (any(!is.finite(filtered_matrix))) {
    stop(accession, "：输出矩阵仍包含NA/Inf。")
  }
  if (!identical(colnames(filtered_matrix), colnames(expression))) {
    stop(accession, "：输出矩阵样本顺序发生变化。")
  }

  output_file <- paste0(accession, "_clean_expression_matrix.csv")
  output_path <- file.path(clean_dir, output_file)
  write.csv(
    filtered_matrix,
    file = output_path,
    row.names = TRUE,
    quote = FALSE,
    na = ""
  )

  final_gene_n <- nrow(filtered_matrix)
  probe_exclusion_rate <- 1 - annotated_probe_n / raw_probe_n
  total_row_reduction_rate <- 1 - final_gene_n / raw_probe_n

  cat("原始维度：", raw_probe_n, " × ", raw_sample_n, "\n", sep = "")
  cat("原始非有限值数量：", raw_na_n, "\n", sep = "")
  cat("缺失比例>20%而移除的探针：", removed_missing_n, "\n", sep = "")
  cat("成功单一注释探针：", annotated_probe_n, "\n", sep = "")
  cat("注释后、合并探针前基因数：", genes_before_filter_n, "\n", sep = "")
  cat("最终维度：", final_gene_n, " × ", raw_sample_n, "\n", sep = "")
  cat("输出：", output_path, "\n\n", sep = "")

  data.frame(
    dataset = accession,
    platform = spec$platform,
    matrix_file = spec$matrix_file,
    sample_info_location = "Series Matrix内部!Sample_*字段",
    annotation_file = spec$annotation_file,
    id_type = "Affymetrix芯片探针ID",
    source_preprocessing = spec$source_preprocessing,
    raw_probes = raw_probe_n,
    samples = raw_sample_n,
    raw_nonfinite_values = raw_na_n,
    probes_removed_missing_gt20pct = removed_missing_n,
    probes_unannotated = unannotated_n,
    probes_ambiguous_annotation = ambiguous_n,
    probes_successfully_annotated = annotated_probe_n,
    annotation_rate_pct = round(100 * annotation_rate, 4),
    genes_after_probe_averaging = genes_before_filter_n,
    zero_variance_genes_removed = removed_constant_n,
    low_expression_method = low_expression_method,
    low_expression_cutoff = low_expression_cutoff,
    low_expression_genes_removed = removed_low_expression_n,
    final_genes = final_gene_n,
    final_samples = ncol(filtered_matrix),
    probe_exclusion_before_aggregation_pct = round(
      100 * probe_exclusion_rate, 4
    ),
    total_row_reduction_pct = round(100 * total_row_reduction_rate, 4),
    output_file = output_file,
    output_path = normalizePath(output_path, winslash = "/"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

metrics_list <- lapply(seq_len(nrow(specs)), function(i) {
  process_one(as.list(specs[i, , drop = FALSE]))
})
metrics <- rbindlist(metrics_list, fill = TRUE)

metrics_path <- file.path(report_dir, "GEO_cleaning_quality_metrics.csv")
fwrite(metrics, metrics_path, bom = TRUE)

report_path <- file.path(report_dir, "GEO_data_cleaning_report.md")
report_lines <- c(
  "# GEO表达芯片数据清洗报告",
  "",
  paste0("- 处理时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- R版本：", R.version.string),
  "- 处理范围：探针注释、缺失值处理、同基因多探针取均值、低表达/非信息性过滤。",
  "- 未执行差异表达分析，也未跨数据集合并或批次校正。",
  "",
  "## raw_data文件识别",
  "",
  "- `GSE*_series_matrix.txt.gz`：GEO Series Matrix；同一文件内包含探针表达矩阵和`!Sample_*`样本元数据。",
  "- `GPL570.annot.gz`：GPL570官方平台注释，用于GSE79973。",
  "- `GPL96.soft.gz`：GPL96官方提交者平台注释，用于GSE29272。",
  "- `download_log.txt`：下载过程记录，不是表达数据或注释表。",
  "",
  "## 通用处理规则",
  "",
  "1. 使用GEOquery读取本地Series Matrix和平台注释。",
  "2. 两个数据集的行ID均判定为Affymetrix芯片探针ID，不是Gene Symbol或Ensembl ID。",
  "3. 移除缺失/非有限值比例大于20%的探针；其余少量缺失值按探针中位数填补。",
  "4. 仅保留能唯一映射到一个合法Gene Symbol的探针；无注释和多基因歧义探针均移除。",
  "5. 同一Gene Symbol对应多个探针时，对各样本的探针表达值取算术平均。",
  "6. 移除零方差基因。对具有绝对log2样强度的矩阵，进一步移除基因中位表达最低10%；对以0为中心的矩阵不进行绝对低表达过滤。",
  "",
  "## 分数据集结果",
  ""
)

for (i in seq_len(nrow(metrics))) {
  x <- metrics[i]
  report_lines <- c(
    report_lines,
    paste0("### ", x$dataset),
    "",
    paste0("- 平台：", x$platform),
    paste0("- 表达矩阵与样本信息：`raw_data/", x$matrix_file, "`"),
    paste0("- 平台注释：`raw_data/", x$annotation_file, "`"),
    paste0("- 数据来源预处理：", x$source_preprocessing),
    paste0("- 原始维度：", x$raw_probes, "个探针 × ", x$samples, "个样本"),
    paste0("- 原始NA/Inf数量：", x$raw_nonfinite_values),
    paste0("- 因缺失比例>20%移除的探针：", x$probes_removed_missing_gt20pct),
    paste0("- 无有效注释探针：", x$probes_unannotated),
    paste0("- 多基因或格式歧义探针：", x$probes_ambiguous_annotation),
    paste0(
      "- 成功唯一注释探针：", x$probes_successfully_annotated,
      "；成功注释率：", format(x$annotation_rate_pct, nsmall = 4), "%"
    ),
    paste0("- 多探针取均值后基因数：", x$genes_after_probe_averaging),
    paste0("- 零方差基因移除数：", x$zero_variance_genes_removed),
    paste0("- 低表达处理：", x$low_expression_method),
    paste0("- 低表达基因移除数：", x$low_expression_genes_removed),
    paste0("- 最终维度：", x$final_genes, "个基因 × ", x$final_samples, "个样本"),
    paste0(
      "- 注释前探针排除比例：",
      format(x$probe_exclusion_before_aggregation_pct, nsmall = 4), "%"
    ),
    paste0(
      "- 总行数缩减比例（含多探针合并）：",
      format(x$total_row_reduction_pct, nsmall = 4), "%"
    ),
    paste0("- 输出：`clean_data/", x$output_file, "`"),
    ""
  )
}

report_lines <- c(
  report_lines,
  "## 解释限制",
  "",
  "- 两个数据集均为表达芯片数据，不是RNA-seq计数矩阵。",
  "- 最低10%过滤是对GSE79973和GSE29272采用的保守操作性阈值；后续如有预先规定的研究方案，可在脚本中调整。",
  "- Gene Symbol来自当前下载的GEO平台注释，未执行跨版本基因命名更新。",
  "",
  "## 可复现文件",
  "",
  "- `scripts/02_clean_GEO_microarray.R`：本报告及清洗矩阵的生成脚本。",
  "- `report/GEO_cleaning_quality_metrics.csv`：全部质量指标的机器可读汇总。",
  "- `report/02_clean_GEO_microarray_run.log`：R运行日志。",
  "- `report/R_sessionInfo.txt`：R及依赖包版本。"
)

writeLines(report_lines, report_path, useBytes = TRUE)
writeLines(
  capture.output(sessionInfo()),
  file.path(report_dir, "R_sessionInfo.txt"),
  useBytes = TRUE
)

cat("质量指标：", metrics_path, "\n", sep = "")
cat("处理报告：", report_path, "\n", sep = "")
cat("完成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
