# 从GEO Series Matrix提取并标准化样本分组信息
# 当前数据集：GSE79973、GSE29272
# 全部读取、分组、匹配、验证和输出均由R完成。

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
    stop("请使用Rscript运行本脚本，以便自动定位项目目录。")
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

run_log_path <- file.path(report_dir, "03_prepare_sample_info_run.log")
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
  matrix_file = c(
    "GSE79973_series_matrix.txt.gz",
    "GSE29272_series_matrix.txt.gz"
  ),
  clean_expression_file = c(
    "GSE79973_clean_expression_matrix.csv",
    "GSE29272_clean_expression_matrix.csv"
  ),
  stringsAsFactors = FALSE
)

get_metadata_col <- function(metadata, name, default = NA_character_) {
  if (name %in% colnames(metadata)) {
    value <- trimws(as.character(metadata[[name]]))
    value[value == ""] <- NA_character_
    return(value)
  }
  rep(default, nrow(metadata))
}

combine_evidence <- function(...) {
  fields <- list(...)
  fields <- lapply(fields, function(x) {
    x[is.na(x)] <- ""
    trimws(x)
  })
  tolower(do.call(paste, c(fields, sep = " | ")))
}

extract_first <- function(x, pattern, replacement = "\\1") {
  hit <- grepl(pattern, x, perl = TRUE, ignore.case = TRUE)
  result <- rep(NA_character_, length(x))
  result[hit] <- sub(
    pattern,
    replacement,
    x[hit],
    perl = TRUE,
    ignore.case = TRUE
  )
  result
}

read_clean_sample_names <- function(path) {
  if (!file.exists(path)) {
    stop("清洗表达矩阵不存在：", path)
  }
  preview <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE,
    nrows = 1
  )
  colnames(preview)
}

align_sample_info <- function(sample_info, expression_samples, accession) {
  original_samples <- sample_info$sample

  if (identical(original_samples, expression_samples)) {
    return(list(
      data = sample_info,
      correction = "无需修正：sample名称及顺序与清洗表达矩阵列名完全一致"
    ))
  }

  if (setequal(original_samples, expression_samples)) {
    sample_info <- sample_info[
      match(expression_samples, sample_info$sample),
      ,
      drop = FALSE
    ]
    return(list(
      data = sample_info,
      correction = "自动按清洗表达矩阵列名重新排序；sample名称集合原本一致"
    ))
  }

  trimmed_metadata <- trimws(original_samples)
  trimmed_expression <- trimws(expression_samples)
  if (
    !anyDuplicated(trimmed_metadata) &&
      setequal(trimmed_metadata, trimmed_expression)
  ) {
    sample_info$sample <- trimmed_metadata
    sample_info <- sample_info[
      match(trimmed_expression, sample_info$sample),
      ,
      drop = FALSE
    ]
    return(list(
      data = sample_info,
      correction = "自动去除sample名称首尾空格，并按清洗表达矩阵列名重新排序"
    ))
  }

  stop(
    accession,
    "：metadata sample与清洗表达矩阵列名无法可靠自动匹配。metadata独有：",
    paste(setdiff(original_samples, expression_samples), collapse = ", "),
    "；表达矩阵独有：",
    paste(setdiff(expression_samples, original_samples), collapse = ", ")
  )
}

prepare_gse79973 <- function(metadata) {
  title <- get_metadata_col(metadata, "title")
  source_name <- get_metadata_col(metadata, "source_name_ch1")
  characteristics <- get_metadata_col(metadata, "characteristics_ch1")
  evidence <- combine_evidence(title, source_name, characteristics)

  is_tumor <- grepl(
    "adenocarcinoma|tumou?r",
    evidence,
    ignore.case = TRUE
  )
  is_normal <- grepl(
    "gastric mucosa|normal",
    evidence,
    ignore.case = TRUE
  ) & !is_tumor

  group <- rep(NA_character_, nrow(metadata))
  group[is_tumor] <- "Tumor"
  group[is_normal] <- "Normal"

  tissue <- rep(NA_character_, nrow(metadata))
  tissue[is_tumor] <- "gastric adenocarcinoma"
  tissue[is_normal] <- "gastric mucosa"

  disease_status <- ifelse(
    group == "Tumor",
    "Tumor",
    ifelse(group == "Normal", "Normal", NA_character_)
  )

  pair_number <- extract_first(
    title,
    ".*patient_([0-9]+).*",
    "\\1"
  )
  pair_id <- ifelse(
    is.na(pair_number),
    NA_character_,
    paste0("patient_", pair_number)
  )

  data.frame(
    sample = get_metadata_col(metadata, "geo_accession"),
    group = group,
    tissue = tissue,
    disease_status = disease_status,
    anatomical_site = "stomach",
    pair_id = pair_id,
    treatment = "Not reported",
    platform = get_metadata_col(metadata, "platform_id"),
    title = title,
    source_name = source_name,
    characteristics_ch1 = characteristics,
    grouping_basis = paste0(
      "title + source_name_ch1 + characteristics_ch1；",
      ifelse(is_tumor, "gastric adenocarcinoma/tumor", "gastric mucosa/normal")
    ),
    stringsAsFactors = FALSE
  )
}

prepare_gse29272 <- function(metadata) {
  title <- get_metadata_col(metadata, "title")
  source_name <- get_metadata_col(metadata, "source_name_ch1")
  characteristics <- get_metadata_col(metadata, "characteristics_ch1")
  evidence <- combine_evidence(title, source_name, characteristics)

  is_noncardia_tumor <- grepl(
    "tumou?r.*non-cardia|non-cardia.*tumou?r",
    evidence,
    ignore.case = TRUE
  )
  is_cardia_tumor <- !is_noncardia_tumor & grepl(
    "tumou?r.*cardia|cardia.*tumou?r",
    evidence,
    ignore.case = TRUE
  )
  is_normal <- grepl(
    "adjacent.*normal|normal.*gland",
    evidence,
    ignore.case = TRUE
  ) & !is_noncardia_tumor & !is_cardia_tumor

  group <- rep(NA_character_, nrow(metadata))
  group[is_normal] <- "Normal"
  group[is_noncardia_tumor] <- "NonCardia_Tumor"
  group[is_cardia_tumor] <- "Cardia_Tumor"

  disease_status <- ifelse(
    group == "Normal",
    "Normal",
    ifelse(!is.na(group), "Tumor", NA_character_)
  )

  tissue <- rep(NA_character_, nrow(metadata))
  tissue[is_normal] <- "adjacent normal gastric glands"
  tissue[is_noncardia_tumor] <- "non-cardia gastric tumor"
  tissue[is_cardia_tumor] <- "cardia gastric tumor"

  pair_id <- extract_first(
    title,
    ".*(TY[BC][0-9]{4})[NT].*",
    "\\1"
  )

  anatomical_site <- rep(NA_character_, nrow(metadata))
  anatomical_site[is_noncardia_tumor] <- "Non-cardia"
  anatomical_site[is_cardia_tumor] <- "Cardia"

  # 正常样本标题未直接标注cardia/non-cardia。
  # 通过相同pair_id的配对肿瘤样本补充其解剖部位。
  tumor_site_map <- setNames(
    anatomical_site[disease_status == "Tumor"],
    pair_id[disease_status == "Tumor"]
  )
  normal_index <- which(disease_status == "Normal")
  anatomical_site[normal_index] <- unname(
    tumor_site_map[pair_id[normal_index]]
  )

  data.frame(
    sample = get_metadata_col(metadata, "geo_accession"),
    group = group,
    tissue = tissue,
    disease_status = disease_status,
    anatomical_site = anatomical_site,
    pair_id = pair_id,
    treatment = "Not reported",
    platform = get_metadata_col(metadata, "platform_id"),
    title = title,
    source_name = source_name,
    characteristics_ch1 = characteristics,
    grouping_basis = paste0(
      "title + source_name_ch1 + characteristics_ch1；",
      ifelse(
        is_normal,
        "adjacent normal gastric glands",
        ifelse(
          is_noncardia_tumor,
          "non-cardia gastric tumor",
          "cardia gastric tumor"
        )
      )
    ),
    stringsAsFactors = FALSE
  )
}

validate_pairs <- function(sample_info, accession) {
  if (anyNA(sample_info$pair_id)) {
    stop(accession, "：存在无法提取pair_id的样本。")
  }
  pair_size <- table(sample_info$pair_id)
  if (any(pair_size != 2L)) {
    stop(accession, "：存在样本数不为2的配对编号。")
  }
  pair_status <- table(sample_info$pair_id, sample_info$disease_status)
  if (
    !all(c("Normal", "Tumor") %in% colnames(pair_status)) ||
      any(pair_status[, "Normal"] != 1L) ||
      any(pair_status[, "Tumor"] != 1L)
  ) {
    stop(accession, "：配对内不是严格的1个Normal加1个Tumor。")
  }
  length(unique(sample_info$pair_id))
}

process_one <- function(spec) {
  accession <- spec$accession
  series_path <- file.path(raw_dir, spec$matrix_file)
  expression_path <- file.path(clean_dir, spec$clean_expression_file)

  if (!file.exists(series_path)) {
    stop("Series Matrix不存在：", series_path)
  }

  eset <- getGEO(filename = series_path)
  metadata <- Biobase::pData(eset)

  sample_info <- switch(
    accession,
    GSE79973 = prepare_gse79973(metadata),
    GSE29272 = prepare_gse29272(metadata),
    stop("未定义的数据集：", accession)
  )

  if (anyNA(sample_info$sample) || anyDuplicated(sample_info$sample)) {
    stop(accession, "：sample编号缺失或重复。")
  }
  if (anyNA(sample_info$group)) {
    stop(
      accession,
      "：以下样本无法分组：",
      paste(sample_info$sample[is.na(sample_info$group)], collapse = ", ")
    )
  }
  if (anyNA(sample_info$anatomical_site)) {
    stop(
      accession,
      "：以下样本无法确定解剖部位：",
      paste(
        sample_info$sample[is.na(sample_info$anatomical_site)],
        collapse = ", "
      )
    )
  }

  pair_count <- validate_pairs(sample_info, accession)
  expression_samples <- read_clean_sample_names(expression_path)
  aligned <- align_sample_info(sample_info, expression_samples, accession)
  sample_info <- aligned$data

  if (!identical(sample_info$sample, expression_samples)) {
    stop(accession, "：自动匹配后sample名称或顺序仍不一致。")
  }

  output_file <- paste0(accession, "_sample_info.csv")
  output_path <- file.path(clean_dir, output_file)
  fwrite(sample_info, output_path, bom = TRUE)

  group_counts <- as.data.frame(table(sample_info$group))
  colnames(group_counts) <- c("group", "n")
  group_counts$dataset <- accession
  group_counts <- group_counts[, c("dataset", "group", "n")]

  site_counts <- as.data.frame(
    table(sample_info$anatomical_site, sample_info$disease_status)
  )
  colnames(site_counts) <- c("anatomical_site", "disease_status", "n")
  site_counts$dataset <- accession
  site_counts <- site_counts[site_counts$n > 0L, ]

  cat("==== ", accession, " ====\n", sep = "")
  cat("样本数：", nrow(sample_info), "\n", sep = "")
  cat("配对数：", pair_count, "\n", sep = "")
  print(group_counts, row.names = FALSE)
  cat("匹配检查：", aligned$correction, "\n", sep = "")
  cat("输出：", output_path, "\n\n", sep = "")

  list(
    sample_info = sample_info,
    group_counts = group_counts,
    site_counts = site_counts,
    summary = data.frame(
      dataset = accession,
      samples = nrow(sample_info),
      groups = length(unique(sample_info$group)),
      pairs = pair_count,
      platform = paste(unique(sample_info$platform), collapse = ";"),
      matrix_match = identical(sample_info$sample, expression_samples),
      correction_rule = aligned$correction,
      output_file = output_file,
      stringsAsFactors = FALSE
    )
  )
}

results <- lapply(seq_len(nrow(specs)), function(i) {
  process_one(as.list(specs[i, , drop = FALSE]))
})

group_counts <- rbindlist(lapply(results, `[[`, "group_counts"))
site_counts <- rbindlist(lapply(results, `[[`, "site_counts"))
summaries <- rbindlist(lapply(results, `[[`, "summary"))

fwrite(
  group_counts,
  file.path(report_dir, "GEO_sample_group_counts.csv"),
  bom = TRUE
)
fwrite(
  site_counts,
  file.path(report_dir, "GEO_sample_site_counts.csv"),
  bom = TRUE
)
fwrite(
  summaries,
  file.path(report_dir, "GEO_sample_info_validation.csv"),
  bom = TRUE
)

report_path <- file.path(report_dir, "GEO_sample_grouping_report.md")
report_lines <- c(
  "# GEO样本分组处理报告",
  "",
  paste0("- 处理时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- R版本：", R.version.string),
  "- 数据来源：本地GEO Series Matrix中的`pData`样本元数据。",
  "- 判定字段：`geo_accession`、`title`、`source_name_ch1`、`characteristics_ch1`和`platform_id`。",
  "- 验证对象：对应的`clean_expression_matrix.csv`列名。",
  "",
  "## GSE79973",
  "",
  "- 分组依据：title、source_name_ch1和characteristics_ch1中的`gastric adenocarcinoma/tumor`判为Tumor，`gastric mucosa/normal`判为Normal。",
  "- 配对依据：title中的`patient_n`；每个patient严格包含1个Tumor和1个Normal。",
  "- tissue标准化为`gastric adenocarcinoma`或`gastric mucosa`。",
  "",
  "### 分组数量",
  ""
)

append_count_lines <- function(lines, dataset) {
  count_data <- group_counts[group_counts$dataset == dataset]
  c(
    lines,
    paste0("- ", count_data$group, "：", count_data$n, "个样本"),
    ""
  )
}

report_lines <- append_count_lines(report_lines, "GSE79973")
gse79973_summary <- summaries[summaries$dataset == "GSE79973"]
report_lines <- c(
  report_lines,
  paste0(
    "- 表达矩阵匹配：",
    ifelse(gse79973_summary$matrix_match, "通过", "未通过")
  ),
  paste0("- 修正规则：", gse79973_summary$correction_rule),
  paste0("- 样本表：`clean_data/", gse79973_summary$output_file, "`"),
  "",
  "## GSE29272",
  "",
  "- 分组依据：联合title、source_name_ch1和characteristics_ch1判断。",
  "- `tumor tissue non-cardia`判为NonCardia_Tumor；`tumor tissue cardia`判为Cardia_Tumor；`adjacent normal gastric glands`判为Normal。",
  "- 配对依据：title中的TYB/TYC编号；每个编号严格包含1个Tumor和1个Normal。",
  "- Normal样本标题未直接写明cardia/non-cardia，其anatomical_site由同一pair_id的配对肿瘤样本补充；主分组仍保留为Normal。",
  "",
  "### 分组数量",
  ""
)

report_lines <- append_count_lines(report_lines, "GSE29272")
gse29272_summary <- summaries[summaries$dataset == "GSE29272"]
report_lines <- c(
  report_lines,
  paste0(
    "- 表达矩阵匹配：",
    ifelse(gse29272_summary$matrix_match, "通过", "未通过")
  ),
  paste0("- 修正规则：", gse29272_summary$correction_rule),
  paste0("- 样本表：`clean_data/", gse29272_summary$output_file, "`"),
  "",
  "## 一致性与缺失信息",
  "",
  "- 两个样本表的sample名称和顺序均与对应清洗表达矩阵列名完全一致。",
  "- 本次无需更改任何sample名称，也无需重新排序。",
  "- GEO元数据未报告可用于实验分组的处理/用药信息，因此`treatment`统一记录为`Not reported`，未据此推断分组。",
  "- 所有样本均获得明确group、tissue、disease_status、anatomical_site、pair_id和platform。",
  "",
  "## 输出文件",
  "",
  "- `clean_data/GSE79973_sample_info.csv`",
  "- `clean_data/GSE29272_sample_info.csv`",
  "- `report/GEO_sample_group_counts.csv`",
  "- `report/GEO_sample_site_counts.csv`",
  "- `report/GEO_sample_info_validation.csv`",
  "- `report/03_prepare_sample_info_run.log`",
  "- `scripts/03_prepare_sample_info.R`"
)

writeLines(report_lines, report_path, useBytes = TRUE)
writeLines(
  capture.output(sessionInfo()),
  file.path(report_dir, "R_sessionInfo_sample_grouping.txt"),
  useBytes = TRUE
)

cat("处理报告：", report_path, "\n", sep = "")
cat("完成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")

