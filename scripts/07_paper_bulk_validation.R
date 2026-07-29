# GSE29272胃上皮身份丢失模块：独立bulk芯片队列验证
# 验证队列：GSE79973（10对）和GSE19826（12对，另有3份未配对正常组织不进入配对主分析）
# 统计单位：患者；Tumor - Normal为效应方向

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(20260729)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) stop("请使用Rscript运行本脚本。")
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "raw_data")
clean_dir <- file.path(project_dir, "clean_data")
result_dir <- file.path(project_dir, "results", "PaperValidation")
report_dir <- file.path(project_dir, "report", "manuscript")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_dir, "07_paper_bulk_validation_run.log")
log_con <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("开始时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
cat("R版本：", R.version.string, "\n")

candidate_path <- file.path(
  project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"
)
candidate <- fread(candidate_path)
stopifnot(nrow(candidate) == 73L, uniqueN(candidate$gene) == 73L)
candidate_genes <- candidate$gene
discovery_lfc <- setNames(candidate$DEG_log2FoldChange, candidate$gene)

find_col <- function(tab, candidates) {
  idx <- match(tolower(candidates), tolower(colnames(tab)))
  idx <- idx[!is.na(idx)]
  if (!length(idx)) stop("平台注释缺少基因Symbol列。")
  colnames(tab)[idx[1L]]
}

classify_symbols <- function(x) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | x == "" | x == "---" |
    toupper(x) %in% c("NA", "N/A", "NULL") |
    grepl("///|//|;|,", x) |
    !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", x)
  x[bad] <- NA_character_
  x
}

read_clean_matrix <- function(path) {
  x <- fread(path)
  genes <- x[[1L]]
  x[[1L]] <- NULL
  mat <- as.matrix(x)
  rownames(mat) <- genes
  storage.mode(mat) <- "double"
  mat
}

prepare_gse19826 <- function() {
  matrix_path <- file.path(raw_dir, "GSE19826_series_matrix.txt.gz")
  gpl_path <- file.path(raw_dir, "GPL570.annot.gz")
  if (!file.exists(matrix_path) || !file.exists(gpl_path)) {
    stop("GSE19826矩阵或GPL570注释文件缺失。")
  }

  eset <- getGEO(filename = matrix_path)
  expr <- Biobase::exprs(eset)
  storage.mode(expr) <- "double"
  # GSE19826 的series matrix为未取对数的正值信号强度。
  # 在探针合并前转换为log2(signal + 1)，使效应量和后续去卷积处于常用芯片尺度。
  q99 <- unname(quantile(expr, 0.99, na.rm = TRUE))
  if (is.finite(q99) && q99 > 100) {
    expr <- log2(pmax(expr, 0) + 1)
  }
  pd <- Biobase::pData(eset)
  stopifnot(identical(colnames(expr), rownames(pd)))

  gpl <- getGEO(filename = gpl_path)
  ann <- GEOquery::Table(gpl)
  id_col <- find_col(ann, c("ID", "ID_REF"))
  symbol_col <- find_col(ann, c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol"))
  symbol <- classify_symbols(
    ann[[symbol_col]][match(rownames(expr), as.character(ann[[id_col]]))]
  )
  keep <- !is.na(symbol) & rowSums(is.finite(expr)) == ncol(expr)
  expr <- expr[keep, , drop = FALSE]
  symbol <- symbol[keep]

  dt <- as.data.table(expr)
  dt[, gene := symbol]
  collapsed <- dt[, lapply(.SD, mean), by = gene, .SDcols = colnames(expr)]
  clean <- as.matrix(collapsed[, -1L])
  rownames(clean) <- collapsed$gene
  storage.mode(clean) <- "double"
  clean <- clean[apply(clean, 1L, sd) > 0, , drop = FALSE]
  clean <- clean[order(rownames(clean)), , drop = FALSE]

  title <- as.character(pd$title)
  source <- as.character(pd$source_name_ch1)
  group <- ifelse(grepl("gastric cancer", source, ignore.case = TRUE), "Tumor", "Normal")
  pair_id <- ifelse(
    grepl("[0-9]+[NT]$", title),
    sub("[NT]$", "", title),
    NA_character_
  )
  meta <- data.frame(
    sample = rownames(pd),
    group = group,
    pair_id = pair_id,
    title = title,
    source_name = source,
    platform = "GPL570",
    stringsAsFactors = FALSE
  )
  stopifnot(identical(colnames(clean), meta$sample))

  fwrite(
    data.table(gene = rownames(clean), as.data.table(clean)),
    file.path(clean_dir, "GSE19826_clean_expression_matrix.csv")
  )
  fwrite(meta, file.path(clean_dir, "GSE19826_sample_info.csv"))
  list(expr = clean, meta = meta)
}

paired_validation <- function(dataset, expr, meta) {
  stopifnot(identical(colnames(expr), meta$sample))
  genes <- intersect(candidate_genes, rownames(expr))
  if (length(genes) < 50L) stop(dataset, "可用候选基因不足50个。")

  # 每个基因在本队列内标准化后等权计算模块分数，避免高表达基因支配分数
  z <- t(scale(t(expr[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  score <- colMeans(z)
  score_dt <- data.table(
    dataset = dataset,
    sample = colnames(expr),
    group = meta$group,
    pair_id = meta$pair_id,
    module_score = score
  )

  pair_counts <- score_dt[!is.na(pair_id), uniqueN(group), by = pair_id]
  good_pairs <- pair_counts[V1 == 2L, pair_id]
  paired <- score_dt[pair_id %in% good_pairs]
  wide <- dcast(paired, pair_id ~ group, value.var = "module_score")
  stopifnot(all(c("Normal", "Tumor") %in% colnames(wide)))
  diff <- wide$Tumor - wide$Normal

  tt <- t.test(wide$Tumor, wide$Normal, paired = TRUE)
  wt <- suppressWarnings(wilcox.test(
    wide$Tumor, wide$Normal, paired = TRUE, exact = FALSE
  ))

  paired_samples <- meta$sample[meta$pair_id %in% good_pairs]
  paired_expr <- expr[genes, paired_samples, drop = FALSE]
  paired_meta <- meta[match(paired_samples, meta$sample), ]
  tumor_mat <- paired_expr[, paired_meta$group == "Tumor", drop = FALSE]
  normal_mat <- paired_expr[, paired_meta$group == "Normal", drop = FALSE]
  tumor_pair <- paired_meta$pair_id[paired_meta$group == "Tumor"]
  normal_pair <- paired_meta$pair_id[paired_meta$group == "Normal"]
  normal_mat <- normal_mat[, match(tumor_pair, normal_pair), drop = FALSE]
  gene_lfc <- rowMeans(tumor_mat - normal_mat)

  gene_dt <- data.table(
    dataset = dataset,
    gene = genes,
    discovery_log2FC = discovery_lfc[genes],
    validation_mean_paired_difference = gene_lfc,
    direction_concordant = gene_lfc < 0
  )
  rho <- suppressWarnings(cor.test(
    gene_dt$discovery_log2FC,
    gene_dt$validation_mean_paired_difference,
    method = "spearman",
    exact = FALSE
  ))
  bt <- binom.test(sum(gene_dt$direction_concordant), nrow(gene_dt), p = 0.5)

  summary <- data.table(
    dataset = dataset,
    platform = unique(meta$platform),
    available_module_genes = length(genes),
    paired_patients = nrow(wide),
    mean_tumor_minus_normal_score = mean(diff),
    score_difference_ci_low = unname(tt$conf.int[1L]),
    score_difference_ci_high = unname(tt$conf.int[2L]),
    paired_t_pvalue = tt$p.value,
    paired_wilcoxon_pvalue = wt$p.value,
    down_direction_genes = sum(gene_dt$direction_concordant),
    direction_concordance = mean(gene_dt$direction_concordant),
    binomial_pvalue = bt$p.value,
    lfc_spearman_rho = unname(rho$estimate),
    lfc_spearman_pvalue = rho$p.value
  )
  list(score = score_dt, paired = paired, gene = gene_dt, summary = summary)
}

# GSE79973
gse79973_expr <- read_clean_matrix(
  file.path(clean_dir, "GSE79973_clean_expression_matrix.csv")
)
gse79973_meta <- fread(file.path(clean_dir, "GSE79973_sample_info.csv"))
gse79973_meta <- as.data.frame(gse79973_meta)
v79973 <- paired_validation("GSE79973", gse79973_expr, gse79973_meta)

# GSE19826
gse19826 <- prepare_gse19826()
v19826 <- paired_validation("GSE19826", gse19826$expr, gse19826$meta)

score_all <- rbindlist(list(v79973$score, v19826$score), fill = TRUE)
paired_all <- rbindlist(list(v79973$paired, v19826$paired), fill = TRUE)
gene_all <- rbindlist(list(v79973$gene, v19826$gene), fill = TRUE)
summary_all <- rbindlist(list(v79973$summary, v19826$summary), fill = TRUE)

fwrite(score_all, file.path(result_dir, "external_bulk_module_scores.csv"))
fwrite(paired_all, file.path(result_dir, "external_bulk_paired_scores.csv"))
fwrite(gene_all, file.path(result_dir, "external_bulk_gene_direction_validation.csv"))
fwrite(summary_all, file.path(result_dir, "external_bulk_validation_summary.csv"))

palette <- c(Normal = "#8FA6A0", Tumor = "#B56F70")
plot_dt <- copy(paired_all)
plot_dt[, group := factor(group, levels = c("Normal", "Tumor"))]
plot_dt[, pair_key := paste(dataset, pair_id, sep = "_")]

p <- ggplot(plot_dt, aes(group, module_score, group = pair_key)) +
  geom_line(colour = "#B8B8B8", linewidth = 0.35, alpha = 0.75) +
  geom_point(aes(colour = group), size = 1.8, alpha = 0.95) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_colour_manual(values = palette) +
  labs(
    x = NULL,
    y = "73-gene module score (within-cohort z score)"
  ) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35)
  )

ggsave(
  file.path(result_dir, "external_bulk_paired_module_validation.pdf"),
  p, width = 183 / 25.4, height = 85 / 25.4, device = cairo_pdf
)
ggsave(
  file.path(result_dir, "external_bulk_paired_module_validation.png"),
  p, width = 183 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white"
)

report <- c(
  "# 外部bulk队列验证报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "- 锁定模块：GSE29272的73个DEG-WGCNA交集基因；未在验证队列重新筛基因。",
  "- 主统计单位：患者；仅完整Tumor-Normal配对进入主检验。",
  "- 模块分数：各队列内逐基因Z标准化后，对可用模块基因等权取平均。",
  "",
  "## 结果",
  "",
  paste(capture.output(print(summary_all)), collapse = "\n"),
  "",
  "## 解释边界",
  "",
  "- 两个验证队列均为Affymetrix GPL570芯片，与发现队列GPL96平台不同。",
  "- GSE19826的3份额外正常组织没有配对肿瘤，不进入配对主分析。",
  "- 队列内Z分数适合检验方向和患者内差异，不用于跨平台比较绝对表达量。",
  "- 外部复现支持模块稳定性，但不能单独区分细胞组成变化和细胞内在调控。"
)
writeLines(report, file.path(result_dir, "external_bulk_validation_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "07_sessionInfo.txt"), useBytes = TRUE)

cat("\n验证摘要：\n")
print(summary_all)
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
