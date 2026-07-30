# GSE270680 第二个独立单细胞队列交叉验证
# 目的：使用作者公开的 Seurat 对象及原始细胞注释，重复检验锁定的
#       73 基因模块在患者配对上皮 pseudobulk 中的肿瘤-邻癌差异。
# 统计原则：
#   1) 患者是独立统计单位，细胞不作为生物学重复。
#   2) 每个患者-组织至少 20 个上皮细胞才进入配对检验。
#   3) 直接使用作者 majorCluster/subCluster，不用 73 基因重新注释细胞。
#   4) 与 GSE206785 分队列分析，仅并列展示标准化配对效应，不计算合并效应。

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(sandwich)
  library(lmtest)
})

options(stringsAsFactors = FALSE)
set.seed(20260729)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) stop("请使用 Rscript 运行本脚本。")
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_dir <- normalizePath(file.path(dirname(script_path), ".."),
                             winslash = "/", mustWork = TRUE)
raw_file <- file.path(project_dir, "raw_data", "GSE270680", "sc.rds")
result_dir <- file.path(
  project_dir, "results", "PaperValidation", "SingleCellReplication"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_dir, "14_GSE270680_singlecell_replication_run.log")
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
if (!file.exists(raw_file)) stop("缺少作者公开的 sc.rds：", raw_file)

candidate <- fread(file.path(
  project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"
))
stopifnot("gene" %in% names(candidate))
candidate_genes <- unique(trimws(candidate$gene))
candidate_genes <- candidate_genes[nzchar(candidate_genes)]
stopifnot(length(candidate_genes) == 73L)

cat("读取作者公开 Seurat 对象（约 14.4 GB 内存对象）……\n")
sc <- readRDS(raw_file)
if (!inherits(sc, "Seurat")) stop("sc.rds 不是 Seurat 对象。")

required_meta <- c(
  "patient", "library", "loc", "majorCluster", "subCluster"
)
if (!all(required_meta %in% colnames(sc@meta.data))) {
  stop("Seurat metadata 缺少字段：",
       paste(setdiff(required_meta, colnames(sc@meta.data)), collapse = ", "))
}
if (!all(c("counts", "data") %in% Layers(sc[["RNA"]]))) {
  stop("RNA assay 必须同时包含 counts 和 data layer。")
}

meta <- as.data.table(sc@meta.data, keep.rownames = "cell")
meta[, `:=`(
  patient = as.character(patient),
  library = as.character(library),
  loc = as.character(loc),
  majorCluster = as.character(majorCluster),
  subCluster = as.character(subCluster)
)]
stopifnot(identical(meta$cell, colnames(sc)))

candidate_available <- intersect(candidate_genes, rownames(sc))
candidate_missing <- setdiff(candidate_genes, candidate_available)
if (length(candidate_available) < 65L) {
  stop("候选基因覆盖不足 65 个。")
}
cat("候选基因覆盖：", length(candidate_available), "/73\n", sep = "")

audit <- data.table(
  metric = c(
    "cells", "genes", "patients", "libraries", "locations",
    "candidate_genes_available", "candidate_genes_missing"
  ),
  value = c(
    ncol(sc), nrow(sc), uniqueN(meta$patient), uniqueN(meta$library),
    uniqueN(meta$loc), length(candidate_available), length(candidate_missing)
  )
)
fwrite(audit, file.path(result_dir, "GSE270680_object_audit.csv"))
fwrite(
  data.table(gene = candidate_missing),
  file.path(result_dir, "GSE270680_missing_candidate_genes.csv")
)

# 仅提取 73 个候选基因的 log-normalized data layer，避免复制完整表达矩阵。
expr_candidate <- LayerData(
  sc, assay = "RNA", layer = "data",
  features = candidate_available, cells = meta$cell
)
stopifnot(nrow(expr_candidate) == length(candidate_available))

aggregate_sparse <- function(expr, group, group_columns) {
  group <- factor(group, levels = unique(group))
  design <- sparse.model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  sums <- expr %*% design
  n_cells <- as.numeric(Matrix::colSums(design))
  means <- sweep(as.matrix(sums), 2L, n_cells, "/")
  out <- as.data.table(t(means), keep.rownames = "group_id")
  out[, n_cells := n_cells]
  parts <- tstrsplit(out$group_id, "__", fixed = TRUE)
  if (length(parts) != length(group_columns)) {
    stop("分组 ID 无法按预期拆分。")
  }
  for (j in seq_along(group_columns)) {
    out[, (group_columns[j]) := parts[[j]]]
  }
  setcolorder(out, c(group_columns, "n_cells", candidate_available))
  out
}

score_pseudobulk <- function(dt, eligible, genes = candidate_available) {
  m <- as.matrix(dt[, ..genes])
  mu <- colMeans(m[eligible, , drop = FALSE])
  sig <- apply(m[eligible, , drop = FALSE], 2L, sd)
  sig[!is.finite(sig) | sig == 0] <- 1
  z <- sweep(sweep(m, 2L, mu, "-"), 2L, sig, "/")
  rowMeans(z)
}

# 各 major cell type 的患者-组织 pseudobulk，仅用于细胞定位描述。
major_group <- paste(meta$patient, meta$loc, meta$majorCluster, sep = "__")
major_pseudo <- aggregate_sparse(
  expr_candidate, major_group, c("patient", "loc", "majorCluster")
)
major_keep <- major_pseudo$n_cells >= 20L
major_pseudo[, module_score := score_pseudobulk(major_pseudo, major_keep)]
major_pseudo[, eligible := major_keep]
fwrite(
  major_pseudo,
  file.path(result_dir, "GSE270680_major_celltype_pseudobulk_scores.csv")
)

# 上皮患者-组织 pseudobulk。
epi_idx <- which(meta$majorCluster == "Epithelium" & meta$loc %in% c("N", "T"))
if (!length(epi_idx)) stop("未找到 N/T 来源的作者注释 Epithelium 细胞。")
epi_meta <- meta[epi_idx]
epi_expr <- expr_candidate[, epi_idx, drop = FALSE]
epi_group <- paste(epi_meta$patient, epi_meta$loc, sep = "__")
epi_pseudo <- aggregate_sparse(
  epi_expr, epi_group, c("patient", "loc")
)
epi_pseudo[, tissue := fifelse(loc == "N", "Normal", "Tumor")]
epi_pseudo[, eligible := n_cells >= 20L]
epi_pseudo[, module_score := score_pseudobulk(
  epi_pseudo, eligible, candidate_available
)]

eligible_epi <- epi_pseudo[eligible == TRUE]
complete_patients <- eligible_epi[
  , .(n_tissues = uniqueN(tissue)), by = patient
][n_tissues == 2L, patient]
epi_pairs <- eligible_epi[patient %in% complete_patients]
epi_wide <- dcast(
  epi_pairs, patient ~ tissue,
  value.var = c("module_score", "n_cells")
)
stopifnot(all(c("module_score_Normal", "module_score_Tumor") %in% names(epi_wide)))
epi_wide[, difference := module_score_Tumor - module_score_Normal]

epi_t <- t.test(
  epi_wide$module_score_Tumor,
  epi_wide$module_score_Normal,
  paired = TRUE
)
epi_w <- suppressWarnings(wilcox.test(
  epi_wide$module_score_Tumor,
  epi_wide$module_score_Normal,
  paired = TRUE, exact = FALSE, conf.int = TRUE
))

fwrite(
  epi_pseudo,
  file.path(result_dir, "GSE270680_epithelial_patient_pseudobulk_scores.csv")
)
fwrite(
  epi_wide,
  file.path(result_dir, "GSE270680_epithelial_complete_pair_scores.csv")
)

# 作者定义上皮亚群内的复核；亚群至少需要 20 个细胞/患者/组织。
lineage_group <- paste(
  epi_meta$patient, epi_meta$loc, epi_meta$subCluster, sep = "__"
)
lineage_pseudo <- aggregate_sparse(
  epi_expr, lineage_group, c("patient", "loc", "subCluster")
)
lineage_pseudo[, tissue := fifelse(loc == "N", "Normal", "Tumor")]
lineage_pseudo[, eligible := n_cells >= 20L]
lineage_pseudo[, module_score := score_pseudobulk(
  lineage_pseudo, eligible, candidate_available
)]
lineage_use <- lineage_pseudo[eligible == TRUE]
lineage_use[, `:=`(
  patient = factor(patient),
  tissue = factor(tissue, levels = c("Normal", "Tumor")),
  subCluster = factor(subCluster)
)]

lineage_fit <- lm(
  module_score ~ tissue + patient + subCluster,
  data = lineage_use,
  weights = sqrt(n_cells)
)
lineage_vcov <- sandwich::vcovCL(
  lineage_fit, cluster = ~ patient, type = "HC1"
)
lineage_coefs <- lmtest::coeftest(lineage_fit, vcov. = lineage_vcov)
lineage_tissue <- lineage_coefs["tissueTumor", ]

lineage_tests <- rbindlist(lapply(
  sort(unique(as.character(lineage_use$subCluster))),
  function(lin) {
    d <- copy(lineage_use[as.character(subCluster) == lin])
    counts <- d[, .(n_tissues = uniqueN(tissue)), by = patient]
    paired_ids <- counts[n_tissues == 2L, patient]
    d <- d[patient %in% paired_ids]
    if (uniqueN(d$patient) < 3L) {
      return(data.table(
        subCluster = lin, paired_patients = uniqueN(d$patient),
        mean_difference = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        paired_t_p = NA_real_, wilcoxon_p = NA_real_
      ))
    }
    w <- dcast(d, patient ~ tissue, value.var = "module_score")
    tt <- t.test(w$Tumor, w$Normal, paired = TRUE)
    ww <- suppressWarnings(wilcox.test(
      w$Tumor, w$Normal, paired = TRUE, exact = FALSE
    ))
    data.table(
      subCluster = lin,
      paired_patients = nrow(w),
      mean_difference = mean(w$Tumor - w$Normal),
      ci_low = unname(tt$conf.int[1]),
      ci_high = unname(tt$conf.int[2]),
      paired_t_p = tt$p.value,
      wilcoxon_p = ww$p.value
    )
  }
), fill = TRUE)
lineage_tests[, paired_t_fdr := p.adjust(paired_t_p, method = "BH")]
fwrite(
  lineage_pseudo,
  file.path(result_dir, "GSE270680_epithelial_lineage_pseudobulk_scores.csv")
)
fwrite(
  lineage_tests,
  file.path(result_dir, "GSE270680_epithelial_lineage_paired_tests.csv")
)

# 作者定义上皮亚群组成。患者-组织是独立单位，细胞比例只作捕获组成描述。
lineage_counts <- epi_meta[
  , .(n_cells = .N), by = .(patient, loc, subCluster)
]
epi_totals <- epi_meta[, .(total_epithelial = .N), by = .(patient, loc)]
lineage_counts <- merge(
  lineage_counts, epi_totals, by = c("patient", "loc"), all.x = TRUE
)
lineage_counts[, proportion := n_cells / total_epithelial]
lineage_counts[, tissue := fifelse(loc == "N", "Normal", "Tumor")]
lineage_counts[, eligible := total_epithelial >= 20L]
fwrite(
  lineage_counts,
  file.path(result_dir, "GSE270680_epithelial_lineage_proportions.csv")
)

# 两个独立单细胞队列的标准化配对效应汇总。
paired_effect <- function(diff, cohort) {
  n <- length(diff)
  dz <- mean(diff) / sd(diff)
  correction <- 1 - 3 / (4 * n - 5)
  g <- correction * dz
  var_g <- correction^2 * (1 / n + dz^2 / (2 * n))
  data.table(
    cohort = cohort,
    n_pairs = n,
    raw_mean_difference = mean(diff),
    raw_ci_low = unname(t.test(diff)$conf.int[1]),
    raw_ci_high = unname(t.test(diff)$conf.int[2]),
    hedges_gz = g,
    variance = var_g,
    se = sqrt(var_g)
  )
}

old_pair_file <- file.path(
  project_dir, "results", "PaperValidation", "SingleCell",
  "epithelial_complete_pair_scores.csv"
)
old_pairs <- fread(old_pair_file)
old_wide <- dcast(
  old_pairs, Patient ~ Tissue, value.var = "module_score"
)
old_diff <- old_wide$Tumor - old_wide$Normal

effects <- rbindlist(list(
  paired_effect(old_diff, "GSE206785"),
  paired_effect(epi_wide$difference, "GSE270680")
))
effects[, `:=`(
  ci_low = hedges_gz - qnorm(0.975) * se,
  ci_high = hedges_gz + qnorm(0.975) * se
)]

effects[, `:=`(
  p_value = 2 * pnorm(-abs(hedges_gz / se)),
  analysis_type = "Descriptive cohort-specific standardized effect; no pooling"
)]
fwrite(
  effects,
  file.path(result_dir, "singlecell_cross_cohort_standardized_effects.csv")
)
old_meta_file <- file.path(
  result_dir, "singlecell_cross_cohort_meta_analysis.csv"
)
if (file.exists(old_meta_file)) file.remove(old_meta_file)

# 所有大对象在绘图前释放。
rm(sc, expr_candidate, epi_expr)
invisible(gc())

morandi <- c("Normal" = "#8FA6A1", "Tumor" = "#B98585")
theme_nature <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

major_plot <- copy(major_pseudo[eligible == TRUE & loc %in% c("N", "T")])
major_plot[, tissue := factor(
  fifelse(loc == "N", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)]
p_type <- ggplot(
  major_plot,
  aes(module_score, reorder(majorCluster, module_score, FUN = median),
      colour = tissue)
) +
  geom_boxplot(
    aes(group = interaction(majorCluster, tissue)),
    outlier.shape = NA, width = 0.52, linewidth = 0.4, alpha = 0.15
  ) +
  geom_jitter(height = 0.12, width = 0, size = 1.15, alpha = 0.68) +
  scale_colour_manual(values = morandi) +
  labs(
    x = "73-gene patient–tissue pseudobulk score",
    y = NULL,
    colour = NULL,
    title = "GSE270680 author-defined major cell types"
  ) +
  theme_nature +
  theme(legend.position = "top")

ggsave(
  file.path(result_dir, "GSE270680_major_celltype_localization.pdf"),
  p_type, width = 6.8, height = 4.9, units = "in", device = cairo_pdf
)
ggsave(
  file.path(result_dir, "GSE270680_major_celltype_localization.png"),
  p_type, width = 6.8, height = 4.9, units = "in", dpi = 600
)

pair_plot <- copy(epi_pairs)
pair_plot[, tissue := factor(tissue, levels = c("Normal", "Tumor"))]
p_pair <- ggplot(
  pair_plot,
  aes(tissue, module_score, group = patient)
) +
  geom_line(colour = "grey70", linewidth = 0.35) +
  geom_point(aes(colour = tissue), size = 1.7) +
  scale_colour_manual(values = morandi) +
  labs(
    x = NULL,
    y = "Epithelial pseudobulk module score",
    title = paste0("GSE270680 paired epithelial test (n = ", nrow(epi_wide), ")")
  ) +
  theme_nature +
  theme(legend.position = "none")

forest_dt <- copy(effects)
forest_dt[, cohort := factor(cohort, levels = rev(cohort))]
p_forest <- ggplot(
  forest_dt,
  aes(hedges_gz, cohort)
) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.16, linewidth = 0.5, colour = "#6E7F80"
  ) +
  geom_point(size = 2.1, colour = "#8C5F62", fill = "#8C5F62") +
  labs(
    x = "Standardized paired effect (Hedges' gz)\nTumor minus normal",
    y = NULL,
    title = "Cross-cohort single-cell evidence"
  ) +
  theme_nature

lineage_plot <- rbindlist(list(
  data.table(
    analysis = paste0("All epithelium, paired (n=", nrow(epi_wide), ")"),
    estimate = mean(epi_wide$difference),
    ci_low = unname(epi_t$conf.int[1]),
    ci_high = unname(epi_t$conf.int[2])
  ),
  data.table(
    analysis = "Lineage-adjusted coefficient",
    estimate = unname(lineage_tissue[1]),
    ci_low = unname(lineage_tissue[1] - qnorm(0.975) * lineage_tissue[2]),
    ci_high = unname(lineage_tissue[1] + qnorm(0.975) * lineage_tissue[2])
  ),
  lineage_tests[
    subCluster == "Mucous_MUC5AC" & is.finite(mean_difference),
    .(
      analysis = paste0("Mucous_MUC5AC, paired (n=", paired_patients, ")"),
      estimate = mean_difference,
      ci_low = ci_low,
      ci_high = ci_high
    )
  ]
))
p_lineage <- ggplot(
  lineage_plot,
  aes(estimate, reorder(analysis, estimate))
) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.16, linewidth = 0.45, colour = "#7C8FA3"
  ) +
  geom_point(size = 1.8, colour = "#9B6D70") +
  labs(
    x = "Module-score effect\nTumor minus normal",
    y = NULL,
    title = "Within-epithelium sensitivity analyses"
  ) +
  theme_nature

p_main <- (p_pair | p_forest) / p_lineage +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(result_dir, "GSE270680_singlecell_replication_figure.pdf"),
  p_main, width = 8.2, height = 7.2, units = "in", device = cairo_pdf
)
ggsave(
  file.path(result_dir, "GSE270680_singlecell_replication_figure.png"),
  p_main, width = 8.2, height = 7.2, units = "in", dpi = 600
)
ggsave(
  file.path(result_dir, "singlecell_cross_cohort_forest.pdf"),
  p_forest, width = 6.5, height = 3.6, units = "in", device = cairo_pdf
)
ggsave(
  file.path(result_dir, "singlecell_cross_cohort_forest.png"),
  p_forest, width = 6.5, height = 3.6, units = "in", dpi = 600
)

# 同步更新主稿图5及其逐面板源数据，避免主稿继续引用旧的合并效应图。
main_figure_dir <- file.path(project_dir, "results", "ManuscriptFigures")
main_source_dir <- file.path(main_figure_dir, "SourceData")
dir.create(main_source_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(
  file.path(result_dir, "GSE270680_singlecell_replication_figure.pdf"),
  file.path(main_figure_dir, "Figure5_independent_single_cell_validation.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(result_dir, "GSE270680_singlecell_replication_figure.png"),
  file.path(main_figure_dir, "Figure5_independent_single_cell_validation.png"),
  overwrite = TRUE
)
fwrite(epi_wide, file.path(main_source_dir, "Figure5a_source_data.csv"))
fwrite(effects, file.path(main_source_dir, "Figure5b_source_data.csv"))
fwrite(lineage_plot, file.path(main_source_dir, "Figure5c_source_data.csv"))

summary_table <- data.table(
  metric = c(
    "cells_total", "patients_total", "libraries_total",
    "normal_tumor_patients_before_cell_threshold",
    "epithelial_cells_normal", "epithelial_cells_tumor",
    "eligible_complete_epithelial_pairs", "candidate_genes_available",
    "mean_tumor_minus_normal", "paired_t_ci_low", "paired_t_ci_high",
    "paired_t_pvalue", "paired_wilcoxon_pvalue",
    "lineage_adjusted_tumor_coefficient",
    "lineage_adjusted_cluster_robust_se",
    "lineage_adjusted_cluster_robust_pvalue"
  ),
  value = c(
    audit[metric == "cells", value],
    audit[metric == "patients", value],
    audit[metric == "libraries", value],
    meta[loc %in% c("N", "T"), .(n_loc = uniqueN(loc)), by = patient][
      n_loc == 2L, .N
    ],
    sum(epi_meta$loc == "N"),
    sum(epi_meta$loc == "T"),
    nrow(epi_wide),
    length(candidate_available),
    mean(epi_wide$difference),
    unname(epi_t$conf.int[1]),
    unname(epi_t$conf.int[2]),
    epi_t$p.value,
    epi_w$p.value,
    lineage_tissue[1],
    lineage_tissue[2],
    lineage_tissue[4]
  )
)
fwrite(
  summary_table,
  file.path(result_dir, "GSE270680_singlecell_replication_summary.csv")
)

fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "fg")
effect206 <- effects[cohort == "GSE206785"]
effect270 <- effects[cohort == "GSE270680"]
report <- c(
  "# GSE270680 第二个独立单细胞队列交叉验证报告",
  "",
  paste0("- 分析日期：", format(Sys.Date(), "%Y-%m-%d")),
  "- 物种：Homo sapiens（人）。",
  "- 数据来源：作者公开的 sc.rds（Mendeley Data DOI 10.17632/559mchb37p.1）。",
  paste0("- 原始对象：", format(audit[metric == "cells", value], big.mark = ","),
         " 个细胞，", audit[metric == "patients", value], " 例供者，",
         audit[metric == "libraries", value], " 个样本库。"),
  paste0("- 73 基因覆盖：", length(candidate_available), "/73。"),
  "- 使用作者 majorCluster/subCluster；未重新聚类，也未用候选基因定义细胞类型。",
  "",
  "## 患者配对上皮 pseudobulk",
  "",
  "- 纳入标准：同一患者的 N 和 T 组织均至少包含 20 个作者注释上皮细胞。",
  paste0("- 合格完整配对：", nrow(epi_wide), " 对。"),
  paste0(
    "- 平均肿瘤－邻癌差值：", fmt(mean(epi_wide$difference)),
    "（95% CI ", fmt(epi_t$conf.int[1]), " 至 ",
    fmt(epi_t$conf.int[2]), "；配对 t 检验 P = ",
    fmt(epi_t$p.value), "；Wilcoxon P = ", fmt(epi_w$p.value), "）。"
  ),
  paste0(
    "- 作者定义上皮亚群校正后的肿瘤系数：", fmt(lineage_tissue[1]),
    "（患者聚类稳健 SE = ", fmt(lineage_tissue[2]),
    "，P = ", fmt(lineage_tissue[4]), "）。"
  ),
  "",
  "## 两个单细胞队列的交叉验证",
  "",
  paste0(
    "- GSE206785 标准化配对效应 Hedges' gz = ",
    fmt(effect206$hedges_gz), "（95% CI ",
    fmt(effect206$ci_low), " 至 ", fmt(effect206$ci_high), "）。"
  ),
  paste0(
    "- GSE270680 标准化配对效应 Hedges' gz = ",
    fmt(effect270$hedges_gz), "（95% CI ",
    fmt(effect270$ci_low), " 至 ", fmt(effect270$ci_high), "）。"
  ),
  "- 仅有两个队列，因此不计算或解释固定效应、随机效应、tau²或I²；两项队列估计仅作描述性并列。",
  "- 未观察到显著差异不能解释为等效；置信区间决定仍不能排除的效应范围。",
  "",
  "## 解释边界",
  "",
  "该分析检验的是作者注释上皮细胞的患者级 pseudobulk，不把肿瘤来源上皮细胞自动等同为恶性细胞。",
  "结果可用于判断 bulk 模块是否在上皮区室内稳定下降，并与组织组成解释交叉核对；",
  "它不能替代基因扰动实验、蛋白验证或因果推断。"
)
writeLines(
  report,
  file.path(result_dir, "GSE270680_singlecell_replication_report.md"),
  useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "14_sessionInfo.txt"),
  useBytes = TRUE
)
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
