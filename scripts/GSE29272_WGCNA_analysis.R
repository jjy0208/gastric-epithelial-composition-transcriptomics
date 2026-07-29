#!/usr/bin/env Rscript

# ==============================================================================
# GSE29272 WGCNA 加权基因共表达网络分析
# 输入：清洗后的完整表达矩阵、样本分组信息
# 输出：results/WGCNA 下的统计表、PDF/PNG 图片、中文报告和运行日志
# 设计原则：
#   1) 使用 MAD 最大的前 5000 个基因建网；
#   2) 样本离群值采用预先定义的连接度 Z.k < -2.5 判定；
#   3) 软阈值选择为有符号无标度拟合 R² 首次达到 0.85；
#   4) 使用 signed bicor 网络，minModuleSize=30，mergeCutHeight=0.25；
#   5) 关键模块定义为非灰色模块中与疾病状态绝对相关性最高者，
#      若并列则优先 P 值更小者。
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) != 1L) stop("请使用Rscript运行本脚本。")
  normalizePath(sub("^--file=", "", f), winslash = "/", mustWork = TRUE)
}

project_dir <- normalizePath(file.path(dirname(get_script_path()), ".."),
                             winslash = "/", mustWork = TRUE)
expr_file <- file.path(project_dir, "clean_data", "GSE29272_clean_expression_matrix.csv")
meta_file <- file.path(project_dir, "clean_data", "GSE29272_sample_info.csv")
deg_file <- file.path(project_dir, "results", "DEG", "GSE29272_DEG_all.csv")
out_dir <- file.path(project_dir, "results", "WGCNA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "WGCNA", "data.table", "matrixStats", "ggplot2", "patchwork", "ggrepel", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("缺少必要 R 包：", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

allowWGCNAThreads(nThreads = max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L)))
set.seed(20260729)

log_file <- file.path(out_dir, "GSE29272_WGCNA_run.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("GSE29272 WGCNA analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Project:", project_dir, "\n\n")

# ------------------------------------------------------------------------------
# 通用工具函数
# ------------------------------------------------------------------------------
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
                                formatC(x, format = "f", digits = 3)))
}

write_csv_utf8 <- function(x, filename) {
  data.table::fwrite(x, filename, bom = TRUE, na = "NA")
}

save_gg_dual <- function(plot_obj, stem, width, height) {
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot_obj, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot_obj, width = width, height = height,
    units = "in", dpi = 600, device = ragg::agg_png, bg = "white"
  )
}

open_pdf <- function(filename, width, height) {
  grDevices::cairo_pdf(filename, width = width, height = height, family = "sans")
}

open_png <- function(filename, width, height) {
  ragg::agg_png(
    filename, width = width, height = height, units = "in",
    res = 600, background = "white"
  )
}

theme_nature <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.title = element_text(color = "#2B2B2B", size = base_size),
      axis.text = element_text(color = "#2B2B2B", size = base_size - 1),
      plot.title = element_text(face = "bold", color = "#222222", size = base_size + 1),
      plot.subtitle = element_text(color = "#555555", size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# 低饱和度莫兰迪色：蓝灰—米灰—陶土红
morandi_blue <- "#7089A6"
morandi_beige <- "#D8D2C8"
morandi_red <- "#B36F63"
morandi_dark <- "#4A4A4A"
module_palette <- c(
  "#78909C", "#A1887F", "#7E8F78", "#9A86A4", "#B08D73",
  "#6F8F9D", "#9C7C83", "#8D9A73", "#7F86A3", "#A98972",
  "#698B78", "#978A75", "#8796A5", "#A77E72", "#798C9C",
  "#93869B", "#8C9975", "#B09A83", "#738F8A", "#A4808E",
  "#8691A7", "#9B8F78", "#7E998F", "#A18470", "#7D8298",
  "#8E827A", "#7896A1", "#9A838A", "#819077", "#AA927A"
)

# ------------------------------------------------------------------------------
# 1. 读取数据并核对样本
# ------------------------------------------------------------------------------
if (!file.exists(expr_file)) stop("表达矩阵不存在：", expr_file)
if (!file.exists(meta_file)) stop("样本信息表不存在：", meta_file)

expr_dt <- data.table::fread(expr_file, check.names = FALSE)
if (ncol(expr_dt) < 3) stop("表达矩阵列数异常。")
gene_col <- names(expr_dt)[1]
gene_ids <- trimws(as.character(expr_dt[[gene_col]]))
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- gene_ids

meta <- data.table::fread(meta_file, check.names = FALSE)
required_meta <- c("sample", "disease_status")
if (!all(required_meta %in% names(meta))) {
  stop("样本信息表缺少列：", paste(setdiff(required_meta, names(meta)), collapse = ", "))
}
meta$sample <- trimws(as.character(meta$sample))
meta$disease_status <- trimws(as.character(meta$disease_status))

if (anyDuplicated(rownames(expr_mat))) stop("表达矩阵存在重复基因名。")
if (anyDuplicated(colnames(expr_mat))) stop("表达矩阵存在重复样本名。")
if (anyDuplicated(meta$sample)) stop("样本信息表存在重复样本名。")
if (any(!is.finite(expr_mat))) stop("表达矩阵包含 NA、Inf 或非数值。")

missing_in_meta <- setdiff(colnames(expr_mat), meta$sample)
missing_in_expr <- setdiff(meta$sample, colnames(expr_mat))
if (length(missing_in_meta) > 0 || length(missing_in_expr) > 0) {
  stop(
    "表达矩阵与样本表不完全一致。表达矩阵独有：",
    paste(missing_in_meta, collapse = ";"),
    "；样本表独有：", paste(missing_in_expr, collapse = ";")
  )
}
meta <- meta[match(colnames(expr_mat), meta$sample), ]

if (!all(meta$disease_status %in% c("Tumor", "Normal"))) {
  stop("disease_status 仅允许 Tumor/Normal；实际值：",
       paste(unique(meta$disease_status), collapse = ", "))
}
meta$trait_disease <- ifelse(meta$disease_status == "Tumor", 1, 0)

n_genes_raw <- nrow(expr_mat)
n_samples_raw <- ncol(expr_mat)
group_counts <- table(meta$disease_status)
sample_size_ok <- n_samples_raw >= 15

cat("Input dimension:", n_genes_raw, "genes x", n_samples_raw, "samples\n")
cat("Groups:", paste(names(group_counts), group_counts, sep = "=", collapse = "; "), "\n")
cat("Sample size >= 15:", sample_size_ok, "\n")
cat("Sample names aligned exactly: TRUE\n\n")

# ------------------------------------------------------------------------------
# 2. MAD 筛选前 5000 个高变基因
# ------------------------------------------------------------------------------
gene_mad <- matrixStats::rowMads(expr_mat, na.rm = TRUE)
mad_table <- data.frame(
  gene = rownames(expr_mat),
  MAD = gene_mad,
  stringsAsFactors = FALSE
)
mad_table <- mad_table[order(-mad_table$MAD, mad_table$gene), ]
n_select <- min(5000L, nrow(mad_table))
selected_genes <- mad_table$gene[seq_len(n_select)]
datExpr0 <- t(expr_mat[selected_genes, , drop = FALSE])

gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
selected_mad <- mad_table[match(colnames(datExpr0), mad_table$gene), ]
write_csv_utf8(selected_mad, file.path(out_dir, "GSE29272_WGCNA_top5000_MAD_genes.csv"))

cat("Selected by MAD:", n_select, "\n")
cat("After goodSamplesGenes:", nrow(datExpr0), "samples x", ncol(datExpr0), "genes\n\n")

# ------------------------------------------------------------------------------
# 3. 样本聚类与离群样本判定
#    预设客观规则：稳健相关网络中的样本连接度 Z.k < -2.5。
# ------------------------------------------------------------------------------
sample_cor <- WGCNA::bicor(
  t(datExpr0), use = "pairwise.complete.obs", maxPOutliers = 0.1
)
diag(sample_cor) <- 1
sample_adj <- ((1 + sample_cor) / 2)^2
sample_k <- rowSums(sample_adj) - 1
sample_zk <- as.numeric(scale(sample_k))
names(sample_zk) <- rownames(datExpr0)
outlier_cutoff <- -2.5
outlier_samples <- names(sample_zk)[sample_zk < outlier_cutoff]

sample_tree <- stats::hclust(stats::as.dist(1 - sample_cor), method = "average")
qc_table <- data.frame(
  sample = rownames(datExpr0),
  disease_status = meta$disease_status[match(rownames(datExpr0), meta$sample)],
  pair_id = if ("pair_id" %in% names(meta)) meta$pair_id[match(rownames(datExpr0), meta$sample)] else NA,
  connectivity = sample_k[rownames(datExpr0)],
  Z_k = sample_zk[rownames(datExpr0)],
  outlier_rule = paste0("Z_k < ", outlier_cutoff),
  removed = rownames(datExpr0) %in% outlier_samples,
  stringsAsFactors = FALSE
)
write_csv_utf8(qc_table, file.path(out_dir, "GSE29272_WGCNA_sample_QC.csv"))

sample_status_colors <- ifelse(qc_table$removed, morandi_red, morandi_blue)
group_colors <- ifelse(qc_table$disease_status == "Tumor", "#A77E72", "#78909C")
names(sample_status_colors) <- qc_table$sample
names(group_colors) <- qc_table$sample
dendro_colors <- cbind(
  "QC status" = sample_status_colors[sample_tree$labels],
  "Disease status" = group_colors[sample_tree$labels]
)

plot_sample_tree <- function() {
  par(mar = c(2, 4, 3, 1), family = "sans")
  WGCNA::plotDendroAndColors(
    sample_tree, dendro_colors,
    groupLabels = c("QC status", "Disease status"),
    main = "GSE29272 sample clustering and outlier assessment",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
  )
  mtext(
    paste0("Outlier rule: sample connectivity Z.k < ", outlier_cutoff,
           "; removed n = ", length(outlier_samples)),
    side = 1, line = 0.2, cex = 0.7, col = morandi_dark
  )
}
open_pdf(file.path(out_dir, "GSE29272_WGCNA_sample_clustering.pdf"), 8.2, 5.3)
plot_sample_tree()
dev.off()
open_png(file.path(out_dir, "GSE29272_WGCNA_sample_clustering.png"), 8.2, 5.3)
plot_sample_tree()
dev.off()

keep_samples <- setdiff(rownames(datExpr0), outlier_samples)
datExpr <- datExpr0[keep_samples, , drop = FALSE]
meta_net <- meta[match(rownames(datExpr), meta$sample), ]
if (nrow(datExpr) < 15) stop("去除离群样本后少于 15 个样本，不适合继续 WGCNA。")

cat("Outlier cutoff:", outlier_cutoff, "\n")
cat("Outliers removed:", ifelse(length(outlier_samples) == 0, "None",
                                paste(outlier_samples, collapse = ", ")), "\n")
cat("Network input:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n\n")

# ------------------------------------------------------------------------------
# 4. 软阈值选择：signed R² 首次达到 0.85；若未达到则透明记录回退规则
# ------------------------------------------------------------------------------
powers <- c(1:10, seq(12, 30, by = 2))
sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = powers,
  RsquaredCut = 0.85,
  corFnc = "bicor",
  corOptions = list(use = "p", maxPOutliers = 0.1),
  networkType = "signed",
  verbose = 5
)
fit_indices <- as.data.frame(sft$fitIndices)
fit_indices$signed_R2 <- -sign(fit_indices$slope) * fit_indices$SFT.R.sq
eligible <- which(fit_indices$signed_R2 >= 0.85)
if (length(eligible) > 0) {
  soft_power <- fit_indices$Power[min(eligible)]
  soft_rule <- "有符号无标度拟合R²首次达到0.85"
  soft_reached <- TRUE
} else {
  best_idx <- which.max(fit_indices$signed_R2)
  soft_power <- fit_indices$Power[best_idx]
  soft_rule <- "候选范围内未达到0.85，回退为有符号R²最大值对应的β"
  soft_reached <- FALSE
}
write_csv_utf8(fit_indices, file.path(out_dir, "GSE29272_WGCNA_soft_threshold_metrics.csv"))

sft_plot_data <- fit_indices
p_r2 <- ggplot(sft_plot_data, aes(x = Power, y = signed_R2)) +
  geom_hline(yintercept = 0.85, linetype = 2, linewidth = 0.45, color = morandi_red) +
  geom_line(linewidth = 0.55, color = morandi_blue) +
  geom_point(size = 1.8, color = morandi_blue) +
  geom_text(aes(label = Power), vjust = -0.7, size = 2.4, color = morandi_dark) +
  geom_point(
    data = subset(sft_plot_data, Power == soft_power),
    size = 3.2, shape = 21, fill = morandi_red, color = "white", stroke = 0.6
  ) +
  scale_x_continuous(breaks = powers) +
  coord_cartesian(ylim = c(min(-0.1, min(sft_plot_data$signed_R2, na.rm = TRUE)),
                           max(0.9, max(sft_plot_data$signed_R2, na.rm = TRUE) + 0.08))) +
  labs(
    title = "A  Scale-free topology fit",
    subtitle = paste0("Selected β = ", soft_power),
    x = "Soft-threshold power (β)", y = "Signed scale-free fit, R²"
  ) +
  theme_nature(9)

p_k <- ggplot(sft_plot_data, aes(x = Power, y = mean.k.)) +
  geom_line(linewidth = 0.55, color = "#7E8F78") +
  geom_point(size = 1.8, color = "#7E8F78") +
  geom_text(aes(label = Power), vjust = -0.7, size = 2.4, color = morandi_dark) +
  geom_point(
    data = subset(sft_plot_data, Power == soft_power),
    size = 3.2, shape = 21, fill = morandi_red, color = "white", stroke = 0.6
  ) +
  scale_x_continuous(breaks = powers) +
  labs(
    title = "B  Mean connectivity",
    x = "Soft-threshold power (β)", y = "Mean connectivity"
  ) +
  theme_nature(9)

p_soft <- p_r2 + p_k + patchwork::plot_layout(ncol = 2)
save_gg_dual(
  p_soft, file.path(out_dir, "GSE29272_WGCNA_soft_threshold"),
  width = 8.4, height = 3.6
)
cat("Soft power:", soft_power, "\n")
cat("Criterion reached:", soft_reached, "\n")
cat("Selection rule:", soft_rule, "\n\n")

# ------------------------------------------------------------------------------
# 5. 一步法构建 signed 共表达网络并识别模块
# ------------------------------------------------------------------------------
net <- WGCNA::blockwiseModules(
  datExpr,
  maxBlockSize = 5000,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.1,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  randomSeed = 20260729,
  verbose = 3
)

numeric_labels <- net$colors
non_grey_labels <- sort(setdiff(unique(numeric_labels), 0))
if (length(non_grey_labels) > length(module_palette)) {
  module_palette <- grDevices::hcl.colors(length(non_grey_labels), palette = "Muted")
}
module_name_map <- setNames(sprintf("M%02d", seq_along(non_grey_labels)), non_grey_labels)
module_color_map <- setNames(module_palette[seq_along(non_grey_labels)], non_grey_labels)
module_ids <- ifelse(numeric_labels == 0, "grey", module_name_map[as.character(numeric_labels)])
module_colors <- ifelse(numeric_labels == 0, "#C8C8C8",
                        module_color_map[as.character(numeric_labels)])
names(module_ids) <- colnames(datExpr)
names(module_colors) <- colnames(datExpr)

module_assign <- data.frame(
  gene = colnames(datExpr),
  module = unname(module_ids[colnames(datExpr)]),
  module_color = unname(module_colors[colnames(datExpr)]),
  MAD = selected_mad$MAD[match(colnames(datExpr), selected_mad$gene)],
  stringsAsFactors = FALSE
)
module_assign <- module_assign[order(module_assign$module, -module_assign$MAD), ]
module_sizes <- as.data.frame(table(module_assign$module), stringsAsFactors = FALSE)
names(module_sizes) <- c("module", "gene_count")
module_sizes$module_color <- vapply(
  module_sizes$module,
  function(m) unique(module_assign$module_color[module_assign$module == m])[1],
  character(1)
)
module_sizes <- module_sizes[order(module_sizes$module == "grey", -module_sizes$gene_count), ]
write_csv_utf8(module_assign, file.path(out_dir, "GSE29272_WGCNA_module_assignment.csv"))
write_csv_utf8(module_sizes, file.path(out_dir, "GSE29272_WGCNA_module_sizes.csv"))

plot_module_tree <- function() {
  par(mar = c(3, 4, 3, 1), family = "sans")
  WGCNA::plotDendroAndColors(
    net$dendrograms[[1]],
    module_colors[net$blockGenes[[1]]],
    groupLabels = "Merged modules",
    main = "GSE29272 gene dendrogram and co-expression modules",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
  )
  mtext(
    "signed bicor network; minModuleSize = 30; mergeCutHeight = 0.25",
    side = 1, line = 0.3, cex = 0.7, col = morandi_dark
  )
}
open_pdf(file.path(out_dir, "GSE29272_WGCNA_module_dendrogram.pdf"), 9.0, 5.0)
plot_module_tree()
dev.off()
open_png(file.path(out_dir, "GSE29272_WGCNA_module_dendrogram.png"), 9.0, 5.0)
plot_module_tree()
dev.off()

# ------------------------------------------------------------------------------
# 6. 模块—性状相关分析
# ------------------------------------------------------------------------------
MEs <- WGCNA::moduleEigengenes(datExpr, colors = module_ids)$eigengenes
MEs <- WGCNA::orderMEs(MEs)
trait_df <- data.frame(Disease = meta_net$trait_disease)
rownames(trait_df) <- meta_net$sample
# 二分类 0/1 性状的 MAD 为 0，不适合直接用 bicor；Pearson 在此等价于点二列相关。
cor_res <- WGCNA::corAndPvalue(
  MEs, trait_df, use = "pairwise.complete.obs"
)
module_trait_cor <- cor_res$cor[, 1]
module_trait_p <- cor_res$p[, 1]
module_names <- sub("^ME", "", names(module_trait_cor))

module_trait <- data.frame(
  module = module_names,
  correlation = as.numeric(module_trait_cor),
  p_value = as.numeric(module_trait_p),
  gene_count = module_sizes$gene_count[match(module_names, module_sizes$module)],
  module_color = module_sizes$module_color[match(module_names, module_sizes$module)],
  stringsAsFactors = FALSE
)
module_trait$FDR <- p.adjust(module_trait$p_value, method = "BH")
module_trait <- module_trait[order(-abs(module_trait$correlation), module_trait$p_value), ]

key_candidates <- subset(module_trait, module != "grey")
if (nrow(key_candidates) == 0) stop("未识别到非灰色模块，无法定义关键模块。")
key_module <- key_candidates$module[1]
key_cor <- key_candidates$correlation[1]
key_p <- key_candidates$p_value[1]
key_fdr <- key_candidates$FDR[1]
key_color <- key_candidates$module_color[1]
write_csv_utf8(module_trait, file.path(out_dir, "GSE29272_WGCNA_module_trait_correlations.csv"))

me_export <- data.frame(sample = rownames(MEs), MEs, check.names = FALSE)
me_export <- merge(
  meta_net[, intersect(c("sample", "disease_status", "group", "pair_id"), names(meta_net)),
           with = FALSE],
  me_export, by = "sample", sort = FALSE
)
write_csv_utf8(me_export, file.path(out_dir, "GSE29272_WGCNA_module_eigengenes.csv"))

heat_data <- module_trait
heat_data$module <- factor(
  heat_data$module,
  levels = rev(module_trait$module[order(module_trait$correlation)])
)
heat_data$trait <- "Disease vs control"
heat_data$label <- paste0(
  "r = ", sprintf("%.2f", heat_data$correlation), "\nP = ", fmt_p(heat_data$p_value)
)

p_heat <- ggplot(heat_data, aes(x = trait, y = module, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 2.7, lineheight = 0.95, color = "#252525") +
  scale_fill_gradient2(
    low = morandi_blue, mid = morandi_beige, high = morandi_red,
    midpoint = 0, limits = c(-1, 1), name = "Pearson r"
  ) +
  labs(
    title = "Module–trait relationships",
    subtitle = "Disease coded as 1; control coded as 0",
    x = NULL, y = "Module"
  ) +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "#2B2B2B"),
    axis.text.y = element_text(color = "#2B2B2B"),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(6, 8, 6, 6)
  )
heat_height <- max(4.0, min(9.0, 1.6 + 0.42 * nrow(heat_data)))
save_gg_dual(
  p_heat, file.path(out_dir, "GSE29272_WGCNA_module_trait_heatmap"),
  width = 5.1, height = heat_height
)

cat("Non-grey modules:", length(non_grey_labels), "\n")
cat("Key module:", key_module, "\n")
cat("Key module correlation:", key_cor, "P:", key_p, "FDR:", key_fdr, "\n\n")

# ------------------------------------------------------------------------------
# 7. 关键模块：MM、GS、候选枢纽基因及配对敏感性检验
# ------------------------------------------------------------------------------
key_genes <- names(module_ids)[module_ids == key_module]
key_me_name <- paste0("ME", key_module)
key_me <- MEs[[key_me_name]]

mm_res <- WGCNA::bicorAndPvalue(
  datExpr[, key_genes, drop = FALSE], key_me,
  use = "pairwise.complete.obs", maxPOutliers = 0.1
)
gs_res <- WGCNA::corAndPvalue(
  datExpr[, key_genes, drop = FALSE], meta_net$trait_disease,
  use = "pairwise.complete.obs"
)
hub_genes <- data.frame(
  gene = key_genes,
  module = key_module,
  module_color = key_color,
  MAD = selected_mad$MAD[match(key_genes, selected_mad$gene)],
  MM = as.numeric(mm_res$bicor[, 1]),
  MM_pvalue = as.numeric(mm_res$p[, 1]),
  GS = as.numeric(gs_res$cor[, 1]),
  GS_pvalue = as.numeric(gs_res$p[, 1]),
  stringsAsFactors = FALSE
)
hub_genes$hub_candidate <- (
  abs(hub_genes$MM) >= 0.80 &
    abs(hub_genes$GS) >= 0.20 &
    hub_genes$MM_pvalue < 0.05 &
    hub_genes$GS_pvalue < 0.05
)

if (file.exists(deg_file)) {
  deg <- data.table::fread(deg_file, check.names = FALSE)
  deg_gene_col <- intersect(c("gene", "Gene", "symbol", "GeneSymbol"), names(deg))[1]
  if (!is.na(deg_gene_col)) {
    keep_deg <- intersect(
      c(deg_gene_col, "log2FoldChange", "P.Value", "PValue", "pvalue",
        "adj.P.Val", "padj", "FDR"),
      names(deg)
    )
    deg_small <- as.data.frame(deg[, ..keep_deg])
    names(deg_small)[names(deg_small) == deg_gene_col] <- "gene"
    hub_genes <- merge(hub_genes, deg_small, by = "gene", all.x = TRUE, sort = FALSE)
  }
}
hub_genes <- hub_genes[order(-abs(hub_genes$MM), -abs(hub_genes$GS), hub_genes$gene), ]
write_csv_utf8(
  hub_genes,
  file.path(out_dir, "GSE29272_WGCNA_hubmodule_genes.csv")
)

label_n <- min(10L, nrow(hub_genes))
label_genes <- hub_genes[order(-(abs(hub_genes$MM) * abs(hub_genes$GS))), ][seq_len(label_n), ]
p_scatter <- ggplot(hub_genes, aes(x = MM, y = GS)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "#BEBEBE") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "#BEBEBE") +
  geom_point(
    aes(fill = hub_candidate), shape = 21, size = 2.1,
    color = "white", stroke = 0.25, alpha = 0.88
  ) +
  ggrepel::geom_text_repel(
    data = label_genes, aes(label = gene),
    size = 2.6, max.overlaps = Inf, box.padding = 0.3,
    point.padding = 0.2, min.segment.length = 0,
    segment.color = "#888888", color = "#303030", seed = 20260729
  ) +
  scale_fill_manual(
    values = c("FALSE" = key_color, "TRUE" = morandi_red),
    labels = c("FALSE" = "Other module genes", "TRUE" = "Hub candidates"),
    name = NULL
  ) +
  labs(
    title = paste0("Key module ", key_module, ": module membership vs gene significance"),
    subtitle = paste0("Module–disease Pearson r = ", sprintf("%.2f", key_cor),
                      "; P = ", fmt_p(key_p)),
    x = paste0("Module membership (MM, ", key_module, ")"),
    y = "Gene significance (GS, disease status)"
  ) +
  theme_nature(9) +
  theme(legend.position = "top")
save_gg_dual(
  p_scatter, file.path(out_dir, "GSE29272_WGCNA_keymodule_GS_MM"),
  width = 6.2, height = 5.1
)

# 配对敏感性检验：仅用于确认关键模块 ME 在同一患者肿瘤与正常间的方向。
paired_test_available <- FALSE
paired_n <- NA_integer_
paired_mean_delta <- NA_real_
paired_p <- NA_real_
if ("pair_id" %in% names(meta_net)) {
  paired_df <- data.frame(
    sample = meta_net$sample,
    pair_id = meta_net$pair_id,
    disease_status = meta_net$disease_status,
    ME = key_me,
    stringsAsFactors = FALSE
  )
  tumor_me <- paired_df[paired_df$disease_status == "Tumor", c("pair_id", "ME")]
  normal_me <- paired_df[paired_df$disease_status == "Normal", c("pair_id", "ME")]
  names(tumor_me)[2] <- "Tumor"
  names(normal_me)[2] <- "Normal"
  paired_wide <- merge(tumor_me, normal_me, by = "pair_id")
  if (nrow(paired_wide) >= 3) {
    paired_test <- stats::t.test(paired_wide$Tumor, paired_wide$Normal, paired = TRUE)
    paired_test_available <- TRUE
    paired_n <- nrow(paired_wide)
    paired_mean_delta <- mean(paired_wide$Tumor - paired_wide$Normal)
    paired_p <- paired_test$p.value
    write_csv_utf8(
      paired_wide,
      file.path(out_dir, "GSE29272_WGCNA_keymodule_paired_ME.csv")
    )
  }
}

# ------------------------------------------------------------------------------
# 8. 自动生成可审计的中文报告
# ------------------------------------------------------------------------------
top_hub <- head(hub_genes$gene, 15)
hub_candidate_n <- sum(hub_genes$hub_candidate, na.rm = TRUE)
grey_n <- if ("grey" %in% module_sizes$module) {
  module_sizes$gene_count[module_sizes$module == "grey"]
} else {
  0L
}
key_gene_n <- nrow(hub_genes)

# 基于关键模块基因构成给出“线索级”解释，不把关键词命中冒充正式富集分析。
marker_sets <- list(
  "细胞外基质/基质重塑" = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2",
                            "FN1", "SPP1", "THBS1", "MMP2", "MMP9", "FAP", "ACTA2"),
  "细胞周期/增殖" = c("MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "UBE2C", "BUB1"),
  "胃上皮分化与分泌功能" = c("ATP4A", "ATP4B", "GIF", "GKN1", "GKN2",
                              "TFF1", "TFF2", "MUC5AC", "PGC"),
  "免疫炎症" = c("PTPRC", "CD3D", "CD3E", "CD68", "HLA-DRA", "CXCL9",
                  "CXCL10", "CCL5", "LST1"),
  "金属离子/应激反应" = c("MT1A", "MT1E", "MT1F", "MT1G", "MT1H", "MT1M", "MT2A")
)
marker_hits <- lapply(marker_sets, function(x) intersect(x, hub_genes$gene))
marker_text <- vapply(
  names(marker_hits),
  function(nm) {
    hits <- marker_hits[[nm]]
    if (length(hits) > 0) paste0(nm, "（", paste(hits, collapse = "、"), "）") else NA_character_
  },
  character(1)
)
marker_text <- marker_text[!is.na(marker_text)]
if (length(marker_text) == 0) {
  bio_hint <- paste0(
    "预设标志基因集合未出现明确集中命中。关键模块的高 MM/GS 基因包括：",
    paste(head(top_hub, 10), collapse = "、"),
    "。这些仅能作为后续 GO/KEGG 或细胞组成分析的候选线索，不能据此直接命名通路。"
  )
} else {
  bio_hint <- paste0(
    "关键模块包含与", paste(marker_text, collapse = "；"), "相关的标志基因。",
    "这提示该模块可能反映相应的组织过程或细胞组成差异；本结论为基因构成线索，",
    "不是独立富集检验，仍需结合正式功能富集和外部队列验证。"
  )
}

outlier_sentence <- if (length(outlier_samples) == 0) {
  paste0("按预设规则 Z.k < ", outlier_cutoff, " 未检出明显离群样本，因此 268 个样本全部保留。")
} else {
  paste0(
    "按预设规则 Z.k < ", outlier_cutoff, " 剔除 ", length(outlier_samples),
    " 个明显离群样本：", paste(outlier_samples, collapse = "、"),
    "；最终 ", nrow(datExpr), " 个样本进入网络分析。"
  )
}
soft_sentence <- if (soft_reached) {
  paste0("有符号无标度拟合 R² 首次达到 0.85 的候选功率为 β=", soft_power, "，据此建网。")
} else {
  paste0(
    "在预设候选功率 1–30 中，有符号无标度拟合 R² 未达到 0.85；",
    "按预先声明的回退规则选择 R² 最大时的 β=", soft_power,
    "。因此网络近似无标度的证据有限，解释时需保守。"
  )
}
paired_sentence <- if (paired_test_available) {
  paste0(
    "考虑到该队列为配对设计，另对关键模块 ME 做了配对敏感性检验：完整配对 n=",
    paired_n, "，肿瘤−正常平均 ME 差=", sprintf("%.3f", paired_mean_delta),
    "，配对 t 检验 P=", fmt_p(paired_p), "。该检验用于验证方向，不替代 WGCNA 模块相关分析。"
  )
} else {
  "未能形成足够完整配对，故未进行关键模块 ME 的配对敏感性检验。"
}

report_lines <- c(
  "# GSE29272 WGCNA 分析报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- R 版本：", R.version.string),
  paste0("- WGCNA 版本：", as.character(packageVersion("WGCNA"))),
  "",
  "## 1. 数据与样本量检查",
  "",
  paste0(
    "输入表达矩阵包含 ", n_genes_raw, " 个基因和 ", n_samples_raw,
    " 个样本；样本表包含 ", unname(group_counts["Tumor"]), " 个疾病样本和 ",
    unname(group_counts["Normal"]), " 个对照样本。表达矩阵列名与样本信息表的 sample 列完全一致。"
  ),
  paste0(
    "样本量为 ", n_samples_raw,
    "，高于常用的 WGCNA 建议下限 15，因此从样本量角度满足网络分析要求。"
  ),
  "表达值来自已清洗、已注释的 GPL96 芯片矩阵；本步骤不重复归一化。",
  "",
  "## 2. 基因筛选",
  "",
  paste0(
    "对全部 ", n_genes_raw, " 个基因逐行计算 MAD，并按 MAD 从大到小选择前 ",
    ncol(datExpr0), " 个基因。MAD 对极端值较稳健，适合保留跨样本变化最明显、",
    "对共表达结构信息量较高的基因，同时控制网络计算规模。"
  ),
  "",
  "## 3. 样本聚类与离群值",
  "",
  outlier_sentence,
  "样本聚类图中的 QC 色条用于标记保留/剔除状态，疾病状态色条仅用于辅助观察，不参与离群值判定。",
  "",
  "## 4. 软阈值选择",
  "",
  soft_sentence,
  "相关性采用 biweight midcorrelation（bicor），网络类型为 signed；软阈值图同时展示有符号无标度拟合 R² 与平均连接度。",
  "",
  "## 5. 网络构建与模块识别",
  "",
  paste0(
    "使用 blockwiseModules 一步法建立 signed 网络，TOMType=signed，",
    "minModuleSize=30，mergeCutHeight=0.25，deepSplit=2。共识别 ",
    length(non_grey_labels), " 个非灰色共表达模块；另有 ", grey_n,
    " 个基因未能归入稳定模块而标记为 grey。"
  ),
  "",
  "## 6. 模块—疾病性状关联",
  "",
  paste0(
    "疾病状态编码为肿瘤=1、正常=0。非灰色模块中，关键模块为 **", key_module,
    "**（", key_gene_n, " 个基因，图示颜色 ", key_color, "），",
    "其模块特征基因与疾病状态的 Pearson 相关（点二列相关）为 ", sprintf("%.3f", key_cor),
    "，P=", fmt_p(key_p), "，模块层面 BH-FDR=", fmt_p(key_fdr), "。"
  ),
  if (key_cor > 0) {
    "正相关表示该模块整体表达模式在疾病组相对更高。"
  } else {
    "负相关表示该模块整体表达模式在疾病组相对更低、在对照组相对更高。"
  },
  paste0(
    "关键模块按 |MM|≥0.80、|GS|≥0.20 且两者 P<0.05 标记出 ",
    hub_candidate_n, " 个候选枢纽基因。该阈值用于候选排序，不等同于功能验证。"
  ),
  paired_sentence,
  "",
  "## 7. 关键模块的生物学线索",
  "",
  bio_hint,
  paste0("按 |MM| 和 |GS| 综合排序靠前的基因包括：", paste(top_hub, collapse = "、"), "。"),
  "",
  "## 8. 结果解释边界",
  "",
  "- WGCNA 识别的是同一队列中的共表达关系及其与分组的关联，不能证明因果关系。",
  "- GSE29272 为组织水平芯片数据，模块可能同时反映肿瘤细胞状态、基质/免疫细胞比例和取材差异。",
  "- 模块—性状 Pearson 相关的逐样本 P 值未显式建模患者配对，因此报告中同时提供关键模块 ME 的配对敏感性检验。",
  "- 关键模块与候选枢纽基因应在独立队列中验证，并可进一步进行 GO/KEGG、细胞类型去卷积或蛋白水平验证。",
  "",
  "## 9. 主要输出文件",
  "",
  "- `GSE29272_WGCNA_hubmodule_genes.csv`：关键模块全部基因及 MM、GS、P 值和候选枢纽标记。",
  "- `GSE29272_WGCNA_module_trait_correlations.csv`：各模块与疾病状态的相关系数、P 值及 FDR。",
  "- `GSE29272_WGCNA_module_assignment.csv`：前 5000 个高 MAD 基因的模块归属。",
  "- `GSE29272_WGCNA_sample_QC.csv`：样本连接度、Z.k 和离群判定。",
  "- `GSE29272_WGCNA_soft_threshold_metrics.csv`：软阈值候选功率的拟合指标。",
  "- 所有核心图均同时输出矢量 PDF 与 600 dpi PNG。",
  "",
  "## 10. 可复现性",
  "",
  "完整 R 脚本位于 `scripts/GSE29272_WGCNA_analysis.R`；运行日志、参数摘要与 sessionInfo 均保存在 `results/WGCNA`。"
)
writeLines(
  report_lines,
  file.path(out_dir, "GSE29272_WGCNA_report.md"),
  useBytes = TRUE
)

summary_table <- data.frame(
  metric = c(
    "input_genes", "input_samples", "tumor_samples", "normal_samples",
    "MAD_selected_genes", "outlier_samples_removed", "network_samples",
    "soft_power", "R2_criterion_reached", "non_grey_modules",
    "grey_genes", "key_module", "key_module_genes", "key_module_correlation",
    "key_module_p_value", "key_module_FDR", "hub_candidates",
    "paired_samples", "paired_mean_ME_delta", "paired_t_p_value"
  ),
  value = as.character(c(
    n_genes_raw, n_samples_raw, unname(group_counts["Tumor"]),
    unname(group_counts["Normal"]), ncol(datExpr0), length(outlier_samples),
    nrow(datExpr), soft_power, soft_reached, length(non_grey_labels), grey_n,
    key_module, key_gene_n, key_cor, key_p, key_fdr, hub_candidate_n,
    paired_n, paired_mean_delta, paired_p
  )),
  stringsAsFactors = FALSE
)
write_csv_utf8(summary_table, file.path(out_dir, "GSE29272_WGCNA_summary.csv"))

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "GSE29272_WGCNA_sessionInfo.txt")
)

cat("\nAnalysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Output directory:", out_dir, "\n")
