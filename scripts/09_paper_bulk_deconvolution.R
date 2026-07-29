# 胃上皮身份丢失模块：全组织单细胞组成拆分与独立bulk标志基因去卷积
# 说明：
# 1) 单细胞部分使用GSE206785作者注释，统计单位为患者。
# 2) bulk芯片为已标准化log2强度，不适合直接输入需要原始计数的BayesPrism。
# 3) 因此bulk采用BisqueRNA MarkerBasedDecomposition（PCA型相对丰度分数），
#    并使用与73基因模块不重叠的胃上皮谱系标志基因，避免循环论证。

suppressPackageStartupMessages({
  library(data.table)
  library(Biobase)
  library(BisqueRNA)
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
  if (length(file_arg) != 1L) stop("请使用Rscript运行本脚本。")
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
clean_dir <- file.path(project_dir, "clean_data")
result_dir <- file.path(project_dir, "results", "PaperValidation", "Deconvolution")
sc_dir <- file.path(project_dir, "results", "PaperValidation", "SingleCell")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_dir, "09_paper_bulk_deconvolution_run.log")
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
cat("BisqueRNA版本：", as.character(packageVersion("BisqueRNA")), "\n")

candidate <- fread(file.path(
  project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"
))
candidate_genes <- candidate$gene
stopifnot(length(candidate_genes) == 73L)

cache <- readRDS(file.path(sc_dir, "GSE206785_selected_expression_and_metadata.rds"))
meta <- as.data.table(cache$metadata)
expr_sc <- cache$expression
candidate_available_sc <- intersect(candidate_genes, colnames(expr_sc))
marker_sets <- cache$marker_sets
lineage_markers <- cache$lineage_markers

read_clean_matrix <- function(path) {
  x <- fread(path)
  genes <- x[[1L]]
  x[[1L]] <- NULL
  mat <- as.matrix(x)
  rownames(mat) <- genes
  storage.mode(mat) <- "double"
  mat
}

standardized_module_score <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  z <- t(scale(t(expr[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

extract_tissue_status <- function(meta) {
  if ("disease_status" %in% colnames(meta)) {
    return(ifelse(meta$disease_status == "Tumor", "Tumor", "Normal"))
  }
  if ("group" %in% colnames(meta)) {
    return(ifelse(grepl("Tumor", meta$group), "Tumor", "Normal"))
  }
  stop("样本表缺少分组字段。")
}

# ---------------------------
# A. GSE206785全组织组成拆分
# ---------------------------
cand_sc <- expr_sc[, candidate_available_sc, drop = FALSE]

allcell_dt <- data.table(
  Patient = meta$Patient,
  Sample = meta$Sample,
  Tissue = meta$Tissue,
  Type = meta$Type,
  as.data.table(cand_sc)
)

sample_gene_means <- allcell_dt[
  ,
  c(list(n_cells = .N), lapply(.SD, mean)),
  by = .(Patient, Sample, Tissue),
  .SDcols = candidate_available_sc
]

# scRNA组织组成容易受消化与捕获影响；预先要求每个样本至少100个细胞
sample_gene_means[, eligible := n_cells >= 100L]
eligible_samples <- sample_gene_means[eligible == TRUE]
sgm <- as.matrix(eligible_samples[, ..candidate_available_sc])
sgm_mu <- colMeans(sgm)
sgm_sd <- apply(sgm, 2L, sd)
sgm_sd[!is.finite(sgm_sd) | sgm_sd == 0] <- 1
sgm_z <- sweep(sweep(sgm, 2L, sgm_mu, "-"), 2L, sgm_sd, "/")
eligible_samples[, observed_module_score := rowMeans(sgm_z)]

type_counts <- allcell_dt[, .N, by = .(Patient, Sample, Tissue, Type)]
type_totals <- type_counts[, .(total_cells = sum(N)), by = .(Patient, Sample, Tissue)]
type_props <- merge(type_counts, type_totals, by = c("Patient", "Sample", "Tissue"))
type_props[, proportion := N / total_cells]
type_props <- type_props[total_cells >= 100L]

all_types <- sort(unique(meta$Type))
sample_index <- unique(type_props[, .(Patient, Sample, Tissue, total_cells)])
prop_grid <- CJ(
  Patient = unique(sample_index$Patient),
  Tissue = c("Normal", "Tumor"),
  Type = all_types,
  unique = TRUE
)
prop_grid <- merge(prop_grid, sample_index, by = c("Patient", "Tissue"), all.x = TRUE)
prop_grid <- merge(
  prop_grid,
  type_props[, .(Patient, Tissue, Type, N, proportion)],
  by = c("Patient", "Tissue", "Type"),
  all.x = TRUE
)
prop_grid[!is.na(total_cells) & is.na(N), `:=`(N = 0L, proportion = 0)]
prop_grid <- prop_grid[!is.na(total_cells)]

epi_props <- prop_grid[Type == "Epithelial", .(
  Patient, Sample, Tissue, epithelial_proportion = proportion, total_cells
)]

# 使用正常组织各大类细胞的候选基因均值作为固定参考
normal_ref <- allcell_dt[Tissue == "Normal", lapply(.SD, mean), by = Type, .SDcols = candidate_available_sc]
ref_mat <- as.matrix(normal_ref[, ..candidate_available_sc])
rownames(ref_mat) <- normal_ref$Type

expected_mat <- matrix(
  NA_real_,
  nrow = nrow(eligible_samples),
  ncol = length(candidate_available_sc),
  dimnames = list(NULL, candidate_available_sc)
)
for (i in seq_len(nrow(eligible_samples))) {
  p <- prop_grid[
    Patient == eligible_samples$Patient[i] &
      Tissue == eligible_samples$Tissue[i] &
      Type %in% rownames(ref_mat),
    .(Type, proportion)
  ]
  w <- setNames(p$proportion, p$Type)[rownames(ref_mat)]
  w[is.na(w)] <- 0
  if (sum(w) > 0) w <- w / sum(w)
  expected_mat[i, ] <- colSums(ref_mat * w)
}
expected_z <- sweep(sweep(expected_mat, 2L, sgm_mu, "-"), 2L, sgm_sd, "/")
eligible_samples[, composition_expected_score := rowMeans(expected_z)]
eligible_samples <- merge(
  eligible_samples,
  epi_props[, .(Patient, Tissue, epithelial_proportion)],
  by = c("Patient", "Tissue"),
  all.x = TRUE
)

complete_sc <- eligible_samples[
  Patient %in% eligible_samples[, uniqueN(Tissue), by = Patient][V1 == 2L, Patient]
]
obs_sc <- dcast(complete_sc, Patient ~ Tissue, value.var = "observed_module_score")
exp_sc <- dcast(complete_sc, Patient ~ Tissue, value.var = "composition_expected_score")
epi_sc <- dcast(complete_sc, Patient ~ Tissue, value.var = "epithelial_proportion")
common <- Reduce(intersect, list(obs_sc$Patient, exp_sc$Patient, epi_sc$Patient))
obs_sc <- obs_sc[match(common, Patient)]
exp_sc <- exp_sc[match(common, Patient)]
epi_sc <- epi_sc[match(common, Patient)]

obs_diff <- obs_sc$Tumor - obs_sc$Normal
exp_diff <- exp_sc$Tumor - exp_sc$Normal
epi_diff <- epi_sc$Tumor - epi_sc$Normal
obs_test <- t.test(obs_sc$Tumor, obs_sc$Normal, paired = TRUE)
exp_test <- t.test(exp_sc$Tumor, exp_sc$Normal, paired = TRUE)
epi_test <- t.test(epi_sc$Tumor, epi_sc$Normal, paired = TRUE)

complete_sc[, Tissue := factor(Tissue, levels = c("Normal", "Tumor"))]
complete_sc[, Patient := factor(Patient)]
fit_sc_unadj <- lm(observed_module_score ~ Patient + Tissue, data = complete_sc)
fit_sc_adj <- lm(
  observed_module_score ~ Patient + epithelial_proportion + Tissue,
  data = complete_sc
)
sc_unadj <- coef(summary(fit_sc_unadj))["TissueTumor", ]
sc_adj <- coef(summary(fit_sc_adj))["TissueTumor", ]

sc_summary <- data.table(
  paired_patients = length(common),
  observed_difference = mean(obs_diff),
  observed_pvalue = obs_test$p.value,
  composition_expected_difference = mean(exp_diff),
  composition_expected_pvalue = exp_test$p.value,
  residual_difference = mean(obs_diff) - mean(exp_diff),
  descriptive_fraction_explained = mean(exp_diff) / mean(obs_diff),
  epithelial_proportion_normal = mean(epi_sc$Normal),
  epithelial_proportion_tumor = mean(epi_sc$Tumor),
  epithelial_proportion_difference = mean(epi_diff),
  epithelial_proportion_pvalue = epi_test$p.value,
  tissue_coef_unadjusted = sc_unadj[1L],
  tissue_p_unadjusted = sc_unadj[4L],
  tissue_coef_adjusted_for_epithelial_prop = sc_adj[1L],
  tissue_p_adjusted_for_epithelial_prop = sc_adj[4L],
  tissue_coef_attenuation = 1 - sc_adj[1L] / sc_unadj[1L]
)

fwrite(prop_grid, file.path(result_dir, "GSE206785_major_celltype_proportions.csv"))
fwrite(eligible_samples, file.path(result_dir, "GSE206785_whole_tissue_module_scores.csv"))
fwrite(sc_summary, file.path(result_dir, "GSE206785_whole_tissue_composition_summary.csv"))

# -----------------------------------
# B. bulk标志基因去卷积和组成校正敏感性
# -----------------------------------
lineage_marker_table <- rbindlist(lapply(names(lineage_markers), function(ct) {
  data.table(cluster = ct, gene = lineage_markers[[ct]])
}))
lineage_marker_table <- lineage_marker_table[!gene %in% candidate_genes]

major_marker_table <- rbindlist(lapply(names(marker_sets), function(ct) {
  data.table(cluster = ct, gene = marker_sets[[ct]])
}))
major_marker_table <- major_marker_table[!gene %in% candidate_genes]

run_bisque_marker <- function(expr, markers) {
  markers <- markers[gene %in% rownames(expr)]
  counts <- markers[, .N, by = cluster]
  keep_types <- counts[N >= 3L, cluster]
  markers <- markers[cluster %in% keep_types]
  if (uniqueN(markers$cluster) < 2L) stop("可用标志基因不足以去卷积。")

  bulk_eset <- ExpressionSet(assayData = expr)
  out <- BisqueRNA::MarkerBasedDecomposition(
    bulk.eset = bulk_eset,
    markers = as.data.frame(markers),
    ct_col = "cluster",
    gene_col = "gene",
    min_gene = 3,
    max_gene = 50,
    weighted = FALSE,
    unique_markers = FALSE,
    verbose = FALSE
  )
  scores <- out$bulk.props

  # PCA第一主成分符号任意；以对应标志基因平均表达的相关方向校正
  for (ct in rownames(scores)) {
    genes <- markers[cluster == ct & gene %in% rownames(expr), gene]
    marker_mean <- colMeans(expr[genes, , drop = FALSE])
    if (cor(scores[ct, ], marker_mean, use = "pairwise.complete.obs") < 0) {
      scores[ct, ] <- -scores[ct, ]
    }
  }
  scores <- t(scale(t(scores)))
  scores[!is.finite(scores)] <- 0
  list(scores = scores, genes = out$genes.used, marker_counts = counts)
}

validate_bulk_deconv <- function(dataset, expr, meta) {
  stopifnot(identical(colnames(expr), meta$sample))
  status <- extract_tissue_status(meta)
  pair_id <- meta$pair_id
  module_score <- standardized_module_score(expr, candidate_genes)

  lin <- run_bisque_marker(expr, lineage_marker_table)
  maj <- run_bisque_marker(expr, major_marker_table)
  # Proliferative是状态程序，不是成熟胃上皮谱系；不纳入胃身份组成复合分数
  identity_rows <- setdiff(rownames(lin$scores), "Proliferative")
  lineage_composite <- colMeans(lin$scores[identity_rows, , drop = FALSE])
  epithelial_score <- if ("Epithelial" %in% rownames(maj$scores)) {
    maj$scores["Epithelial", ]
  } else {
    rep(NA_real_, ncol(expr))
  }

  sample_dt <- data.table(
    dataset = dataset,
    sample = colnames(expr),
    pair_id = pair_id,
    Tissue = status,
    module_score = module_score,
    gastric_lineage_deconv_score = lineage_composite,
    broad_epithelial_deconv_score = epithelial_score
  )
  pair_ok <- sample_dt[
    !is.na(pair_id),
    uniqueN(Tissue),
    by = pair_id
  ][V1 == 2L, pair_id]
  paired <- sample_dt[pair_id %in% pair_ok]
  paired[, Tissue := factor(Tissue, levels = c("Normal", "Tumor"))]
  paired[, pair_id := factor(pair_id)]

  unadj <- lm(module_score ~ pair_id + Tissue, data = paired)
  adj_lineage <- lm(
    module_score ~ pair_id + gastric_lineage_deconv_score + Tissue,
    data = paired
  )
  un <- coef(summary(unadj))["TissueTumor", ]
  ad <- coef(summary(adj_lineage))["TissueTumor", ]
  cor_lineage <- cor.test(
    paired$module_score,
    paired$gastric_lineage_deconv_score,
    method = "spearman",
    exact = FALSE
  )

  summary <- data.table(
    dataset = dataset,
    paired_patients = uniqueN(paired$pair_id),
    module_genes = sum(candidate_genes %in% rownames(expr)),
    lineage_celltypes_estimated = nrow(lin$scores),
    identity_lineages_in_composite = length(identity_rows),
    lineage_marker_genes_used = length(unique(unlist(lin$genes))),
    unadjusted_tissue_coefficient = un[1L],
    unadjusted_tissue_pvalue = un[4L],
    adjusted_tissue_coefficient = ad[1L],
    adjusted_tissue_pvalue = ad[4L],
    coefficient_attenuation = 1 - ad[1L] / un[1L],
    module_lineage_score_spearman = unname(cor_lineage$estimate),
    module_lineage_score_pvalue = cor_lineage$p.value
  )

  lin_long <- as.data.table(t(lin$scores), keep.rownames = "sample")
  lin_long <- melt(lin_long, id.vars = "sample", variable.name = "celltype", value.name = "score")
  lin_long[, dataset := dataset]
  maj_long <- as.data.table(t(maj$scores), keep.rownames = "sample")
  maj_long <- melt(maj_long, id.vars = "sample", variable.name = "celltype", value.name = "score")
  maj_long[, dataset := dataset]
  list(sample = sample_dt, paired = paired, summary = summary, lineage = lin_long, major = maj_long)
}

datasets <- list(
  GSE29272 = list(
    expr = read_clean_matrix(file.path(clean_dir, "GSE29272_clean_expression_matrix.csv")),
    meta = as.data.frame(fread(file.path(clean_dir, "GSE29272_sample_info.csv")))
  ),
  GSE79973 = list(
    expr = read_clean_matrix(file.path(clean_dir, "GSE79973_clean_expression_matrix.csv")),
    meta = as.data.frame(fread(file.path(clean_dir, "GSE79973_sample_info.csv")))
  ),
  GSE19826 = list(
    expr = read_clean_matrix(file.path(clean_dir, "GSE19826_clean_expression_matrix.csv")),
    meta = as.data.frame(fread(file.path(clean_dir, "GSE19826_sample_info.csv")))
  )
)

bulk_results <- lapply(names(datasets), function(ds) {
  validate_bulk_deconv(ds, datasets[[ds]]$expr, datasets[[ds]]$meta)
})
names(bulk_results) <- names(datasets)

bulk_summary <- rbindlist(lapply(bulk_results, `[[`, "summary"))
bulk_samples <- rbindlist(lapply(bulk_results, `[[`, "sample"))
bulk_lineage <- rbindlist(lapply(bulk_results, `[[`, "lineage"))
bulk_major <- rbindlist(lapply(bulk_results, `[[`, "major"))
fwrite(bulk_summary, file.path(result_dir, "bulk_deconvolution_adjustment_summary.csv"))
fwrite(bulk_samples, file.path(result_dir, "bulk_deconvolution_sample_scores.csv"))
fwrite(bulk_lineage, file.path(result_dir, "bulk_lineage_deconvolution_scores.csv"))
fwrite(bulk_major, file.path(result_dir, "bulk_major_celltype_deconvolution_scores.csv"))
fwrite(lineage_marker_table, file.path(result_dir, "bulk_lineage_marker_panel.csv"))
fwrite(major_marker_table, file.path(result_dir, "bulk_major_celltype_marker_panel.csv"))

# ---------------------------
# C. 出版级结果图
# ---------------------------
sc_plot_dt <- complete_sc[, .(
  Patient = as.character(Patient),
  Tissue = as.character(Tissue),
  observed_module_score,
  composition_expected_score,
  epithelial_proportion
)]

p_sc_prop <- ggplot(sc_plot_dt, aes(Tissue, epithelial_proportion, group = Patient)) +
  geom_line(colour = "#B8B8B8", linewidth = 0.35, alpha = 0.8) +
  geom_point(aes(colour = Tissue), size = 1.7) +
  scale_colour_manual(values = c(Normal = "#8FA6A0", Tumor = "#B56F70")) +
  labs(x = NULL, y = "Epithelial-cell proportion") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "none")

sc_decomp_plot <- data.table(
  component = factor(
    c("Observed", "Composition expected", "Residual"),
    levels = c("Observed", "Composition expected", "Residual")
  ),
  difference = c(
    mean(obs_diff),
    mean(exp_diff),
    mean(obs_diff) - mean(exp_diff)
  )
)
p_sc_decomp <- ggplot(sc_decomp_plot, aes(component, difference, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#666666") +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(
    "Observed" = "#6F8795",
    "Composition expected" = "#B0A58B",
    "Residual" = "#9A86A4"
  )) +
  labs(x = NULL, y = "Tumor - Normal module score") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

scatter_dt <- bulk_samples
p_bulk_scatter <- ggplot(
  scatter_dt,
  aes(gastric_lineage_deconv_score, module_score, colour = Tissue)
) +
  geom_point(size = 0.9, alpha = 0.65) +
  geom_smooth(method = "lm", se = FALSE, colour = "#444444", linewidth = 0.45) +
  facet_wrap(~dataset, scales = "free") +
  scale_colour_manual(values = c(Normal = "#8FA6A0", Tumor = "#B56F70")) +
  labs(
    x = "Independent gastric-lineage deconvolution score",
    y = "73-gene module score",
    colour = NULL
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "top", strip.background = element_blank())

attenuation_dt <- bulk_summary[, .(
  dataset,
  unadjusted = unadjusted_tissue_coefficient,
  adjusted = adjusted_tissue_coefficient
)]
attenuation_long <- melt(
  attenuation_dt,
  id.vars = "dataset",
  variable.name = "model",
  value.name = "tissue_coefficient"
)
p_attenuation <- ggplot(
  attenuation_long,
  aes(model, tissue_coefficient, group = dataset, colour = dataset)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#888888") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c(
    GSE29272 = "#6F8795",
    GSE79973 = "#9A86A4",
    GSE19826 = "#B56F70"
  )) +
  labs(x = NULL, y = "Tumor coefficient", colour = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "top")

fig <- (p_sc_prop | p_sc_decomp) / (p_bulk_scatter | p_attenuation) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

ggsave(
  file.path(result_dir, "composition_aware_validation_figure.pdf"),
  fig, width = 183 / 25.4, height = 150 / 25.4, device = cairo_pdf
)
ggsave(
  file.path(result_dir, "composition_aware_validation_figure.png"),
  fig, width = 183 / 25.4, height = 150 / 25.4, dpi = 600, bg = "white"
)
fwrite(sc_plot_dt, file.path(result_dir, "figure_sc_proportion_source_data.csv"))
fwrite(sc_decomp_plot, file.path(result_dir, "figure_sc_decomposition_source_data.csv"))
fwrite(scatter_dt, file.path(result_dir, "figure_bulk_scatter_source_data.csv"))
fwrite(attenuation_long, file.path(result_dir, "figure_bulk_attenuation_source_data.csv"))

report <- c(
  "# 细胞组成校正与bulk去卷积报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- BisqueRNA版本：", as.character(packageVersion("BisqueRNA"))),
  "- 单细胞主验证：GSE206785，作者注释，24位配对患者。",
  "- bulk队列：GSE29272、GSE79973和GSE19826。",
  "- bulk输入为已标准化Affymetrix log2强度，因此使用BisqueRNA标志基因PCA相对丰度分数，而非要求原始计数的BayesPrism。",
  "- 去卷积标志基因预先按胃上皮谱系与主要细胞类型定义，并移除所有73模块基因。",
  "",
  "## GSE206785全组织组成拆分",
  "",
  paste(capture.output(print(sc_summary)), collapse = "\n"),
  "",
  "## bulk组成校正敏感性分析",
  "",
  paste(capture.output(print(bulk_summary)), collapse = "\n"),
  "",
  "## 解释边界",
  "",
  "- 单细胞捕获比例受组织消化、细胞活性和平台影响，不等同于组织学面积比例。",
  "- BisqueRNA标志基因结果为队列内相对丰度分数，不是绝对细胞百分比。",
  "- 模块分数与独立胃上皮谱系分数的相关及系数衰减支持组成混杂，但不构成因果证明。",
  "- 若校正后Tumor系数消失或显著减弱，应把bulk Hub信号解释为组织组成敏感标志，而非肿瘤细胞内在抑制程序。"
)
writeLines(report, file.path(result_dir, "composition_aware_validation_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "09_sessionInfo.txt"), useBytes = TRUE)

cat("\n单细胞全组织组成摘要：\n")
print(sc_summary)
cat("\nbulk去卷积摘要：\n")
print(bulk_summary)
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
