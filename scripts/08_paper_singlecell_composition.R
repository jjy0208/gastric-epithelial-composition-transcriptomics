# GSE206785单细胞组成感知分析
# 目的：确定73基因模块的细胞类型来源，并区分上皮谱系比例变化与谱系内表达变化
# 统计单位：患者；不把单个细胞作为独立生物学重复

suppressPackageStartupMessages({
  library(data.table)
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
raw_dir <- file.path(project_dir, "raw_data")
result_dir <- file.path(project_dir, "results", "PaperValidation", "SingleCell")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_dir, "08_paper_singlecell_composition_run.log")
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

candidate <- fread(file.path(
  project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"
))
stopifnot(nrow(candidate) == 73L, uniqueN(candidate$gene) == 73L)
candidate_genes <- candidate$gene

# 用于细胞类型核对和后续去卷积的独立标志基因。
# 后续会自动移除与73基因模块重叠的基因，避免循环论证。
marker_sets <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "TACSTD2", "CLDN4", "CLDN7"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "PDGFRA", "CFD", "C7"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "ESAM", "RAMP2"),
  Myeloid = c("LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "LYZ", "CSF1R"),
  T_NK = c("CD3D", "CD3E", "TRAC", "CD247", "NKG7", "GNLY", "KLRD1", "TRDC"),
  B_Plasma = c("MS4A1", "CD79A", "CD37", "CD74", "JCHAIN", "MZB1", "SDC1", "IGHG1"),
  Mast = c("KIT", "TPSAB1", "TPSB2", "CPA3"),
  Mural = c("RGS5", "CSPG4", "MCAM", "ACTA2", "NOTCH3"),
  Glial = c("S100B", "SOX10", "S100A1", "PLP1")
)

lineage_markers <- list(
  Parietal = c("CKB", "SLC26A7", "SLC4A2", "KCNQ1", "CLIC6", "SLC9A4"),
  Chief = c("PGA3", "PGA4", "PGA5", "PNLIPRP2", "PRSS1"),
  Pit_mucous = c("GKN2", "CAPN8", "VSIG1", "SULT1C2"),
  Gland_mucous = c("AQP5", "CLDN18", "PRR4", "LYZ", "WFDC2", "AGR2"),
  Enterocyte = c("APOA1", "APOA4", "FABP1", "ALPI", "KRT20"),
  Proliferative = c("MKI67", "TOP2A", "STMN1", "HMGB2")
)

marker_sets <- lapply(marker_sets, setdiff, y = candidate_genes)
lineage_markers <- lapply(lineage_markers, setdiff, y = candidate_genes)

meta_path <- file.path(raw_dir, "GSE206785_metadata.txt.gz")
expr_path <- file.path(raw_dir, "GSE206785_scgex.txt.gz")
if (!file.exists(meta_path) || !file.exists(expr_path)) {
  stop("缺少GSE206785单细胞表达或metadata文件。")
}

meta <- fread(meta_path)
required_meta <- c("Sample", "Patient", "Tissue", "Platform", "Subtype", "Type", "Annotation")
stopifnot(all(required_meta %in% colnames(meta)))
stopifnot(uniqueN(meta$Patient) == 24L, uniqueN(meta$Sample) == 48L)

# 只提取锁定模块和预先定义的细胞标志基因，不载入约2万基因的完整宽矩阵
con <- gzfile(expr_path, open = "rt")
header <- readLines(con, n = 1L, warn = FALSE)
close(con)
all_genes <- strsplit(header, ",", fixed = TRUE)[[1L]]
all_genes <- gsub('^"|"$', "", all_genes)
requested <- unique(c(
  "Cell", candidate_genes,
  unlist(marker_sets, use.names = FALSE),
  unlist(lineage_markers, use.names = FALSE)
))
available <- intersect(requested, all_genes)
missing_candidate <- setdiff(candidate_genes, available)
if (length(missing_candidate)) {
  warning("单细胞矩阵缺少候选基因：", paste(missing_candidate, collapse = ", "))
}
if (sum(candidate_genes %in% available) < 65L) {
  stop("单细胞矩阵中的候选基因覆盖不足。")
}

expr <- fread(expr_path, select = available, showProgress = TRUE)
stopifnot(nrow(expr) == nrow(meta), expr[[1L]][1L] == "AAACCTGAGAGTACCG-1-171012N")
cell_id <- expr$Cell
expr[, Cell := NULL]
expr_mat <- as.matrix(expr)
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- cell_id

sample_from_cell <- sub("^.*-1-", "", cell_id)
sample_mismatch <- sum(sample_from_cell != meta$Sample)
if (sample_mismatch != 0L) {
  stop("表达矩阵与metadata行顺序不一致，错配细胞数：", sample_mismatch)
}

candidate_available <- intersect(candidate_genes, colnames(expr_mat))
marker_available <- setdiff(colnames(expr_mat), candidate_available)

# 数据值审计：判断是否与log1p整数计数一致；只作为文件语义记录，不反推缺失全基因库大小
audit_values <- as.numeric(expr_mat[
  seq_len(min(5000L, nrow(expr_mat))),
  seq_len(min(30L, ncol(expr_mat))),
  drop = FALSE
])
audit_values <- audit_values[is.finite(audit_values) & audit_values > 0]
log1p_integer_max_error <- if (length(audit_values)) {
  max(abs(expm1(audit_values) - round(expm1(audit_values))))
} else {
  NA_real_
}

lineage_from_annotation <- function(x) {
  fifelse(grepl("Parietal", x), "Parietal",
  fifelse(grepl("Chief", x), "Chief",
  fifelse(grepl("_PMC_", x), "Pit_mucous",
  fifelse(grepl("_GMC_", x), "Gland_mucous",
  fifelse(grepl("Enterocyte", x), "Enterocyte",
  fifelse(grepl("_MSC_", x), "Metaplastic_stem", "Unidentified"))))))
}

meta[, EpithelialLineage := ifelse(
  Type == "Epithelial",
  lineage_from_annotation(Annotation),
  NA_character_
)]

# 逐基因标准化后等权计算每个细胞的模块分数；细胞只用于描述，患者才是推断统计单位
cand_mat <- expr_mat[, candidate_available, drop = FALSE]
cand_center <- colMeans(cand_mat)
cand_scale <- apply(cand_mat, 2L, sd)
cand_scale[!is.finite(cand_scale) | cand_scale == 0] <- 1
cand_z <- sweep(sweep(cand_mat, 2L, cand_center, "-"), 2L, cand_scale, "/")
cell_score <- rowMeans(cand_z)
rm(cand_z)

cell_meta <- copy(meta)
cell_meta[, Cell := cell_id]
cell_meta[, module_score := cell_score]

# 各作者定义细胞类型的患者-组织pseudobulk分数
type_pseudo <- cell_meta[
  ,
  .(
    n_cells = .N,
    module_score = mean(module_score),
    module_score_median = median(module_score)
  ),
  by = .(Patient, Sample, Tissue, Type)
][n_cells >= 20L]
fwrite(type_pseudo, file.path(result_dir, "singlecell_type_pseudobulk_scores.csv"))

# 上皮患者-组织pseudobulk：先对每个候选基因按细胞求平均，再在pseudobulk层面逐基因标准化
epi_idx <- which(meta$Type == "Epithelial")
epi_dt <- data.table(
  Patient = meta$Patient[epi_idx],
  Sample = meta$Sample[epi_idx],
  Tissue = meta$Tissue[epi_idx],
  as.data.table(cand_mat[epi_idx, , drop = FALSE])
)
epi_pseudo <- epi_dt[
  ,
  c(list(n_cells = .N), lapply(.SD, mean)),
  by = .(Patient, Sample, Tissue),
  .SDcols = candidate_available
]
epi_keep <- epi_pseudo$n_cells >= 20L
epi_gene_mat <- as.matrix(epi_pseudo[, ..candidate_available])
gene_mu <- colMeans(epi_gene_mat[epi_keep, , drop = FALSE])
gene_sd <- apply(epi_gene_mat[epi_keep, , drop = FALSE], 2L, sd)
gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- 1
epi_gene_z <- sweep(sweep(epi_gene_mat, 2L, gene_mu, "-"), 2L, gene_sd, "/")
epi_pseudo[, module_score := rowMeans(epi_gene_z)]
epi_pseudo[, eligible := epi_keep]

eligible_epi <- epi_pseudo[eligible == TRUE]
pair_counts <- eligible_epi[, uniqueN(Tissue), by = Patient]
complete_patients <- pair_counts[V1 == 2L, Patient]
epi_paired <- eligible_epi[Patient %in% complete_patients]
epi_wide <- dcast(epi_paired, Patient ~ Tissue, value.var = "module_score")
stopifnot(all(c("Normal", "Tumor") %in% colnames(epi_wide)))
epi_diff <- epi_wide$Tumor - epi_wide$Normal
epi_t <- t.test(epi_wide$Tumor, epi_wide$Normal, paired = TRUE)
epi_w <- suppressWarnings(wilcox.test(
  epi_wide$Tumor, epi_wide$Normal, paired = TRUE, exact = FALSE
))

fwrite(epi_pseudo, file.path(result_dir, "epithelial_patient_pseudobulk_scores.csv"))
fwrite(epi_paired, file.path(result_dir, "epithelial_complete_pair_scores.csv"))

# 上皮谱系内pseudobulk与患者固定效应模型
epi_meta <- meta[epi_idx]
lineage_dt <- data.table(
  Patient = epi_meta$Patient,
  Sample = epi_meta$Sample,
  Tissue = epi_meta$Tissue,
  Lineage = epi_meta$EpithelialLineage,
  as.data.table(cand_mat[epi_idx, , drop = FALSE])
)
lineage_pseudo <- lineage_dt[
  ,
  c(list(n_cells = .N), lapply(.SD, mean)),
  by = .(Patient, Sample, Tissue, Lineage),
  .SDcols = candidate_available
][n_cells >= 20L]

lin_gene_mat <- as.matrix(lineage_pseudo[, ..candidate_available])
lin_mu <- colMeans(lin_gene_mat)
lin_sd <- apply(lin_gene_mat, 2L, sd)
lin_sd[!is.finite(lin_sd) | lin_sd == 0] <- 1
lin_z <- sweep(sweep(lin_gene_mat, 2L, lin_mu, "-"), 2L, lin_sd, "/")
lineage_pseudo[, module_score := rowMeans(lin_z)]
lineage_pseudo[, Tissue := factor(Tissue, levels = c("Normal", "Tumor"))]
lineage_pseudo[, Patient := factor(Patient)]
lineage_pseudo[, Lineage := factor(Lineage)]

fit <- lm(
  module_score ~ Tissue + Patient + Lineage,
  data = lineage_pseudo,
  weights = sqrt(n_cells)
)
robust_vcov <- sandwich::vcovCL(fit, cluster = ~Patient, type = "HC1")
robust_test <- lmtest::coeftest(fit, vcov. = robust_vcov)
tissue_row <- robust_test["TissueTumor", ]

lineage_tests <- rbindlist(lapply(
  levels(lineage_pseudo$Lineage),
  function(lin) {
    d <- lineage_pseudo[Lineage == lin]
    cc <- d[, uniqueN(Tissue), by = Patient][V1 == 2L, Patient]
    d <- d[Patient %in% cc]
    if (uniqueN(d$Patient) < 5L) return(NULL)
    w <- dcast(d, Patient ~ Tissue, value.var = "module_score")
    tt <- t.test(w$Tumor, w$Normal, paired = TRUE)
    data.table(
      Lineage = as.character(lin),
      paired_patients = nrow(w),
      mean_tumor_minus_normal = mean(w$Tumor - w$Normal),
      ci_low = unname(tt$conf.int[1L]),
      ci_high = unname(tt$conf.int[2L]),
      pvalue = tt$p.value
    )
  }
), fill = TRUE)
if (nrow(lineage_tests)) {
  lineage_tests[, padj := p.adjust(pvalue, method = "BH")]
}
fwrite(lineage_pseudo, file.path(result_dir, "epithelial_lineage_pseudobulk_scores.csv"))
fwrite(lineage_tests, file.path(result_dir, "epithelial_lineage_paired_tests.csv"))

# 组成变化：每位患者每种组织的上皮谱系比例
lineage_counts <- epi_meta[
  ,
  .N,
  by = .(Patient, Sample, Tissue, EpithelialLineage)
]
sample_totals <- lineage_counts[, .(total_epithelial = sum(N)), by = .(Patient, Sample, Tissue)]
lineage_counts <- merge(
  lineage_counts, sample_totals,
  by = c("Patient", "Sample", "Tissue")
)
lineage_counts[, proportion := N / total_epithelial]
lineage_counts <- lineage_counts[total_epithelial >= 20L]

all_lin <- sort(unique(epi_meta$EpithelialLineage))
all_sample <- unique(lineage_counts[, .(Patient, Sample, Tissue, total_epithelial)])
lineage_grid <- CJ(
  Patient = unique(all_sample$Patient),
  Tissue = c("Normal", "Tumor"),
  EpithelialLineage = all_lin,
  unique = TRUE
)
lineage_grid <- merge(
  lineage_grid,
  all_sample,
  by = c("Patient", "Tissue"),
  all.x = TRUE
)
lineage_grid <- merge(
  lineage_grid,
  lineage_counts[, .(Patient, Tissue, EpithelialLineage, N, proportion)],
  by = c("Patient", "Tissue", "EpithelialLineage"),
  all.x = TRUE
)
lineage_grid[!is.na(total_epithelial) & is.na(N), `:=`(N = 0L, proportion = 0)]
lineage_grid <- lineage_grid[!is.na(total_epithelial)]
fwrite(lineage_grid, file.path(result_dir, "epithelial_lineage_proportions.csv"))

composition_tests <- rbindlist(lapply(all_lin, function(lin) {
  d <- lineage_grid[EpithelialLineage == lin]
  cc <- d[, uniqueN(Tissue), by = Patient][V1 == 2L, Patient]
  d <- d[Patient %in% cc]
  if (uniqueN(d$Patient) < 5L) return(NULL)
  w <- dcast(d, Patient ~ Tissue, value.var = "proportion")
  tt <- t.test(w$Tumor, w$Normal, paired = TRUE)
  data.table(
    Lineage = lin,
    paired_patients = nrow(w),
    mean_normal = mean(w$Normal),
    mean_tumor = mean(w$Tumor),
    mean_difference = mean(w$Tumor - w$Normal),
    ci_low = unname(tt$conf.int[1L]),
    ci_high = unname(tt$conf.int[2L]),
    pvalue = tt$p.value
  )
}), fill = TRUE)
if (nrow(composition_tests)) composition_tests[, padj := p.adjust(pvalue, method = "BH")]
fwrite(composition_tests, file.path(result_dir, "epithelial_lineage_composition_tests.csv"))

# 组成预期拆分：
# 使用正常组织各谱系的候选基因均值作为固定参考，
# 将每个样本的谱系比例加权得到“若只有组成变化”的预期表达。
normal_idx <- epi_idx[meta$Tissue[epi_idx] == "Normal"]
normal_lineage_ref <- data.table(
  Lineage = meta$EpithelialLineage[normal_idx],
  as.data.table(cand_mat[normal_idx, , drop = FALSE])
)[
  ,
  lapply(.SD, mean),
  by = Lineage,
  .SDcols = candidate_available
]
ref_lineages <- normal_lineage_ref$Lineage
ref_mat <- as.matrix(normal_lineage_ref[, ..candidate_available])
rownames(ref_mat) <- ref_lineages

decomp_samples <- eligible_epi[, .(Patient, Sample, Tissue, module_score)]
expected_mat <- matrix(
  NA_real_,
  nrow = nrow(decomp_samples),
  ncol = length(candidate_available),
  dimnames = list(NULL, candidate_available)
)
for (i in seq_len(nrow(decomp_samples))) {
  prop <- lineage_grid[
    Patient == decomp_samples$Patient[i] &
      Tissue == decomp_samples$Tissue[i] &
      EpithelialLineage %in% ref_lineages,
    .(EpithelialLineage, proportion)
  ]
  weights <- setNames(prop$proportion, prop$EpithelialLineage)
  weights <- weights[ref_lineages]
  weights[is.na(weights)] <- 0
  if (sum(weights) > 0) weights <- weights / sum(weights)
  expected_mat[i, ] <- colSums(ref_mat * weights)
}

obs_idx <- match(
  paste(decomp_samples$Patient, decomp_samples$Tissue),
  paste(epi_pseudo$Patient, epi_pseudo$Tissue)
)
observed_mat <- as.matrix(epi_pseudo[obs_idx, ..candidate_available])
observed_z <- sweep(sweep(observed_mat, 2L, gene_mu, "-"), 2L, gene_sd, "/")
expected_z <- sweep(sweep(expected_mat, 2L, gene_mu, "-"), 2L, gene_sd, "/")
decomp_samples[, observed_score := rowMeans(observed_z)]
decomp_samples[, composition_expected_score := rowMeans(expected_z)]

decomp_complete <- decomp_samples[
  Patient %in% decomp_samples[, uniqueN(Tissue), by = Patient][V1 == 2L, Patient]
]
obs_wide <- dcast(decomp_complete, Patient ~ Tissue, value.var = "observed_score")
exp_wide <- dcast(
  decomp_complete, Patient ~ Tissue,
  value.var = "composition_expected_score"
)
common_patients <- intersect(obs_wide$Patient, exp_wide$Patient)
obs_wide <- obs_wide[match(common_patients, Patient)]
exp_wide <- exp_wide[match(common_patients, Patient)]
observed_difference <- mean(obs_wide$Tumor - obs_wide$Normal)
composition_difference <- mean(exp_wide$Tumor - exp_wide$Normal)
residual_difference <- observed_difference - composition_difference
fraction_explained <- composition_difference / observed_difference

decomposition_summary <- data.table(
  paired_patients = length(common_patients),
  observed_tumor_minus_normal = observed_difference,
  composition_expected_difference = composition_difference,
  residual_after_composition = residual_difference,
  descriptive_fraction_explained = fraction_explained
)
fwrite(decomp_samples, file.path(result_dir, "epithelial_composition_expected_scores.csv"))
fwrite(decomposition_summary, file.path(result_dir, "epithelial_composition_decomposition_summary.csv"))

analysis_summary <- data.table(
  metric = c(
    "cells", "patients", "samples", "candidate_genes_available",
    "expression_metadata_mismatches", "log1p_integer_max_error",
    "normal_epithelial_cells", "tumor_epithelial_cells",
    "eligible_complete_epithelial_pairs",
    "epithelial_score_mean_difference",
    "epithelial_score_paired_t_pvalue",
    "epithelial_score_paired_wilcoxon_pvalue",
    "lineage_adjusted_tissue_coefficient",
    "lineage_adjusted_cluster_robust_se",
    "lineage_adjusted_cluster_robust_pvalue"
  ),
  value = as.character(c(
    nrow(meta), uniqueN(meta$Patient), uniqueN(meta$Sample),
    length(candidate_available), sample_mismatch,
    log1p_integer_max_error,
    sum(meta$Tissue == "Normal" & meta$Type == "Epithelial"),
    sum(meta$Tissue == "Tumor" & meta$Type == "Epithelial"),
    nrow(epi_wide), mean(epi_diff), epi_t$p.value, epi_w$p.value,
    tissue_row[1L], tissue_row[2L], tissue_row[4L]
  ))
)
fwrite(analysis_summary, file.path(result_dir, "singlecell_analysis_summary.csv"))

# 图1：细胞类型定位；图中每个点为患者-组织-细胞类型pseudobulk
type_order <- type_pseudo[
  ,
  .(median_score = median(module_score)),
  by = Type
][order(median_score), Type]
type_pseudo[, Type := factor(Type, levels = type_order)]
p_type <- ggplot(type_pseudo, aes(module_score, Type, colour = Tissue)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#BDBDBD") +
  geom_point(
    position = position_jitter(height = 0.12, width = 0),
    size = 1.1, alpha = 0.65
  ) +
  scale_colour_manual(values = c(Normal = "#8FA6A0", Tumor = "#B56F70")) +
  labs(x = "73-gene module score", y = NULL, colour = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "top")

# 图2：上皮pseudobulk患者配对验证（主面板）
epi_plot <- copy(epi_paired)
epi_plot[, Tissue := factor(Tissue, levels = c("Normal", "Tumor"))]
p_pair <- ggplot(epi_plot, aes(Tissue, module_score, group = Patient)) +
  geom_line(colour = "#B8B8B8", linewidth = 0.35, alpha = 0.8) +
  geom_point(aes(colour = Tissue), size = 1.8) +
  scale_colour_manual(values = c(Normal = "#8FA6A0", Tumor = "#B56F70")) +
  labs(x = NULL, y = "Epithelial pseudobulk module score") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "none")

# 图3：谱系内效应
if (nrow(lineage_tests)) {
  lineage_tests[, Lineage := factor(Lineage, levels = Lineage[order(mean_tumor_minus_normal)])]
  p_lineage <- ggplot(
    lineage_tests,
    aes(mean_tumor_minus_normal, Lineage)
  ) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#BDBDBD") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.4) +
    geom_point(aes(colour = padj < 0.05), size = 1.8) +
    scale_colour_manual(values = c(`TRUE` = "#B56F70", `FALSE` = "#8D8D8D")) +
    labs(x = "Tumor - Normal module score", y = NULL, colour = "BH FDR < 0.05") +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(legend.position = "top")
} else {
  p_lineage <- ggplot() + theme_void() + labs(title = "No lineage met the paired threshold")
}

# 图4：组成预期与剩余变化
decomp_plot <- data.table(
  component = factor(
    c("Observed", "Composition expected", "Residual"),
    levels = c("Observed", "Composition expected", "Residual")
  ),
  difference = c(observed_difference, composition_difference, residual_difference)
)
p_decomp <- ggplot(decomp_plot, aes(component, difference, fill = component)) +
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

fig <- (p_pair | p_decomp) / (p_type | p_lineage) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

ggsave(
  file.path(result_dir, "singlecell_composition_figure.pdf"),
  fig, width = 183 / 25.4, height = 150 / 25.4, device = cairo_pdf
)
ggsave(
  file.path(result_dir, "singlecell_composition_figure.png"),
  fig, width = 183 / 25.4, height = 150 / 25.4, dpi = 600, bg = "white"
)
fwrite(decomp_plot, file.path(result_dir, "singlecell_decomposition_plot_source_data.csv"))

report <- c(
  "# GSE206785单细胞组成感知分析报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "- 队列：24位患者、48份配对胃癌/正常胃组织、111,140个作者质控并注释的细胞。",
  "- 统计单位：患者；单个细胞不作为独立重复。",
  "- 上皮pseudobulk最低细胞数：每个患者-组织20个。",
  "- 谱系内pseudobulk最低细胞数：每个患者-组织-谱系20个。",
  "- 肿瘤组织来源的Epithelial细胞表述为tumor-derived epithelial cells；未凭组织来源直接命名为恶性细胞。",
  "",
  "## 文件一致性",
  "",
  paste0("- 表达矩阵与metadata行数：", nrow(expr_mat), " / ", nrow(meta)),
  paste0("- 细胞ID推导样本与metadata错配：", sample_mismatch),
  paste0("- 可用候选基因：", length(candidate_available), " / 73"),
  paste0("- log1p整数一致性最大误差：", format(log1p_integer_max_error, scientific = TRUE)),
  "",
  "## 患者级结果",
  "",
  paste0("- 可用完整上皮配对：", nrow(epi_wide)),
  paste0("- 上皮模块分数Tumor-Normal均值：", signif(mean(epi_diff), 4)),
  paste0("- 配对t检验P值：", format(epi_t$p.value, scientific = TRUE)),
  paste0("- 配对Wilcoxon P值：", format(epi_w$p.value, scientific = TRUE)),
  paste0("- 谱系和患者校正后的Tumor系数：", signif(tissue_row[1L], 4)),
  paste0("- 患者聚类稳健P值：", format(tissue_row[4L], scientific = TRUE)),
  "",
  "## 组成拆分",
  "",
  paste(capture.output(print(decomposition_summary)), collapse = "\n"),
  "",
  "## 解释边界",
  "",
  "- 组成预期值使用正常组织各上皮谱系的平均表达作为固定参考，是描述性反事实而非因果分解。",
  "- 作者提供的是处理后单细胞表达矩阵；本研究没有重新进行细胞过滤、聚类或CNV推断。",
  "- 因缺少CNV或体细胞突变证据，肿瘤来源上皮细胞可能混有非恶性上皮细胞。",
  "- 谱系校正后仍存在的组织效应可支持细胞内在表达改变，但不能证明调控机制。"
)
writeLines(report, file.path(result_dir, "singlecell_composition_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "08_sessionInfo.txt"), useBytes = TRUE)

saveRDS(
  list(
    metadata = meta,
    cell_id = cell_id,
    expression = expr_mat,
    candidate_genes = candidate_available,
    marker_sets = marker_sets,
    lineage_markers = lineage_markers
  ),
  file.path(result_dir, "GSE206785_selected_expression_and_metadata.rds"),
  compress = "xz"
)

cat("\n分析摘要：\n")
print(analysis_summary)
cat("\n组成拆分：\n")
print(decomposition_summary)
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
