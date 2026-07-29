# 投稿论文主图：胃上皮转录模块的发现、复现与细胞组成校正
# 所有图形仅使用项目中已经计算并保存的结果；不在本脚本中修改统计结果。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
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
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
deg_dir <- file.path(project_dir, "results", "DEG")
wgcna_dir <- file.path(project_dir, "results", "WGCNA")
hub_dir <- file.path(project_dir, "results", "HubGene")
pv_dir <- file.path(project_dir, "results", "PaperValidation")
sc_dir <- file.path(pv_dir, "SingleCell")
dc_dir <- file.path(pv_dir, "Deconvolution")
out_dir <- file.path(project_dir, "results", "ManuscriptFigures")
source_dir <- file.path(out_dir, "SourceData")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, "10_paper_figures_run.log")
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

pal <- c(
  normal = "#8FA8A1", tumor = "#BC706D", blue = "#728C9A",
  purple = "#9A86A4", olive = "#A9A087", green = "#7E8F78",
  grey = "#B8B8B8", dark = "#333333", light = "#E8E4DF"
)

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.title = element_text(colour = pal["dark"]),
      axis.text = element_text(colour = pal["dark"]),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, colour = "#555555"),
      legend.title = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(7, 7, 7, 7)
    )
}

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 600, bg = "white")
}

# -------------------------
# Figure 1：发现队列
# -------------------------
deg_all <- fread(file.path(deg_dir, "GSE29272_DEG_all.csv"))
deg_sig <- fread(file.path(deg_dir, "GSE29272_DEG_significant.csv"))
module_cor <- fread(file.path(wgcna_dir, "GSE29272_WGCNA_module_trait_correlations.csv"))
candidate <- fread(file.path(hub_dir, "GSE29272_candidate_genes.csv"))

stopifnot(nrow(deg_sig) == 122L, nrow(candidate) == 73L)

cohort_df <- data.table(
  order = 1:4,
  stage = c("Discovery", "Bulk validation", "Single-cell test", "Composition adjustment"),
  dataset = c("GSE29272\n134 paired patients",
              "GSE79973 + GSE19826\n10 + 12 paired patients",
              "GSE206785\n24 paired patients; 111,140 cells",
              "GSE206785 + 3 bulk cohorts\npatient-level models")
)
f1a <- ggplot(cohort_df, aes(order, 1)) +
  geom_segment(aes(x = 1, xend = 4, y = 1, yend = 1),
               linewidth = 0.8, colour = pal["grey"]) +
  geom_point(size = 5, colour = unname(c(pal["blue"], pal["purple"], pal["green"], pal["olive"]))) +
  geom_label(aes(label = stage), vjust = -1.5, size = 3.0, linewidth = 0,
             fill = "white", fontface = "bold") +
  geom_text(aes(label = dataset), vjust = 2.0, size = 2.7, lineheight = 0.95) +
  scale_x_continuous(limits = c(0.7, 4.3)) +
  coord_cartesian(ylim = c(0.55, 1.45), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(18, 10, 18, 10))

deg_all[, status := fifelse(padj < 0.05 & log2FoldChange > 1, "Up",
                            fifelse(padj < 0.05 & log2FoldChange < -1, "Down", "Not significant"))]
deg_all[, neglog10 := -log10(pmax(padj, .Machine$double.xmin))]
label_genes <- unique(c(head(deg_all[status == "Down"][order(padj)]$gene, 6),
                        head(deg_all[status == "Up"][order(padj)]$gene, 4)))
f1b <- ggplot(deg_all, aes(log2FoldChange, neglog10, colour = status)) +
  geom_point(size = 0.8, alpha = 0.55) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35, colour = "#777777") +
  ggrepel::geom_text_repel(
    data = deg_all[gene %in% label_genes], aes(label = gene),
    size = 2.5, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.25,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Down" = pal[["blue"]], "Not significant" = "#D3D3D3",
                                 "Up" = pal[["tumor"]])) +
  labs(x = "Tumor–normal log2 fold change", y = expression(-log[10]("FDR")),
       subtitle = "Paired limma-TREAT; |log2FC| > 1 and FDR < 0.05") +
  theme_nature()

module_cor[, module := factor(module, levels = module[order(correlation)])]
f1c <- ggplot(module_cor, aes(module, correlation, fill = correlation)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_gradient2(low = pal[["blue"]], mid = "white", high = pal[["tumor"]], midpoint = 0) +
  coord_flip() +
  labs(x = NULL, y = "Module–disease correlation",
       subtitle = "M04: r = −0.756, FDR = 3.24 × 10⁻49") +
  theme_nature() +
  theme(legend.position = "none")

set_counts <- data.table(
  set = factor(c("Significant DEGs", "M04 genes", "Intersection"),
               levels = c("Significant DEGs", "M04 genes", "Intersection")),
  n = c(nrow(deg_sig), 192L, nrow(candidate))
)
f1d <- ggplot(set_counts, aes(set, n, fill = set)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.35, fontface = "bold", size = 3.2) +
  scale_fill_manual(values = unname(c(pal["tumor"], pal["purple"], pal["green"]))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of genes",
       subtitle = "All 73 intersecting genes were lower in tumor") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")

fig1 <- f1a / (f1b | f1c | f1d) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 13)))
save_figure(fig1, "Figure1_discovery", 12.0, 7.2)

fwrite(deg_all[, .(gene, log2FoldChange, padj, status)], file.path(source_dir, "Figure1b_source_data.csv"))
fwrite(module_cor, file.path(source_dir, "Figure1c_source_data.csv"))
fwrite(set_counts, file.path(source_dir, "Figure1d_source_data.csv"))

# -------------------------
# Figure 2：跨队列复现
# -------------------------
bulk_scores <- fread(file.path(dc_dir, "bulk_deconvolution_sample_scores.csv"))
ext_scores <- fread(file.path(pv_dir, "external_bulk_paired_scores.csv"))
direction <- fread(file.path(pv_dir, "external_bulk_gene_direction_validation.csv"))
validation_summary <- fread(file.path(pv_dir, "external_bulk_validation_summary.csv"))

disc_scores <- bulk_scores[dataset == "GSE29272",
                           .(dataset, sample, group = Tissue, pair_id, module_score)]
score_all <- rbindlist(list(disc_scores, ext_scores), use.names = TRUE)
score_all[, group := factor(group, levels = c("Normal", "Tumor"))]
score_all[, dataset := factor(dataset, levels = c("GSE29272", "GSE79973", "GSE19826"))]

f2a <- ggplot(score_all, aes(group, module_score, group = pair_id)) +
  geom_line(colour = "#C6C6C6", linewidth = 0.3, alpha = 0.7) +
  geom_point(aes(colour = group), size = 1.2, alpha = 0.8) +
  facet_wrap(~dataset, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("Normal" = pal[["normal"]], "Tumor" = pal[["tumor"]])) +
  labs(x = NULL, y = "73-gene module score",
       subtitle = "Patient-paired bulk cohorts") +
  theme_nature() +
  theme(legend.position = "top")

f2b <- ggplot(direction, aes(discovery_log2FC, validation_mean_paired_difference)) +
  geom_hline(yintercept = 0, colour = "#888888", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.3) +
  geom_point(colour = pal["purple"], alpha = 0.75, size = 1.3) +
  geom_smooth(method = "lm", se = FALSE, colour = pal["dark"], linewidth = 0.55) +
  facet_wrap(~dataset, nrow = 1) +
  labs(x = "Discovery log2 fold change",
       y = "Validation paired difference",
       subtitle = "72/72 assayed genes were directionally concordant in each cohort") +
  theme_nature()

validation_summary[, dataset := factor(dataset, levels = c("GSE79973", "GSE19826"))]
f2c <- ggplot(validation_summary,
              aes(dataset, mean_tumor_minus_normal_score, colour = dataset)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.35) +
  geom_errorbar(aes(ymin = score_difference_ci_low, ymax = score_difference_ci_high),
                width = 0.12, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("GSE79973" = pal[["purple"]], "GSE19826" = pal[["tumor"]])) +
  labs(x = NULL, y = "Mean paired score difference\n(tumor − normal)",
       subtitle = "Points show mean; bars show 95% CI") +
  theme_nature() +
  theme(legend.position = "none")

top_gene <- direction[, .(mean_abs = mean(abs(discovery_log2FC))), by = gene][order(-mean_abs)][1:12, gene]
heat_dt <- dcast(direction[gene %in% top_gene], gene ~ dataset,
                 value.var = "validation_mean_paired_difference")
heat_dt <- merge(deg_all[, .(gene, GSE29272 = log2FoldChange)], heat_dt, by = "gene")
heat_long <- melt(heat_dt, id.vars = "gene", variable.name = "dataset", value.name = "effect")
heat_long[, gene := factor(gene, levels = rev(top_gene))]
f2d <- ggplot(heat_long, aes(dataset, gene, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = pal[["blue"]], mid = "white", high = pal[["tumor"]], midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Effect",
       subtitle = "Largest discovery effects") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

fig2 <- (f2a | f2c) / (f2b | f2d) +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 13)))
save_figure(fig2, "Figure2_bulk_replication", 12.0, 8.2)

fwrite(score_all, file.path(source_dir, "Figure2a_source_data.csv"))
fwrite(validation_summary, file.path(source_dir, "Figure2c_source_data.csv"))
fwrite(direction, file.path(source_dir, "Figure2b_source_data.csv"))
fwrite(heat_long, file.path(source_dir, "Figure2d_source_data.csv"))

# -------------------------
# Figure 3：单细胞定位与患者级检验
# -------------------------
type_scores <- fread(file.path(sc_dir, "singlecell_type_pseudobulk_scores.csv"))
epi_scores <- fread(file.path(sc_dir, "epithelial_patient_pseudobulk_scores.csv"))
lineage_tests <- fread(file.path(sc_dir, "epithelial_lineage_paired_tests.csv"))
sc_summary <- fread(file.path(sc_dir, "singlecell_analysis_summary.csv"))
sc_prop <- fread(file.path(dc_dir, "figure_sc_proportion_source_data.csv"))

type_order <- type_scores[, .(med = median(module_score, na.rm = TRUE)), by = Type][order(med)]$Type
type_scores[, Type := factor(Type, levels = type_order)]
f3a <- ggplot(type_scores, aes(module_score, Type, colour = Tissue)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35, width = 0.6) +
  geom_jitter(size = 0.65, alpha = 0.45, height = 0.13, width = 0) +
  scale_colour_manual(values = c("Normal" = pal[["normal"]], "Tumor" = pal[["tumor"]])) +
  labs(x = "73-gene pseudobulk score", y = NULL,
       subtitle = "Author-annotated cell types; patient–sample pseudobulks") +
  theme_nature() +
  theme(legend.position = "top")

epi_pair <- epi_scores[eligible == TRUE]
f3b <- ggplot(epi_pair, aes(Tissue, module_score, group = Patient)) +
  geom_line(colour = "#BEBEBE", linewidth = 0.45) +
  geom_point(aes(colour = Tissue), size = 2) +
  scale_colour_manual(values = c("Normal" = pal[["normal"]], "Tumor" = pal[["tumor"]])) +
  labs(x = NULL, y = "Epithelial pseudobulk score",
       subtitle = "10 complete pairs; paired t-test P = 0.907") +
  theme_nature() +
  theme(legend.position = "none")

epi_wide <- dcast(epi_pair, Patient ~ Tissue, value.var = "module_score")
epi_tt <- t.test(epi_wide$Tumor, epi_wide$Normal, paired = TRUE)
get_metric <- function(x) as.numeric(sc_summary[metric == x, value])
model_dt <- data.table(
  model = c("Paired epithelial pseudobulk", "Patient + lineage adjusted"),
  estimate = c(mean(epi_wide$Tumor - epi_wide$Normal),
               get_metric("lineage_adjusted_tissue_coefficient")),
  ci_low = c(unname(epi_tt$conf.int[1]),
             get_metric("lineage_adjusted_tissue_coefficient") -
               1.96 * get_metric("lineage_adjusted_cluster_robust_se")),
  ci_high = c(unname(epi_tt$conf.int[2]),
              get_metric("lineage_adjusted_tissue_coefficient") +
                1.96 * get_metric("lineage_adjusted_cluster_robust_se")),
  pvalue = c(epi_tt$p.value, get_metric("lineage_adjusted_cluster_robust_pvalue"))
)
model_dt[, model := factor(model, levels = rev(model))]
f3c <- ggplot(model_dt, aes(estimate, model)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.14,
                colour = pal["blue"], linewidth = 0.7) +
  geom_point(colour = pal["blue"], size = 2.8) +
  geom_text(aes(label = sprintf("P = %.3f", pvalue)), nudge_y = 0.18,
            hjust = 0, size = 2.8) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.35))) +
  labs(x = "Tumor coefficient for epithelial module score", y = NULL,
       subtitle = "Estimates with 95% confidence intervals") +
  theme_nature() +
  theme(legend.position = "none")

prop_pair <- sc_prop[, .(Patient, Tissue, epithelial_proportion)]
f3d <- ggplot(prop_pair, aes(Tissue, epithelial_proportion, group = Patient)) +
  geom_line(colour = "#BEBEBE", linewidth = 0.4) +
  geom_point(aes(colour = Tissue), size = 1.8) +
  scale_colour_manual(values = c("Normal" = pal[["normal"]], "Tumor" = pal[["tumor"]])) +
  labs(x = NULL, y = "Epithelial-cell proportion",
       subtitle = "20 complete pairs; paired t-test P = 0.0177") +
  theme_nature() +
  theme(legend.position = "none")

fig3 <- (f3a | (f3b / f3d)) / f3c +
  plot_layout(heights = c(1.7, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 13)))
save_figure(fig3, "Figure3_single_cell_localization", 12.0, 9.0)

fwrite(type_scores, file.path(source_dir, "Figure3a_source_data.csv"))
fwrite(epi_pair, file.path(source_dir, "Figure3b_source_data.csv"))
fwrite(model_dt, file.path(source_dir, "Figure3c_source_data.csv"))
fwrite(lineage_tests, file.path(source_dir, "Supplementary_lineage_tests.csv"))
fwrite(prop_pair, file.path(source_dir, "Figure3d_source_data.csv"))

# -------------------------
# Figure 4：组成拆分与bulk校正
# -------------------------
decomp <- fread(file.path(dc_dir, "figure_sc_decomposition_source_data.csv"))
bulk_scatter <- fread(file.path(dc_dir, "figure_bulk_scatter_source_data.csv"))
atten <- fread(file.path(dc_dir, "figure_bulk_attenuation_source_data.csv"))
dc_summary <- fread(file.path(dc_dir, "bulk_deconvolution_adjustment_summary.csv"))

f4a <- ggplot(sc_prop, aes(Tissue, epithelial_proportion, group = Patient)) +
  geom_line(colour = "#C0C0C0", linewidth = 0.4) +
  geom_point(aes(colour = Tissue), size = 1.8) +
  scale_colour_manual(values = c("Normal" = pal[["normal"]], "Tumor" = pal[["tumor"]])) +
  labs(x = NULL, y = "Epithelial-cell proportion",
       subtitle = "GSE206785 whole-tissue cell composition") +
  theme_nature() +
  theme(legend.position = "none")

decomp[, component := factor(component, levels = c("Observed", "Composition expected", "Residual"))]
f4b <- ggplot(decomp, aes(component, difference, fill = component)) +
  geom_hline(yintercept = 0, colour = "#666666", linewidth = 0.35) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = unname(c(pal["blue"], pal["olive"], pal["purple"]))) +
  labs(x = NULL, y = "Tumor − normal module score",
       subtitle = "Observed, composition-expected and residual differences") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

f4c <- ggplot(bulk_scatter, aes(gastric_lineage_deconv_score, module_score, colour = Tissue)) +
  geom_point(size = 0.9, alpha = 0.65) +
  geom_smooth(method = "lm", se = FALSE, colour = pal["dark"], linewidth = 0.5) +
  facet_wrap(~dataset, scales = "free", nrow = 1) +
  scale_colour_manual(values = c("Normal" = pal[["normal"]], "Tumor" = pal[["tumor"]])) +
  labs(x = "Independent gastric-lineage score", y = "73-gene module score",
       subtitle = "Candidate genes were excluded from all marker panels") +
  theme_nature() +
  theme(legend.position = "top")

atten[, model := factor(model, levels = c("unadjusted", "adjusted"))]
f4d <- ggplot(atten, aes(model, tissue_coefficient, group = dataset, colour = dataset)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.35) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.7) +
  scale_colour_manual(values = c("GSE29272" = pal[["blue"]], "GSE79973" = pal[["purple"]],
                                 "GSE19826" = pal[["tumor"]])) +
  labs(x = NULL, y = "Tumor coefficient",
       subtitle = "Adjustment for independent mature gastric-lineage scores") +
  theme_nature() +
  theme(legend.position = "top")

fig4 <- (f4a | f4b) / (f4c | f4d) +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 13)))
save_figure(fig4, "Figure4_composition_adjustment", 12.0, 8.2)

fwrite(sc_prop, file.path(source_dir, "Figure4a_source_data.csv"))
fwrite(decomp, file.path(source_dir, "Figure4b_source_data.csv"))
fwrite(bulk_scatter, file.path(source_dir, "Figure4c_source_data.csv"))
fwrite(atten, file.path(source_dir, "Figure4d_source_data.csv"))
fwrite(dc_summary, file.path(source_dir, "Figure4_adjustment_statistics.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "10_sessionInfo.txt"), useBytes = TRUE)
cat("完成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
cat("已生成4张主图（PDF + 600 dpi PNG）及逐面板源数据。\n")
