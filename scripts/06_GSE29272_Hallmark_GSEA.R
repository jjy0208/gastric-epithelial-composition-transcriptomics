# GSE29272 MSigDB Hallmark preranked GSEA
# 物种：Homo sapiens（人类）
# 排序指标：Disease vs Control 的 log2FoldChange，从大到小排序
# 方法：clusterProfiler::GSEA，fgseaMultilevel 引擎，BH 方法控制 FDR
# 说明：正 NES 表示疾病组方向富集；负 NES 表示对照方向富集/疾病组相对抑制

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(fgsea)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
  library(ragg)
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
deg_path <- file.path(
  project_dir, "results", "DEG", "GSE29272_DEG_all.csv"
)
gsea_dir <- file.path(project_dir, "results", "GSEA")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(gsea_dir, "GSE29272_GSEA_run.log")
session_path <- file.path(gsea_dir, "GSE29272_GSEA_sessionInfo.txt")
report_path <- file.path(gsea_dir, "GSE29272_GSEA_report.md")

log_con <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("开始时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
cat("数据集：GSE29272\n")
cat("比较：Disease vs Control\n")
cat("物种：Homo sapiens（人类）\n")
cat("基因集：MSigDB Hallmark collection H\n")
cat("排序指标：log2FoldChange（从大到小）\n")
cat("GSEA引擎：multilevel；exponent=1；P值数值下限=1e-10；RNG seed=20260729\n\n")

required_packages <- c(
  "clusterProfiler", "fgsea", "msigdbr", "enrichplot",
  "ggplot2", "patchwork", "dplyr", "stringr", "ragg"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少必需 R 包：", pkg)
  }
}
if (!file.exists(deg_path)) {
  stop("全部基因差异分析结果不存在：", deg_path)
}

de_all <- read.csv(deg_path, check.names = FALSE, encoding = "UTF-8")
required_columns <- c("gene", "log2FoldChange", "PValue", "padj")
if (!all(required_columns %in% colnames(de_all))) {
  stop(
    "差异分析结果缺少字段：",
    paste(setdiff(required_columns, colnames(de_all)), collapse = ", ")
  )
}
if (nrow(de_all) == 0L) stop("全部基因差异分析结果为空。")
if (anyDuplicated(de_all$gene) > 0L) stop("全部基因表中存在重复 Gene Symbol。")

valid <- (
  !is.na(de_all$gene) & nzchar(de_all$gene) &
    is.finite(de_all$log2FoldChange)
)
ranking_df <- de_all[valid, c("gene", "log2FoldChange"), drop = FALSE]
if (nrow(ranking_df) < 1000L) {
  stop("有效排序基因少于1000个，不适合本次全基因Hallmark GSEA。")
}
if (anyDuplicated(ranking_df$gene) > 0L) {
  stop("有效排序列表中存在重复 Gene Symbol。")
}

gene_list <- ranking_df$log2FoldChange
names(gene_list) <- ranking_df$gene
gene_list <- sort(gene_list, decreasing = TRUE)
if (!identical(gene_list, sort(gene_list, decreasing = TRUE))) {
  stop("基因排序向量未严格按log2FoldChange降序排列。")
}

ranking_output <- data.frame(
  rank = seq_along(gene_list),
  gene_symbol = names(gene_list),
  log2FoldChange = unname(gene_list)
)
write.csv(
  ranking_output,
  file.path(gsea_dir, "GSE29272_GSEA_ranked_gene_list.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 首次运行时从 msigdbr 获取人类 Hallmark H 集合并保存快照；
# 以后优先读取快照，避免MSigDB版本更新导致结果漂移。
hallmark_snapshot_path <- file.path(
  gsea_dir, "MSigDB_Hallmark_Homo_sapiens_snapshot.csv"
)
if (file.exists(hallmark_snapshot_path)) {
  hallmark <- read.csv(
    hallmark_snapshot_path,
    check.names = FALSE,
    encoding = "UTF-8"
  )
  hallmark_source <- "读取项目中已有的Hallmark快照"
} else {
  hallmark <- msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "H"
  )
  hallmark <- as.data.frame(hallmark)
  write.csv(
    hallmark,
    hallmark_snapshot_path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  hallmark_source <- "通过msigdbr获取并保存Hallmark快照"
}

required_hallmark_columns <- c(
  "gs_name", "gene_symbol", "gs_collection",
  "gs_source_species", "db_version", "db_target_species"
)
if (!all(required_hallmark_columns %in% colnames(hallmark))) {
  stop(
    "Hallmark快照缺少字段：",
    paste(
      setdiff(required_hallmark_columns, colnames(hallmark)),
      collapse = ", "
    )
  )
}
if (!all(hallmark$gs_collection == "H")) {
  stop("基因集快照包含非Hallmark集合。")
}

msigdb_version <- paste(unique(hallmark$db_version), collapse = ", ")
source_species <- paste(unique(hallmark$gs_source_species), collapse = ", ")
target_species <- paste(unique(hallmark$db_target_species), collapse = ", ")
hallmark_sets <- length(unique(hallmark$gs_name))
hallmark_genes <- length(unique(hallmark$gene_symbol))

term2gene <- unique(
  hallmark[, c("gs_name", "gene_symbol"), drop = FALSE]
)
colnames(term2gene) <- c("term", "gene")
term2name <- unique(
  hallmark[, c("gs_name", "gs_description"), drop = FALSE]
)
colnames(term2name) <- c("term", "name")

overlap_genes <- intersect(names(gene_list), unique(term2gene$gene))
ranking_overlap_rate <- length(overlap_genes) / length(gene_list)
if (length(overlap_genes) < 1000L) {
  stop("排序基因与Hallmark成员的重叠基因少于1000个。")
}

mapping_summary <- data.frame(
  input_all_genes = nrow(de_all),
  valid_unique_ranked_symbols = length(gene_list),
  hallmark_gene_sets = hallmark_sets,
  unique_hallmark_symbols = hallmark_genes,
  ranked_symbols_in_hallmark = length(overlap_genes),
  ranked_symbol_overlap_rate = ranking_overlap_rate,
  msigdb_version = msigdb_version,
  source_species = source_species,
  target_species = target_species
)
write.csv(
  mapping_summary,
  file.path(gsea_dir, "GSE29272_GSEA_input_mapping_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("全部基因行数：", nrow(de_all), "\n", sep = "")
cat("有效唯一排序基因：", length(gene_list), "\n", sep = "")
cat("排序并列数量：", sum(duplicated(unname(gene_list))), "\n", sep = "")
cat("Hallmark基因集数量：", hallmark_sets, "\n", sep = "")
cat("Hallmark唯一Gene Symbol：", hallmark_genes, "\n", sep = "")
cat("排序列表与Hallmark重叠基因：", length(overlap_genes), "\n", sep = "")
cat("MSigDB版本：", msigdb_version, "\n", sep = "")
cat("Hallmark来源：", hallmark_source, "\n\n", sep = "")

set.seed(20260729)
gsea_object <- GSEA(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  nPerm = 1000,
  method = "multilevel",
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)

gsea_all <- as.data.frame(gsea_object)
if (nrow(gsea_all) == 0L) {
  stop("Hallmark GSEA未返回任何可检验基因集。")
}

clean_hallmark_name <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  stringr::str_to_title(tolower(x))
}

leading_edge_size <- function(x) {
  if (is.na(x) || !nzchar(x)) return(0L)
  length(unique(strsplit(x, "/", fixed = TRUE)[[1]]))
}

gsea_all$pathway <- clean_hallmark_name(gsea_all$ID)
gsea_all$direction <- ifelse(
  gsea_all$NES > 0,
  "Activated in disease",
  "Suppressed in disease"
)
gsea_all$leading_edge_size <- vapply(
  gsea_all$core_enrichment,
  leading_edge_size,
  integer(1)
)
gsea_all <- gsea_all[
  order(gsea_all$p.adjust, -abs(gsea_all$NES), gsea_all$ID),
  ,
  drop = FALSE
]
gsea_all <- gsea_all[
  c(
    "pathway", "ID", "Description", "setSize", "enrichmentScore",
    "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge",
    "leading_edge_size", "core_enrichment", "direction"
  )
]
rownames(gsea_all) <- NULL

gsea_sig <- gsea_all[
  !is.na(gsea_all$p.adjust) & gsea_all$p.adjust < 0.05,
  ,
  drop = FALSE
]
gsea_activated <- gsea_sig[gsea_sig$NES > 0, , drop = FALSE]
gsea_suppressed <- gsea_sig[gsea_sig$NES < 0, , drop = FALSE]

write.csv(
  gsea_sig,
  file.path(gsea_dir, "GSE29272_GSEA_result.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  gsea_all,
  file.path(gsea_dir, "GSE29272_GSEA_all_tested.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  gsea_activated,
  file.path(gsea_dir, "GSE29272_GSEA_activated.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  gsea_suppressed,
  file.path(gsea_dir, "GSE29272_GSEA_suppressed.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary_table <- data.frame(
  dataset = "GSE29272",
  species = "Homo sapiens",
  collection = "MSigDB Hallmark H",
  msigdb_version = msigdb_version,
  ranking_metric = "log2FoldChange",
  ranked_genes = length(gene_list),
  tested_gene_sets = nrow(gsea_all),
  significant_gene_sets_padj_lt_0.05 = nrow(gsea_sig),
  activated_NES_positive = nrow(gsea_activated),
  suppressed_NES_negative = nrow(gsea_suppressed),
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalue_numerical_floor = 1e-10,
  seed = 20260729,
  engine = "clusterProfiler/enrichit multilevel",
  p_adjust_method = "BH"
)
write.csv(
  summary_table,
  file.path(gsea_dir, "GSE29272_GSEA_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 图形契约：
# 核心结论：Hallmark通路沿疾病vs对照的全基因log2FC排序呈现方向性富集。
# 图形类型：quantitative grid；总览图为主图，前5条running-score图为验证图。
# 颜色：低饱和度蓝色表示疾病组抑制，低饱和度陶红表示疾病组激活。
color_activated <- "#B36F63"
color_suppressed <- "#7089A6"
color_neutral <- "#D8D2C8"

theme_gsea <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2E2E2E"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#2E2E2E"),
      axis.text = element_text(colour = "#333333"),
      axis.title = element_text(colour = "#222222"),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "#EAE6E0"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 1.2),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#5A5652"),
      plot.caption = element_text(
        size = base_size - 1.2,
        colour = "#6B6661",
        hjust = 0
      ),
      plot.margin = margin(6, 8, 6, 6)
    )
}

save_plot <- function(plot, stem, width_mm, height_mm, dpi = 600L) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "sans"
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    paste0(stem, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

top_overview <- head(gsea_sig, 5L)
if (nrow(top_overview) > 0L) {
  top_overview$minus_log10_padj <- -log10(
    pmax(top_overview$p.adjust, .Machine$double.xmin)
  )
  top_overview$pathway_plot <- factor(
    top_overview$pathway,
    levels = rev(top_overview$pathway)
  )
  overview_plot <- ggplot(
    top_overview,
    aes(
      x = NES,
      y = pathway_plot,
      size = minus_log10_padj,
      colour = direction
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      colour = "#77736E"
    ) +
    geom_point(alpha = 0.95) +
    scale_colour_manual(
      values = c(
        "Activated in disease" = color_activated,
        "Suppressed in disease" = color_suppressed
      ),
      breaks = c("Activated in disease", "Suppressed in disease"),
      labels = c("Activated in disease", "Suppressed in disease"),
      name = "Direction"
    ) +
    scale_size_continuous(
      range = c(3.2, 8),
      name = expression(-log[10]("adjusted P"))
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 34)) +
    labs(
      title = "Hallmark GSEA overview",
      subtitle = "Five most significant pathways ordered by BH-adjusted P value",
      x = "Normalized enrichment score (NES)",
      y = NULL,
      caption = paste0(
        "Positive NES: enriched toward disease. Negative NES: enriched toward control.\n",
        "Human MSigDB Hallmark ", msigdb_version,
        "; ranking metric: log2FoldChange; BH-adjusted P < 0.05."
      )
    ) +
    theme_gsea() +
    theme(legend.position = "right")
} else {
  overview_plot <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "BH-adjusted P < 0.05 时无显著Hallmark通路",
      colour = "#655D58", size = 3.2
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = "Hallmark GSEA overview") +
    theme_void(base_family = "sans")
}

save_plot(
  overview_plot,
  file.path(gsea_dir, "GSE29272_GSEA_overview_dotplot"),
  width_mm = 175,
  height_mm = 105
)
write.csv(
  top_overview,
  file.path(gsea_dir, "GSE29272_GSEA_overview_source_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 为最显著的前5条通路分别输出经典running-score图。
safe_file_component <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  substr(x, 1, 80)
}

curve_index <- data.frame()
top_curve_rows <- head(gsea_sig, 5L)
if (nrow(top_curve_rows) > 0L) {
  for (i in seq_len(nrow(top_curve_rows))) {
    pathway_id <- top_curve_rows$ID[i]
    direction_label <- top_curve_rows$direction[i]
    curve_color <- if (
      direction_label == "Activated in disease"
    ) color_activated else color_suppressed
    curve_title <- paste0(
      top_curve_rows$pathway[i],
      "  |  NES=", sprintf("%.2f", top_curve_rows$NES[i]),
      ", padj=", formatC(
        top_curve_rows$p.adjust[i],
        format = "e",
        digits = 2
      )
    )
    stem_name <- paste0(
      "GSE29272_GSEA_curve_",
      sprintf("%02d", i),
      "_",
      safe_file_component(pathway_id)
    )

    # 按经典加权GSEA公式计算逐基因running score。
    # 命中基因按|log2FC|^exponent归一化加分，非命中基因等权扣分。
    pathway_genes <- unique(
      term2gene$gene[term2gene$term == pathway_id]
    )
    hit <- names(gene_list) %in% pathway_genes
    n_ranked <- length(gene_list)
    n_hits <- sum(hit)
    if (n_hits == 0L || n_hits == n_ranked) {
      stop("无法为通路计算running score：", pathway_id)
    }
    hit_weights <- abs(gene_list)^1
    hit_norm <- sum(hit_weights[hit])
    running_increment <- ifelse(
      hit,
      hit_weights / hit_norm,
      -1 / (n_ranked - n_hits)
    )
    running_score <- cumsum(running_increment)
    peak_rank <- if (
      top_curve_rows$NES[i] > 0
    ) which.max(running_score) else which.min(running_score)

    curve_data <- data.frame(
      rank = seq_along(gene_list),
      gene_symbol = names(gene_list),
      log2FoldChange = unname(gene_list),
      pathway_hit = hit,
      running_score = running_score
    )
    source_file <- paste0(stem_name, "_source_data.csv")
    write.csv(
      curve_data,
      file.path(gsea_dir, source_file),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    p_running <- ggplot(
      curve_data,
      aes(x = rank, y = running_score)
    ) +
      geom_hline(
        yintercept = 0,
        linewidth = 0.35,
        colour = "#77736E"
      ) +
      geom_vline(
        xintercept = peak_rank,
        linetype = "dashed",
        linewidth = 0.35,
        colour = "#9A958E"
      ) +
      geom_line(
        linewidth = 0.75,
        colour = curve_color,
        lineend = "round"
      ) +
      labs(
        title = stringr::str_wrap(curve_title, width = 75),
        x = NULL,
        y = "Running enrichment score"
      ) +
      theme_gsea(base_size = 8) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(5, 6, 0, 6)
      )

    p_hits <- ggplot(
      curve_data[curve_data$pathway_hit, , drop = FALSE],
      aes(x = rank)
    ) +
      geom_segment(
        aes(xend = rank, y = 0, yend = 1),
        linewidth = 0.25,
        colour = "#252525"
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      theme_classic(base_size = 8, base_family = "sans") +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0, 6, 0, 6)
      )

    p_metric <- ggplot(
      curve_data,
      aes(x = rank, y = log2FoldChange)
    ) +
      geom_hline(
        yintercept = 0,
        linewidth = 0.3,
        colour = "#77736E"
      ) +
      geom_area(
        data = curve_data[curve_data$log2FoldChange >= 0, , drop = FALSE],
        fill = color_activated,
        alpha = 0.45
      ) +
      geom_area(
        data = curve_data[curve_data$log2FoldChange < 0, , drop = FALSE],
        fill = color_suppressed,
        alpha = 0.45
      ) +
      labs(
        x = "Rank in ordered dataset",
        y = "log2FoldChange"
      ) +
      theme_gsea(base_size = 8) +
      theme(
        panel.grid.major.x = element_blank(),
        plot.margin = margin(0, 6, 5, 6)
      )

    curve_plot <- p_running / p_hits / p_metric +
      patchwork::plot_layout(heights = c(1.7, 0.35, 1))

    save_plot(
      curve_plot,
      file.path(gsea_dir, stem_name),
      width_mm = 165,
      height_mm = 115
    )
    curve_index <- bind_rows(
      curve_index,
      data.frame(
        rank = i,
        pathway = top_curve_rows$pathway[i],
        ID = pathway_id,
        NES = top_curve_rows$NES[i],
        pvalue = top_curve_rows$pvalue[i],
        padj = top_curve_rows$p.adjust[i],
        direction = direction_label,
        pdf_file = paste0(stem_name, ".pdf"),
        png_file = paste0(stem_name, ".png"),
        source_data_file = source_file
      )
    )
  }
}
write.csv(
  curve_index,
  file.path(gsea_dir, "GSE29272_GSEA_curve_file_index.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

format_p <- function(x) formatC(x, format = "e", digits = 2)
markdown_top_table <- function(df, n = 5L) {
  if (nrow(df) == 0L) return("未检出符合条件的通路。\n")
  z <- head(df, n)
  lines <- c(
    "| 排名 | 通路 | NES | P值 | padj | 核心基因数 |",
    "|---:|---|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(z))) {
    lines <- c(
      lines,
      paste0(
        "| ", i,
        " | ", gsub("\\|", "/", z$pathway[i]),
        " | ", sprintf("%.3f", z$NES[i]),
        " | ", format_p(z$pvalue[i]),
        " | ", format_p(z$p.adjust[i]),
        " | ", z$leading_edge_size[i],
        " |"
      )
    )
  }
  paste(lines, collapse = "\n")
}

sig_ids <- toupper(gsea_sig$ID)
is_positive <- gsea_sig$NES > 0
is_negative <- gsea_sig$NES < 0
interpretation <- character(0)

if (any(grepl("EPITHELIAL_MESENCHYMAL_TRANSITION|ANGIOGENESIS", sig_ids) & is_positive)) {
  interpretation <- c(
    interpretation,
    "- **间质重塑与侵袭相关程序：** EMT或血管生成Hallmark在疾病方向呈正NES，支持胃癌组织中基质重塑、迁移和肿瘤间质反应增强。该信号也可能受到成纤维细胞和内皮细胞比例变化影响。"
  )
}
if (any(grepl("TNFA_SIGNALING_VIA_NFKB|INFLAMMATORY_RESPONSE|IL6_JAK_STAT3_SIGNALING|COMPLEMENT", sig_ids) & is_positive)) {
  interpretation <- c(
    interpretation,
    "- **炎症与细胞因子信号：** NF-κB、炎症反应、IL6-JAK-STAT3或补体相关Hallmark若为正NES，提示胃癌组织的炎症和免疫微环境整体向疾病方向偏移；混合组织数据不能区分肿瘤细胞内信号与免疫浸润。"
  )
}
if (any(grepl("E2F_TARGETS|G2M_CHECKPOINT|MYC_TARGETS|MITOTIC_SPINDLE", sig_ids) & is_positive)) {
  interpretation <- c(
    interpretation,
    "- **增殖与细胞周期：** E2F、G2M、MYC或有丝分裂相关Hallmark的正NES提示肿瘤组织具有更强的细胞周期和增殖转录程序。"
  )
}
if (any(grepl("OXIDATIVE_PHOSPHORYLATION|FATTY_ACID_METABOLISM|XENOBIOTIC_METABOLISM|BILE_ACID_METABOLISM", sig_ids) & is_negative)) {
  interpretation <- c(
    interpretation,
    "- **分化与代谢功能下降：** 氧化磷酸化、脂肪酸、异生物或胆汁酸代谢相关Hallmark的负NES表示这些程序更偏向对照胃腺体，可能反映胃癌中正常上皮分化和代谢功能丧失。"
  )
}
if (any(grepl("KRAS_SIGNALING_DN", sig_ids) & is_negative)) {
  interpretation <- c(
    interpretation,
    "- **KRAS命名注意：** `KRAS signaling DN` 是在KRAS激活实验中下调的基因集合；其负NES不能直接等价为“KRAS被抑制”，必须结合基因集定义和其他证据解释。"
  )
}
if (length(interpretation) == 0L && nrow(gsea_sig) > 0L) {
  interpretation <- paste0(
    "- 最显著通路包括：",
    paste(head(gsea_sig$pathway, 5), collapse = "、"),
    "。这些结果反映Hallmark成员在全基因排序两端的聚集，需要结合核心富集基因和独立队列验证。"
  )
}
if (nrow(gsea_sig) == 0L) {
  interpretation <- "- 在BH校正padj < 0.05时未检出显著Hallmark通路，不能据此声称存在整体激活或抑制。"
}

report_lines <- c(
  "# GSE29272 MSigDB Hallmark GSEA分析报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "- 数据集：GSE29272",
  "- 研究比较：胃癌肿瘤组织（贲门与非贲门合并）vs 配对癌旁正常胃腺体",
  "- 物种：Homo sapiens（人类）",
  paste0("- MSigDB版本：", msigdb_version),
  "- 基因集：Hallmark collection H",
  "",
  "## 1. 输入与排序",
  "",
  paste0("- 读取全部差异分析结果：", nrow(de_all), " 个基因。"),
  paste0("- 有效唯一Gene Symbol：", length(gene_list), " 个。"),
  "- 排序依据：Disease vs Control的log2FoldChange，从大到小。",
  "- 正值位于排序顶部，代表疾病组表达更高；负值位于底部，代表对照组表达更高。",
  paste0("- 排序列表与Hallmark成员重叠：", length(overlap_genes), " 个Gene Symbol。"),
  paste0("- Hallmark快照来源：", hallmark_source, "。"),
  "",
  "## 2. GSEA参数",
  "",
  "- 方法：clusterProfiler::GSEA，multilevel preranked引擎。",
  "- 排序权重指数：exponent = 1。",
  "- 基因集大小：10–500。",
  "- R随机种子：20260729。当前clusterProfiler 4.20接口由内部multilevel引擎读取该随机状态。",
  "- 当前接口的P值数值下限为1×10⁻¹⁰；达到该下限的通路应解释为P≤1×10⁻¹⁰，而不是精确等于该数值。",
  "- 多重检验：Benjamini-Hochberg；显著标准：padj < 0.05。",
  "- preranked GSEA使用基因置换近似，不保留基因间相关结构。",
  "",
  "## 3. 结果概览",
  "",
  paste0("- 实际检验Hallmark通路：", nrow(gsea_all), " 条。"),
  paste0("- 显著通路：", nrow(gsea_sig), " 条。"),
  paste0("- 疾病组方向富集（NES > 0）：", nrow(gsea_activated), " 条。"),
  paste0("- 对照方向富集/疾病组相对抑制（NES < 0）：", nrow(gsea_suppressed), " 条。"),
  "",
  "## 4. 疾病组显著激活通路（NES > 0）前5条",
  "",
  markdown_top_table(gsea_activated),
  "",
  "## 5. 疾病组显著抑制通路（NES < 0）前5条",
  "",
  markdown_top_table(gsea_suppressed),
  "",
  "## 6. 生物学解读",
  "",
  interpretation,
  "",
  "## 7. 结果边界与注意事项",
  "",
  "- 用户指定以log2FoldChange排序，本分析已严格执行。log2FoldChange不是方差校准统计量，理论上limma moderated t更适合作为GSEA排序指标；极端但不稳定的倍数变化可能影响leading edge。",
  "- 正NES表示通路成员集中于疾病组高表达端，负NES表示集中于对照高表达端；NES本身不等同于直接测量的通路活性。",
  "- preranked fgsea采用基因置换，可能低估共表达基因集的零假设方差。若用于论文核心结论，建议用表达矩阵和配对设计补充CAMERA敏感性分析。",
  "- 癌旁正常胃腺体不等同于健康志愿者胃组织，肿瘤邻近效应可能影响排序。",
  "- 芯片来自混合组织，GSEA信号可能同时反映细胞内调控和细胞组成变化。",
  "- 结果是探索性功能假设，需要独立队列或实验验证。",
  "",
  "## 8. 输出文件",
  "",
  "- `GSE29272_GSEA_result.csv`：padj < 0.05的显著Hallmark结果。",
  "- `GSE29272_GSEA_all_tested.csv`：全部受检Hallmark通路。",
  "- `GSE29272_GSEA_activated.csv`、`GSE29272_GSEA_suppressed.csv`：按NES方向拆分的显著结果。",
  "- `GSE29272_GSEA_ranked_gene_list.csv`：实际使用的全基因排序。",
  "- `GSE29272_GSEA_overview_dotplot.pdf/png`：最显著5条通路总览。",
  "- `GSE29272_GSEA_curve_*.pdf/png`：最显著前5条通路的经典running-score图。",
  "- `GSE29272_GSEA_curve_*_source_data.csv`：每张running-score图的逐基因源数据。",
  "- `GSE29272_GSEA_curve_file_index.csv`：曲线图与通路对应关系。",
  "- `MSigDB_Hallmark_Homo_sapiens_snapshot.csv`：本次Hallmark基因集快照。",
  "- `GSE29272_GSEA_run.log`：运行日志。",
  "- `GSE29272_GSEA_sessionInfo.txt`：R及软件包环境。",
  "- `scripts/06_GSE29272_Hallmark_GSEA.R`：完整可复现R脚本。"
)
writeLines(report_lines, report_path, useBytes = TRUE)

capture.output(sessionInfo(), file = session_path)

cat("检验Hallmark通路：", nrow(gsea_all), "\n", sep = "")
cat("显著通路：", nrow(gsea_sig), "\n", sep = "")
cat("NES正向：", nrow(gsea_activated), "\n", sep = "")
cat("NES负向：", nrow(gsea_suppressed), "\n", sep = "")
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
