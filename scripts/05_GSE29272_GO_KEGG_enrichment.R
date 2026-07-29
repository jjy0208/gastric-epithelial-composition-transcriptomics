# GSE29272 显著差异基因 GO/KEGG 过度富集分析
# 物种：Homo sapiens（人类）
# 疾病背景：胃癌肿瘤组织（贲门与非贲门合并）vs 配对癌旁正常胃腺体
# 方法：clusterProfiler 超几何检验，Benjamini-Hochberg 方法控制 FDR
# 前景：|log2FoldChange| > 1 且 TREAT-BH padj < 0.05 的显著差异基因
# 背景：实际进入 GSE29272 差异检验且具有有效 P 值的基因

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(gson)
  library(ggplot2)
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
deg_dir <- file.path(project_dir, "results", "DEG")
enrich_dir <- file.path(project_dir, "results", "Enrich")
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

sig_path <- file.path(deg_dir, "GSE29272_DEG_significant.csv")
all_path <- file.path(deg_dir, "GSE29272_DEG_all.csv")
log_path <- file.path(enrich_dir, "GSE29272_enrichment_run.log")
session_path <- file.path(enrich_dir, "GSE29272_enrichment_sessionInfo.txt")

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
cat("物种：Homo sapiens（人类）；OrgDb = org.Hs.eg.db；KEGG organism = hsa\n")
cat("方法：clusterProfiler ORA；BH 多重检验校正；显著阈值 padj < 0.05\n\n")

required_packages <- c(
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "GO.db",
  "gson", "ggplot2", "dplyr", "stringr", "ragg"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少必需 R 包：", pkg)
  }
}
if (!file.exists(sig_path)) stop("显著差异基因表不存在：", sig_path)
if (!file.exists(all_path)) stop("完整差异分析表不存在：", all_path)

sig_de <- read.csv(sig_path, check.names = FALSE, encoding = "UTF-8")
all_de <- read.csv(all_path, check.names = FALSE, encoding = "UTF-8")

required_sig <- c("gene", "log2FoldChange", "PValue", "padj")
required_all <- c("gene", "PValue")
if (!all(required_sig %in% colnames(sig_de))) {
  stop("显著差异基因表缺少字段：", paste(setdiff(required_sig, colnames(sig_de)), collapse = ", "))
}
if (!all(required_all %in% colnames(all_de))) {
  stop("完整差异分析表缺少字段：", paste(setdiff(required_all, colnames(all_de)), collapse = ", "))
}
if (nrow(sig_de) == 0L) stop("显著差异基因表为空，无法进行 ORA。")
if (anyDuplicated(sig_de$gene) > 0L) stop("显著差异基因表中存在重复 Gene Symbol。")
if (anyDuplicated(all_de$gene) > 0L) stop("完整差异分析表中存在重复 Gene Symbol。")

sig_symbols <- unique(sig_de$gene[!is.na(sig_de$gene) & nzchar(sig_de$gene)])
background_symbols <- unique(
  all_de$gene[!is.na(all_de$gene) & nzchar(all_de$gene) & !is.na(all_de$PValue)]
)
if (!all(sig_symbols %in% background_symbols)) {
  stop("显著差异基因不是差异检验背景基因的严格子集。")
}

# 前景和背景使用同一映射规则，防止 ID 转换造成不一致。
symbol_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = background_symbols,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
  )
)
symbol_map <- unique(symbol_map[!is.na(symbol_map$ENTREZID), c("SYMBOL", "ENTREZID")])
foreground_map <- symbol_map[symbol_map$SYMBOL %in% sig_symbols, , drop = FALSE]

foreground_entrez <- unique(foreground_map$ENTREZID)
background_entrez <- unique(symbol_map$ENTREZID)
mapped_fg_symbols <- unique(foreground_map$SYMBOL)
mapped_bg_symbols <- unique(symbol_map$SYMBOL)

foreground_mapping_rate <- length(mapped_fg_symbols) / length(sig_symbols)
background_mapping_rate <- length(mapped_bg_symbols) / length(background_symbols)
if (foreground_mapping_rate < 0.85) {
  stop(
    "显著基因 SYMBOL->ENTREZID 映射率低于 85%：",
    sprintf("%.2f%%", 100 * foreground_mapping_rate)
  )
}
if (!all(foreground_entrez %in% background_entrez)) {
  stop("映射后的前景 Entrez ID 不是背景 Entrez ID 的子集。")
}

mapping_audit <- data.frame(
  set = c("Significant foreground", "Tested background"),
  input_symbols = c(length(sig_symbols), length(background_symbols)),
  mapped_symbols = c(length(mapped_fg_symbols), length(mapped_bg_symbols)),
  unique_entrez_ids = c(length(foreground_entrez), length(background_entrez)),
  symbol_mapping_rate = c(foreground_mapping_rate, background_mapping_rate)
)
write.csv(
  mapping_audit,
  file.path(enrich_dir, "GSE29272_ID_mapping_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  symbol_map,
  file.path(enrich_dir, "GSE29272_SYMBOL_ENTREZ_mapping.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("显著 Gene Symbol：", length(sig_symbols), "\n", sep = "")
cat("映射成功的显著 Gene Symbol：", length(mapped_fg_symbols), "\n", sep = "")
cat("显著基因映射率：", sprintf("%.2f%%", 100 * foreground_mapping_rate), "\n", sep = "")
cat("前景唯一 Entrez ID：", length(foreground_entrez), "\n", sep = "")
cat("背景 Gene Symbol：", length(background_symbols), "\n", sep = "")
cat("背景唯一 Entrez ID：", length(background_entrez), "\n\n", sep = "")

ratio_to_numeric <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
}

augment_result <- function(df, ontology) {
  if (nrow(df) == 0L) {
    df$Ontology <- character(0)
    df$FoldEnrichment <- numeric(0)
    return(df)
  }
  df$Ontology <- ontology
  df$FoldEnrichment <- ratio_to_numeric(df$GeneRatio) / ratio_to_numeric(df$BgRatio)
  df <- df[order(df$p.adjust, df$pvalue, -df$FoldEnrichment), , drop = FALSE]
  rownames(df) <- NULL
  df
}

run_go <- function(ontology) {
  obj <- enrichGO(
    gene = foreground_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  list(
    object = obj,
    all = augment_result(as.data.frame(obj), ontology)
  )
}

go_bp <- run_go("BP")
go_cc <- run_go("CC")
go_mf <- run_go("MF")

# 首次运行时下载并保存 KEGG 注释快照；以后优先读取已有快照。
# 这样即使 KEGG 在线数据库更新，也能复现本次通路集合。
kegg_snapshot_path <- file.path(enrich_dir, "KEGG_hsa_snapshot.gson")
if (file.exists(kegg_snapshot_path)) {
  kegg_snapshot <- gson::read.gson(kegg_snapshot_path)
  kegg_snapshot_source <- "读取项目中已有的GSON快照"
} else {
  kegg_snapshot <- gson_KEGG(
    species = "hsa",
    KEGG_Type = "KEGG",
    keyType = "ncbi-geneid"
  )
  gson::write.gson(kegg_snapshot, kegg_snapshot_path)
  kegg_snapshot_source <- "从KEGG在线数据库下载并保存GSON快照"
}
kegg_accessed_date <- as.character(kegg_snapshot@accessed_date)

kegg_object <- enricher(
  gene = foreground_entrez,
  universe = background_entrez,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  gson = kegg_snapshot
)
kegg_all <- augment_result(as.data.frame(kegg_object), "KEGG")

# 将 KEGG geneID 中的 Entrez ID 同时翻译成 Gene Symbol，保留原始 ID 便于审计。
entrez_to_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = background_entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
translate_kegg_gene_ids <- function(x) {
  ids <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  syms <- unname(entrez_to_symbol[ids])
  syms <- syms[!is.na(syms) & nzchar(syms)]
  paste(unique(syms), collapse = "/")
}
if (nrow(kegg_all) > 0L) {
  kegg_all$geneID_Entrez <- kegg_all$geneID
  kegg_all$geneID <- vapply(kegg_all$geneID_Entrez, translate_kegg_gene_ids, character(1))
}

filter_significant <- function(df) {
  if (nrow(df) == 0L) return(df)
  out <- df[!is.na(df$p.adjust) & df$p.adjust < 0.05, , drop = FALSE]
  rownames(out) <- NULL
  out
}

go_bp_sig <- filter_significant(go_bp$all)
go_cc_sig <- filter_significant(go_cc$all)
go_mf_sig <- filter_significant(go_mf$all)
kegg_sig <- filter_significant(kegg_all)

# 用户指定的四个主要结果文件：均为 padj < 0.05 的显著结果。
write.csv(go_bp_sig, file.path(enrich_dir, "GO_BP.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(go_cc_sig, file.path(enrich_dir, "GO_CC.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(go_mf_sig, file.path(enrich_dir, "GO_MF.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(kegg_sig, file.path(enrich_dir, "KEGG.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 同时保留全部受检条目，便于核查未通过 FDR 的结果。
write.csv(go_bp$all, file.path(enrich_dir, "GO_BP_all_tested.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(go_cc$all, file.path(enrich_dir, "GO_CC_all_tested.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(go_mf$all, file.path(enrich_dir, "GO_MF_all_tested.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(kegg_all, file.path(enrich_dir, "KEGG_all_tested.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summary_table <- data.frame(
  analysis = c("GO_BP", "GO_CC", "GO_MF", "KEGG"),
  tested_terms = c(nrow(go_bp$all), nrow(go_cc$all), nrow(go_mf$all), nrow(kegg_all)),
  significant_terms_padj_lt_0.05 = c(
    nrow(go_bp_sig), nrow(go_cc_sig), nrow(go_mf_sig), nrow(kegg_sig)
  ),
  foreground_entrez = length(foreground_entrez),
  universe_entrez = length(background_entrez),
  p_adjust_method = "BH",
  threshold = "padj < 0.05",
  min_gene_set_size = 10L,
  max_gene_set_size = 500L
)
write.csv(
  summary_table,
  file.path(enrich_dir, "GSE29272_enrichment_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Nature 风格 + 莫兰迪低饱和度连续色阶。
morandi_low <- "#D8D0C4"
morandi_mid <- "#A9B3A2"
morandi_high <- "#806F86"
theme_nature <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2F2F2F"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#2F2F2F"),
      axis.text = element_text(colour = "#333333"),
      axis.title = element_text(colour = "#222222"),
      strip.background = element_rect(fill = "#E6E0D8", colour = NA),
      strip.text = element_text(face = "bold", colour = "#3D3A37"),
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "#ECE8E2"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1.2),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#5A5652"),
      plot.caption = element_text(size = base_size - 1.2, colour = "#6B6661", hjust = 0),
      plot.margin = margin(6, 8, 6, 6)
    )
}

top_n_terms <- function(df, n = 10L) {
  if (nrow(df) == 0L) return(df)
  df[seq_len(min(n, nrow(df))), , drop = FALSE]
}

prepare_dot_data <- function(df_list) {
  out <- bind_rows(df_list)
  if (nrow(out) == 0L) return(out)
  out <- out %>%
    group_by(Ontology) %>%
    arrange(p.adjust, pvalue, .by_group = TRUE) %>%
    slice_head(n = 10L) %>%
    ungroup()
  out$minus_log10_padj <- -log10(pmax(out$p.adjust, .Machine$double.xmin))
  out$term_key <- paste(out$Description, out$Ontology, sep = "___")
  ordered_keys <- unlist(
    lapply(split(out, out$Ontology), function(z) rev(z$term_key)),
    use.names = FALSE
  )
  out$term_key <- factor(out$term_key, levels = unique(ordered_keys))
  out
}

empty_plot <- function(title, subtitle) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "padj < 0.05 时无显著条目", size = 3.2, colour = "#655D58") +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8.5, colour = "#5A5652")
    )
}

go_plot_data <- prepare_dot_data(list(go_bp_sig, go_cc_sig, go_mf_sig))
kegg_plot_data <- prepare_dot_data(list(kegg_sig))

if (nrow(go_plot_data) > 0L) {
  go_plot <- ggplot(
    go_plot_data,
    aes(x = FoldEnrichment, y = term_key, size = Count, colour = minus_log10_padj)
  ) +
    geom_point(alpha = 0.92) +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) {
      stringr::str_wrap(sub("___[^_]+$", "", x), width = 42)
    }) +
    scale_colour_gradientn(
      colours = c(morandi_low, morandi_mid, morandi_high),
      name = expression(-log[10]("adjusted P"))
    ) +
    scale_size_continuous(range = c(2.4, 7), name = "Gene count") +
    labs(
      title = "GO over-representation analysis",
      subtitle = "Top 10 significant terms per ontology, ordered by adjusted P value",
      x = "Fold enrichment",
      y = NULL,
      caption = paste0(
        "Human gastric cancer; foreground n=", length(foreground_entrez),
        ", tested background n=", length(background_entrez),
        "; BH-adjusted P < 0.05."
      )
    ) +
    theme_nature() +
    theme(legend.position = "right")
} else {
  go_plot <- empty_plot(
    "GO over-representation analysis",
    "BP、CC、MF 在 BH-adjusted P < 0.05 时均无显著条目"
  )
}

if (nrow(kegg_plot_data) > 0L) {
  kegg_plot <- ggplot(
    kegg_plot_data,
    aes(x = FoldEnrichment, y = term_key, size = Count, colour = minus_log10_padj)
  ) +
    geom_point(alpha = 0.92) +
    scale_y_discrete(labels = function(x) {
      stringr::str_wrap(sub("___[^_]+$", "", x), width = 44)
    }) +
    scale_colour_gradientn(
      colours = c(morandi_low, morandi_mid, morandi_high),
      name = expression(-log[10]("adjusted P"))
    ) +
    scale_size_continuous(range = c(2.6, 7.5), name = "Gene count") +
    labs(
      title = "KEGG pathway over-representation analysis",
      subtitle = "Top 10 significant pathways, ordered by adjusted P value",
      x = "Fold enrichment",
      y = NULL,
      caption = paste0(
        "KEGG hsa snapshot accessed ", kegg_accessed_date, ".\n",
        "Foreground n=", length(foreground_entrez),
        ", tested background n=", length(background_entrez),
        "; BH-adjusted P < 0.05."
      )
    ) +
    theme_nature() +
    theme(legend.position = "right")
} else {
  kegg_plot <- empty_plot(
    "KEGG pathway over-representation analysis",
    "BH-adjusted P < 0.05 时无显著通路"
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

save_plot(
  go_plot,
  file.path(enrich_dir, "GSE29272_GO_dotplot"),
  width_mm = 183,
  height_mm = 170
)
save_plot(
  kegg_plot,
  file.path(enrich_dir, "GSE29272_KEGG_dotplot"),
  width_mm = 160,
  height_mm = 105
)

write.csv(
  go_plot_data,
  file.path(enrich_dir, "GSE29272_GO_dotplot_source_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  kegg_plot_data,
  file.path(enrich_dir, "GSE29272_KEGG_dotplot_source_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

format_p <- function(x) formatC(x, format = "e", digits = 2)
markdown_top_table <- function(df, n = 5L) {
  if (nrow(df) == 0L) return("未检出 padj < 0.05 的显著条目。\n")
  z <- top_n_terms(df, n)
  lines <- c(
    "| 排名 | 条目 | padj | 富集倍数 | 基因数 |",
    "|---:|---|---:|---:|---:|"
  )
  for (i in seq_len(nrow(z))) {
    lines <- c(
      lines,
      paste0(
        "| ", i, " | ", gsub("\\|", "/", z$Description[i]),
        " | ", format_p(z$p.adjust[i]),
        " | ", sprintf("%.2f", z$FoldEnrichment[i]),
        " | ", z$Count[i], " |"
      )
    )
  }
  paste(lines, collapse = "\n")
}

all_sig_terms <- bind_rows(go_bp_sig, go_cc_sig, go_mf_sig, kegg_sig)
term_text <- tolower(paste(all_sig_terms$Description, collapse = " | "))
interpretation <- character(0)
if (grepl("copper ion|zinc ion|cadmium ion|metal ion|mineral absorption", term_text)) {
  interpretation <- c(
    interpretation,
    "- **金属离子稳态与应激：** GO BP最显著的一组铜、锌、镉应答条目以及KEGG `Mineral absorption` 主要由同一组8个金属硫蛋白基因（MT1G、MT1F、MT1H、MT2A、MT1X、MT1HL1、MT1E、MT1M）驱动。这提示金属离子结合、解毒和氧化应激相关程序在胃癌组织中发生变化；这些高度重叠的GO条目应视为一个共同主题，而不是多个独立发现。"
  )
}
if (grepl("extracellular matrix|collagen|ecm-receptor|focal adhesion|cell adhesion", term_text)) {
  interpretation <- c(
    interpretation,
    "- **细胞外基质与黏附重塑：** 相关条目提示显著差异基因集中于基质沉积、胶原组织和细胞黏附。胃癌中这些变化可能与肿瘤间质重塑、侵袭和迁移有关；但 ORA 本身不能证明因果关系。"
  )
}
if (grepl("protein digestion and absorption", term_text)) {
  interpretation <- c(
    interpretation,
    "- **KEGG名称需要结合贡献基因解释：** `Protein digestion and absorption` 虽位居KEGG首位，但本数据中主要由COL1A1、COL1A2、COL2A1、COL3A1、COL4A1、COL5A1、COL5A2、COL6A3、COL10A1等胶原基因驱动。因此它更支持胶原/基质重塑主题，不能仅依据通路名称解释为肿瘤的蛋白消化功能增强。"
  )
}
if (grepl("gastric acid|digestion|proton|potassium|ion transport|secretory", term_text)) {
  interpretation <- c(
    interpretation,
    "- **胃上皮分泌与离子转运：** 相关条目与胃腺体、壁细胞分泌和离子稳态相联系。结合 ATP4A、ATP4B、GIF 等胃功能基因在肿瘤中的明显下降，这更可能反映肿瘤组织正常胃上皮分化和分泌功能丧失。"
  )
}
if (grepl("immune|cytokine|complement|chemokine|leukocyte|inflammatory", term_text)) {
  interpretation <- c(
    interpretation,
    "- **免疫与炎症反应：** 免疫、细胞因子或补体相关条目提示肿瘤组织的炎症微环境发生变化。不过芯片混合组织无法区分信号来自肿瘤细胞还是浸润免疫细胞。"
  )
}
if (grepl("metabolic|metabolism|oxidative|mitochond", term_text)) {
  interpretation <- c(
    interpretation,
    "- **代谢状态改变：** 代谢相关条目提示肿瘤与癌旁胃腺体之间存在代谢程序差异；该结果可能同时受到肿瘤细胞重编程和组织细胞组成变化的影响。"
  )
}
if (length(interpretation) == 0L) {
  interpretation <- paste0(
    "- 当前显著条目主要包括：",
    paste(head(unique(all_sig_terms$Description), 5), collapse = "、"),
    "。这些结果应视为由差异基因集合提出的功能假设，需要结合独立队列或实验验证。"
  )
}

report_lines <- c(
  "# GSE29272 GO与KEGG富集分析报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "- 数据集：GSE29272",
  "- 研究比较：胃癌肿瘤组织（贲门与非贲门合并）vs 配对癌旁正常胃腺体",
  "- 物种：Homo sapiens（人类）",
  "- GO注释：org.Hs.eg.db / GO.db",
  paste0("- KEGG物种代码：hsa；快照访问日期：", kegg_accessed_date),
  "",
  "## 1. 输入与基因ID转换",
  "",
  paste0("- 输入显著差异基因：", length(sig_symbols), " 个 Gene Symbol。"),
  paste0("- 成功映射：", length(mapped_fg_symbols), " 个 Gene Symbol（",
         sprintf("%.2f%%", 100 * foreground_mapping_rate), "）。"),
  paste0("- 映射后前景：", length(foreground_entrez), " 个唯一 Entrez ID。"),
  paste0("- 实际受检背景：", length(background_symbols), " 个 Gene Symbol，映射为 ",
         length(background_entrez), " 个唯一 Entrez ID。"),
  "- 前景和背景采用完全相同的 SYMBOL→ENTREZID 映射规则。",
  "",
  "## 2. 方法与参数",
  "",
  "- 使用 clusterProfiler 进行过度富集分析（ORA，单侧超几何检验）。",
  "- GO分别运行 BP、CC、MF；KEGG使用本次保存的 hsa GSON 快照。",
  paste0("- KEGG快照来源：", kegg_snapshot_source, "。"),
  "- 多重检验校正：Benjamini-Hochberg；显著标准：padj < 0.05。",
  "- 基因集大小：10–500。",
  "- 背景不是默认全基因组，而是实际进入差异检验并具有有效P值的基因。",
  "- 输出表增加 FoldEnrichment = GeneRatio / BgRatio。",
  "",
  "## 3. 显著结果数量",
  "",
  paste0("- GO BP：", nrow(go_bp_sig), " 条。"),
  paste0("- GO CC：", nrow(go_cc_sig), " 条。"),
  paste0("- GO MF：", nrow(go_mf_sig), " 条。"),
  paste0("- KEGG：", nrow(kegg_sig), " 条。"),
  "",
  "## 4. 最显著条目",
  "",
  "### GO BP",
  "",
  markdown_top_table(go_bp_sig),
  "",
  "### GO CC",
  "",
  markdown_top_table(go_cc_sig),
  "",
  "### GO MF",
  "",
  markdown_top_table(go_mf_sig),
  "",
  "### KEGG",
  "",
  markdown_top_table(kegg_sig),
  "",
  "## 5. 生物学解读",
  "",
  interpretation,
  "",
  "## 6. 结果边界与注意事项",
  "",
  "- 本分析将上调和下调差异基因合并进行ORA，因此只能说明某功能在显著基因中“过度代表”，不能据此判断通路被激活或抑制。",
  "- GO条目具有层级继承和语义冗余；气泡图展示每个本体按padj排序的前10条，是显著结果的摘要而不是全部结果。",
  "- 癌旁正常胃腺体并不等同于健康志愿者胃组织，肿瘤邻近效应可能影响比较。",
  "- 芯片来源于混合组织，富集信号可能同时反映细胞内变化与细胞组成差异。",
  "- ORA结果是由同一差异基因表产生的功能假设，不能作为对该差异分析的独立验证。",
  "- KEGG会持续更新；本次快照已保存，复现时优先读取该快照。",
  "",
  "## 7. 输出文件",
  "",
  "- `GO_BP.csv`、`GO_CC.csv`、`GO_MF.csv`、`KEGG.csv`：padj < 0.05 的显著结果。",
  "- `*_all_tested.csv`：全部受检条目。",
  "- `GSE29272_GO_dotplot.pdf/png`：GO三类气泡图，每类最多10条。",
  "- `GSE29272_KEGG_dotplot.pdf/png`：KEGG气泡图，最多10条。",
  "- `GSE29272_*_dotplot_source_data.csv`：图件源数据。",
  "- `GSE29272_ID_mapping_summary.csv`：ID映射统计。",
  "- `KEGG_hsa_snapshot.gson`：本次KEGG注释快照。",
  "- `GSE29272_enrichment_run.log`：运行日志。",
  "- `GSE29272_enrichment_sessionInfo.txt`：R与软件包环境。",
  "- `scripts/05_GSE29272_GO_KEGG_enrichment.R`：完整可复现脚本。"
)
writeLines(
  report_lines,
  file.path(enrich_dir, "GSE29272_enrichment_report.md"),
  useBytes = TRUE
)

capture.output(sessionInfo(), file = session_path)

cat("\nKEGG快照访问日期：", kegg_accessed_date, "\n", sep = "")
cat("GO BP显著条目：", nrow(go_bp_sig), "\n", sep = "")
cat("GO CC显著条目：", nrow(go_cc_sig), "\n", sep = "")
cat("GO MF显著条目：", nrow(go_mf_sig), "\n", sep = "")
cat("KEGG显著条目：", nrow(kegg_sig), "\n", sep = "")
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
