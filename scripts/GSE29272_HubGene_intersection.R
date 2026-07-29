#!/usr/bin/env Rscript

# ==============================================================================
# GSE29272 差异基因与 WGCNA 关键模块基因交集分析
#
# 输入：
#   1. results/DEG/GSE29272_DEG_significant.csv
#   2. results/WGCNA/GSE29272_WGCNA_hubmodule_genes.csv
#
# 输出：
#   results/HubGene/GSE29272_candidate_genes.csv
#   results/HubGene/GSE29272_DEG_WGCNA_Venn.pdf
#   results/HubGene/GSE29272_DEG_WGCNA_Venn.png
#   以及质量检查、摘要、日志、sessionInfo 和中文报告
#
# 图形契约：
#   核心结论：展示同时满足显著差异表达和关键模块归属的候选基因数量。
#   证据形式：双集合定量韦恩图，直接标注集合总数和三个互斥区域计数。
#   输出规格：白色背景、莫兰迪配色、矢量 PDF、600 dpi PNG。
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
deg_file <- file.path(
  project_dir, "results", "DEG", "GSE29272_DEG_significant.csv"
)
wgcna_file <- file.path(
  project_dir, "results", "WGCNA", "GSE29272_WGCNA_hubmodule_genes.csv"
)
out_dir <- file.path(project_dir, "results", "HubGene")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("data.table", "ggplot2", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("缺少必要 R 包：", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

log_file <- file.path(out_dir, "GSE29272_HubGene_run.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("GSE29272 DEG-WGCNA intersection analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

write_csv_utf8 <- function(x, filename) {
  data.table::fwrite(x, filename, bom = TRUE, na = "NA")
}

detect_id_type <- function(ids) {
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0) return("Unknown")
  prop_ensembl <- mean(grepl("^ENS[A-Z]*G[0-9]+(?:\\.[0-9]+)?$", ids))
  prop_entrez <- mean(grepl("^[0-9]+$", ids))
  prop_affy <- mean(grepl("(_at|_s_at|_x_at)$", ids, ignore.case = TRUE))
  if (prop_ensembl >= 0.80) return("Ensembl Gene ID")
  if (prop_entrez >= 0.80) return("Entrez Gene ID")
  if (prop_affy >= 0.80) return("Affymetrix probe ID")
  "Gene Symbol"
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
           formatC(x, format = "f", digits = 3))
  )
}

# ------------------------------------------------------------------------------
# 1. 读取输入并检查基因列
# ------------------------------------------------------------------------------
if (!file.exists(deg_file)) stop("显著差异基因表不存在：", deg_file)
if (!file.exists(wgcna_file)) stop("WGCNA 关键模块基因表不存在：", wgcna_file)

deg <- data.table::fread(deg_file, check.names = FALSE)
wgcna <- data.table::fread(wgcna_file, check.names = FALSE)

if (!"gene" %in% names(deg)) stop("差异基因表缺少 gene 列。")
if (!"gene" %in% names(wgcna)) stop("WGCNA 基因表缺少 gene 列。")

deg$gene_original <- as.character(deg$gene)
wgcna$gene_original <- as.character(wgcna$gene)
deg$gene <- trimws(deg$gene_original)
wgcna$gene <- trimws(wgcna$gene_original)

deg_empty <- sum(is.na(deg$gene) | !nzchar(deg$gene))
wgcna_empty <- sum(is.na(wgcna$gene) | !nzchar(wgcna$gene))
deg <- deg[!is.na(gene) & nzchar(gene)]
wgcna <- wgcna[!is.na(gene) & nzchar(gene)]

deg_duplicate_rows <- sum(duplicated(deg$gene))
wgcna_duplicate_rows <- sum(duplicated(wgcna$gene))

# 如果输入中意外存在重复基因，仅保留第一条记录，并在 QC 表中记录。
deg <- deg[!duplicated(gene)]
wgcna <- wgcna[!duplicated(gene)]

deg_id_type <- detect_id_type(deg$gene)
wgcna_id_type <- detect_id_type(wgcna$gene)

# 人类 Gene Symbol 采用大写键进行稳健匹配，同时保留原始显示符号。
deg$gene_key <- toupper(deg$gene)
wgcna$gene_key <- toupper(wgcna$gene)

deg_case_changed <- sum(deg$gene_key != deg$gene)
wgcna_case_changed <- sum(wgcna$gene_key != wgcna$gene)
exact_overlap_n <- length(intersect(deg$gene, wgcna$gene))
case_overlap_n <- length(intersect(deg$gene_key, wgcna$gene_key))

if (deg_id_type != wgcna_id_type) {
  stop(
    "两张表的基因ID类型不一致：DEG=", deg_id_type,
    "；WGCNA=", wgcna_id_type, "。请先进行显式ID转换。"
  )
}
if (deg_id_type != "Gene Symbol") {
  stop(
    "当前脚本要求两张表均为 Gene Symbol；检测到：", deg_id_type,
    "。为避免未经版本固定的在线映射，分析已停止。"
  )
}

cat("DEG unique genes:", nrow(deg), "\n")
cat("WGCNA unique genes:", nrow(wgcna), "\n")
cat("DEG ID type:", deg_id_type, "\n")
cat("WGCNA ID type:", wgcna_id_type, "\n")
cat("Exact overlap:", exact_overlap_n, "\n")
cat("Case-normalized overlap:", case_overlap_n, "\n")

# ------------------------------------------------------------------------------
# 2. 计算交集并整合两侧统计字段
# ------------------------------------------------------------------------------
deg_keys <- unique(deg$gene_key)
wgcna_keys <- unique(wgcna$gene_key)
candidate_keys <- intersect(deg_keys, wgcna_keys)

n_deg <- length(deg_keys)
n_wgcna <- length(wgcna_keys)
n_overlap <- length(candidate_keys)
n_deg_only <- length(setdiff(deg_keys, wgcna_keys))
n_wgcna_only <- length(setdiff(wgcna_keys, deg_keys))
overlap_pct_deg <- if (n_deg > 0) 100 * n_overlap / n_deg else NA_real_
overlap_pct_wgcna <- if (n_wgcna > 0) 100 * n_overlap / n_wgcna else NA_real_

# 给来源字段增加前缀，避免 log2FoldChange、padj 等同名字段覆盖。
deg_join <- as.data.frame(deg)
wgcna_join <- as.data.frame(wgcna)
deg_cols <- setdiff(names(deg_join), c("gene", "gene_key", "gene_original"))
wgcna_cols <- setdiff(names(wgcna_join), c("gene", "gene_key", "gene_original"))
names(deg_join)[match(deg_cols, names(deg_join))] <- paste0("DEG_", deg_cols)
names(wgcna_join)[match(wgcna_cols, names(wgcna_join))] <- paste0("WGCNA_", wgcna_cols)

candidate <- merge(
  deg_join[, c("gene", "gene_key", paste0("DEG_", deg_cols)), drop = FALSE],
  wgcna_join[, c("gene_key", paste0("WGCNA_", wgcna_cols)), drop = FALSE],
  by = "gene_key",
  all = FALSE,
  sort = FALSE
)
candidate <- candidate[candidate$gene_key %in% candidate_keys, , drop = FALSE]

# 输出列首位固定为 gene；gene_key 仅用于审计，不作为主显示字段。
candidate <- candidate[, c("gene", "gene_key", setdiff(names(candidate), c("gene", "gene_key")))]
if ("DEG_padj" %in% names(candidate)) {
  candidate <- candidate[order(candidate$DEG_padj, -abs(candidate$DEG_log2FoldChange)), ]
} else {
  candidate <- candidate[order(candidate$gene), ]
}
rownames(candidate) <- NULL

candidate_file <- file.path(out_dir, "GSE29272_candidate_genes.csv")
write_csv_utf8(candidate, candidate_file)

region_counts <- data.frame(
  region = c("DEG_only", "Overlap", "WGCNA_only"),
  count = c(n_deg_only, n_overlap, n_wgcna_only),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  region_counts,
  file.path(out_dir, "GSE29272_Venn_source_data.csv")
)

id_qc <- data.frame(
  dataset = c("DEG_significant", "WGCNA_key_module"),
  input_rows = c(
    nrow(deg) + deg_duplicate_rows + deg_empty,
    nrow(wgcna) + wgcna_duplicate_rows + wgcna_empty
  ),
  unique_nonempty_genes = c(nrow(deg), nrow(wgcna)),
  empty_removed = c(deg_empty, wgcna_empty),
  duplicate_rows_removed = c(deg_duplicate_rows, wgcna_duplicate_rows),
  detected_id_type = c(deg_id_type, wgcna_id_type),
  case_normalized_entries = c(deg_case_changed, wgcna_case_changed),
  stringsAsFactors = FALSE
)
write_csv_utf8(id_qc, file.path(out_dir, "GSE29272_gene_ID_QC.csv"))

# ------------------------------------------------------------------------------
# 3. 绘制双集合韦恩图
# ------------------------------------------------------------------------------
# 使用参数方程绘制两个等面积圆；区域数值来自真实集合运算。
circle_points <- function(cx, cy, radius, set_name, n = 720) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = cx + radius * cos(theta),
    y = cy + radius * sin(theta),
    set = set_name,
    stringsAsFactors = FALSE
  )
}

circle_radius <- 1.18
left_center <- -0.63
right_center <- 0.63
venn_df <- rbind(
  circle_points(left_center, 0, circle_radius, "Significant DEGs"),
  circle_points(right_center, 0, circle_radius, "WGCNA key module")
)

# 莫兰迪蓝灰与陶土色；透明叠加自然形成重叠区。
morandi_colors <- c(
  "Significant DEGs" = "#7890A8",
  "WGCNA key module" = "#B4877E"
)

venn_plot <- ggplot() +
  geom_polygon(
    data = venn_df,
    aes(x = x, y = y, group = set, fill = set),
    color = "#4A4A4A", linewidth = 0.45, alpha = 0.58
  ) +
  annotate(
    "text", x = left_center - 0.52, y = 0,
    label = n_deg_only, size = 6.0, fontface = "bold", color = "#222222"
  ) +
  annotate(
    "text", x = 0, y = 0,
    label = n_overlap, size = 6.4, fontface = "bold", color = "#222222"
  ) +
  annotate(
    "text", x = right_center + 0.52, y = 0,
    label = n_wgcna_only, size = 6.0, fontface = "bold", color = "#222222"
  ) +
  annotate(
    "text", x = left_center - 0.28, y = 1.38,
    label = paste0("Significant DEGs\nn = ", n_deg),
    size = 3.4, lineheight = 0.95, fontface = "bold", color = "#394A59"
  ) +
  annotate(
    "text", x = right_center + 0.28, y = 1.38,
    label = paste0("WGCNA key module\nn = ", n_wgcna),
    size = 3.4, lineheight = 0.95, fontface = "bold", color = "#6D4742"
  ) +
  annotate(
    "text", x = 0, y = -1.45,
    label = paste0("Candidate genes in overlap: n = ", n_overlap),
    size = 3.0, color = "#4A4A4A"
  ) +
  scale_fill_manual(values = morandi_colors) +
  coord_fixed(
    xlim = c(-2.05, 2.05), ylim = c(-1.70, 1.72),
    clip = "off", expand = FALSE
  ) +
  labs(title = "GSE29272 DEG–WGCNA intersection") +
  theme_void(base_family = "sans", base_size = 7) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      size = 9, face = "bold", hjust = 0.5, color = "#222222",
      margin = margin(b = 4)
    ),
    plot.margin = margin(7, 7, 7, 7),
    plot.background = element_rect(fill = "white", color = NA)
  )

figure_stem <- file.path(out_dir, "GSE29272_DEG_WGCNA_Venn")
ggsave(
  paste0(figure_stem, ".pdf"), venn_plot,
  width = 100 / 25.4, height = 88 / 25.4,
  device = grDevices::cairo_pdf, units = "in", bg = "white"
)
ggsave(
  paste0(figure_stem, ".png"), venn_plot,
  width = 100 / 25.4, height = 88 / 25.4,
  device = ragg::agg_png, units = "in", dpi = 600, bg = "white"
)

# ------------------------------------------------------------------------------
# 4. 生成摘要和中文报告
# ------------------------------------------------------------------------------
hub_candidate_overlap_n <- if ("WGCNA_hub_candidate" %in% names(candidate)) {
  sum(candidate$WGCNA_hub_candidate %in% TRUE, na.rm = TRUE)
} else {
  NA_integer_
}
up_n <- if ("DEG_direction" %in% names(candidate)) {
  sum(tolower(candidate$DEG_direction) == "up", na.rm = TRUE)
} else if ("DEG_log2FoldChange" %in% names(candidate)) {
  sum(candidate$DEG_log2FoldChange > 0, na.rm = TRUE)
} else {
  NA_integer_
}
down_n <- if ("DEG_direction" %in% names(candidate)) {
  sum(tolower(candidate$DEG_direction) == "down", na.rm = TRUE)
} else if ("DEG_log2FoldChange" %in% names(candidate)) {
  sum(candidate$DEG_log2FoldChange < 0, na.rm = TRUE)
} else {
  NA_integer_
}

summary_table <- data.frame(
  metric = c(
    "DEG_unique_genes", "WGCNA_key_module_unique_genes",
    "DEG_only", "WGCNA_only", "candidate_overlap_genes",
    "overlap_percent_of_DEG", "overlap_percent_of_WGCNA",
    "candidate_upregulated", "candidate_downregulated",
    "candidate_meeting_WGCNA_hub_threshold",
    "exact_overlap", "case_normalized_overlap"
  ),
  value = as.character(c(
    n_deg, n_wgcna, n_deg_only, n_wgcna_only, n_overlap,
    overlap_pct_deg, overlap_pct_wgcna, up_n, down_n,
    hub_candidate_overlap_n, exact_overlap_n, case_overlap_n
  )),
  stringsAsFactors = FALSE
)
write_csv_utf8(summary_table, file.path(out_dir, "GSE29272_HubGene_summary.csv"))

top_candidate_text <- if (nrow(candidate) > 0) {
  paste(head(candidate$gene, 15), collapse = "、")
} else {
  "无"
}

report_lines <- c(
  "# GSE29272 差异基因–WGCNA关键模块交集分析报告",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- R版本：", R.version.string),
  "",
  "## 1. 输入文件与基因格式检查",
  "",
  paste0(
    "读取现有显著差异基因表和WGCNA关键模块基因表。本步骤直接使用已有筛选结果，",
    "没有重新设定差异分析阈值，也没有只保留WGCNA表中 `hub_candidate=TRUE` 的子集。"
  ),
  paste0(
    "显著差异基因表包含 ", n_deg, " 个唯一基因；WGCNA关键模块表包含 ",
    n_wgcna, " 个唯一基因。两张表的基因列均识别为 ", deg_id_type, "。"
  ),
  paste0(
    "去除首尾空格并转为大写匹配键后，精确匹配交集和大小写标准化交集均为 ",
    n_overlap, "，说明没有因大小写差异额外获得或丢失基因。"
  ),
  paste0(
    "空基因名移除数：DEG=", deg_empty, "，WGCNA=", wgcna_empty,
    "；重复行移除数：DEG=", deg_duplicate_rows, "，WGCNA=", wgcna_duplicate_rows, "。"
  ),
  "",
  "## 2. 交集结果",
  "",
  paste0(
    "两个集合取交集后得到 **", n_overlap, " 个候选基因**。其中仅属于显著差异基因集合的有 ",
    n_deg_only, " 个，仅属于WGCNA关键模块的有 ", n_wgcna_only, " 个。"
  ),
  paste0(
    "交集基因占全部显著差异基因的 ", sprintf("%.1f%%", overlap_pct_deg),
    "，占WGCNA关键模块基因的 ", sprintf("%.1f%%", overlap_pct_wgcna), "。"
  ),
  paste0(
    "候选基因中上调 ", up_n, " 个、下调 ", down_n, " 个；其中 ",
    hub_candidate_overlap_n, " 个同时满足上一阶段预设的WGCNA候选枢纽阈值。"
  ),
  paste0("按差异分析校正后P值排序靠前的候选基因包括：", top_candidate_text, "。"),
  "",
  "## 3. 生物学意义",
  "",
  paste0(
    "交集基因同时具有两类证据：一方面，它们在疾病组与对照组之间达到既定差异表达标准；",
    "另一方面，它们属于与疾病状态关联最强的WGCNA共表达模块。因此，交集能够从",
    "“单基因表达改变”和“共表达网络归属”两个角度缩小后续验证范围。"
  ),
  paste0(
    "但两类结果均来自同一个GSE29272队列，二者并不是相互独立的外部验证证据；",
    "取交集会提高候选优先级，同时也可能漏掉效应较小但具有调控作用的基因。",
    "这些候选基因应继续结合独立队列、功能富集、蛋白互作网络及实验验证进行筛选。"
  ),
  "",
  "## 4. 输出文件",
  "",
  paste0(
    "- `GSE29272_candidate_genes.csv`：", n_overlap,
    "个交集候选基因及DEG、WGCNA两侧统计字段。"
  ),
  "- `GSE29272_DEG_WGCNA_Venn.pdf/png`：莫兰迪配色韦恩图，分别为矢量PDF和600 dpi PNG。",
  "- `GSE29272_Venn_source_data.csv`：韦恩图三个互斥区域的源数据。",
  "- `GSE29272_gene_ID_QC.csv`：基因ID类型、空值、重复和大小写检查。",
  "- `GSE29272_HubGene_summary.csv`：主要数量统计。",
  "- `GSE29272_HubGene_run.log` 和 `GSE29272_HubGene_sessionInfo.txt`：运行记录与环境信息。",
  "",
  "## 5. 可复现性",
  "",
  "完整R代码保存在 `scripts/GSE29272_HubGene_intersection.R`，代码包含中文注释。"
)
writeLines(
  report_lines,
  file.path(out_dir, "GSE29272_HubGene_report.md"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "GSE29272_HubGene_sessionInfo.txt")
)

cat("Candidate overlap genes:", n_overlap, "\n")
cat("DEG only:", n_deg_only, "\n")
cat("WGCNA only:", n_wgcna_only, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
