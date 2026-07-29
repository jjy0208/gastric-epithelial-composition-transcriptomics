# GSE29272 配对差异表达分析
# 数据类型：Affymetrix GPL96芯片，GEO提交者已进行RMA标准化
# 比较：Disease（Cardia_Tumor + NonCardia_Tumor）vs Control（Normal）
# 方法：每个TYB/TYC配对内计算Tumor - Normal，再使用limma经验贝叶斯模型
# 主检验：limma TREAT，检验 |log2FC| > 1；BH方法控制FDR
# 注意：本数据不是RNA-seq整数计数，因此不使用DESeq2或edgeR。

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ragg)
})

options(stringsAsFactors = FALSE)
options(ggrepel.max.overlaps = Inf)

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
clean_dir <- file.path(project_dir, "clean_data")
deg_dir <- file.path(project_dir, "results", "DEG")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

accession <- "GSE29272"
expression_path <- file.path(
  clean_dir,
  "GSE29272_clean_expression_matrix.csv"
)
sample_info_path <- file.path(
  clean_dir,
  "GSE29272_sample_info.csv"
)

run_log_path <- file.path(deg_dir, "GSE29272_DEG_run.log")
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
cat("数据集：", accession, "\n", sep = "")
cat("比较：Disease (Tumor) vs Control (Normal)\n\n")

if (!file.exists(expression_path)) {
  stop("表达矩阵不存在：", expression_path)
}
if (!file.exists(sample_info_path)) {
  stop("样本信息表不存在：", sample_info_path)
}

expression <- as.matrix(
  read.csv(
    expression_path,
    row.names = 1,
    check.names = FALSE
  )
)
storage.mode(expression) <- "double"
sample_info <- fread(sample_info_path, encoding = "UTF-8", data.table = FALSE)

required_sample_columns <- c(
  "sample",
  "group",
  "disease_status",
  "anatomical_site",
  "pair_id",
  "platform"
)
missing_columns <- setdiff(required_sample_columns, colnames(sample_info))
if (length(missing_columns) > 0L) {
  stop("样本信息缺少字段：", paste(missing_columns, collapse = ", "))
}
if (!identical(colnames(expression), sample_info$sample)) {
  stop("样本信息sample名称或顺序与表达矩阵列名不完全一致。")
}
if (any(!is.finite(expression))) {
  stop("表达矩阵包含NA/Inf。")
}
if (anyDuplicated(rownames(expression)) > 0L) {
  stop("表达矩阵Gene Symbol重复。")
}

sample_info$analysis_group <- ifelse(
  sample_info$disease_status == "Tumor",
  "Disease",
  ifelse(sample_info$disease_status == "Normal", "Control", NA_character_)
)
if (anyNA(sample_info$analysis_group)) {
  stop("存在无法映射为Disease/Control的样本。")
}

pair_table <- table(sample_info$pair_id, sample_info$analysis_group)
if (
  !all(c("Control", "Disease") %in% colnames(pair_table)) ||
    any(pair_table[, "Control"] != 1L) ||
    any(pair_table[, "Disease"] != 1L)
) {
  stop("每个pair_id必须严格包含1个Disease和1个Control样本。")
}

pair_ids <- sort(unique(sample_info$pair_id))
disease_samples <- vapply(
  pair_ids,
  function(id) {
    sample_info$sample[
      sample_info$pair_id == id &
        sample_info$analysis_group == "Disease"
    ]
  },
  character(1)
)
control_samples <- vapply(
  pair_ids,
  function(id) {
    sample_info$sample[
      sample_info$pair_id == id &
        sample_info$analysis_group == "Control"
    ]
  },
  character(1)
)

# 配对差值为Tumor - Normal，因此正log2FC表示疾病组上调。
paired_difference <- (
  expression[, disease_samples, drop = FALSE] -
    expression[, control_samples, drop = FALSE]
)
colnames(paired_difference) <- pair_ids

design <- matrix(
  1,
  nrow = length(pair_ids),
  ncol = 1,
  dimnames = list(pair_ids, "Disease_vs_Control")
)

fit_base <- lmFit(paired_difference, design)
fit_standard <- eBayes(
  fit_base,
  trend = TRUE,
  robust = TRUE
)
fit_treat <- treat(
  fit_base,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

standard_table <- topTable(
  fit_standard,
  coef = "Disease_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
treat_table <- topTreat(
  fit_treat,
  coef = "Disease_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

if (
  !identical(rownames(standard_table), rownames(expression)) ||
    !identical(rownames(treat_table), rownames(expression))
) {
  stop("limma结果行顺序与输入基因顺序不一致。")
}

deg_all <- data.frame(
  gene = rownames(expression),
  log2FoldChange = unname(treat_table$logFC),
  meanExpression = rowMeans(expression),
  t = unname(treat_table$t),
  PValue = unname(treat_table$P.Value),
  padj = unname(treat_table$adj.P.Val),
  standard_t = unname(standard_table$t),
  standard_PValue = unname(standard_table$P.Value),
  standard_padj = unname(standard_table$adj.P.Val),
  B = unname(standard_table$B),
  stringsAsFactors = FALSE
)
deg_all$direction <- ifelse(
  deg_all$log2FoldChange > 1 & deg_all$padj < 0.05,
  "Up",
  ifelse(
    deg_all$log2FoldChange < -1 & deg_all$padj < 0.05,
    "Down",
    "Not_significant"
  )
)
deg_all <- deg_all[
  order(deg_all$padj, -abs(deg_all$log2FoldChange)),
  ,
  drop = FALSE
]

deg_significant <- deg_all[
  is.finite(deg_all$padj) &
    deg_all$padj < 0.05 &
    abs(deg_all$log2FoldChange) > 1,
  ,
  drop = FALSE
]

threshold_rule <- "|log2FoldChange| > 1 且 TREAT-BH padj < 0.05"
threshold_relaxed <- FALSE
if (nrow(deg_significant) == 0L) {
  # 不进行未经预先规定的数据驱动放宽；保留空表并在报告中说明。
  threshold_rule <- paste0(
    threshold_rule,
    "；未检出基因，未自动放宽阈值以避免选择性报告"
  )
}

all_path <- file.path(deg_dir, "GSE29272_DEG_all.csv")
significant_path <- file.path(
  deg_dir,
  "GSE29272_DEG_significant.csv"
)
fwrite(deg_all, all_path, bom = TRUE)
fwrite(deg_significant, significant_path, bom = TRUE)

n_up <- sum(deg_significant$direction == "Up")
n_down <- sum(deg_significant$direction == "Down")
n_significant <- nrow(deg_significant)
n_standard_bh <- sum(
  is.finite(deg_all$standard_padj) &
    deg_all$standard_padj < 0.05
)

summary_table <- data.frame(
  dataset = accession,
  method = "Paired differences + limma TREAT",
  comparison = "Disease_vs_Control",
  disease_samples = sum(sample_info$analysis_group == "Disease"),
  control_samples = sum(sample_info$analysis_group == "Control"),
  pairs = length(pair_ids),
  tested_genes = nrow(deg_all),
  standard_BH_padj_lt_0.05 = n_standard_bh,
  threshold = threshold_rule,
  significant_genes = n_significant,
  upregulated = n_up,
  downregulated = n_down,
  threshold_relaxed = threshold_relaxed,
  stringsAsFactors = FALSE
)
fwrite(
  summary_table,
  file.path(deg_dir, "GSE29272_DEG_summary.csv"),
  bom = TRUE
)

# ------------------------------
# Nature风格火山图
# ------------------------------

palette <- c(
  Up = "#B2182B",
  Down = "#2166AC",
  Not_significant = "#BDBDBD"
)

plot_data <- deg_all
positive_padj <- plot_data$padj[
  is.finite(plot_data$padj) & plot_data$padj > 0
]
minimum_plot_padj <- if (
  length(positive_padj) > 0L
) {
  min(positive_padj) / 10
} else {
  .Machine$double.xmin
}
plot_data$plot_padj <- pmax(
  ifelse(
    is.finite(plot_data$padj),
    plot_data$padj,
    1
  ),
  minimum_plot_padj
)
plot_data$minus_log10_padj <- -log10(plot_data$plot_padj)

label_up <- head(
  plot_data[plot_data$direction == "Up", ],
  4
)
label_down <- head(
  plot_data[plot_data$direction == "Down", ],
  4
)
label_data <- rbind(label_up, label_down)

x_limit <- max(
  2,
  ceiling(
    max(abs(plot_data$log2FoldChange), na.rm = TRUE)
  )
)

theme_nature <- function(base_size = 7.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.7),
      plot.title = element_text(
        size = base_size + 0.8,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4D4D4D"),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
}

volcano_plot <- ggplot(
  plot_data,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    colour = direction
  )
) +
  geom_point(
    alpha = 0.72,
    size = 0.8,
    stroke = 0
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.35,
    colour = "#636363"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.35,
    colour = "#636363"
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 2.25,
    family = "sans",
    box.padding = 0.25,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.colour = "#737373",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = palette,
    breaks = c("Up", "Down", "Not_significant"),
    labels = c("Up", "Down", "Not significant")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  coord_cartesian(xlim = c(-x_limit, x_limit), clip = "off") +
  labs(
    title = "GSE29272 paired differential expression",
    subtitle = paste0(
      "Tumor vs paired normal; limma TREAT, |log2FC| > 1; ",
      n_up,
      " up, ",
      n_down,
      " down"
    ),
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * " adjusted P value")
  ) +
  theme_nature() +
  theme(
    legend.position = "top",
    legend.justification = "left"
  )

save_ggplot_pdf_png <- function(
  plot,
  basename,
  width_mm,
  height_mm,
  dpi = 600
) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    paste0(basename, ".pdf"),
    width = width_in,
    height = height_in,
    family = "sans",
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(basename, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

volcano_base <- file.path(deg_dir, "GSE29272_volcano")
save_ggplot_pdf_png(
  volcano_plot,
  volcano_base,
  width_mm = 90,
  height_mm = 85,
  dpi = 600
)
fwrite(
  plot_data[
    ,
    c(
      "gene",
      "log2FoldChange",
      "PValue",
      "padj",
      "direction",
      "minus_log10_padj"
    )
  ],
  file.path(deg_dir, "GSE29272_volcano_source_data.csv"),
  bom = TRUE
)

# ------------------------------
# P值分布诊断
# ------------------------------

pvalue_data <- data.frame(
  PValue = deg_all$standard_PValue
)
pvalue_plot <- ggplot(pvalue_data, aes(x = PValue)) +
  geom_histogram(
    bins = 50,
    boundary = 0,
    fill = "#4C78A8",
    colour = "white",
    linewidth = 0.18
  ) +
  labs(
    title = "P-value distribution",
    subtitle = "Standard paired limma eBayes test",
    x = "Raw P value",
    y = "Number of genes"
  ) +
  theme_nature()

pvalue_base <- file.path(deg_dir, "GSE29272_pvalue_histogram")
save_ggplot_pdf_png(
  pvalue_plot,
  pvalue_base,
  width_mm = 90,
  height_mm = 70,
  dpi = 600
)

# ------------------------------
# Nature风格热图
# ------------------------------

heatmap_gene_n <- min(50L, nrow(deg_significant))
if (heatmap_gene_n == 0L) {
  # 仅在严格阈值无结果时用TREAT排序前50基因作探索性热图。
  heatmap_genes <- head(deg_all$gene, 50L)
  heatmap_selection_note <- paste0(
    "严格阈值无显著基因；热图展示TREAT padj排序前",
    length(heatmap_genes),
    "个基因，仅供探索"
  )
} else {
  heatmap_genes <- head(deg_significant$gene, heatmap_gene_n)
  heatmap_selection_note <- paste0(
    "热图展示严格阈值显著基因中TREAT padj最小的前",
    heatmap_gene_n,
    "个基因"
  )
}

group_order <- c("Normal", "NonCardia_Tumor", "Cardia_Tumor")
sample_order <- sample_info$sample[
  order(
    match(sample_info$group, group_order),
    sample_info$pair_id
  )
]
heatmap_meta <- sample_info[
  match(sample_order, sample_info$sample),
  ,
  drop = FALSE
]

heatmap_raw <- expression[
  heatmap_genes,
  sample_order,
  drop = FALSE
]
heatmap_z <- t(scale(t(heatmap_raw)))
if (any(!is.finite(heatmap_z))) {
  stop("热图行Z分数包含NA/Inf。")
}
heatmap_z[heatmap_z > 2.5] <- 2.5
heatmap_z[heatmap_z < -2.5] <- -2.5

direction_map <- setNames(
  deg_all$direction,
  deg_all$gene
)
heatmap_direction <- factor(
  direction_map[rownames(heatmap_z)],
  levels = c("Down", "Up", "Not_significant")
)
heatmap_group <- factor(
  heatmap_meta$group,
  levels = group_order
)

group_colours <- c(
  Normal = "#BDBDBD",
  NonCardia_Tumor = "#4C78A8",
  Cardia_Tumor = "#B2182B"
)
site_colours <- c(
  `Non-cardia` = "#4C78A8",
  Cardia = "#E69F00"
)

column_annotation <- HeatmapAnnotation(
  Group = heatmap_group,
  Site = heatmap_meta$anatomical_site,
  col = list(
    Group = group_colours,
    Site = site_colours
  ),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(
    fontfamily = "sans",
    fontsize = 6
  ),
  annotation_legend_param = list(
    title_gp = gpar(fontfamily = "sans", fontsize = 6.5, fontface = "bold"),
    labels_gp = gpar(fontfamily = "sans", fontsize = 6)
  )
)

heatmap_object <- Heatmap(
  heatmap_z,
  name = "Row Z-score",
  col = colorRamp2(
    c(-2.5, 0, 2.5),
    c("#2166AC", "#F7F7F7", "#B2182B")
  ),
  top_annotation = column_annotation,
  row_split = heatmap_direction,
  column_split = heatmap_group,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  show_row_names = TRUE,
  row_names_gp = gpar(
    fontfamily = "sans",
    fontsize = 5.4
  ),
  show_column_names = FALSE,
  row_title_gp = gpar(
    fontfamily = "sans",
    fontsize = 6.5,
    fontface = "bold"
  ),
  column_title_gp = gpar(
    fontfamily = "sans",
    fontsize = 6.5,
    fontface = "bold"
  ),
  heatmap_legend_param = list(
    title_gp = gpar(fontfamily = "sans", fontsize = 6.5, fontface = "bold"),
    labels_gp = gpar(fontfamily = "sans", fontsize = 6),
    at = c(-2.5, 0, 2.5),
    labels = c("-2.5", "0", "2.5"),
    legend_height = unit(30, "mm")
  ),
  border = FALSE,
  use_raster = TRUE,
  raster_quality = 3
)

draw_heatmap <- function() {
  draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    padding = unit(c(3, 3, 3, 3), "mm")
  )
}

heatmap_width_in <- 183 / 25.4
heatmap_height_in <- 130 / 25.4
heatmap_pdf <- file.path(deg_dir, "GSE29272_DEG_heatmap.pdf")
heatmap_png <- file.path(deg_dir, "GSE29272_DEG_heatmap.png")

grDevices::cairo_pdf(
  heatmap_pdf,
  width = heatmap_width_in,
  height = heatmap_height_in,
  family = "sans",
  bg = "white"
)
draw_heatmap()
grDevices::dev.off()

ragg::agg_png(
  heatmap_png,
  width = heatmap_width_in,
  height = heatmap_height_in,
  units = "in",
  res = 600,
  background = "white"
)
draw_heatmap()
grDevices::dev.off()

fwrite(
  data.table(
    gene = rownames(heatmap_raw),
    as.data.table(heatmap_raw, keep.rownames = FALSE)
  ),
  file.path(deg_dir, "GSE29272_heatmap_log2_expression.csv"),
  bom = TRUE
)
fwrite(
  data.table(
    gene = rownames(heatmap_z),
    as.data.table(heatmap_z, keep.rownames = FALSE)
  ),
  file.path(deg_dir, "GSE29272_heatmap_row_zscore.csv"),
  bom = TRUE
)

# ------------------------------
# 报告与可复现信息
# ------------------------------

top_gene_n <- min(10L, nrow(deg_significant))
if (top_gene_n > 0L) {
  top_genes <- deg_significant[seq_len(top_gene_n), ]
  top_gene_lines <- paste0(
    "- ",
    top_genes$gene,
    "：log2FC=",
    format(round(top_genes$log2FoldChange, 4), nsmall = 4),
    "，padj=",
    format(top_genes$padj, scientific = TRUE, digits = 3)
  )
} else {
  top_gene_lines <- "- 严格阈值下无显著基因。"
}

report_path <- file.path(deg_dir, "GSE29272_DEG_report.md")
report_lines <- c(
  "# GSE29272配对差异表达分析报告",
  "",
  paste0("- 分析时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- R版本：", R.version.string),
  paste0("- limma版本：", as.character(packageVersion("limma"))),
  "- 数据类型：Affymetrix GPL96表达芯片；GEO提交者已进行RMA标准化。",
  "- 统计方法：配对差值（Tumor - Normal）加limma经验贝叶斯模型。",
  "- 主效应检验：limma TREAT直接检验绝对log2FoldChange是否大于1，BH方法校正多重检验。",
  "- 普通eBayes非零效应检验同时保存在完整结果表的standard_*列中。",
  "- 未使用DESeq2或edgeR，因为输入不是RNA-seq整数计数。",
  "",
  "## 输入与分组",
  "",
  paste0("- 输入基因数：", nrow(expression)),
  paste0("- Disease样本数：", sum(sample_info$analysis_group == "Disease")),
  paste0("- Control样本数：", sum(sample_info$analysis_group == "Control")),
  paste0("- 有效配对数：", length(pair_ids)),
  "- Disease包含Cardia_Tumor和NonCardia_Tumor；Control为相应的配对Normal样本。",
  "- 每个TYB/TYC pair_id均严格包含1个Tumor和1个Normal，样本名与表达矩阵列名完全一致。",
  "",
  "## 模型与参数",
  "",
  "- 对每个基因、每个pair_id计算Tumor表达减Normal表达。",
  "- 使用`limma::lmFit`拟合配对差值的截距模型，系数即平均log2FoldChange。",
  "- 使用`trend = TRUE`建模平均表达相关的方差趋势。",
  "- 使用`robust = TRUE`降低高变异基因对经验贝叶斯方差估计的影响。",
  "- 主阈值：|log2FoldChange| > 1且TREAT-BH padj < 0.05。",
  "- 阈值未放宽；若严格阈值无结果，脚本保留空显著表而不会数据驱动地降低标准。",
  "",
  "## 差异结果",
  "",
  paste0("- 检验基因数：", nrow(deg_all)),
  paste0(
    "- 普通eBayes非零差异检验BH padj < 0.05：",
    n_standard_bh,
    "；该数量不包含两倍效应要求，不作为本任务主显著基因数。"
  ),
  paste0("- 显著差异基因：", n_significant),
  paste0("- 上调基因：", n_up),
  paste0("- 下调基因：", n_down),
  paste0("- 实际筛选规则：", threshold_rule),
  "",
  "### padj最小的前10个显著基因",
  "",
  top_gene_lines,
  "",
  "## 图形",
  "",
  "- 火山图横轴为Tumor相对Normal的log2FoldChange，纵轴为-log10(TREAT BH padj)。虚线表示|log2FC|=1和padj=0.05。",
  "- 普通eBayes原始P值直方图在0附近明显富集，与大样本配对肿瘤-正常比较中存在广泛非零差异一致；主结果进一步由TREAT限制为效应超过两倍。",
  paste0("- ", heatmap_selection_note, "。"),
  "- 热图数值为每个基因在全部样本内计算的行Z分数，并截断至[-2.5, 2.5]；它展示相对表达模式，不代表绝对表达量。",
  "- 热图列按Normal、NonCardia_Tumor、Cardia_Tumor分区，并在区内聚类；该分区属于预先定义的展示结构，不是无监督发现。",
  "- PDF保留矢量文字与线条；PNG以600 dpi导出。",
  "",
  "## 结果解释限制",
  "",
  "- 本分析是总体Tumor vs Normal比较，不能把结果解释为cardia特异或non-cardia特异效应。",
  "- 当前Gene Symbol来自GEO平台注释；未执行HGNC新旧命名更新。",
  "- TREAT padj控制的是效应量超过|log2FC|=1这一主张的FDR；standard_padj对应普通非零差异检验，两者不可互换。",
  "- 热图选取显著性最高的基因，不能作为独立验证证据。",
  "",
  "## 输出文件",
  "",
  "- `GSE29272_DEG_all.csv`：所有基因的完整结果。",
  "- `GSE29272_DEG_significant.csv`：严格阈值显著结果。",
  "- `GSE29272_DEG_summary.csv`：样本数、阈值及上下调计数。",
  "- `GSE29272_volcano.pdf/png`：Nature风格火山图。",
  "- `GSE29272_DEG_heatmap.pdf/png`：Nature风格热图。",
  "- `GSE29272_pvalue_histogram.pdf/png`：普通eBayes原始P值分布诊断。",
  "- `GSE29272_volcano_source_data.csv`：火山图源数据。",
  "- `GSE29272_heatmap_log2_expression.csv`：热图原始log2表达值。",
  "- `GSE29272_heatmap_row_zscore.csv`：热图实际绘制的行Z分数。",
  "- `GSE29272_DEG_run.log`：R运行日志。",
  "- `GSE29272_R_sessionInfo.txt`：R及软件包版本。",
  "- `scripts/04_GSE29272_paired_limma_DEG.R`：完整可复现脚本。"
)
writeLines(report_lines, report_path, useBytes = TRUE)

writeLines(
  capture.output(sessionInfo()),
  file.path(deg_dir, "GSE29272_R_sessionInfo.txt"),
  useBytes = TRUE
)

required_outputs <- c(
  all_path,
  significant_path,
  file.path(deg_dir, "GSE29272_DEG_summary.csv"),
  paste0(volcano_base, ".pdf"),
  paste0(volcano_base, ".png"),
  heatmap_pdf,
  heatmap_png,
  paste0(pvalue_base, ".pdf"),
  paste0(pvalue_base, ".png"),
  report_path
)
if (!all(file.exists(required_outputs))) {
  stop("部分必需输出文件未生成。")
}
if (any(file.info(required_outputs)$size <= 0L)) {
  stop("部分必需输出文件为空。")
}

cat("输入维度：", nrow(expression), " × ", ncol(expression), "\n", sep = "")
cat("配对数：", length(pair_ids), "\n", sep = "")
cat("显著DEG：", n_significant, "\n", sep = "")
cat("上调：", n_up, "；下调：", n_down, "\n", sep = "")
cat("结果目录：", deg_dir, "\n", sep = "")
cat("完成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
