# GSE270678 空间转录组正交验证
# 目的：在独立胃癌空间转录组队列中，定位锁定的 73 基因模块，
#       并检验其空间分布是否与独立定义的胃上皮/分化谱系标志一致。
# 统计原则：
#   1) 只使用 33 张主研究切片；39 张技术重复不作为独立样本。
#   2) spot 仅用于切片内计算；正式推断以供者为统计单位，避免 spot 伪重复。
#   3) 细胞谱系标志集合自动剔除 73 个候选基因，避免循环论证。

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(UCell)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(jsonlite)
  library(scales)
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
raw_dir <- file.path(project_dir, "raw_data", "GSE270678", "extracted")
result_dir <- file.path(project_dir, "results", "PaperValidation", "Spatial")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_dir, "15_GSE270678_spatial_validation_run.log")
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

candidate_file <- file.path(
  project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"
)
candidate <- fread(candidate_file)
stopifnot("gene" %in% names(candidate))
candidate_genes <- unique(trimws(candidate$gene))
candidate_genes <- candidate_genes[nzchar(candidate_genes)]
stopifnot(length(candidate_genes) == 73L)

# 这些标志基因在分析前定义，并自动去除与候选模块重叠的基因。
marker_sets <- list(
  Epithelial = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "TACSTD2", "CLDN4", "CLDN7"
  ),
  Gastric_differentiated = c(
    "GKN1", "GKN2", "TFF1", "TFF2", "MUC5AC", "MUC6", "PGC", "GIF",
    "ATP4A", "ATP4B", "PGA3", "PGA4", "PGA5", "VSIG1", "CAPN8"
  ),
  Cancer_epithelial = c(
    "CEACAM5", "CEACAM6", "KRT17", "KRT23", "TACSTD2", "CLDN4",
    "CLDN7", "EPCAM", "KRT19"
  ),
  Intestinal_metaplasia = c(
    "KRT20", "FABP1", "APOA1", "APOA4", "ALPI", "MUC13", "TFF3"
  ),
  Proliferative = c("MKI67", "TOP2A", "STMN1", "HMGB2", "TUBA1B"),
  Stromal = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1", "COL6A2", "SPARC", "VIM"
  ),
  Immune = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "CD3D", "CD3E",
    "MS4A1", "CD79A", "NKG7"
  )
)
marker_sets <- lapply(marker_sets, setdiff, y = candidate_genes)
if (any(lengths(marker_sets) < 3L)) {
  stop("候选基因剔除后，至少一个独立标志集合少于 3 个基因。")
}
signatures <- c(list(Candidate_73 = candidate_genes), marker_sets)

h5_files_all <- sort(list.files(
  raw_dir,
  pattern = "_filtered_feature_bc_matrix\\.h5$",
  full.names = TRUE
))
if (length(h5_files_all) != 72L) {
  stop("预期找到 72 个空间矩阵，实际为：", length(h5_files_all))
}

sample_table <- data.table(
  h5_path = h5_files_all,
  file = basename(h5_files_all)
)
sample_table[, gsm := sub("_.*$", "", file)]
sample_table[, sample_id := sub("_filtered_feature_bc_matrix\\.h5$", "", file)]
sample_table[, donor := sub("^.*_(D[0-9]+)T.*$", "\\1", sample_id)]
sample_table[, technical_replication := grepl("_[123]_filtered_feature", file)]
sample_table[, primary_section := !technical_replication]
sample_table[, section_label := sub("^GSM[0-9]+_", "", sample_id)]
fwrite(sample_table, file.path(result_dir, "GSE270678_section_manifest.csv"))

primary <- sample_table[primary_section == TRUE]
stopifnot(nrow(primary) == 33L, uniqueN(primary$donor) == 19L)

read_gzip_raw <- function(path) {
  con <- gzfile(path, open = "rb")
  on.exit(close(con))
  readBin(con, what = "raw", n = file.info(path)$size * 30 + 1e6)
}

read_scalefactor <- function(path) {
  raw <- read_gzip_raw(path)
  fromJSON(rawToChar(raw))
}

read_positions <- function(path) {
  pos <- fread(path, header = FALSE)
  if (ncol(pos) != 6L) stop("空间坐标文件应包含 6 列：", path)
  setnames(pos, c(
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row_fullres", "pxl_col_fullres"
  ))
  pos
}

serial_param <- SerialParam(progressbar = FALSE)
spot_cache <- file.path(result_dir, "GSE270678_primary_spot_UCell_scores.csv.gz")
coverage_cache <- file.path(result_dir, "GSE270678_signature_gene_coverage.csv")

if (file.exists(spot_cache) && file.exists(coverage_cache)) {
  cat("读取已完成的 33 张主切片 UCell 缓存。\n")
  spot_scores <- fread(spot_cache)
  coverage <- fread(coverage_cache)
} else {
  spot_results <- vector("list", nrow(primary))
  coverage_results <- vector("list", nrow(primary))

  for (i in seq_len(nrow(primary))) {
    rec <- primary[i]
    cat(sprintf("[%02d/%02d] %s\n", i, nrow(primary), rec$sample_id))

    counts <- Read10X_h5(rec$h5_path, use.names = TRUE, unique.features = TRUE)
    if (is.list(counts)) {
      if (!"Gene Expression" %in% names(counts)) {
        stop("H5 中未找到 Gene Expression：", rec$h5_path)
      }
      counts <- counts[["Gene Expression"]]
    }

    present_signatures <- lapply(signatures, intersect, y = rownames(counts))
    if (length(present_signatures$Candidate_73) < 65L) {
      stop(rec$sample_id, " 中候选基因覆盖不足 65 个。")
    }
    if (any(lengths(present_signatures[-1L]) < 3L)) {
      stop(rec$sample_id, " 中至少一个独立标志集合覆盖少于 3 个基因。")
    }

    scores <- as.data.table(ScoreSignatures_UCell(
      matrix = counts,
      features = present_signatures,
      maxRank = 1500,
      missing_genes = "skip",
      BPPARAM = serial_param,
      force.gc = TRUE
    ), keep.rownames = "barcode")

    prefix <- sub("_filtered_feature_bc_matrix\\.h5$", "", rec$h5_path)
    pos_path <- paste0(prefix, "_tissue_positions_list.csv.gz")
    pos <- read_positions(pos_path)
    dat <- merge(pos[in_tissue == 1L], scores, by = "barcode", all = FALSE)
    if (nrow(dat) < 100L) stop(rec$sample_id, " 的组织内有效 spot 少于 100。")

    dat[, `:=`(
      gsm = rec$gsm,
      section = rec$sample_id,
      section_label = rec$section_label,
      donor = rec$donor
    )]
    setcolorder(dat, c(
      "gsm", "section", "section_label", "donor", "barcode",
      "in_tissue", "array_row", "array_col", "pxl_row_fullres",
      "pxl_col_fullres"
    ))
    spot_results[[i]] <- dat

    coverage_results[[i]] <- data.table(
      section = rec$sample_id,
      signature = names(present_signatures),
      genes_requested = lengths(signatures),
      genes_present = lengths(present_signatures),
      genes_present_names = vapply(
        present_signatures, paste, collapse = ";", FUN.VALUE = character(1)
      )
    )
    rm(counts, scores, dat)
    invisible(gc())
  }

  spot_scores <- rbindlist(spot_results, use.names = TRUE, fill = TRUE)
  coverage <- rbindlist(coverage_results)
  fwrite(spot_scores, spot_cache)
  fwrite(coverage, coverage_cache)
}

score_cols <- grep("_UCell$", names(spot_scores), value = TRUE)
candidate_col <- "Candidate_73_UCell"
comparison_cols <- setdiff(score_cols, candidate_col)

section_summary_long <- melt(
  spot_scores,
  id.vars = c("gsm", "section", "section_label", "donor", "barcode"),
  measure.vars = score_cols,
  variable.name = "signature",
  value.name = "score"
)[, .(
  n_spots = .N,
  median_score = median(score, na.rm = TRUE),
  mean_score = mean(score, na.rm = TRUE)
), by = .(gsm, section, section_label, donor, signature)]
fwrite(section_summary_long, file.path(
  result_dir, "GSE270678_section_score_summary.csv"
))

section_cor <- rbindlist(lapply(comparison_cols, function(comp) {
  spot_scores[, {
    ok <- is.finite(get(candidate_col)) & is.finite(get(comp))
    rho <- if (sum(ok) >= 50L) {
      suppressWarnings(cor(
        get(candidate_col)[ok], get(comp)[ok], method = "spearman"
      ))
    } else {
      NA_real_
    }
    .(n_spots = sum(ok), spearman_rho = rho)
  }, by = .(gsm, section, section_label, donor)][
    , comparison := sub("_UCell$", "", comp)
  ]
}))
fwrite(section_cor, file.path(result_dir, "GSE270678_section_spatial_correlations.csv"))

# 同一供者可有多张空间切片，先在供者内平均相关系数，再进行组水平检验。
donor_cor <- section_cor[
  is.finite(spearman_rho),
  .(
    n_sections = .N,
    mean_spearman_rho = mean(spearman_rho),
    median_spearman_rho = median(spearman_rho)
  ),
  by = .(donor, comparison)
]
fwrite(donor_cor, file.path(result_dir, "GSE270678_donor_spatial_correlations.csv"))

cor_tests <- donor_cor[, {
  x <- mean_spearman_rho[is.finite(mean_spearman_rho)]
  wt <- suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE, conf.int = TRUE))
  tt <- t.test(x, mu = 0)
  .(
    n_donors = length(x),
    median_rho = median(x),
    q1_rho = unname(quantile(x, 0.25)),
    q3_rho = unname(quantile(x, 0.75)),
    mean_rho = mean(x),
    mean_rho_ci_low = unname(tt$conf.int[1]),
    mean_rho_ci_high = unname(tt$conf.int[2]),
    wilcoxon_p = wt$p.value,
    one_sample_t_p = tt$p.value
  )
}, by = comparison]
cor_tests[, wilcoxon_fdr := p.adjust(wilcoxon_p, method = "BH")]
cor_tests[, abs_median_rho := abs(median_rho)]
setorder(cor_tests, wilcoxon_fdr, -abs_median_rho)
cor_tests[, abs_median_rho := NULL]
fwrite(cor_tests, file.path(result_dir, "GSE270678_spatial_concordance_tests.csv"))

morandi <- c(
  "Epithelial" = "#8FA6A1",
  "Gastric_differentiated" = "#B6A38A",
  "Cancer_epithelial" = "#B98585",
  "Intestinal_metaplasia" = "#9C91A7",
  "Proliferative" = "#C0A36E",
  "Stromal" = "#8395A7",
  "Immune" = "#8FA17B"
)
theme_nature <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    plot.title = element_text(face = "bold", size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

p_cor <- ggplot(
  donor_cor,
  aes(comparison, mean_spearman_rho, colour = comparison)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, linewidth = 0.4,
    alpha = 0.18, aes(fill = comparison)
  ) +
  geom_jitter(width = 0.12, height = 0, size = 1.4, alpha = 0.8) +
  scale_colour_manual(values = morandi) +
  scale_fill_manual(values = morandi) +
  labs(
    x = NULL,
    y = "Within-section Spearman rho\n(donor-level mean)",
    title = "Spatial concordance of the 73-gene module"
  ) +
  theme_nature +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(result_dir, "GSE270678_spatial_correlations.pdf"),
  p_cor, width = 7.2, height = 4.6, units = "in", device = cairo_pdf
)
ggsave(
  file.path(result_dir, "GSE270678_spatial_correlations.png"),
  p_cor, width = 7.2, height = 4.6, units = "in", dpi = 600
)

plot_spatial_score <- function(section_id, score_col, title_text) {
  rec <- primary[sample_id == section_id]
  if (nrow(rec) != 1L) stop("无法唯一定位切片：", section_id)
  dat <- spot_scores[section == rec$sample_id]
  prefix <- sub("_filtered_feature_bc_matrix\\.h5$", "", rec$h5_path)

  img_path <- paste0(prefix, "_tissue_hires_image.png.gz")
  img <- readPNG(read_gzip_raw(img_path))
  sf <- read_scalefactor(paste0(prefix, "_scalefactors_json.json.gz"))
  hires_scale <- as.numeric(sf$tissue_hires_scalef)
  dat[, x_hires := pxl_col_fullres * hires_scale]
  dat[, y_hires := pxl_row_fullres * hires_scale]
  width <- dim(img)[2]
  height <- dim(img)[1]

  ggplot(dat, aes(x_hires, y_hires, colour = get(score_col))) +
    annotation_raster(img, xmin = 0, xmax = width, ymin = 0, ymax = height) +
    geom_point(size = 0.42, alpha = 0.78, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#D8D4CC", "#A9B5B1", "#7C8FA3", "#9B6D70"),
      limits = range(dat[[score_col]], finite = TRUE),
      oob = squish
    ) +
    coord_fixed(xlim = c(0, width), ylim = c(height, 0), expand = FALSE) +
    labs(title = title_text, colour = "UCell") +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(12, "mm")
    )
}

candidate_sections <- section_summary_long[
  signature == candidate_col
][order(median_score)]
target_quantiles <- quantile(candidate_sections$median_score, c(0.25, 0.75))
representative <- unique(vapply(target_quantiles, function(q) {
  candidate_sections$section[which.min(abs(candidate_sections$median_score - q))]
}, FUN.VALUE = character(1)))
if (length(representative) < 2L) {
  representative <- candidate_sections$section[c(
    max(1L, floor(.25 * nrow(candidate_sections))),
    min(nrow(candidate_sections), ceiling(.75 * nrow(candidate_sections)))
  )]
}

map_plots <- lapply(representative[1:2], function(sec) {
  label <- primary[sample_id == sec, section_label]
  plot_spatial_score(sec, candidate_col, paste0(label, ": 73-gene module"))
})

p_main <- (map_plots[[1]] | map_plots[[2]]) / p_cor +
  plot_layout(heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(result_dir, "GSE270678_spatial_validation_figure.pdf"),
  p_main, width = 8.2, height = 8.0, units = "in", device = cairo_pdf
)
ggsave(
  file.path(result_dir, "GSE270678_spatial_validation_figure.png"),
  p_main, width = 8.2, height = 8.0, units = "in", dpi = 600
)

fwrite(
  data.table(section = representative[1:2]),
  file.path(result_dir, "GSE270678_representative_sections.csv")
)

top_positive <- cor_tests[order(wilcoxon_fdr, -median_rho)][median_rho > 0]
top_negative <- cor_tests[order(wilcoxon_fdr, median_rho)][median_rho < 0]

fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "fg")
report <- c(
  "# GSE270678 空间转录组正交验证报告",
  "",
  paste0("- 分析日期：", format(Sys.Date(), "%Y-%m-%d")),
  "- 物种：Homo sapiens（人）。",
  "- 数据集：GSE270678；19 例胃癌患者的 33 张主研究空间切片。",
  "- GEO 中另有 39 张技术重复切片；已登记在清单中，但未作为独立推断单位。",
  paste0("- 组织内有效 spot 总数：", format(nrow(spot_scores), big.mark = ","), "。"),
  "- 打分方法：UCell（maxRank = 1500），直接基于各切片原始 UMI 矩阵。",
  "- 73 基因模块在分析前锁定；全部独立谱系标志均剔除了与 73 基因的重叠。",
  "",
  "## 统计设计",
  "",
  "先在每张切片内计算 73 基因模块与独立标志模块的 Spearman 相关系数；",
  "同一供者的多张切片先取平均，再以供者（n = 19）进行单样本 Wilcoxon 检验。",
  "因此 spot 仅用于切片内空间相关计算，不被当作独立生物学重复。",
  "",
  "## 主要结果",
  "",
  paste0(
    "- 最强正向空间相关：",
    if (nrow(top_positive)) {
      paste0(
        top_positive$comparison[1], "（供者中位 rho = ",
        fmt(top_positive$median_rho[1]), "，FDR = ",
        fmt(top_positive$wilcoxon_fdr[1]), "）。"
      )
    } else {
      "未观察到方向为正的比较。"
    }
  ),
  paste0(
    "- 最强负向空间相关：",
    if (nrow(top_negative)) {
      paste0(
        top_negative$comparison[1], "（供者中位 rho = ",
        fmt(top_negative$median_rho[1]), "，FDR = ",
        fmt(top_negative$wilcoxon_fdr[1]), "）。"
      )
    } else {
      "未观察到方向为负的比较。"
    }
  ),
  "",
  "## 解释边界",
  "",
  "该队列只包含肿瘤空间切片，不能替代新的肿瘤-邻癌配对检验，也不能证明因果关系。",
  "H&E 叠加图和独立标志的空间一致性用于判断该模块在组织中的定位，",
  "支持或反驳“组织组成/谱系代表性影响 bulk 信号”的解释；它不是蛋白水平 IHC 验证。",
  "",
  "## 代表切片",
  "",
  paste0("- ", paste(primary[sample_id %in% representative[1:2], section_label],
                     collapse = "、"),
         "；按切片 73 基因模块中位数的第 25 和第 75 百分位附近客观选取。")
)
writeLines(report, file.path(result_dir, "GSE270678_spatial_validation_report.md"),
           useBytes = TRUE)

writeLines(capture.output(sessionInfo()),
           file.path(result_dir, "15_sessionInfo.txt"), useBytes = TRUE)
cat("结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
