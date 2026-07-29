#!/usr/bin/env Rscript

# ==============================================================================
# GSE29272 交集候选基因 STRING 蛋白互作网络分析
#
# 数据库：STRING 12.0
# 物种：Homo sapiens（NCBI Taxonomy ID: 9606）
# 置信度：0.4（STRING combined score >= 400）
# 网络扩展：0，即不添加输入列表之外的节点
#
# 输入：
#   results/HubGene/GSE29272_candidate_genes.csv
#
# 输出：
#   results/PPI 下的映射表、未映射基因、孤立节点、边表、节点表、
#   证据通道汇总、网络图、中文报告、运行日志和R环境信息。
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
set.seed(20260729)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) != 1L) stop("请使用Rscript运行本脚本。")
  normalizePath(sub("^--file=", "", f), winslash = "/", mustWork = TRUE)
}

project_dir <- normalizePath(file.path(dirname(get_script_path()), ".."),
                             winslash = "/", mustWork = TRUE)
candidate_file <- file.path(
  project_dir, "results", "HubGene", "GSE29272_candidate_genes.csv"
)
out_dir <- file.path(project_dir, "results", "PPI")
cache_dir <- file.path(out_dir, "string_cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "STRINGdb", "data.table", "igraph", "ggraph", "tidygraph",
  "ggplot2", "ggrepel", "ragg", "httr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("缺少必要R包：", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(STRINGdb)
  library(data.table)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(ggrepel)
  library(httr)
})

log_file <- file.path(out_dir, "GSE29272_STRING_PPI_run.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("GSE29272 STRING PPI analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

write_csv_utf8 <- function(x, filename) {
  data.table::fwrite(x, filename, bom = TRUE, na = "NA")
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
           formatC(x, format = "f", digits = 3))
  )
}

# ------------------------------------------------------------------------------
# 1. 读取并清理交集基因
# ------------------------------------------------------------------------------
if (!file.exists(candidate_file)) {
  stop("交集候选基因文件不存在：", candidate_file)
}
candidate <- data.table::fread(candidate_file, check.names = FALSE)
if (!"gene" %in% names(candidate)) stop("候选基因表缺少 gene 列。")

input_genes_raw <- as.character(candidate$gene)
clean_genes <- toupper(trimws(input_genes_raw))
valid <- !is.na(clean_genes) & nzchar(clean_genes)
cleaning_table <- data.frame(
  original_gene = input_genes_raw,
  cleaned_gene = clean_genes,
  removed_as_empty = !valid,
  stringsAsFactors = FALSE
)

genes <- unique(clean_genes[valid])
duplicate_count <- sum(duplicated(clean_genes[valid]))
if (length(genes) == 0) stop("清理后没有可用于STRING查询的基因。")

write_csv_utf8(
  cleaning_table,
  file.path(out_dir, "GSE29272_STRING_gene_name_cleaning.csv")
)
cat("Input rows:", length(input_genes_raw), "\n")
cat("Unique cleaned genes:", length(genes), "\n")
cat("Duplicates removed:", duplicate_count, "\n")

# ------------------------------------------------------------------------------
# 2. 初始化 STRINGdb 并完成人类基因—蛋白映射
# ------------------------------------------------------------------------------
string_db_version_requested <- "12.0"
taxonomy_id <- 9606
confidence <- 0.4
score_threshold <- as.integer(confidence * 1000)
query_date <- as.character(Sys.Date())
random_seed <- 20260729

string_db <- STRINGdb$new(
  version = string_db_version_requested,
  species = taxonomy_id,
  score_threshold = score_threshold,
  input_directory = cache_dir
)

mapping_input <- data.frame(gene = genes, stringsAsFactors = FALSE)
mapped_raw <- string_db$map(
  mapping_input,
  "gene",
  removeUnmappedRows = FALSE,
  takeFirst = TRUE,
  quiet = TRUE
)

if (!"STRING_id" %in% names(mapped_raw)) {
  stop("STRINGdb映射结果缺少STRING_id列，无法继续。")
}

protein_info <- string_db$get_proteins()
protein_id_col <- intersect(
  c("protein_external_id", "STRING_id", "string_external_id"),
  names(protein_info)
)[1]
preferred_col <- intersect(
  c("preferred_name", "preferredName", "protein_name"),
  names(protein_info)
)[1]
annotation_col <- intersect(c("annotation", "description"), names(protein_info))[1]

if (is.na(protein_id_col) || is.na(preferred_col)) {
  stop(
    "无法识别STRING蛋白信息字段。实际列：",
    paste(names(protein_info), collapse = ", ")
  )
}

protein_keep <- unique(c(protein_id_col, preferred_col, annotation_col))
protein_keep <- protein_keep[!is.na(protein_keep)]
protein_small <- as.data.frame(protein_info[, protein_keep, drop = FALSE])
names(protein_small)[names(protein_small) == protein_id_col] <- "STRING_id"
names(protein_small)[names(protein_small) == preferred_col] <- "STRING_preferred_name"
if (!is.na(annotation_col)) {
  names(protein_small)[names(protein_small) == annotation_col] <- "STRING_annotation"
}

mapping_table <- merge(
  as.data.frame(mapped_raw),
  protein_small,
  by = "STRING_id",
  all.x = TRUE,
  sort = FALSE
)
mapping_table <- mapping_table[
  match(mapping_input$gene, mapping_table$gene),
  ,
  drop = FALSE
]
mapping_table$mapped <- !is.na(mapping_table$STRING_id) &
  nzchar(mapping_table$STRING_id)
mapping_table$taxonomy_id <- taxonomy_id
mapping_table$mapping_method <- "STRINGdb::map(takeFirst=TRUE)"

write_csv_utf8(
  mapping_table,
  file.path(out_dir, "GSE29272_STRING_gene_protein_mapping.csv")
)

unmapped <- mapping_table[
  !mapping_table$mapped,
  intersect(c("gene", "STRING_id", "mapped", "mapping_method"), names(mapping_table)),
  drop = FALSE
]
if (nrow(unmapped) == 0) {
  unmapped <- data.frame(
    gene = character(),
    STRING_id = character(),
    mapped = logical(),
    mapping_method = character(),
    stringsAsFactors = FALSE
  )
}
write_csv_utf8(
  unmapped,
  file.path(out_dir, "GSE29272_STRING_unmapped_genes.csv")
)

mapped_table <- mapping_table[mapping_table$mapped, , drop = FALSE]
mapped_string_ids <- unique(mapped_table$STRING_id)
mapped_gene_n <- length(unique(mapped_table$gene))
mapped_protein_n <- length(mapped_string_ids)

cat("Mapped genes:", mapped_gene_n, "\n")
cat("Mapped STRING proteins:", mapped_protein_n, "\n")
cat("Unmapped genes:", nrow(unmapped), "\n")

# ------------------------------------------------------------------------------
# 3. 获取输入蛋白内部的 STRING 网络；不添加任何额外节点
# ------------------------------------------------------------------------------
interactions_raw <- string_db$get_interactions(mapped_string_ids)
if (!all(c("from", "to", "combined_score") %in% names(interactions_raw))) {
  stop(
    "STRINGdb交互表缺少必要字段。实际列：",
    paste(names(interactions_raw), collapse = ", ")
  )
}

interactions <- as.data.frame(interactions_raw)
interactions <- interactions[
  interactions$from %in% mapped_string_ids &
    interactions$to %in% mapped_string_ids &
    interactions$combined_score >= score_threshold,
  ,
  drop = FALSE
]

# 无向网络去重：同一对蛋白如有多条记录，保留最高combined score。
interactions$edge_key <- vapply(
  seq_len(nrow(interactions)),
  function(i) paste(sort(c(interactions$from[i], interactions$to[i])), collapse = "|"),
  character(1)
)
if (nrow(interactions) > 0) {
  interactions <- interactions[
    order(interactions$edge_key, -interactions$combined_score),
    ,
    drop = FALSE
  ]
  interactions <- interactions[!duplicated(interactions$edge_key), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# 4. 查询同一STRING版本的边级证据通道
#    网络边仍由STRINGdb构建；REST返回仅补充每条边的证据通道分数。
# ------------------------------------------------------------------------------
evidence_channels <- c(
  nscore = "neighborhood",
  fscore = "gene_fusion",
  pscore = "cooccurrence",
  ascore = "coexpression",
  escore = "experiments",
  dscore = "database",
  tscore = "textmining"
)

api_url <- paste0(
  "https://version-", gsub("\\.", "-", string_db_version_requested),
  ".string-db.org/api/tsv/network"
)
api_response <- httr::POST(
  api_url,
  body = list(
    identifiers = paste(mapped_table$gene, collapse = "\r"),
    species = taxonomy_id,
    required_score = score_threshold,
    network_type = "functional",
    caller_identity = "GSE29272_R_STRINGdb"
  ),
  encode = "form",
  httr::timeout(600)
)
httr::stop_for_status(api_response)
api_text <- httr::content(api_response, as = "text", encoding = "UTF-8")

if (nzchar(trimws(api_text))) {
  evidence_raw <- utils::read.delim(
    text = api_text, sep = "\t", header = TRUE,
    check.names = FALSE, stringsAsFactors = FALSE
  )
} else {
  evidence_raw <- data.frame()
}

required_evidence_cols <- c(
  "stringId_A", "stringId_B", "score", names(evidence_channels)
)
if (nrow(evidence_raw) > 0 &&
    !all(required_evidence_cols %in% names(evidence_raw))) {
  stop(
    "STRING证据通道返回字段不完整。缺少：",
    paste(setdiff(required_evidence_cols, names(evidence_raw)), collapse = ", ")
  )
}

if (nrow(evidence_raw) > 0) {
  evidence_raw$edge_key <- vapply(
    seq_len(nrow(evidence_raw)),
    function(i) {
      paste(
        sort(c(evidence_raw$stringId_A[i], evidence_raw$stringId_B[i])),
        collapse = "|"
      )
    },
    character(1)
  )
  evidence_raw <- evidence_raw[
    order(evidence_raw$edge_key, -evidence_raw$score),
    ,
    drop = FALSE
  ]
  evidence_raw <- evidence_raw[!duplicated(evidence_raw$edge_key), , drop = FALSE]
  evidence_keep <- evidence_raw[
    ,
    c(
      "edge_key", "score", names(evidence_channels),
      "preferredName_A", "preferredName_B"
    ),
    drop = FALSE
  ]
  names(evidence_keep)[names(evidence_keep) == "score"] <- "api_combined_score"
  interactions <- merge(
    interactions,
    evidence_keep,
    by = "edge_key",
    all.x = TRUE,
    sort = FALSE
  )
} else {
  for (nm in c("api_combined_score", names(evidence_channels))) {
    interactions[[nm]] <- numeric(nrow(interactions))
  }
}

# 为边表补充基因显示名。
id_to_gene <- tapply(
  mapped_table$gene,
  mapped_table$STRING_id,
  function(x) paste(sort(unique(x)), collapse = "/")
)
id_to_preferred <- setNames(
  mapped_table$STRING_preferred_name,
  mapped_table$STRING_id
)
interactions$gene_A <- unname(id_to_gene[interactions$from])
interactions$gene_B <- unname(id_to_gene[interactions$to])
interactions$preferred_name_A <- unname(id_to_preferred[interactions$from])
interactions$preferred_name_B <- unname(id_to_preferred[interactions$to])
interactions$confidence <- interactions$combined_score / 1000

edge_output_columns <- c(
  "gene_A", "gene_B", "preferred_name_A", "preferred_name_B",
  "from", "to", "combined_score", "confidence",
  "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore",
  "edge_key"
)
edge_output_columns <- intersect(edge_output_columns, names(interactions))
write_csv_utf8(
  interactions[, edge_output_columns, drop = FALSE],
  file.path(out_dir, "GSE29272_STRING_PPI_edges.csv")
)

# 证据通道统计：分数>0表示该通道对该边提供了证据。
channel_summary <- do.call(
  rbind,
  lapply(names(evidence_channels), function(code) {
    values <- if (code %in% names(interactions)) interactions[[code]] else numeric()
    data.frame(
      channel_code = code,
      evidence_channel = unname(evidence_channels[code]),
      edges_with_channel_score_gt_0 = sum(values > 0, na.rm = TRUE),
      percent_of_edges = if (nrow(interactions) > 0) {
        100 * sum(values > 0, na.rm = TRUE) / nrow(interactions)
      } else {
        NA_real_
      },
      mean_channel_score = if (length(values) > 0) mean(values, na.rm = TRUE) else NA_real_,
      max_channel_score = if (length(values) > 0) max(values, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
)
write_csv_utf8(
  channel_summary,
  file.path(out_dir, "GSE29272_STRING_evidence_channels.csv")
)

# ------------------------------------------------------------------------------
# 5. 构建无向网络，计算Degree并识别孤立节点
# ------------------------------------------------------------------------------
node_table <- data.frame(
  STRING_id = mapped_string_ids,
  gene = unname(id_to_gene[mapped_string_ids]),
  STRING_preferred_name = unname(id_to_preferred[mapped_string_ids]),
  stringsAsFactors = FALSE
)
node_table$label <- ifelse(
  !is.na(node_table$gene) & nzchar(node_table$gene),
  node_table$gene,
  node_table$STRING_preferred_name
)

graph_edges <- interactions[, c("from", "to", "combined_score", "confidence"), drop = FALSE]
ppi_graph <- igraph::graph_from_data_frame(
  graph_edges,
  directed = FALSE,
  vertices = data.frame(
    name = node_table$STRING_id,
    label = node_table$label,
    gene = node_table$gene,
    stringsAsFactors = FALSE
  )
)
ppi_graph <- igraph::simplify(
  ppi_graph,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = list(
    combined_score = "max",
    confidence = "max",
    "ignore"
  )
)

node_table$Degree <- igraph::degree(ppi_graph, v = node_table$STRING_id)
node_table$is_isolated <- node_table$Degree == 0
node_table$node_class <- ifelse(
  node_table$is_isolated, "Isolated node", "Connected node"
)
node_table <- node_table[order(-node_table$Degree, node_table$gene), ]

write_csv_utf8(
  node_table,
  file.path(out_dir, "GSE29272_STRING_PPI_nodes.csv")
)
isolated_nodes <- node_table[
  node_table$is_isolated,
  c("gene", "STRING_id", "STRING_preferred_name", "Degree", "node_class"),
  drop = FALSE
]
if (nrow(isolated_nodes) == 0) {
  isolated_nodes <- data.frame(
    gene = character(),
    STRING_id = character(),
    STRING_preferred_name = character(),
    Degree = integer(),
    node_class = character(),
    stringsAsFactors = FALSE
  )
}
write_csv_utf8(
  isolated_nodes,
  file.path(out_dir, "GSE29272_STRING_isolated_nodes.csv")
)

network_nodes <- igraph::vcount(ppi_graph)
network_edges <- igraph::ecount(ppi_graph)
isolated_n <- sum(node_table$is_isolated)
connected_n <- network_nodes - isolated_n
component_n <- igraph::components(ppi_graph)$no
network_density <- igraph::edge_density(ppi_graph, loops = FALSE)
mean_degree <- mean(node_table$Degree)
max_degree <- max(node_table$Degree)

if (network_nodes != mapped_protein_n) {
  stop("网络节点数与映射蛋白数不一致，可能意外添加或遗漏节点。")
}
if (!all(igraph::V(ppi_graph)$name %in% mapped_string_ids)) {
  stop("网络含有输入映射集合之外的额外节点。")
}

cat("Network nodes:", network_nodes, "\n")
cat("Network edges:", network_edges, "\n")
cat("Isolated nodes:", isolated_n, "\n")
cat("Connected components:", component_n, "\n")

# ------------------------------------------------------------------------------
# 6. 绘制适合发表的PPI网络图
# ------------------------------------------------------------------------------
# 核心网络采用Fruchterman-Reingold力导向布局；
# 孤立节点没有边，因此将其均匀放在外围，避免相互遮挡。
set.seed(random_seed)
layout_obj <- ggraph::create_layout(
  ppi_graph,
  layout = "fr",
  weights = igraph::E(ppi_graph)$confidence,
  niter = 3000
)

layout_df <- as.data.frame(layout_obj)
layout_df$Degree <- igraph::degree(ppi_graph, v = layout_df$name)
layout_df$is_isolated <- layout_df$Degree == 0
layout_df$node_class <- ifelse(
  layout_df$is_isolated, "Isolated node", "Connected node"
)

if (any(layout_df$is_isolated)) {
  connected_xy <- layout_df[!layout_df$is_isolated, c("x", "y"), drop = FALSE]
  core_radius <- if (nrow(connected_xy) > 0) {
    max(sqrt(connected_xy$x^2 + connected_xy$y^2), na.rm = TRUE)
  } else {
    1
  }
  outer_radius <- max(core_radius * 1.35, 1)
  iso_idx <- which(layout_df$is_isolated)
  iso_angles <- seq(0, 2 * pi, length.out = length(iso_idx) + 1)[
    seq_along(iso_idx)
  ]
  layout_df$x[iso_idx] <- outer_radius * cos(iso_angles)
  layout_df$y[iso_idx] <- outer_radius * sin(iso_angles)
}

# 把调整后的坐标写回layout对象，供ggraph绘边。
layout_obj$x <- layout_df$x
layout_obj$y <- layout_df$y
layout_obj$Degree <- layout_df$Degree
layout_obj$is_isolated <- layout_df$is_isolated
layout_obj$node_class <- layout_df$node_class

morandi_node_colors <- c(
  "Connected node" = "#718CA4",
  "Isolated node" = "#B87970"
)

ppi_plot <- ggraph::ggraph(layout_obj) +
  ggraph::geom_edge_link(
    aes(width = confidence, alpha = confidence),
    color = "#8C8C8C", lineend = "round", show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(size = Degree, fill = node_class),
    shape = 21, color = "white", stroke = 0.35
  ) +
  ggrepel::geom_text_repel(
    data = layout_df,
    aes(x = x, y = y, label = label),
    size = 1.85, family = "sans", color = "#252525",
    box.padding = 0.25, point.padding = 0.18, point.size = 4,
    min.segment.length = 0, segment.color = "#9A9A9A",
    segment.size = 0.20, max.overlaps = Inf,
    force = 6, force_pull = 0.25, max.time = 10, max.iter = 100000,
    seed = random_seed
  ) +
  ggraph::scale_edge_width_continuous(
    range = c(0.25, 1.65),
    limits = c(confidence, 1),
    name = "STRING confidence"
  ) +
  ggraph::scale_edge_alpha_continuous(
    range = c(0.25, 0.75),
    limits = c(confidence, 1),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(2.3, 7.2),
    breaks = scales::pretty_breaks(n = 4),
    name = "Degree"
  ) +
  scale_fill_manual(
    values = morandi_node_colors,
    name = "Node class"
  ) +
  guides(
    edge_width = guide_legend(order = 1),
    size = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21, fill = "#718CA4", color = "white", linetype = 0
      )
    ),
    fill = guide_legend(order = 3, override.aes = list(size = 4))
  ) +
  labs(
    title = "GSE29272 STRING protein association network",
    subtitle = paste0(
      "Homo sapiens (9606) | STRING ", string_db$version,
      " | confidence ≥ ", confidence, " | no added nodes"
    ),
    caption = paste0(
      "Nodes: ", network_nodes, "  Edges: ", network_edges,
      "  Isolates: ", isolated_n,
      "  | Edge width: combined STRING confidence"
    )
  ) +
  theme_void(base_family = "sans", base_size = 7) +
  theme(
    plot.title = element_text(
      size = 10, face = "bold", color = "#222222", hjust = 0
    ),
    plot.subtitle = element_text(
      size = 7, color = "#555555", margin = margin(b = 5)
    ),
    plot.caption = element_text(
      size = 6.2, color = "#555555", hjust = 0, margin = margin(t = 5)
    ),
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    plot.margin = margin(8, 8, 8, 8),
    plot.background = element_rect(fill = "white", color = NA)
  )

figure_stem <- file.path(out_dir, "GSE29272_STRING_PPI_network")
ggsave(
  paste0(figure_stem, ".pdf"), ppi_plot,
  width = 183 / 25.4, height = 165 / 25.4,
  units = "in", device = grDevices::cairo_pdf, bg = "white"
)
ggsave(
  paste0(figure_stem, ".png"), ppi_plot,
  width = 183 / 25.4, height = 165 / 25.4,
  units = "in", dpi = 600, device = ragg::agg_png, bg = "white"
)

# ------------------------------------------------------------------------------
# 7. 参数、摘要和中文报告
# ------------------------------------------------------------------------------
parameters <- data.frame(
  parameter = c(
    "query_date", "species", "taxonomy_id", "STRING_database_version",
    "STRINGdb_R_package_version", "confidence_threshold",
    "STRING_score_threshold", "network_type", "added_nodes",
    "mapping_method", "layout", "random_seed", "evidence_channels"
  ),
  value = c(
    query_date, "Homo sapiens", taxonomy_id, string_db$version,
    as.character(packageVersion("STRINGdb")), confidence,
    score_threshold, "functional association", 0,
    "STRINGdb::map(takeFirst=TRUE)",
    "Fruchterman-Reingold; isolated nodes placed on outer ring",
    random_seed,
    paste(unname(evidence_channels), collapse = "; ")
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  parameters,
  file.path(out_dir, "GSE29272_STRING_PPI_parameters.csv")
)

summary_table <- data.frame(
  metric = c(
    "input_unique_genes", "mapped_genes", "mapped_STRING_proteins",
    "unmapped_genes", "network_nodes", "network_edges",
    "connected_nodes", "isolated_nodes", "connected_components",
    "network_density", "mean_degree", "max_degree",
    "minimum_observed_confidence", "maximum_observed_confidence",
    "extra_nodes_added"
  ),
  value = as.character(c(
    length(genes), mapped_gene_n, mapped_protein_n, nrow(unmapped),
    network_nodes, network_edges, connected_n, isolated_n, component_n,
    network_density, mean_degree, max_degree,
    if (network_edges > 0) min(interactions$confidence) else NA,
    if (network_edges > 0) max(interactions$confidence) else NA,
    0
  )),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  summary_table,
  file.path(out_dir, "GSE29272_STRING_PPI_summary.csv")
)

top_degree <- head(node_table[order(-node_table$Degree, node_table$gene), ], 10)
top_degree_text <- paste0(
  top_degree$gene, "（Degree=", top_degree$Degree, "）",
  collapse = "、"
)
unmapped_text <- if (nrow(unmapped) == 0) {
  "无"
} else {
  paste(unmapped$gene, collapse = "、")
}
isolated_text <- if (nrow(isolated_nodes) == 0) {
  "无"
} else {
  paste(isolated_nodes$gene, collapse = "、")
}

report_lines <- c(
  "# GSE29272 STRING PPI网络分析报告",
  "",
  paste0("- 查询日期：", query_date),
  paste0("- R版本：", R.version.string),
  paste0("- STRINGdb R包版本：", as.character(packageVersion("STRINGdb"))),
  paste0("- STRING数据库版本：", string_db$version),
  "",
  "## 1. 输入与基因—蛋白映射",
  "",
  paste0(
    "读取上一阶段得到的 ", length(genes),
    " 个唯一交集候选基因。基因名经过去除首尾空格并统一为大写的人类Gene Symbol；",
    "未重新筛选基因。"
  ),
  paste0(
    "STRINGdb成功映射 ", mapped_gene_n, " 个基因，对应 ",
    mapped_protein_n, " 个唯一STRING蛋白；未映射基因 ",
    nrow(unmapped), " 个。未映射基因：", unmapped_text, "。"
  ),
  "",
  "## 2. STRING查询参数与证据通道",
  "",
  paste0(
    "物种设为Homo sapiens（Taxonomy ID: ", taxonomy_id,
    "），STRING版本 ", string_db$version,
    "，最低置信度为 ", confidence, "（combined score ≥ ",
    score_threshold, "），不添加额外节点。"
  ),
  paste0(
    "网络类型为STRING functional association。combined score综合以下证据通道：",
    "neighborhood、gene fusion、cooccurrence、coexpression、experiments、",
    "database和textmining。每条边的通道分数与汇总统计已分别输出。"
  ),
  "STRING的combined score代表蛋白间功能关联置信度，不应全部解释为直接物理结合。",
  "",
  "## 3. 网络规模",
  "",
  paste0(
    "最终网络包含 ", network_nodes, " 个节点和 ", network_edges,
    " 条无向边；其中有连接节点 ", connected_n, " 个，孤立节点 ",
    isolated_n, " 个，共 ", component_n, " 个连通分量。"
  ),
  paste0(
    "网络密度为 ", sprintf("%.4f", network_density),
    "，平均Degree为 ", sprintf("%.2f", mean_degree),
    "，最大Degree为 ", max_degree, "。孤立节点：", isolated_text, "。"
  ),
  paste0("Degree最高的节点包括：", top_degree_text, "。"),
  "",
  "## 4. 网络图说明",
  "",
  paste0(
    "网络主体使用Fruchterman–Reingold力导向布局并设置随机种子",
    random_seed, "；没有互作边的孤立节点被均匀放置在外围以减少重叠。"
  ),
  paste0(
    "节点大小表示Degree；莫兰迪蓝灰色表示普通连接节点，陶土色表示孤立节点；",
    "边宽及透明度表示STRING combined confidence；全部节点使用避让标签。"
  ),
  "网络图仅用于展示本次输入基因之间的STRING关联，没有根据Degree或其他指标继续筛选基因。",
  "",
  "## 5. 输出文件",
  "",
  "- `GSE29272_STRING_gene_protein_mapping.csv`：基因—STRING蛋白映射表。",
  "- `GSE29272_STRING_unmapped_genes.csv`：未映射基因。",
  "- `GSE29272_STRING_isolated_nodes.csv`：成功映射但在阈值0.4下无内部互作边的节点。",
  "- `GSE29272_STRING_PPI_nodes.csv`：网络节点、Degree和节点类型。",
  "- `GSE29272_STRING_PPI_edges.csv`：网络边、combined score及7类证据通道分数。",
  "- `GSE29272_STRING_evidence_channels.csv`：证据通道支持边数与分数汇总。",
  "- `GSE29272_STRING_PPI_network.pdf/png`：矢量PDF和600 dpi PNG网络图。",
  "- `GSE29272_STRING_PPI_parameters.csv`：数据库、物种、阈值、日期和随机种子。",
  "- `GSE29272_STRING_PPI_summary.csv`：网络规模与拓扑摘要。",
  "",
  "## 6. 可复现性",
  "",
  "完整R代码位于 `scripts/GSE29272_STRING_PPI_analysis.R`；",
  "STRING下载缓存、运行日志和sessionInfo均保存在 `results/PPI`。"
)
writeLines(
  report_lines,
  file.path(out_dir, "GSE29272_STRING_PPI_report.md"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "GSE29272_STRING_PPI_sessionInfo.txt")
)

cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
