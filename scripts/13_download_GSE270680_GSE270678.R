# 下载第二个独立单细胞队列 GSE270680、作者注释 Seurat 对象，
# 以及配套空间转录组队列 GSE270678。
# 数据来源：NCBI GEO 官方 supplementary files 和作者公开的 Mendeley Data。
# 说明：这些是过滤后表达矩阵、作者注释、空间坐标和组织图像，不是 FASTQ。

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) stop("请使用 Rscript 运行本脚本。")
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_dir <- normalizePath(file.path(dirname(script_path), ".."),
                             winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "raw_data")

datasets <- data.frame(
  accession = c("GSE270680", "GSE270678"),
  modality = c("single-cell RNA-seq", "10x spatial transcriptomics"),
  expected_bytes = c(3336919040, 775833600),
  url = c(
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/",
      "GSE270680/suppl/GSE270680_RAW.tar"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/",
      "GSE270678/suppl/GSE270678_RAW.tar"
    )
  )
)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
log_rows <- vector("list", nrow(datasets))

for (i in seq_len(nrow(datasets))) {
  acc <- datasets$accession[i]
  target_dir <- file.path(raw_dir, acc)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(target_dir, paste0(acc, "_RAW.tar"))

  if (!file.exists(target) ||
      is.na(file.info(target)$size) ||
      file.info(target)$size != datasets$expected_bytes[i]) {
    message("正在下载 ", acc, "：", datasets$url[i])
    download.file(
      datasets$url[i],
      destfile = target,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    )
  }

  actual_bytes <- file.info(target)$size
  if (is.na(actual_bytes) || actual_bytes != datasets$expected_bytes[i]) {
    stop(
      acc, " 文件大小校验失败：expected=", datasets$expected_bytes[i],
      ", actual=", actual_bytes
    )
  }

  log_rows[[i]] <- data.frame(
    accession = acc,
    query_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    method = "NCBI GEO supplementary file; R download.file(method='libcurl')",
    modality = datasets$modality[i],
    file_name = basename(target),
    file_path = normalizePath(target, winslash = "/", mustWork = TRUE),
    bytes = actual_bytes,
    md5 = unname(tools::md5sum(target)),
    sha256 = NA_character_,
    source_url = datasets$url[i]
  )
}

log_df <- do.call(rbind, log_rows)

# 作者公开的完整 Seurat 对象保留了 majorCluster/subCluster 注释，
# 用于避免根据候选基因重新注释细胞。
sc_url <- paste0(
  "https://data.mendeley.com/public-files/datasets/559mchb37p/files/",
  "be10b62d-c0e4-4208-9cb1-8fee71f3c1ab/file_downloaded"
)
sc_target <- file.path(raw_dir, "GSE270680", "sc.rds")
sc_expected_bytes <- 2938913974
sc_expected_sha256 <- paste0(
  "3aa4ec260c9552acbf377a697dd962ac8",
  "287aa0c79c9b435714b1697bf7e6265"
)
if (!file.exists(sc_target) ||
    is.na(file.info(sc_target)$size) ||
    file.info(sc_target)$size != sc_expected_bytes) {
  message("正在下载作者公开的 sc.rds：", sc_url)
  download.file(
    sc_url, destfile = sc_target, mode = "wb",
    method = "libcurl", quiet = FALSE
  )
}
if (file.info(sc_target)$size != sc_expected_bytes) {
  stop("sc.rds 文件大小校验失败。")
}
sc_sha256 <- digest::digest(
  sc_target, file = TRUE, algo = "sha256", serialize = FALSE
)
if (!identical(tolower(sc_sha256), sc_expected_sha256)) {
  stop("sc.rds SHA-256 校验失败。")
}
log_df <- rbind(
  log_df,
  data.frame(
    accession = "GSE270680",
    query_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    method = paste0(
      "Mendeley Data author-annotated Seurat object; ",
      "R download.file(method='libcurl')"
    ),
    modality = "single-cell RNA-seq author annotation",
    file_name = basename(sc_target),
    file_path = normalizePath(sc_target, winslash = "/", mustWork = TRUE),
    bytes = file.info(sc_target)$size,
    md5 = unname(tools::md5sum(sc_target)),
    sha256 = sc_sha256,
    source_url = sc_url
  )
)
write.table(
  log_df,
  file.path(raw_dir, "GSE270680_GSE270678_download_log.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8"
)

message("下载及文件大小校验完成。")
print(log_df)
