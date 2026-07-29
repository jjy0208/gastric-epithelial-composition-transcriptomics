# 下载 GSE79973 和 GSE29272 的 GEO Series Matrix
# Series Matrix 同时包含处理后的表达矩阵和样本元数据。
# 本脚本只负责下载与完整性检查，不执行差异表达分析。

options(timeout = max(1200, getOption("timeout")))

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg), winslash = "/"))
  }
  stop("请使用 Rscript 运行本脚本，以便自动定位项目目录。")
}

script_path <- get_script_path()
project_dir <- normalizePath(
  file.path(dirname(script_path), ".."),
  winslash = "/",
  mustWork = TRUE
)
raw_dir <- file.path(project_dir, "raw_data")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase", ask = FALSE, update = FALSE)
}

datasets <- data.frame(
  accession = c("GSE79973", "GSE29272"),
  prefix = c("GSE79nnn", "GSE29nnn"),
  stringsAsFactors = FALSE
)

download_one <- function(accession, prefix) {
  file_name <- paste0(accession, "_series_matrix.txt.gz")
  file_path <- file.path(raw_dir, file_name)
  method <- "GEOquery::getGEO(GSEMatrix=TRUE)"
  geo_object <- NULL

  geo_object <- tryCatch(
    GEOquery::getGEO(
      accession,
      GSEMatrix = TRUE,
      getGPL = FALSE,
      destdir = raw_dir
    ),
    error = function(e) {
      message(accession, " 的 GEOquery 下载或解析失败：", conditionMessage(e))
      NULL
    }
  )

  valid_object <- FALSE
  if (!is.null(geo_object) && length(geo_object) >= 1L) {
    valid_object <- all(vapply(
      geo_object,
      function(eset) {
        nrow(Biobase::exprs(eset)) > 0L &&
          ncol(Biobase::exprs(eset)) > 0L &&
          nrow(Biobase::pData(eset)) == ncol(Biobase::exprs(eset))
      },
      logical(1)
    ))
  }

  if (!valid_object) {
    # GEO 页面显示这两套数据的 processed data 位于 Sample table；
    # 因此回退下载官方 Series Matrix，而不是下载原始 CEL 压缩包。
    method <- "NCBI GEO 官方 Series Matrix 直接下载（GEOquery 回退）"
    url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s",
      prefix, accession, file_name
    )
    download.file(url, file_path, mode = "wb", method = "libcurl", quiet = FALSE)
  }

  if (!file.exists(file_path) || file.info(file_path)$size <= 0L) {
    stop(accession, " 下载后文件不存在或为空。")
  }

  data.frame(
    数据集编号 = accession,
    下载时间 = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    下载方式 = method,
    下载文件名称 = file_name,
    文件路径 = normalizePath(file_path, winslash = "/"),
    文件大小_字节 = file.info(file_path)$size,
    MD5 = unname(tools::md5sum(file_path)),
    check.names = FALSE
  )
}

records <- Map(
  download_one,
  datasets$accession,
  datasets$prefix
)
download_log <- do.call(rbind, records)

log_path <- file.path(raw_dir, "download_log.txt")
write.table(
  download_log,
  file = log_path,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  fileEncoding = "UTF-8"
)

print(download_log)
message("下载日志：", normalizePath(log_path, winslash = "/"))
