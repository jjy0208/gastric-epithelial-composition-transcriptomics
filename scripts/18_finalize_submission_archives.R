# 为既有 v2.2.0 投稿目录创建上传压缩包，并刷新 SHA-256 清单。
suppressPackageStartupMessages({
  library(digest)
  library(zip)
})

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) != 1L) stop("请使用 Rscript 运行本脚本。")
  normalizePath(sub("^--file=", "", f), winslash = "/", mustWork = TRUE)
}

project_dir <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  winslash = "/",
  mustWork = TRUE
)
package_name <- "gastric-epithelial-composition-transcriptomics-v2.2.0-submission"
package_dir <- file.path(project_dir, "release", package_name)
if (!dir.exists(package_dir)) stop("投稿目录不存在：", package_dir)

source_zip <- file.path(package_dir, "submission", "SourceData.zip")
tables_zip <- file.path(package_dir, "submission", "SupplementaryTables.zip")
zip::zipr(
  source_zip,
  list.files(file.path(package_dir, "source_data"), full.names = TRUE),
  root = file.path(package_dir, "source_data"),
  include_directories = FALSE
)
zip::zipr(
  tables_zip,
  list.files(file.path(package_dir, "supplementary_tables"), full.names = TRUE),
  root = file.path(package_dir, "supplementary_tables"),
  include_directories = FALSE
)

manifest_path <- file.path(package_dir, "SHA256_manifest.csv")
all_files <- list.files(
  package_dir, recursive = TRUE, full.names = TRUE
)
all_files <- setdiff(
  normalizePath(all_files, winslash = "/", mustWork = TRUE),
  normalizePath(manifest_path, winslash = "/", mustWork = FALSE)
)
root_norm <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
manifest <- data.frame(
  relative_path = substring(all_files, nchar(root_norm) + 2L),
  bytes = file.info(all_files)$size,
  sha256 = vapply(
    all_files,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$relative_path), ]
write.csv(
  manifest, manifest_path,
  row.names = FALSE, fileEncoding = "UTF-8"
)

archive_path <- file.path(
  project_dir, "release",
  paste0(package_name, "-full.zip")
)
if (file.exists(archive_path)) {
  stop("总压缩包已存在，为避免覆盖请人工核对：", archive_path)
}
zip::zipr(
  archive_path,
  package_name,
  root = dirname(package_dir),
  include_directories = FALSE
)
archive_sha256 <- digest::digest(
  archive_path, algo = "sha256", file = TRUE, serialize = FALSE
)
writeLines(
  paste(archive_sha256, basename(archive_path)),
  paste0(archive_path, ".sha256.txt"),
  useBytes = TRUE
)

cat("SourceData.zip：", source_zip, "\n")
cat("SupplementaryTables.zip：", tables_zip, "\n")
cat("投稿总包：", archive_path, "\n")
cat("SHA-256：", archive_sha256, "\n")
