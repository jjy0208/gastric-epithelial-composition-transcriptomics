# 审核正文参考文献的近五年比例、未引用条目、编号连续性和 DOI 可解析性。
suppressPackageStartupMessages(library(jsonlite))

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) != 1L) stop("请使用 Rscript 运行本脚本。")
  normalizePath(sub("^--file=", "", f), winslash = "/", mustWork = TRUE)
}

project_dir <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  winslash = "/", mustWork = TRUE
)
manuscript_path <- file.path(
  project_dir, "report", "manuscript", "manuscript_en.md"
)
lines <- readLines(manuscript_path, warn = FALSE, encoding = "UTF-8")
reference_heading <- which(lines == "## References")
if (length(reference_heading) != 1L) stop("未找到唯一的 References 标题。")
body <- paste(lines[seq_len(reference_heading - 1L)], collapse = "\n")
reference_lines <- lines[(reference_heading + 1L):length(lines)]
reference_lines <- reference_lines[nzchar(trimws(reference_lines))]

reference_number <- as.integer(sub("^([0-9]+)\\..*$", "\\1", reference_lines))
if (!identical(reference_number, seq_along(reference_lines))) {
  stop("参考文献编号不连续。")
}
year <- as.integer(sub(".*\\(([12][0-9]{3})\\).*", "\\1", reference_lines))
doi <- sub(
  ".*https://doi\\.org/([^[:space:]]+).*$", "\\1", reference_lines
)
if (any(doi == reference_lines)) stop("存在未提取到 DOI 的参考文献。")

expand_citation <- function(x) {
  x <- gsub("^\\[|\\]$", "", x)
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- integer()
  for (piece in pieces) {
    if (grepl("[–-]", piece)) {
      ends <- as.integer(strsplit(piece, "[–-]")[[1]])
      out <- c(out, seq(ends[1], ends[2]))
    } else {
      out <- c(out, as.integer(piece))
    }
  }
  out
}
citation_tokens <- regmatches(
  body,
  gregexpr("\\[[0-9]+(?:[,–-][0-9]+)*\\]", body, perl = TRUE)
)[[1]]
cited_numbers <- sort(unique(unlist(lapply(citation_tokens, expand_citation))))
uncited <- setdiff(reference_number, cited_numbers)
undefined <- setdiff(cited_numbers, reference_number)

crossref_title <- rep(NA_character_, length(doi))
crossref_year <- rep(NA_integer_, length(doi))
crossref_status <- rep("not_checked", length(doi))
for (i in seq_along(doi)) {
  endpoint <- paste0(
    "https://api.crossref.org/works/",
    URLencode(doi[i], reserved = TRUE),
    "?mailto=364929601%40qq.com"
  )
  record <- tryCatch(
    jsonlite::fromJSON(endpoint, simplifyVector = FALSE),
    error = function(e) e
  )
  if (inherits(record, "error")) {
    crossref_status[i] <- paste0("error: ", conditionMessage(record))
  } else {
    item <- record$message
    crossref_title[i] <- if (length(item$title)) item$title[[1]] else NA_character_
    issued <- item$issued$`date-parts`[[1]]
    crossref_year[i] <- if (length(issued)) as.integer(issued[[1]]) else NA_integer_
    crossref_status[i] <- "resolved"
  }
  Sys.sleep(0.1)
}

cutoff_year <- 2021L
old_rationale <- rep("", length(reference_number))
names(old_rationale) <- reference_number
old_rationale["2"] <- "Foundational TCGA gastric-cancer molecular classification; retained as the original source."
old_rationale["12"] <- "Original BisqueRNA deconvolution method used by the analysis; no equivalent recent primary method citation."
old_rationale["14"] <- "Original source publication for the GSE29272 discovery cohort."
old_rationale["16"] <- "Original GEOquery software paper for the exact download interface used."
old_rationale["17"] <- "Primary limma methods paper for the exact differential-expression framework used."
old_rationale["18"] <- "Original WGCNA software and method paper for the exact network framework used."
old_rationale["20"] <- "Original paper defining the MSigDB Hallmark collection used in GSEA."

audit <- data.frame(
  reference_number = reference_number,
  manuscript_year = year,
  recent_2021_2026 = year >= cutoff_year,
  cited_in_text = reference_number %in% cited_numbers,
  doi = doi,
  crossref_status = crossref_status,
  crossref_year = crossref_year,
  crossref_title = crossref_title,
  old_reference_rationale = unname(old_rationale[as.character(reference_number)]),
  stringsAsFactors = FALSE
)

old_without_rationale <- audit[
  !audit$recent_2021_2026 & !nzchar(audit$old_reference_rationale),
  ,
  drop = FALSE
]
year_mismatch <- audit[
  audit$crossref_status == "resolved" &
    !is.na(audit$crossref_year) &
    abs(audit$manuscript_year - audit$crossref_year) > 1L,
  ,
  drop = FALSE
]

report_dir <- file.path(project_dir, "report", "manuscript")
write.csv(
  audit,
  file.path(report_dir, "reference_recency_and_doi_audit.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
recent_n <- sum(audit$recent_2021_2026)
old_n <- nrow(audit) - recent_n
report <- c(
  "# Reference recency and DOI audit",
  "",
  paste0("- Audit date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- Total references: ", nrow(audit)),
  paste0(
    "- References published in 2021–2026: ",
    recent_n, "/", nrow(audit), " (",
    sprintf("%.1f", 100 * recent_n / nrow(audit)), "%)"
  ),
  paste0("- References older than 2021: ", old_n),
  paste0("- Crossref-resolved DOIs: ", sum(audit$crossref_status == "resolved"),
         "/", nrow(audit)),
  paste0("- Uncited reference numbers: ",
         if (length(uncited)) paste(uncited, collapse = ", ") else "None"),
  paste0("- Undefined in-text citation numbers: ",
         if (length(undefined)) paste(undefined, collapse = ", ") else "None"),
  "",
  "## Retained older references and rationale",
  "",
  paste0(
    "- [", audit$reference_number[!audit$recent_2021_2026], "] ",
    audit$manuscript_year[!audit$recent_2021_2026], ": ",
    audit$old_reference_rationale[!audit$recent_2021_2026]
  ),
  "",
  "## Verification exceptions",
  "",
  paste0(
    "- Old references without a documented rationale: ",
    if (nrow(old_without_rationale)) {
      paste(old_without_rationale$reference_number, collapse = ", ")
    } else {
      "None"
    }
  ),
  paste0(
    "- Crossref online-year deviations greater than one year: ",
    if (nrow(year_mismatch)) {
      paste(year_mismatch$reference_number, collapse = ", ")
    } else {
      "None"
    }
  )
)
writeLines(
  report,
  file.path(report_dir, "reference_recency_and_doi_audit.md"),
  useBytes = TRUE
)

if (length(uncited) || length(undefined) || nrow(old_without_rationale) ||
    nrow(year_mismatch) ||
    any(audit$crossref_status != "resolved")) {
  print(audit[audit$crossref_status != "resolved", , drop = FALSE])
  stop("参考文献审计未通过，请查看 reference_recency_and_doi_audit.md。")
}
cat("参考文献审计通过：", recent_n, "/", nrow(audit),
    "篇为2021–2026年文献；其余均有不可替代性说明。\n")
