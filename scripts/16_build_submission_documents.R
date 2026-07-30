# 生成投稿信和合并版补充材料 Word 文件。
suppressPackageStartupMessages({
  library(officer)
  library(data.table)
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
manuscript_dir <- file.path(project_dir, "report", "manuscript")

base_doc <- function() {
  doc <- read_docx()
  doc <- docx_set_paragraph_style(
    doc, style_id = "SubmissionBody", style_name = "Submission Body",
    base_on = "Normal",
    fp_p = fp_par(
      text.align = "justify", line_spacing = 1.25,
      padding.top = 0, padding.bottom = 7
    ),
    fp_t = fp_text(font.family = "Arial", font.size = 10.5)
  )
  doc <- docx_set_paragraph_style(
    doc, style_id = "SubmissionTitle", style_name = "Submission Title",
    base_on = "Normal",
    fp_p = fp_par(
      text.align = "center", line_spacing = 1.15,
      padding.top = 0, padding.bottom = 12, keep_with_next = TRUE
    ),
    fp_t = fp_text(font.family = "Arial", font.size = 17, bold = TRUE)
  )
  doc <- docx_set_paragraph_style(
    doc, style_id = "SubmissionH1", style_name = "Submission H1",
    base_on = "Normal",
    fp_p = fp_par(
      text.align = "left", line_spacing = 1.15,
      padding.top = 12, padding.bottom = 6, keep_with_next = TRUE
    ),
    fp_t = fp_text(font.family = "Arial", font.size = 13, bold = TRUE)
  )
  doc <- docx_set_paragraph_style(
    doc, style_id = "SubmissionH2", style_name = "Submission H2",
    base_on = "Normal",
    fp_p = fp_par(
      text.align = "left", line_spacing = 1.15,
      padding.top = 8, padding.bottom = 4, keep_with_next = TRUE
    ),
    fp_t = fp_text(font.family = "Arial", font.size = 11.5, bold = TRUE)
  )
  doc
}

clean_inline <- function(x) {
  x <- gsub("`", "", x, fixed = TRUE)
  x <- gsub("\\*\\*", "", x)
  x <- gsub("\\*", "", x)
  x <- gsub("\\^([^\\^]+)\\^", "\\1", x)
  x
}

add_markdown <- function(doc, path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  for (line in lines) {
    if (!nzchar(trimws(line))) next
    if (grepl("^# ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^# ", "", line)),
        style = "Submission Title"
      )
    } else if (grepl("^## ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^## ", "", line)),
        style = "Submission H1"
      )
    } else if (grepl("^### ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^### ", "", line)),
        style = "Submission H2"
      )
    } else if (grepl("^- ", line)) {
      doc <- body_add_par(
        doc, paste0("• ", clean_inline(sub("^- ", "", line))),
        style = "Submission Body"
      )
    } else if (grepl("^[0-9]+\\. ", line)) {
      doc <- body_add_par(doc, clean_inline(line), style = "Submission Body")
    } else {
      doc <- body_add_par(doc, clean_inline(line), style = "Submission Body")
    }
  }
  doc
}

# Cover letter
cover <- base_doc()
cover <- add_markdown(
  cover,
  file.path(manuscript_dir, "cover_letter_scientific_reports.md")
)
print(
  cover,
  target = file.path(manuscript_dir, "cover_letter_scientific_reports.docx")
)

# Supplementary Information: narrative plus embedded figure files.
si <- base_doc()
si_lines <- readLines(
  file.path(manuscript_dir, "supplementary_information.md"),
  warn = FALSE, encoding = "UTF-8"
)
figure_heading_index <- which(grepl("^### Supplementary Figure S[0-9]+", si_lines))
table_heading_index <- which(si_lines == "## Supplementary Tables")

# Add title, Supplementary Methods and introductory text up to the first figure.
tmp_intro <- tempfile(fileext = ".md")
writeLines(si_lines[seq_len(figure_heading_index[1] - 1L)], tmp_intro, useBytes = TRUE)
si <- add_markdown(si, tmp_intro)
unlink(tmp_intro)

figure_specs <- list(
  list(
    title = "Supplementary Figure S1 | Discovery sample quality control",
    legend = "Standardized sample-connectivity clustering identified five samples below Z.k=−2.5; these samples were excluded only from WGCNA network construction.",
    files = file.path(project_dir, "results", "WGCNA",
                      "GSE29272_WGCNA_sample_clustering.png")
  ),
  list(
    title = "Supplementary Figure S2 | Soft-threshold assessment",
    legend = "No tested power reached signed scale-free R²≥0.85. The documented fallback β=28 maximized signed fit over the tested grid.",
    files = file.path(project_dir, "results", "WGCNA",
                      "GSE29272_WGCNA_soft_threshold.png")
  ),
  list(
    title = "Supplementary Figure S3 | Module dendrogram and module–trait relationship",
    legend = "Network modules and their association with tumor status. The high fallback power requires exploratory interpretation.",
    files = file.path(
      project_dir, "results", "WGCNA",
      c(
        "GSE29272_WGCNA_module_dendrogram.png",
        "GSE29272_WGCNA_module_trait_heatmap.png"
      )
    )
  ),
  list(
    title = "Supplementary Figure S4 | Differential-expression heat map",
    legend = "Heat map of thresholded paired differential-expression results; the inferential unit was the patient pair.",
    files = file.path(project_dir, "results", "DEG",
                      "GSE29272_DEG_heatmap.png")
  ),
  list(
    title = "Supplementary Figure S5 | Hallmark GSEA",
    legend = "Overview of significant Hallmark pathways from the complete ranked differential-expression result.",
    files = file.path(project_dir, "results", "GSEA",
                      "GSE29272_GSEA_overview_dotplot.png")
  ),
  list(
    title = "Supplementary Figure S6 | STRING PPI network",
    legend = "STRING v12.0 network at combined confidence score≥0.4 without added nodes; the network was descriptive.",
    files = file.path(project_dir, "results", "PPI",
                      "GSE29272_STRING_PPI_network.png")
  ),
  list(
    title = "Supplementary Figure S7 | Original intersection visualization",
    legend = "Intersection of 122 significant DEGs and 192 M04 genes yielded 73 candidates.",
    files = file.path(project_dir, "results", "HubGene",
                      "GSE29272_DEG_WGCNA_Venn.png")
  ),
  list(
    title = "Supplementary Figure S8 | GSE270680 major-cell-type localization",
    legend = "Each point represents a patient–tissue–major-cell-type pseudobulk with at least 20 cells, using author-defined annotations.",
    files = file.path(
      project_dir, "results", "PaperValidation", "SingleCellReplication",
      "GSE270680_major_celltype_localization.png"
    )
  ),
  list(
    title = "Supplementary Figure S9 | Donor-level spatial concordance",
    legend = "Each point is a donor-level mean of section-specific Spearman correlations after removing candidate-overlapping genes from comparator marker sets.",
    files = file.path(
      project_dir, "results", "PaperValidation", "Spatial",
      "GSE270678_spatial_correlations.png"
    )
  )
)

for (spec in figure_specs) {
  if (any(!file.exists(spec$files))) {
    stop("缺少补充图：", paste(basename(spec$files[!file.exists(spec$files)]),
                               collapse = ", "))
  }
  si <- body_add_break(si)
  si <- body_add_par(si, spec$title, style = "Submission H1")
  si <- body_add_par(si, spec$legend, style = "Submission Body")
  for (img in spec$files) {
    si <- body_add_img(si, src = img, width = 6.2, height = 4.7,
                       style = "centered")
  }
}

# Embed the new cohort-comparison tables so reviewers can inspect them directly.
cohort_characteristics <- fread(
  file.path(
    project_dir, "results", "PaperValidation", "SingleCellReplication",
    "singlecell_cohort_characteristics.csv"
  ),
  encoding = "UTF-8"
)
cohort_tests <- fread(
  file.path(
    project_dir, "results", "PaperValidation", "SingleCellReplication",
    "singlecell_cohort_comparison_tests.csv"
  ),
  encoding = "UTF-8"
)
cohort_tests[, `:=`(
  p_value = ifelse(is.na(p_value), "", formatC(p_value, digits = 3, format = "fg")),
  FDR = ifelse(is.na(FDR), "", formatC(FDR, digits = 3, format = "fg"))
)]
cohort_tests[, interpretation := NULL]

si <- body_add_break(si)
si <- body_add_par(
  si,
  "Supplementary Table S6 | Characteristics of patients eligible for paired epithelial analysis",
  style = "Submission H1"
)
si <- body_add_par(
  si,
  paste(
    "Values are n/available unless otherwise stated.",
    "GSE206785 clinical data were unavailable for one eligible patient.",
    "Gastric position was not reported for GSE270680."
  ),
  style = "Submission Body"
)
si <- body_add_table(
  si, value = as.data.frame(cohort_characteristics),
  first_row = TRUE
)

si <- body_add_break(si)
si <- body_add_par(
  si,
  "Supplementary Table S7 | Exploratory comparison of eligible analysis populations",
  style = "Submission H1"
)
si <- body_add_par(
  si,
  paste(
    "These tests describe small selected analysis subsets and do not identify",
    "the cause of the cross-cohort transcriptomic difference.",
    "FDR was controlled across five comparisons."
  ),
  style = "Submission Body"
)
si <- body_add_table(
  si, value = as.data.frame(cohort_tests),
  first_row = TRUE
)

# Add the tables inventory and limitations checklist from the existing Markdown.
# Insert the section heading explicitly after the page break; this avoids a Word
# pagination artefact that can place the first heading partly outside the page.
tmp_tail <- tempfile(fileext = ".md")
writeLines(si_lines[(table_heading_index + 1L):length(si_lines)], tmp_tail, useBytes = TRUE)
si <- body_add_break(si)
si <- body_add_par(si, "Supplementary Tables", style = "Submission H1")
si <- add_markdown(si, tmp_tail)
unlink(tmp_tail)

print(
  si,
  target = file.path(manuscript_dir, "supplementary_information_v2.docx")
)
cat("投稿信和补充材料 DOCX 已生成。\n")
