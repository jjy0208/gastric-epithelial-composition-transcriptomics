# 将已审计的 Markdown 主文、图注和六张主图组装为可编辑 Word 稿件。
suppressPackageStartupMessages(library(officer))

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
figure_dir <- file.path(project_dir, "results", "ManuscriptFigures")
output_path <- file.path(
  manuscript_dir,
  "gastric_epithelial_composition_manuscript_v2.docx"
)

doc <- read_docx()
doc <- docx_set_paragraph_style(
  doc, style_id = "ManuscriptBody", style_name = "Manuscript Body",
  base_on = "Normal",
  fp_p = fp_par(
    text.align = "justify", line_spacing = 1.333,
    padding.top = 0, padding.bottom = 8
  ),
  fp_t = fp_text(font.family = "Calibri", font.size = 11, color = "#000000")
)
doc <- docx_set_paragraph_style(
  doc, style_id = "ManuscriptTitle", style_name = "Manuscript Title",
  base_on = "Normal",
  fp_p = fp_par(
    text.align = "center", line_spacing = 1.15,
    padding.top = 0, padding.bottom = 14, keep_with_next = TRUE
  ),
  fp_t = fp_text(
    font.family = "Calibri", font.size = 20,
    bold = TRUE, color = "#203748"
  )
)
doc <- docx_set_paragraph_style(
  doc, style_id = "ManuscriptH1", style_name = "Manuscript H1",
  base_on = "Normal",
  fp_p = fp_par(
    text.align = "left", line_spacing = 1.15,
    padding.top = 18, padding.bottom = 10, keep_with_next = TRUE
  ),
  fp_t = fp_text(
    font.family = "Calibri", font.size = 16,
    bold = TRUE, color = "#2E74B5"
  )
)
doc <- docx_set_paragraph_style(
  doc, style_id = "ManuscriptH2", style_name = "Manuscript H2",
  base_on = "Normal",
  fp_p = fp_par(
    text.align = "left", line_spacing = 1.15,
    padding.top = 12, padding.bottom = 6, keep_with_next = TRUE
  ),
  fp_t = fp_text(
    font.family = "Calibri", font.size = 13,
    bold = TRUE, color = "#2E74B5"
  )
)
doc <- docx_set_paragraph_style(
  doc, style_id = "ManuscriptLegend", style_name = "Manuscript Legend",
  base_on = "Normal",
  fp_p = fp_par(
    text.align = "justify", line_spacing = 1.15,
    padding.top = 0, padding.bottom = 5
  ),
  fp_t = fp_text(font.family = "Calibri", font.size = 10, color = "#000000")
)

clean_inline <- function(x) {
  x <- gsub("`", "", x, fixed = TRUE)
  x <- gsub("\\*\\*", "", x)
  x <- gsub("\\*", "", x)
  x <- gsub("\\^([^\\^]+)\\^", "\\1", x)
  x
}

add_markdown <- function(
    doc, path, skip_first_title = FALSE,
    body_style = "Manuscript Body") {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (skip_first_title && length(lines) > 0 && grepl("^# ", lines[1])) {
    abstract_start <- which(lines == "## Abstract")
    if (length(abstract_start) == 1L) {
      lines <- lines[abstract_start:length(lines)]
    } else {
      lines <- lines[-1]
    }
  }
  for (line in lines) {
    if (!nzchar(trimws(line))) next
    if (grepl("^# ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^# ", "", line)),
        style = "Manuscript Title"
      )
    } else if (grepl("^## ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^## ", "", line)),
        style = "Manuscript H1"
      )
    } else if (grepl("^### ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^### ", "", line)),
        style = "Manuscript H2"
      )
    } else if (grepl("^> ", line)) {
      doc <- body_add_fpar(
        doc,
        fpar(
          ftext(
            clean_inline(sub("^> ", "", line)),
            fp_text(
              font.family = "Calibri", font.size = 10.5,
              italic = TRUE, color = "#555555"
            )
          ),
          fp_p = fp_par(
            text.align = "left", line_spacing = 1.2,
            padding = 8, shading.color = "#F4F6F9"
          )
        )
      )
    } else if (grepl("^- ", line)) {
      doc <- body_add_par(
        doc, clean_inline(sub("^- ", "", line)),
        style = "List Bullet"
      )
    } else {
      doc <- body_add_par(doc, clean_inline(line), style = body_style)
    }
  }
  doc
}

# 标题页
doc <- body_add_par(
  doc,
  paste(
    "Independent single-cell and spatial validation reveals mixed compositional",
    "and epithelial contributions to a gastric cancer transcriptomic module"
  ),
  style = "Manuscript Title"
)
doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      "Original Article | Public-data transcriptomic reanalysis",
      fp_text(
        font.family = "Calibri", font.size = 12,
        italic = TRUE, color = "#666666"
      )
    ),
    fp_p = fp_par(text.align = "center", padding.bottom = 18)
  )
)
doc <- body_add_par(
  doc,
  "Junyi Jia¹, Mingming Zhu¹, Mingxiong Zhang¹, Likun Luan¹ & Changlong Yang¹,*",
  style = "Manuscript Body"
)
doc <- body_add_par(
  doc,
  "¹The Third Affiliated Hospital of Kunming Medical University/Yunnan Cancer Hospital",
  style = "Manuscript Body"
)
doc <- body_add_par(
  doc,
  paste(
    "Correspondence: Changlong Yang | 965084451@qq.com |",
    "ORCID: https://orcid.org/0009-0001-6553-4711"
  ),
  style = "Manuscript Body"
)
doc <- body_add_par(
  doc,
  paste0(
    "Non-corresponding author email addresses are intentionally omitted from ",
    "the public repository copy of this script. ",
    "Junyi Jia, ORCID 0009-0002-0277-8521; Mingming Zhu; ",
    "Mingxiong Zhang, ORCID 0000-0002-8673-6507; Likun Luan."
  ),
  style = "Manuscript Body"
)
doc <- body_add_break(doc)

# 主文和图注
doc <- add_markdown(
  doc, file.path(manuscript_dir, "manuscript_en.md"),
  skip_first_title = TRUE
)
doc <- body_add_break(doc)
doc <- add_markdown(
  doc, file.path(manuscript_dir, "figure_legends.md"),
  skip_first_title = TRUE,
  body_style = "Manuscript Legend"
)

# 六张主图附于稿件末尾
doc <- body_add_break(doc)
doc <- body_add_par(doc, "Main figures", style = "Manuscript H1")
figure_specs <- data.frame(
  title = c(
    "Figure 1 | Discovery",
    "Figure 2 | Bulk replication",
    "Figure 3 | First single-cell cohort",
    "Figure 4 | Composition adjustment",
    "Figure 5 | Independent single-cell validation",
    "Figure 6 | Spatial validation"
  ),
  file = file.path(
    figure_dir,
    c(
      "Figure1_discovery.png",
      "Figure2_bulk_replication.png",
      "Figure3_single_cell_localization.png",
      "Figure4_composition_adjustment.png",
      "Figure5_independent_single_cell_validation.png",
      "Figure6_spatial_validation.png"
    )
  ),
  height = c(3.90, 4.44, 4.88, 4.44, 5.71, 6.34),
  stringsAsFactors = FALSE
)
if (any(!file.exists(figure_specs$file))) {
  stop("缺少主图：", paste(basename(figure_specs$file[!file.exists(figure_specs$file)]),
                           collapse = ", "))
}
for (i in seq_len(nrow(figure_specs))) {
  if (i > 1) doc <- body_add_break(doc)
  doc <- body_add_par(doc, figure_specs$title[i], style = "Manuscript H2")
  doc <- body_add_img(
    doc, src = figure_specs$file[i],
    width = 6.5, height = figure_specs$height[i],
    style = "centered"
  )
}

header <- block_list(
  fpar(
    ftext(
      "Gastric epithelial composition-aware transcriptomics",
      fp_text(font.family = "Calibri", font.size = 9, color = "#777777")
    ),
    fp_p = fp_par(text.align = "left", padding.bottom = 2)
  )
)
footer <- block_list(
  fpar(
    ftext(
      "Gastric epithelial transcriptomic module  |  Page ",
      fp_text(font.family = "Calibri", font.size = 9, color = "#777777")
    ),
    run_word_field("PAGE"),
    fp_p = fp_par(text.align = "right")
  )
)
sec <- prop_section(
  page_size = page_size(orient = "portrait", width = 8.5, height = 11),
  page_margins = page_mar(
    top = 1, bottom = 1, left = 1, right = 1,
    header = 0.492, footer = 0.492
  ),
  header_default = header,
  footer_default = footer,
  type = "continuous"
)
doc <- body_set_default_section(doc, sec)
print(doc, target = output_path)
cat("Word 稿件已生成：", output_path, "\n")
