# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
      library(ggplot2)
      library(patchwork)
      library(ggrastr)
      library(scales)
      library(RColorBrewer)
      library(cowplot)
  })
  args <- commandArgs(trailingOnly = FALSE)
  entry <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
  code_root <- normalizePath(file.path(dirname(entry), "../.."), winslash = "/")
  package_root <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(package_root)) stop("Set DBSEQ_OUTPUT_ROOT outside the code package.")
  output_resolved <- normalizePath(package_root, winslash = "/", mustWork = FALSE)
  if (identical(tolower(output_resolved), tolower(code_root)) || startsWith(tolower(output_resolved), paste0(tolower(code_root), 
      "/"))) stop("Outputs must remain outside the public code package")
  workspace_root <- package_root
  data_root <- Sys.getenv("DBSEQ_DATA_ROOT")
  input_root <- Sys.getenv("DBSEQ_INPUT_PACKAGE_ROOT")
  stopifnot(nzchar(data_root), nzchar(workspace_root))
  panel_dir <- function(kind, n) {
      p <- file.path(package_root, kind, paste("extended data figure", n))
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
      p
  }
  save_data <- function(df, n, panel) {
      write.csv(df, file.path(panel_dir("sourcedata", n), paste0("ed", n, panel, ".csv")), row.names = FALSE)
  }
  save_png <- function(plot, n, panel, w, h) {
      png(file.path(panel_dir("figures", n), paste0("ed", n, panel, ".png")), width = w, height = h, units = "mm", 
          res = 600, type = "cairo")
      if (inherits(plot, "grob") || inherits(plot, "gtable")) 
          grid::grid.draw(plot)
      else print(plot)
      dev.off()
  }
  set.seed(42)
  pdf(file = NULL)
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  OUT_DIR <- panel_dir("figures", 9)
  CACHE_FILE <- file.path(workspace_root, "temp/s7_precomputed_data.qs")
  cat(">>> Loading original cohort objects...\n")
  list2env(supplied_inputs[["helpers_2"]], envir = environment())
  qc_df <- s7_data$qc
  dp_all <- s7_data$dotplot
  target_cts <- s7_data$target_cts
  markers <- s7_data$human_markers
  extra_types <- setdiff(target_cts, names(global_ct_colors))
  if (length(extra_types)) global_ct_colors <- c(global_ct_colors, setNames(hcl.colors(length(extra_types), 
      "Dark 3"), extra_types))
  cohort_names <- names(cohort_paths)
  grp_colors <- c(Control = "#18a799", Epilepsy = "#ac00cc")
  set.seed(42)
  save_cairo_pdf <- function(plot_obj, filename, w_mm, h_mm) {
      keys <- c(FigS7_B_UMAP_Celltype_4Cohorts_172mm.pdf = "b", FigS7_D_UMAP_Sample_4Cohorts_172mm.pdf = "c", 
          FigS7_E_Marker_Bubble_Human_172mm.pdf = "d", FigS7_F_QC_Violin_1x4_172mm.pdf = "e")
      if (basename(filename) %in% names(keys)) 
          save_png(plot_obj, 9, keys[[basename(filename)]], w_mm, h_mm)
  }
  X_LIMITS <- c(-14, 18.5)
  Y_LIMITS <- c(-15, 17.5)
  cat(">>> [Fig S7 V3] Rendering Panel E (Marker Bubble strictly 33mm height)...\n")
  dp_all$id <- factor(dp_all$id, levels = rev(target_cts))
  dp_all$features.plot <- factor(dp_all$features.plot, levels = markers)
  dp_all$Cohort <- factor(dp_all$Cohort, levels = cohort_names)
  p_e_bubble <- ggplot(dp_all, aes(x = features.plot, y = id)) + geom_point(aes(size = pct.exp, fill = avg.exp.scaled), 
      shape = 21, color = "black", stroke = 0.25) + facet_grid(. ~ Cohort) + scale_fill_gradientn(colors = c("#f4f9f4", 
      "#7fcdbb", "#2c7fb8", "#253494"), name = "Scaled Mean\nExpression") + scale_size_continuous(range = c(0.3, 
      2.2), name = "Expressed (%)", breaks = c(25, 50, 75)) + labs(x = "Canonical Marker Genes", y = "Cell Lineage") + 
      theme_classic(base_size = 5, base_family = "sans") + theme(axis.text.x = element_text(angle = 65, 
      hjust = 1, vjust = 1, face = "italic", size = 3.8, color = "black"), axis.text.y = element_text(face = "bold", 
      size = 4.2, color = "black"), axis.title.x = element_text(size = 4.8, face = "bold", margin = margin(t = 1, 
      unit = "pt")), axis.title.y = element_text(size = 4.8, face = "bold", margin = margin(r = 1, unit = "pt")), 
      axis.line = element_line(linewidth = 0.35, color = "black"), axis.ticks = element_line(linewidth = 0.35, 
          color = "black"), axis.ticks.length = unit(0.7, "mm"), strip.background = element_rect(fill = "grey93", 
          color = "black", linewidth = 0.35), strip.text = element_text(size = 4.8, face = "bold", color = "black", 
          margin = margin(1.5, 1.5, 1.5, 1.5, "pt")), panel.grid.major = element_line(color = "grey92", 
          linetype = "dashed", linewidth = 0.2), legend.position = "right", legend.title = element_text(size = 4.2, 
          face = "bold"), legend.text = element_text(size = 3.8), legend.key.size = unit(1.8, "mm"), legend.spacing.y = unit(0.6, 
          "mm"), legend.margin = margin(0, 0, 0, 0, unit = "pt"), plot.margin = margin(1, 2, 1, 2, "pt"))
  panel_e_pdf <- file.path(OUT_DIR, "FigS7_E_Marker_Bubble_Human_172mm.pdf")
  save_cairo_pdf(p_e_bubble, panel_e_pdf, 172, 33)
  cat(">>> [Fig S7 V3] Rendering Panel F (1x4 Pure QC Violins, no boxplot)...\n")
  qc_df$Cohort <- factor(qc_df$Cohort, levels = cohort_names)
  cohort_short <- c(`Cohort 1 (Temporal)` = "Cohort 1 (Temporal)", `Cohort 2 (Frontal)` = "Cohort 2 (Frontal)", 
      `Cohort 3 (Hippocampus)` = "Cohort 3 (Hippocampus)", `Cohort 4 (Amygdala)` = "Cohort 4 (Amygdala)")
  qc_df$Cohort_Short <- factor(cohort_short[as.character(qc_df$Cohort)], levels = unname(cohort_short))
  base_pure_vln <- function(df, y_var, y_label, title_text, log_y = FALSE) {
      p <- ggplot(df, aes(x = Cohort_Short, y = .data[[y_var]], fill = Group)) + geom_violin(position = position_dodge(0.75), 
          trim = TRUE, scale = "width", linewidth = 0.35, alpha = 0.88) + scale_fill_manual(values = grp_colors, 
          name = "Condition:") + labs(title = title_text, x = NULL, y = y_label) + theme_classic(base_size = 5, 
          base_family = "sans") + theme(plot.title = element_text(size = 5, face = "bold", hjust = 0.5, 
          margin = margin(b = 2)), axis.title.y = element_text(size = 4.6, face = "bold", margin = margin(r = 1.5, 
          unit = "pt")), axis.text.y = element_text(size = 4.2, color = "black"), axis.text.x = element_text(size = 3.8, 
          face = "bold", color = "black", angle = 28, hjust = 1, vjust = 1), axis.ticks.x = element_line(linewidth = 0.35, 
          color = "black"), axis.line = element_line(linewidth = 0.35, color = "black"), axis.ticks.y = element_line(linewidth = 0.35, 
          color = "black"), axis.ticks.length = unit(0.8, "mm"), legend.position = "none", plot.margin = margin(1, 
          2, 1, 2, "pt"))
      if (log_y) {
          p <- p + scale_y_log10(labels = label_log())
      }
      else {
          p <- p + scale_y_continuous(labels = label_comma())
      }
      p
  }
  pv1 <- base_pure_vln(qc_df, "nFeature_RNA", "Detected Genes", "Detected Genes (nFeature)", log_y = FALSE)
  pv2 <- base_pure_vln(qc_df, "nCount_RNA", "nCount (log10)", "Sequencing Depth (nCount)", log_y = TRUE)
  pv3 <- base_pure_vln(qc_df, "percent.mt", "Mito Reads (%)", "Mitochondrial Ratio (%mito)", log_y = FALSE)
  pv4 <- base_pure_vln(qc_df, "percent.rb", "Ribo Reads (%)", "Ribosomal Ratio (%ribo)", log_y = FALSE)
  p_qc_for_leg <- pv1 + theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 4.8, 
      face = "bold"), legend.text = element_text(size = 4.5), legend.key.size = unit(2.2, "mm"), legend.spacing.x = unit(2.5, 
      "mm"), legend.margin = margin(t = 0, b = 1, unit = "pt")) + guides(fill = guide_legend(override.aes = list(linewidth = 0.35, 
      alpha = 0.9)))
  p_qc_leg <- cowplot::get_legend(p_qc_for_leg)
  p_f_grid <- wrap_elements(full = p_qc_leg)/(pv1 | pv2 | pv3 | pv4) + plot_layout(heights = c(0.08, 1))
  panel_f_pdf <- file.path(OUT_DIR, "FigS7_F_QC_Violin_1x4_172mm.pdf")
  save_cairo_pdf(p_f_grid, panel_f_pdf, 172, 36)
  cat(">>> [Fig S7 V3] All R panels generated successfully!\n")
  save_data(dp_all, 9, "d")
  save_data(qc_df[, c("Cell_ID", "Cohort", "Sample", "Group", "celltype", "nFeature_RNA", "nCount_RNA", 
      "percent.mt", "percent.rb")], 9, "e")
  invisible(as.list(environment()))
}
