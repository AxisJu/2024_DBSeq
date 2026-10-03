# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(Seurat)
      library(dplyr)
      library(tidyr)
      library(ggplot2)
      library(ggdist)
      library(ggrastr)
  })
  base_dir <- data_root
  output_dir_fig <- output_path("figures/figure 3")
  output_dir_data <- output_path("sourcedata/figure 3")
  dir.create(output_dir_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  cohort1_path <- file.path(base_dir, "data/snRNAseq_human_DRE/cohort1_MTLE_brainbenigh_temporalcortex.qs")
  cohort2_path <- file.path(base_dir, "data/snRNAseq_human_DRE/cohort2_FCD2b_selfcontrol_frontalcortex.qs")
  cohort3_path <- file.path(base_dir, "data/snRNAseq_human_DRE/cohort3_MTLE_autopsy_hippocampus.qs")
  cohort4_path <- file.path(base_dir, "data/snRNAseq_human_DRE/cohort4_MTLE_autopsy_amygdala.qs")
  cohort_labels <- c("Cohort 1 (Temporal)", "Cohort 2 (Frontal)", "Cohort 3 (Hippocampus)", "Cohort 4 (Amygdala)")
  set.seed(42)
  if ("--from-source" %in% commandArgs(trailingOnly = TRUE)) {
      obs <- supplied_inputs[["obs_2"]][[as.character(input_path("sourcedata/figure 3/3e.csv"))]]
      df_long <- obs %>% transmute(Cell_ID, Cohort_ID = Cohort, Cohort = factor(cohort_labels[match(Cohort, 
          paste("Cohort", 1:4))], levels = cohort_labels), Region, Pathology, Group = factor(Condition, 
          levels = c("Control", "Epilepsy")), Gene = factor(Gene, levels = c("NLGN1", "NRXN1")), Expression = Expression_logCPM)
  } else {
      extract_gene_data <- function(obj, cohort_id, cohort_label, region, pathology) {
          df <- FetchData(obj, vars = c("NLGN1", "NRXN1", "Group"))
          colnames(df)[3] <- "Group"
          df$Cell_ID <- rownames(df)
          df$Cohort_ID <- cohort_id
          df$Cohort <- cohort_label
          df$Region <- region
          df$Pathology <- pathology
          return(df)
      }
      df1 <- extract_gene_data(supplied_inputs[["df1_3"]][[as.character(cohort1_path)]] %>% subset(celltype_coarse == 
          "Glu.N"), "Cohort 1", "Cohort 1 (Temporal)", "Temporal Ctx.", "MTLE")
      gc()
      gc()
      df2 <- extract_gene_data(supplied_inputs[["df2_4"]][[as.character(cohort2_path)]] %>% subset(celltype_coarse == 
          "Glu.N"), "Cohort 2", "Cohort 2 (Frontal)", "Frontal Ctx.", "FCD2b")
      gc()
      df3 <- extract_gene_data(supplied_inputs[["df3_5"]][[as.character(cohort3_path)]] %>% subset(celltype_coarse == 
          "Glu.N"), "Cohort 3", "Cohort 3 (Hippocampus)", "Hippocampus", "MTLE")
      df4 <- extract_gene_data(supplied_inputs[["df4_6"]][[as.character(cohort4_path)]] %>% subset(celltype_coarse == 
          "Glu.N"), "Cohort 4", "Cohort 4 (Amygdala)", "Amygdala", "MTLE, FCD2a")
      gc()
      df_combined <- bind_rows(df1, df2, df3, df4)
      df_long <- df_combined %>% pivot_longer(cols = c("NLGN1", "NRXN1"), names_to = "Gene", values_to = "Expression") %>% 
          mutate(Group = factor(Group, levels = c("Control", "Epilepsy")), Cohort = factor(Cohort, levels = c("Cohort 1 (Temporal)", 
              "Cohort 2 (Frontal)", "Cohort 3 (Hippocampus)", "Cohort 4 (Amygdala)")), Gene = factor(Gene, 
              levels = c("NLGN1", "NRXN1")))
  }
  sourcedata_e <- df_long %>% dplyr::select(Cell_ID, Cohort = Cohort_ID, Region, Pathology, Condition = Group, 
      Gene, Expression_logCPM = Expression)
  write.csv(sourcedata_e, file.path(output_dir_data, "3e.csv"), row.names = FALSE)
  cat(sprintf("Exported Source Data: %s (Rows: %d)\n", file.path(output_dir_data, "3e.csv"), nrow(sourcedata_e)))
  stat_annotations <- bind_rows(lapply(1:4, function(i) {
      sc_path <- file.path(base_dir, sprintf("results/Table/human_cohort_de/sc/SC_DE_cohort%d_Glu.N.csv", 
          i))
      pb_path <- file.path(base_dir, sprintf("results/Table/human_cohort_de/pb/PB_DE_cohort%d_Glu.N.csv", 
          i))
      sc <- supplied_inputs[["sc_7"]][[as.character(sc_path)]] %>% filter(Gene %in% c("NLGN1", "NRXN1")) %>% 
          transmute(Gene, SC_P = p_val, SC_Adjusted_P = p_val_adj, SC_CohenD_Control_minus_Epilepsy = CohenD)
      pb <- supplied_inputs[["pb_8"]][[as.character(pb_path)]] %>% filter(Gene %in% c("NLGN1", "NRXN1")) %>% 
          transmute(Gene, PB_P = PValue, PB_FDR = p_val_adj, PB_CohenD_Control_minus_Epilepsy = CohenD)
      stopifnot(!anyDuplicated(sc$Gene), !anyDuplicated(pb$Gene))
      full_join(sc, pb, by = "Gene") %>% mutate(Cohort_ID = paste("Cohort", i), Cohort = factor(cohort_labels[i], 
          levels = cohort_labels), Celltype = "Glu.N", CohenD_Epilepsy_minus_Control = -SC_CohenD_Control_minus_Epilepsy, 
          SC_Method = "MAST; max.cells.per.ident=2000; logfc.threshold=0.1; min.pct=0.1; available Age/Gender covariates", 
          SC_Unit = "nucleus", SC_Adjustment = "Bonferroni (Seurat FindMarkers)", PB_Method = "edgeR quasi-likelihood F-test", 
          PB_Unit = "sample pseudobulk", PB_Adjustment = "BH", SC_Source = sc_path, PB_Source = pb_path, 
          Analysis_Code = "src/200_TRG-identification/206_analysis_TRG-common_260413.R", SC_Source_MD5 = unname(tools::md5sum(sc_path)), 
          PB_Source_MD5 = unname(tools::md5sum(pb_path)))
  }))
  stopifnot(nrow(stat_annotations) == 8L)
  stat_annotations <- stat_annotations %>% mutate(Gene = factor(Gene, levels = c("NLGN1", "NRXN1")), Label = paste0(ifelse(Cohort_ID == 
      "Cohort 2", "One patient; descriptive", "Cell distributions; descriptive"), "\nCohen's d = ", sprintf("%.2f", 
      CohenD_Epilepsy_minus_Control)), y_pos = ifelse(Gene == "NLGN1", 6.8, 8.8))
  write.csv(stat_annotations, output_path("statistics/figure 3/3e_statistics.csv"), row.names = FALSE)
  colors_group <- c(Control = "#18a799", Epilepsy = "#ac00cc")
  p_e <- ggplot(df_long, aes(x = Cohort, y = Expression, fill = Group, color = Group)) + geom_point_rast(position = position_jitterdodge(jitter.width = 0.15, 
      dodge.width = 0.8), size = 0.1, alpha = 0.3, stroke = 0, raster.dpi = 300) + geom_boxplot(width = 0.12, 
      position = position_dodge(width = 0.8), fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.45, 
      show.legend = FALSE) + stat_halfeye(adjust = 0.5, justification = -0.25, .width = 0, point_colour = NA, 
      position = position_dodge(width = 0.8), alpha = 0.5, linewidth = 0.45, normalize = "groups", scale = 0.5) + 
      geom_text(data = stat_annotations, aes(x = Cohort, y = y_pos, label = Label), inherit.aes = FALSE, 
          size = 2.7, lineheight = 0.9, color = "black") + facet_wrap(~Gene, ncol = 1, scales = "free_y") + 
      scale_fill_manual(values = colors_group) + scale_color_manual(values = colors_group) + scale_x_discrete(labels = function(x) sub(" (", 
      "\n(", x, fixed = TRUE)) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + theme_bw(base_size = 11) + 
      labs(x = NULL, y = "Log-normalized expression (Seurat RNA data)", fill = "Condition", color = "Condition") + 
      theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face = "bold", 
          color = "black", size = 10), axis.text.y = element_text(color = "black", size = 9), axis.title.y = element_text(face = "bold", 
          size = 10), strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "bold.italic", size = 11, color = "black"), legend.position = "top", 
          legend.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 9))
  fig_pdf <- file.path(output_dir_fig, "3e.pdf")
  fig_png <- file.path(output_dir_fig, "3e.png")
  ggsave(fig_png, p_e, width = 6.5, height = 7.5, dpi = 300)
  invisible(as.list(environment()))
}
