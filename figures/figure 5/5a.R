# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(qs)
  out_csv_path <- output_path("sourcedata/figure 5/5a.csv")
  out_pdf_path <- output_path("figures/figure 5/5a.pdf")
  out_png_path <- output_path("figures/figure 5/5a.png")
  dir.create(dirname(out_csv_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_pdf_path), recursive = TRUE, showWarnings = FALSE)
  detailed_path <- data_path("results/Table/Pathway_Circuit_Integration/7_Gene_CohenD_Splicing_Detailed.csv")
  summary_path <- data_path("results/Table/Pathway_Circuit_Integration/splicing-gene-stat.csv")
  plot_df_path <- data_path("results/Table/CosSim_v2/plot_df.qs")
  detailed <- supplied_inputs[["detailed_2"]][[as.character(detailed_path)]]
  summary_data <- supplied_inputs[["summary_data_3"]][[as.character(summary_path)]]
  plot_df <- supplied_inputs[["plot_df_4"]][[as.character(plot_df_path)]]
  meta_path <- data_path("data/snRNAseq_mouse/processed/metadata/metadata_obs_109.csv")
  meta <- supplied_inputs[["meta_5"]][[as.character(meta_path)]]
  region_to_celltype <- meta %>% mutate(Region = case_when(celltype_level3 == "ExN_ATN" ~ "ATN", !is.na(brainregion_projection) & 
      brainregion_projection != "/" ~ brainregion_projection, TRUE ~ "Unknown")) %>% filter(Region != "Unknown", 
      Sample %in% c("DBS_I_1", "DBS_I_2", "Sham_I_1", "Sham_I_2", "Sham_I_3")) %>% group_by(Region) %>% 
      summarise(Main_Celltype = names(which.max(table(celltype_level3))), .groups = "drop")
  write_csv(region_to_celltype, sub("5a.csv$", "5a_region_celltype.csv", out_csv_path))
  splicing_genes <- summary_data$Gene
  trg_per_region <- plot_df %>% filter(genes %in% splicing_genes, is_TRG) %>% inner_join(region_to_celltype, 
      by = c(Celltype = "Main_Celltype")) %>% distinct(Gene = genes, Region)
  df_merged <- detailed %>% inner_join(summary_data %>% select(Gene, Rho_SC2_All, P_SC2_All), by = "Gene")
  regions <- c("ATN", "CP", "RSPd", "RT", "RE", "GPe", "LH", "PF", "CM", "DG", "LSc", "MH", "TRS")
  gene_order <- summary_data %>% arrange(desc(Rho_SC2_All)) %>% pull(Gene)
  df_plot <- df_merged %>% filter(Region %in% regions) %>% mutate(Region = factor(Region, levels = regions), 
      Gene = factor(Gene, levels = rev(gene_order)))
  df_border <- trg_per_region %>% filter(Region %in% regions) %>% mutate(Region = factor(Region, levels = regions), 
      Gene = factor(Gene, levels = rev(gene_order))) %>% filter(!is.na(Gene))
  out_csv <- df_plot %>% left_join(region_to_celltype, by = "Region") %>% left_join(trg_per_region %>% 
      mutate(is_TRG = TRUE), by = c("Gene", "Region")) %>% mutate(is_TRG = coalesce(is_TRG, FALSE), SC2_significant = P_SC2_All < 
      0.05) %>% arrange(desc(Gene), Region)
  write_csv(out_csv, out_csv_path)
  max_val <- max(abs(df_merged$CohenD_All_v_Sham), na.rm = TRUE)
  p_main <- ggplot(df_plot, aes(x = Region, y = Gene)) + geom_tile(aes(fill = CohenD_All_v_Sham), color = "white", 
      linewidth = 0.3) + geom_tile(data = df_border, color = "black", fill = NA, linewidth = 0.6) + scale_fill_gradientn(name = expression(paste("Cohen's ", 
      italic(d), " (DBS vs Sham)")), colors = rev(brewer.pal(11, "RdBu"))[2:10], values = rescale(c(-max_val, 
      0, max_val)), limits = c(-max_val, max_val), na.value = "grey90") + theme_minimal(base_size = 12) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"), axis.text.y = element_blank(), 
          axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.title = element_text(size = 10), 
          plot.margin = margin(5, 10, 5, 2))
  df_bar <- summary_data %>% mutate(Gene = factor(Gene, levels = rev(gene_order)), x_col = "Rho")
  sig_genes <- summary_data$Gene[!is.na(summary_data$P_SC2_All) & summary_data$P_SC2_All < 0.05]
  p_bar <- ggplot(df_bar, aes(x = x_col, y = Gene)) + geom_tile(aes(fill = Rho_SC2_All), color = "white", 
      linewidth = 0.3) + geom_text(data = df_bar %>% filter(Gene %in% sig_genes), aes(label = "*"), color = "white", 
      size = 6, vjust = 0.75) + scale_fill_gradientn(name = expression(paste("Spearman's ", rho, " with 2"^"nd", 
      " SC")), colors = c("#f4f9f4", "#7fcdbb", "#2c7fb8", "#253494"), limits = range(c(0, summary_data$Rho_SC2_All), 
      na.rm = TRUE), oob = squish, na.value = "grey90") + theme_minimal(base_size = 12) + theme(axis.text.x = element_blank(), 
      axis.text.y = element_text(face = "italic", size = 9, color = "black"), axis.title = element_blank(), 
      axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.title = element_text(size = 10), 
      plot.margin = margin(5, 0, 5, 5))
  p_combined <- (p_bar + p_main + plot_layout(widths = c(1, 13), guides = "collect")) & theme(legend.position = "bottom", 
      legend.title.position = "top") & guides(fill = guide_colourbar(title.position = "top", barwidth = unit(35, 
      "mm")))
  ggsave(out_png_path, p_combined, width = 6, height = 10, units = "in", dpi = 600)
  cat(">>> Figure 5a complete.\n")
  invisible(as.list(environment()))
}
