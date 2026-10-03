# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  proj_dir <- data_root
  out_fig_dir <- output_path("figures/figure 2")
  out_data_dir <- output_path("sourcedata/figure 2")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  plot_df <- supplied_inputs[["plot_df_2"]]
  mouse_non_neuron <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  mouse_neuron <- setdiff(unique(plot_df$Celltype), mouse_non_neuron)
  n_total_neuron <- length(mouse_neuron)
  n_total_nonneuron <- length(mouse_non_neuron)
  gene_stats <- plot_df %>% filter(is_TRG == TRUE) %>% group_by(genes) %>% summarise(n_celltypes = n_distinct(Celltype), 
      trg_in_neuron = sum(Celltype %in% mouse_neuron), trg_in_nonneuron = sum(Celltype %in% mouse_non_neuron), 
      .groups = "drop") %>% mutate(pct_neuron = round((trg_in_neuron/n_total_neuron) * 100, 1), pct_nonneuron = round((trg_in_nonneuron/n_total_nonneuron) * 
      100, 1))
  source_g <- gene_stats %>% select(Gene = genes, Recurrence_Count = n_celltypes, Neuronal_Celltypes = trg_in_neuron, 
      NonNeuronal_Celltypes = trg_in_nonneuron, Pct_Neuronal = pct_neuron, Pct_NonNeuronal = pct_nonneuron)
  write.csv(source_g, file.path(out_data_dir, "2g.csv"), row.names = FALSE)
  cat("Saved source data: sourcedata/figure 2/2g.csv\n")
  trg_counts <- gene_stats %>% group_by(n_celltypes) %>% mutate(count = n()) %>% ungroup()
  common_genes <- c("Rnpc3", "Crebzf", "Ndfip1", "Ckb", "Pura")
  neuron_genes <- c("Ube3a", "Bsg", "Matk", "Timm9", "Myl12b")
  nonneuron_genes <- c("Mertk", "Ccnd3", "Rab31", "Csrp1", "Sox10")
  build_atop <- function(strings) {
      if (length(strings) == 1) 
          return(strings[1])
      if (length(strings) == 2) 
          return(paste0("atop(", strings[1], ", ", strings[2], ")"))
      return(paste0("atop(", strings[1], ", ", build_atop(strings[-1]), ")"))
  }
  label_df <- trg_counts %>% filter(genes %in% c(common_genes, neuron_genes, nonneuron_genes)) %>% mutate(Category = case_when(genes %in% 
      common_genes ~ "Common", genes %in% neuron_genes ~ "Neuron", genes %in% nonneuron_genes ~ "NonNeuron"), 
      Category = factor(Category, levels = c("Common", "Neuron", "NonNeuron")), math_str = sprintf("paste(italic('%s'), ': %.1f%% Neuronal | %.1f%% Non-Neuronal')", 
          genes, pct_neuron, pct_nonneuron)) %>% group_by(n_celltypes, count, Category) %>% summarise(final_label = build_atop(math_str), 
      .groups = "drop")
  threshold_common <- 3
  p <- ggplot(trg_counts, aes(x = n_celltypes)) + geom_bar(aes(fill = factor(n_celltypes > threshold_common)), 
      color = "white", width = 0.95) + scale_fill_manual(values = c(`FALSE` = "#d4d4d4", `TRUE` = "#5e5e5e"), 
      guide = "none") + geom_vline(xintercept = threshold_common + 0.5, linetype = "dashed", color = "black", 
      linewidth = 0.5) + geom_text_repel(data = label_df, aes(x = n_celltypes, y = count, label = final_label, 
      color = Category), parse = TRUE, nudge_y = 60, direction = "both", segment.color = "grey40", segment.size = 0.5, 
      min.segment.length = 0, box.padding = 0.8, max.overlaps = Inf, arrow = arrow(length = unit(0.001, 
          "npc"), type = "closed")) + scale_color_manual(values = c(Common = "#ffd700", Neuron = "#18a799", 
      NonNeuron = "#000000"), name = "Gene Module") + labs(x = "Number of subclasses or supertypes", y = "Count") + 
      theme_classic(base_size = 14) + theme(axis.line = element_line(color = "black", linewidth = 0.5), 
      axis.ticks = element_line(color = "black", linewidth = 0.5), axis.text = element_text(color = "black"), 
      panel.grid = element_blank(), legend.position = c(0.85, 0.85), legend.background = element_rect(fill = "transparent"))
  pdf_path <- file.path(out_fig_dir, "2g.pdf")
  png_path <- file.path(out_fig_dir, "2g.png")
  ggsave(png_path, p, width = 7.3, height = 9, dpi = 300)
  cat("Saved figures: figures/figure 2/2g.pdf and 2g.png\n")
  invisible(as.list(environment()))
}
