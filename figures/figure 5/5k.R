# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(qs)
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(STRINGdb)
  dir.create(output_path("figures/figure 5"), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_path("sourcedata/figure 5"), recursive = TRUE, showWarnings = FALSE)
  plot_df <- supplied_inputs[["plot_df_2"]][[as.character(data_path("results/Table/CosSim_v2/plot_df.qs"))]]
  mouse_non_neuron <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  mouse_neuron <- setdiff(unique(plot_df$Celltype), mouse_non_neuron)
  trg_class <- plot_df %>% filter(is_TRG == TRUE, Celltype %in% mouse_neuron) %>% group_by(genes) %>% summarise(n_neurons = n_distinct(Celltype), 
      is_in_ATN = "ExN_ATN" %in% Celltype, .groups = "drop") %>% filter(is_in_ATN == TRUE) %>% mutate(TRG_Group = ifelse(n_neurons == 
      1, "ATN-Specific", "Neuron-Common"))
  list2env(supplied_inputs[["helpers_3"]], envir = environment())
  publication_universe <- get_publication_universe()
  trg_class <- trg_class %>% filter(toupper(genes) %in% publication_universe)
  union_genes <- unique(trg_class$genes)
  write.csv(data.frame(Gene = sort(publication_universe)), output_path("sourcedata/figure 5/5kl_go_universe.csv"), 
      row.names = FALSE)
  str_dir <- data_path("results/260428 Cytoscape/STRING_data")
  string_db <- STRINGdb$new(version = "12.0", species = 10090, score_threshold = 900, network_type = "full", 
      input_directory = str_dir)
  df_genes <- data.frame(Gene = as.character(union_genes), stringsAsFactors = FALSE)
  df_mapped <- string_db$map(df_genes, "Gene", removeUnmappedRows = TRUE)
  ppi_graph <- string_db$get_graph()
  global_degrees <- igraph::degree(ppi_graph)
  df_mapped$Degree <- as.numeric(global_degrees[as.character(df_mapped$STRING_id)])
  df_degree_clean <- df_mapped %>% mutate(JoinKey = toupper(trimws(as.character(Gene)))) %>% group_by(JoinKey) %>% 
      summarise(Degree = if (all(is.na(Degree))) 0 else max(Degree, na.rm = TRUE), .groups = "drop") %>% 
      mutate(Degree = ifelse(is.infinite(Degree) | is.na(Degree), 0, as.numeric(Degree)))
  df_ppi_stats <- trg_class %>% mutate(JoinKey = toupper(trimws(as.character(genes)))) %>% left_join(df_degree_clean, 
      by = "JoinKey") %>% mutate(Degree = ifelse(is.na(Degree), 0, Degree)) %>% filter(Degree > 0) %>% 
      mutate(Log10_Degree_Plus1 = log10(Degree + 1))
  csv_out <- df_ppi_stats %>% select(Gene = genes, Category = TRG_Group, n_neurons, Degree, Log10_Degree_Plus1)
  write.csv(csv_out, output_path("sourcedata/figure 5/5k.csv"), row.names = FALSE)
  p_value <- wilcox.test(Log10_Degree_Plus1 ~ TRG_Group, df_ppi_stats, exact = FALSE)$p.value
  write.csv(data.frame(n_genes = nrow(df_ppi_stats), test = "Wilcoxon rank-sum", p_value = p_value), output_path("statistics/figure 5/5k_statistics.csv"), 
      row.names = FALSE)
  df_ppi_stats$Category_plot <- factor(df_ppi_stats$TRG_Group, levels = c("ATN-Specific", "Neuron-Common"), 
      labels = c("ATN\nSpecific", "Neuronal\nCommon"))
  p_box <- ggplot(df_ppi_stats, aes(x = Category_plot, y = Log10_Degree_Plus1, fill = Category_plot)) + 
      geom_boxplot(color = "black", alpha = 0.75, outlier.shape = 21, outlier.alpha = 0, width = 0.5) + 
      scale_fill_manual(values = c(`ATN\nSpecific` = "#E64B35", `Neuronal\nCommon` = "#4DBBD5")) + stat_compare_means(comparisons = list(c("ATN\nSpecific", 
      "Neuronal\nCommon")), method = "wilcox.test", label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.1, 
      0.2))) + theme_bw(base_size = 14) + theme(legend.position = "none", plot.title = element_text(size = 11), 
      axis.title.x = element_blank(), axis.text = element_text(face = "bold", color = "black"), panel.grid.minor = element_blank()) + 
      labs(y = "Log10(PPI Degree + 1)", title = "Centrality in\nHighest Conf. PPI Network")
  ggsave(output_path("figures/figure 5/5k.png"), p_box, width = 3.5, height = 4.5, dpi = 600)
  cat(">>> Figure 5k complete.\n")
  all_neuronal_trg <- plot_df %>% filter(is_TRG == TRUE, Celltype %in% mouse_neuron) %>% group_by(genes) %>% 
      summarise(n_neurons = n_distinct(Celltype), is_in_ATN = "ExN_ATN" %in% Celltype, .groups = "drop")
  saveRDS(list(schema_version = 2L, publication_universe = publication_universe, trg_class = trg_class, 
      all_neuronal_trg = all_neuronal_trg, df_degree_clean = df_degree_clean), cache_path("5kl_shared_data.rds"))
  invisible(as.list(environment()))
}
