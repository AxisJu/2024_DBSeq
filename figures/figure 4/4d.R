# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(scales)
  cache_file <- cache_path("fig4_def_master.qs")
  if (!supplied_inputs[["data_2"]] || Sys.getenv("DBSEQ_REBUILD_SOURCE") == "1") {
      list2env(supplied_inputs[["helpers_3"]], envir = environment())
      prepare_4_ora()
  }
  data_master <- supplied_inputs[["data_master_4"]][[as.character(cache_file)]]
  df_pathway_stats <- data_master$df_pathway_stats
  df_pie <- df_pathway_stats %>% dplyr::group_by(Category) %>% dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>% 
      dplyr::mutate(Prop = Count/sum(Count), Category = factor(Category, levels = c("Only ATN-specific", 
          "Only Neuron-Common", "Both Shared"))) %>% arrange(Category)
  df_pie$y_pos <- cumsum(df_pie$Prop) - df_pie$Prop/2
  category_labels <- c(`Only ATN-specific` = "ATN-specific", `Only Neuron-Common` = "Neuronal\nCommon", 
      `Both Shared` = "Shared")
  df_pie$Label <- paste0(category_labels[as.character(df_pie$Category)], "\n", sprintf("%.1f%%", 100 * 
      df_pie$Prop))
  p_pie <- ggplot(df_pie, aes(x = "", y = Prop, fill = Category)) + geom_bar(stat = "identity", width = 1, 
      color = "white", linewidth = 1.2, position = position_stack(reverse = TRUE)) + coord_polar("y", start = 0) + 
      scale_fill_manual(values = c(`Both Shared` = "#1B9E77", `Only ATN-specific` = "#E64B35", `Only Neuron-Common` = "#4DBBD5")) + 
      geom_text(aes(y = y_pos, label = Label), color = "white", fontface = "bold", size = 3.8) + theme_void(base_size = 12) + 
      theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11, hjust = 0.5, 
          lineheight = 1.2, margin = margin(b = 10)), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + 
      labs(title = "ATN-enriched pathways\n(BH FDR < 0.05; TRG hits > 5)")
  out_fig_dir <- output_path("figures/figure 4")
  out_data_dir <- output_path("sourcedata/figure 4")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_fig_dir, "4d.png"), plot = p_pie, width = 4, height = 4, dpi = 600)
  df_source_d <- df_pie %>% select(Category, Count, Prop) %>% rename(Pathway_Category = Category, Pathway_Count = Count, 
      Proportion = Prop)
  write.csv(df_source_d, file.path(out_data_dir, "4d.csv"), row.names = FALSE)
  cat("Figure 4d generated successfully.\n")
  invisible(as.list(environment()))
}
