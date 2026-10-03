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
  bg_rate <- data_master$bg_rate
  p_val <- data_master$p_value
  stopifnot(isTRUE(all.equal(df_pathway_stats$Prop_ATN, df_pathway_stats$Count_ATN/df_pathway_stats$Total_Genes)))
  df_box <- df_pathway_stats %>% mutate(Category = "ATN-Specific")
  p_box <- ggplot(df_box, aes(x = Category, y = Prop_ATN)) + geom_boxplot(fill = "#E64B35", color = "black", 
      width = 0.45, linewidth = 0.7, alpha = 0.85, outlier.shape = 16, outlier.size = 2, outlier.color = "grey60") + 
      geom_hline(yintercept = bg_rate, linetype = "dashed", color = "grey40", linewidth = 0.6) + annotate("text", 
      x = 0.32, y = bg_rate + 0.01, label = sprintf("Conditional null\n%.3f%%", 100 * bg_rate), fontface = "plain", 
      size = 2.9, hjust = 1, lineheight = 0.95, color = "black") + annotate("segment", x = 0.36, xend = 0.48, 
      y = bg_rate + 0.009, yend = bg_rate, arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", 
      linewidth = 0.55) + annotate("segment", x = 0.78, xend = 1.22, y = 0.22, yend = 0.22, linewidth = 0.55, 
      color = "black") + annotate("text", x = 1, y = 0.235, label = paste0("Permutation P = ", format.pval(p_val, 
      digits = 3, eps = .Machine$double.xmin)), parse = FALSE, size = 3.6, color = "black") + scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
      breaks = c(0, 0.1, 0.2), expand = c(0, 0)) + scale_x_discrete(expand = expansion(mult = c(0.7, 0.5))) + 
      coord_cartesian(ylim = c(-0.01, max(0.26, max(df_box$Prop_ATN, na.rm = TRUE) * 1.1)), clip = "off") + 
      theme_bw(base_size = 11) + theme(axis.title = element_blank(), axis.text.x = element_text(face = "bold", 
      color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10), axis.ticks = element_line(color = "black", 
      linewidth = 0.5), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), panel.grid.major.x = element_blank(), 
      panel.grid.minor = element_blank(), panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4), 
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5, lineheight = 1.15, margin = margin(b = 10)), 
      plot.margin = margin(t = 12, r = 10, b = 10, l = 15)) + labs(title = "Proportion of ATN-specific\nTRGs in involved pathways")
  out_fig_dir <- output_path("figures/figure 4")
  out_data_dir <- output_path("sourcedata/figure 4")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_fig_dir, "4f.png"), plot = p_box, width = 3, height = 4.2, dpi = 600)
  df_source_f <- df_pathway_stats %>% select(Term, Clean_Term, Count_ATN, Total_Genes, Prop_ATN) %>% rename(Pathway_Term = Term, 
      Pathway_Name = Clean_Term, ATN_Specific_Gene_Count = Count_ATN, Total_Pathway_Genes = Total_Genes, 
      Proportion_ATN_Specific = Prop_ATN)
  write.csv(df_source_f, file.path(out_data_dir, "4f.csv"), row.names = FALSE)
  cat("Figure 4f generated successfully.\n")
  invisible(as.list(environment()))
}
