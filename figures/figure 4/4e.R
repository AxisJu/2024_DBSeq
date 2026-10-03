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
  grn_counts <- data_master$grn_counts
  df_rect <- grn_counts %>% select(Category, Count) %>% mutate(Category = factor(Category, levels = c("Other Genes", 
      "Neuronal-Common", "ATN-specific"))) %>% arrange(Category) %>% mutate(Prop = Count/sum(Count), ymax = cumsum(Prop), 
      ymin = ymax - Prop, ymid = (ymin + ymax)/2, Label_Pct = sprintf("%.1f%%", 100 * Prop), Label_y = ifelse(Category == 
          "ATN-specific", ymax + 0.025, ymid))
  p_bar <- ggplot(df_rect) + geom_rect(aes(xmin = 0.82, xmax = 1.18, ymin = ymin, ymax = ymax, fill = Category), 
      color = "black", linewidth = 0.65) + scale_fill_manual(values = c(`Other Genes` = "#D9D9D9", `Neuronal-Common` = "#4DBBD5", 
      `ATN-specific` = "#E64B35")) + geom_text(aes(x = 1, y = Label_y, label = Label_Pct), fontface = "bold", 
      size = 4.2) + annotate("text", x = 1, y = 1.075, label = "ATN-specific", color = "#E64B35", fontface = "bold", 
      size = 4.3) + geom_text(data = subset(df_rect, Category != "ATN-specific"), aes(x = 1.28, y = ymid, 
      label = Category), color = "grey40", fontface = "bold", size = 4, angle = -90) + scale_x_continuous(limits = c(0.7, 
      1.45), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.12), expand = c(0, 0)) + theme_void() + 
      theme(legend.position = "none", plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
  out_fig_dir <- output_path("figures/figure 4")
  out_data_dir <- output_path("sourcedata/figure 4")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_fig_dir, "4e.png"), plot = p_bar, width = 2.4, height = 7, dpi = 600)
  df_source_e <- df_rect %>% select(Category, Count, Prop) %>% rename(Gene_Category = Category, Gene_Count = Count, 
      Proportion = Prop)
  write.csv(df_source_e, file.path(out_data_dir, "4e.csv"), row.names = FALSE)
  cat("Figure 4e generated successfully.\n")
  invisible(as.list(environment()))
}
