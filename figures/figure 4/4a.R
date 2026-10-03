# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtern)
  library(RColorBrewer)
  library(scales)
  library(cowplot)
  DATA_COMBINED <- cache_path("df_combined_tern.qs")
  MOUSE_DF_PATH <- data_path("results/Table/CosSim_v2/plot_df.qs")
  OUT_PNG <- output_path("figures/figure 4/4a.png")
  OUT_CSV <- output_path("sourcedata/figure 4/4a.csv")
  list2env(supplied_inputs[["helpers_2"]], envir = environment())
  if (!supplied_inputs[["data_3"]] || Sys.getenv("DBSEQ_REBUILD_SOURCE") == "1") prepare_4a()
  df_combined <- supplied_inputs[["df_combined_4"]][[as.character(DATA_COMBINED)]]
  mouse_df <- supplied_inputs[["mouse_df_5"]]
  mouse_non_neuron <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  mouse_neuron <- setdiff(unique(mouse_df$Celltype), mouse_non_neuron)
  df_heatmap <- mouse_df %>% filter(is_TRG == TRUE & Celltype %in% mouse_neuron) %>% group_by(genes) %>% 
      summarise(n_neurons = n_distinct(Celltype), is_in_ATN = "ExN_ATN" %in% Celltype, ATN_CoD = -mean(CoD_DBSvSham_cohend[Celltype == 
          "ExN_ATN"], na.rm = TRUE), .groups = "drop") %>% filter(is_in_ATN == TRUE) %>% mutate(Direction = ifelse(ATN_CoD > 
      0, "Up", "Down"), Specificity = case_when(n_neurons == 1 ~ "1", n_neurons >= 2 & n_neurons <= 4 ~ 
      "2-4", n_neurons >= 5 & n_neurons <= 7 ~ "5-7", n_neurons >= 8 & n_neurons <= 10 ~ "8-10", n_neurons >= 
      11 & n_neurons <= 13 ~ "11-13", n_neurons >= 14 & n_neurons <= 16 ~ "14-16", n_neurons >= 17 & n_neurons <= 
      19 ~ "17-19", n_neurons >= 20 ~ "20+"))
  df_source <- df_combined %>% left_join(df_heatmap %>% select(genes, n_neurons, Specificity, Direction), 
      by = c(gene = "genes")) %>% select(gene, DBS_I, Sham_I, Saline_I, is_TRG, delta_CohenD, median_cosine, 
      mean_cosine, n_neurons, Specificity_Bin = Specificity, Direction)
  write.csv(df_source, OUT_CSV, row.names = FALSE)
  df_tern_sig <- df_combined %>% filter(is_TRG == TRUE)
  df_tern_nsig <- df_combined %>% filter(is_TRG == FALSE)
  p_tern <- ggtern(mapping = aes(x = DBS_I, y = Sham_I, z = Saline_I)) + geom_point(data = df_tern_nsig, 
      alpha = 0.2, size = 0.7, color = "#d0d0d0") + geom_point(data = df_tern_sig %>% arrange(median_cosine), 
      shape = 16, aes(color = delta_CohenD, alpha = median_cosine, size = median_cosine)) + geom_Tline(Tintercept = 1/3, 
      linetype = "twodash", linewidth = 0.6, color = "#ac00cc") + geom_Lline(Lintercept = 1/3, linetype = "twodash", 
      linewidth = 0.6, color = "#ffd700") + geom_Rline(Rintercept = 1/3, linetype = "twodash", linewidth = 0.6, 
      color = "#18a799") + theme_bw() + theme(axis.line = element_line(linewidth = 0.8), axis.ticks = element_blank(), 
      panel.grid = element_blank(), tern.panel.background = element_rect(fill = "white", color = NA), legend.position = c(0.1, 
          0.15), legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_text(size = 9, 
          face = "bold"), legend.text = element_text(size = 8), axis.title = element_blank(), tern.axis.title.T = element_blank(), 
      tern.axis.title.L = element_blank(), tern.axis.title.R = element_blank()) + guides(alpha = "none", 
      size = "none") + scale_alpha_continuous(range = c(0.4, 1)) + scale_size_continuous(range = c(0.6, 
      2.5)) + scale_color_gradientn(colors = c(brewer.pal(11, "RdBu")[11], "white", brewer.pal(11, "RdBu")[1]), 
      values = scales::rescale(c(-2, 0, 2), from = c(-2, 2)), limits = c(-2, 2), oob = scales::squish, 
      name = "Delta Cohen d: disease minus treatment")
  df_counts <- df_heatmap %>% group_by(Specificity, Direction) %>% summarise(Count = n(), .groups = "drop")
  levels_spec <- c("20+", "17-19", "14-16", "11-13", "8-10", "5-7", "2-4", "1")
  df_plot <- expand_grid(Specificity = levels_spec, Direction = c("Down", "Up")) %>% left_join(df_counts, 
      by = c("Specificity", "Direction")) %>% mutate(Count = replace_na(Count, 0), Plot_Value = ifelse(Direction == 
      "Down", -Count, Count)) %>% mutate(Specificity = factor(Specificity, levels = levels_spec), Direction = factor(Direction, 
      levels = c("Down", "Up")))
  p_table <- ggplot(df_plot, aes(x = Direction, y = Specificity, fill = Plot_Value)) + geom_tile(color = "black", 
      linewidth = 0.5) + geom_text(aes(label = Count), size = 4.2, color = "black", fontface = "bold") + 
      scale_fill_gradientn(colors = c(brewer.pal(11, "RdBu")[10], "white", brewer.pal(11, "RdBu")[2]), 
          values = scales::rescale(c(min(df_plot$Plot_Value), 0, max(df_plot$Plot_Value))), guide = "none") + 
      scale_x_discrete(position = "bottom") + theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank(), 
      axis.text.x = element_text(size = 11, face = "bold", color = "black"), axis.text.y = element_text(size = 11, 
          face = "bold", color = "black", margin = margin(r = 5)), axis.ticks = element_blank(), plot.margin = margin(5, 
          5, 5, 5)) + coord_fixed(ratio = 0.45)
  canvas <- ggdraw() + draw_plot(p_tern, x = 0.12, y = 0.02, width = 0.56, height = 0.94) + draw_text("High in DBS (Treatment)", 
      x = 0.2, y = 0.45, size = 10.5, angle = 60, color = "#BF8F51", fontface = "bold") + draw_text("High in Sham (Epilepsy)", 
      x = 0.59, y = 0.58, size = 10.5, angle = -60, color = "#73355E", fontface = "bold") + draw_text("High in Saline (Healthy)", 
      x = 0.4, y = 0.08, size = 10.5, color = "#18a799", fontface = "bold") + draw_text("ATN TRGs: neuronal commonness", 
      x = 0.84, y = 0.94, size = 9.5, fontface = "bold", hjust = 0.5) + draw_text("1 = ATN-specific; 2+ = shared", 
      x = 0.84, y = 0.89, size = 9, hjust = 0.5) + draw_plot(p_table, x = 0.72, y = 0.02, width = 0.27, 
      height = 0.82) + draw_text("a", x = 0.01, y = 0.97, size = 18, fontface = "bold")
  png(OUT_PNG, width = 12.5, height = 6.2, units = "in", res = 600)
  print(canvas)
  dev.off()
  cat("Figure 4a composite complete.\n")
  invisible(as.list(environment()))
}
