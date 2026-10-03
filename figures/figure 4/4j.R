# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(qs)
  list2env(supplied_inputs[["helpers_2"]], envir = environment())
  set.seed(42)
  region_to_division <- c(CP = "CNU", GPe = "CNU", LSc = "CNU", TRS = "CNU", CM = "TH", PF = "TH", RE = "TH", 
      RT = "TH", MH = "TH", LH = "HY", RSPd = "Isocortex", DG = "HPF")
  division_colors <- c(HY = "#E74538", CNU = "#9BD5F4", Isocortex = "#18A799", HPF = "#84C551", TH = "#E0ABCE")
  metrics_4j <- c("ED", "FC", "SC", "SC2")
  metric_labels_4j <- c(ED = "Euclidean Distance", FC = "Functional Connectivity", SC = "1st Structural Connectivity", 
      SC2 = "2nd Structural Connectivity")
  df_plot_4j <- df_jk %>% dplyr::mutate(Division = region_to_division[Region]) %>% tidyr::pivot_longer(cols = all_of(metrics_4j), 
      names_to = "Metric", values_to = "Value") %>% dplyr::mutate(Metric_Name = factor(metric_labels_4j[Metric], 
      levels = metric_labels_4j))
  stats_4j <- bind_rows(lapply(metrics_4j, function(metric) {
      z <- df_jk %>% dplyr::select(Region, Enriched_Count, all_of(unique(c(metric, "ED")))) %>% filter(if_all(-Region, 
          is.finite))
      stopifnot(nrow(z) >= 4)
      test <- cor.test(z$Enriched_Count, z[[metric]], method = "spearman", exact = FALSE)
      data.frame(Metric = metric, N = nrow(z), Rho = unname(test$estimate), P = test$p.value)
  }))
  stats_4j$FDR <- p.adjust(stats_4j$P, "BH")
  stats_4j$Metric_Name <- factor(metric_labels_4j[stats_4j$Metric], levels = metric_labels_4j)
  stats_layer1_4j <- stats_4j %>% mutate(Label1 = sprintf("Spearman's <U+03C1> = %.2f", Rho))
  stats_layer2_4j <- stats_4j %>% mutate(Label2 = paste0("BH FDR = ", format.pval(FDR, digits = 3, eps = .Machine$double.xmin)))
  write.csv(stats_4j %>% dplyr::select(-Metric_Name), output_path("statistics/figure 4/4j_statistics.csv"), 
      row.names = FALSE)
  p_4j <- ggplot(df_plot_4j, aes(x = Enriched_Count, y = Value)) + geom_smooth(method = "lm", se = TRUE, 
      color = "black", linetype = "dashed", linewidth = 0.8, fill = "grey85", alpha = 0.5) + geom_point(aes(color = Division), 
      size = 3.5, alpha = 0.85) + geom_text_repel(aes(label = Region), fontface = "italic", size = 3.4, 
      color = "grey25", box.padding = 0.35, point.padding = 0.25, max.overlaps = 30) + geom_text(data = stats_layer1_4j, 
      aes(x = -Inf, y = Inf, label = Label1), hjust = -0.06, vjust = 1.4, size = 3.5, color = "black", 
      inherit.aes = FALSE) + geom_text(data = stats_layer2_4j, aes(x = -Inf, y = Inf, label = Label2), 
      hjust = -0.06, vjust = 2.9, size = 3.5, fontface = "italic", color = "black", inherit.aes = FALSE) + 
      facet_wrap(~Metric_Name, scales = "free_y", ncol = 2) + scale_color_manual(values = division_colors) + 
      scale_x_continuous(breaks = c(5, 10, 15, 20, 25), limits = c(4, 27)) + scale_y_continuous(expand = expansion(mult = c(0.08, 
      0.28))) + theme_bw(base_size = 11) + theme(legend.position = "none", strip.background = element_rect(fill = "grey92", 
      color = "black", linewidth = 0.7), strip.text = element_text(face = "bold", size = 10.5, margin = margin(t = 4, 
      b = 4)), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7), panel.grid.major = element_line(color = "grey95", 
      linewidth = 0.3), panel.grid.minor = element_blank(), axis.title.x = element_text(face = "bold", 
      size = 11, margin = margin(t = 8)), axis.title.y = element_blank(), axis.text = element_text(color = "black", 
      size = 9.5), axis.ticks = element_line(color = "black", linewidth = 0.5), plot.margin = margin(t = 8, 
      r = 12, b = 8, l = 8)) + labs(x = "Num. of enriched TRGs in RNA Splicing")
  out_fig_dir <- output_path("figures/figure 4")
  out_data_dir <- output_path("sourcedata/figure 4")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_fig_dir, "4j.png"), plot = p_4j, width = 6.2, height = 5.2, dpi = 600, type = "cairo")
  df_source_j <- df_jk %>% dplyr::mutate(Division = region_to_division[Region]) %>% dplyr::select(Region_Acronym = Region, 
      Brain_Division = Division, Enriched_TRG_Count = Enriched_Count, Euclidean_Distance = ED, Functional_Connectivity = FC, 
      First_Structural_Connectivity = SC, Second_Structural_Connectivity = SC2)
  write.csv(df_source_j, file.path(out_data_dir, "4j.csv"), row.names = FALSE)
  cat("Figure 4j generated successfully.\n")
  invisible(as.list(environment()))
}
