# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ggplot2)
  dir.create(output_path("figures/figure 5"), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_path("sourcedata/figure 5"), recursive = TRUE, showWarnings = FALSE)
  shared <- supplied_inputs[["shared_2"]][[as.character(cache_path("5kl_shared_data.rds"))]]
  stopifnot(identical(shared$schema_version, 2L))
  all_neuronal_trg <- shared$trg_class
  df_degree_clean <- shared$df_degree_clean
  df_scatter_ppi <- all_neuronal_trg %>% mutate(JoinKey = toupper(trimws(as.character(genes)))) %>% inner_join(df_degree_clean %>% 
      select(JoinKey, Degree), by = "JoinKey") %>% filter(Degree > 0) %>% mutate(Log_Degree = log10(Degree + 
      1), Color_Value = (Log_Degree - min(Log_Degree))/(max(Log_Degree) - min(Log_Degree)))
  cor_res <- cor.test(df_scatter_ppi$n_neurons, df_scatter_ppi$Log_Degree, method = "spearman", exact = FALSE)
  rho_label <- sprintf("%.2f", cor_res$estimate)
  p_label <- formatC(cor_res$p.value, format = "e", digits = 2)
  cor_text <- paste0("Spearman's <U+03C1> = ", rho_label, "\n", "P = ", p_label)
  csv_out <- df_scatter_ppi %>% select(Gene = genes, n_neurons, Degree, Log10_Degree_Plus1 = Log_Degree)
  write.csv(csv_out, output_path("sourcedata/figure 5/5l.csv"), row.names = FALSE)
  write.csv(data.frame(n_genes = nrow(df_scatter_ppi), rho = unname(cor_res$estimate), p_value = cor_res$p.value), 
      output_path("statistics/figure 5/5l_statistics.csv"), row.names = FALSE)
  set.seed(42)
  p_scatter <- ggplot(df_scatter_ppi, aes(x = n_neurons, y = Log_Degree)) + geom_boxplot(aes(group = n_neurons), 
      fill = "grey90", color = "grey60", alpha = 0.5, outlier.shape = NA, width = 0.6) + geom_jitter(aes(color = Color_Value), 
      alpha = 0.6, size = 1.5, position = position_jitter(width = 0.25, height = 0, seed = 42)) + scale_color_gradientn(colors = c("#7fcdbb", 
      "#2c7fb8", "#253494")) + geom_smooth(method = "lm", color = "#000000", fill = "#cccccc", alpha = 0.3, 
      linewidth = 1.2, linetype = "dashed") + annotate("text", x = Inf, y = Inf, label = cor_text, hjust = 1.05, 
      vjust = 1.3, size = 4, fontface = "bold") + scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + 
      theme_bw(base_size = 14) + theme(legend.position = "none", axis.title = element_text(face = "bold"), 
      axis.text = element_text(face = "bold", color = "black"), panel.grid.minor = element_blank()) + labs(x = "Commonness (Num. of Neuronal Types)", 
      y = expression(log[10](PPI ~ Degree + 1)))
  ggsave(output_path("figures/figure 5/5l.png"), p_scatter, width = 5, height = 4.5, dpi = 600)
  cat(">>> Figure 5l complete.\n")
  invisible(as.list(environment()))
}
