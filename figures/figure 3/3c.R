# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
      library(tidyr)
      library(stringr)
      library(ggplot2)
      library(patchwork)
  })
  base_dir <- data_root
  output_dir_fig <- output_path("figures/figure 3")
  output_dir_data <- output_path("sourcedata/figure 3")
  dir.create(output_dir_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  dbs_qs_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/NeuronChat/DBS_I.qs")
  sham_qs_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/NeuronChat/Sham_I.qs")
  saline_qs_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/NeuronChat/Saline_I.qs")
  x_DBS_I <- supplied_inputs[["x_DBS_I_2"]][[as.character(dbs_qs_path)]]
  x_Sham_I <- supplied_inputs[["x_Sham_I_3"]][[as.character(sham_qs_path)]]
  x_Saline_I <- supplied_inputs[["x_Saline_I_4"]][[as.character(saline_qs_path)]]
  all_samples <- list(x_Saline_I@net, x_Sham_I@net, x_DBS_I@net)
  target_total <- 500
  all_samples_norm <- lapply(all_samples, function(sample) {
      current_total <- sum(sapply(sample, sum))
      scale_factor <- target_total/current_total
      lapply(sample, function(mat) mat * scale_factor)
  })
  net_saline <- all_samples_norm[[1]]
  net_sham <- all_samples_norm[[2]]
  net_dbs <- all_samples_norm[[3]]
  result_in <- NULL
  for (pair in names(net_saline)) {
      Pair <- pair
      DBSvSham <- sum(net_dbs[[pair]][, "ExN_ATN"] - net_sham[[pair]][, "ExN_ATN"])
      ShamvSaline <- sum(net_sham[[pair]][, "ExN_ATN"] - net_saline[[pair]][, "ExN_ATN"])
      if (abs(DBSvSham) > 0.05 | abs(ShamvSaline) > 0.05) {
          result_in <- rbind(result_in, data.frame(Direction = "Input to ATN", Pair, DBSvSham, ShamvSaline))
      }
  }
  df_in_top5 <- result_in %>% arrange(desc(abs(DBSvSham)), Pair) %>% slice_head(n = 5) %>% mutate(Pair_Clean = str_replace_all(Pair, 
      "_", "-"))
  result_out <- NULL
  for (pair in names(net_saline)) {
      Pair <- pair
      DBSvSham <- sum(net_dbs[[pair]]["ExN_ATN", ] - net_sham[[pair]]["ExN_ATN", ])
      ShamvSaline <- sum(net_sham[[pair]]["ExN_ATN", ] - net_saline[[pair]]["ExN_ATN", ])
      if (abs(DBSvSham) > 0.05 | abs(ShamvSaline) > 0.05) {
          result_out <- rbind(result_out, data.frame(Direction = "Output from ATN", Pair, DBSvSham, ShamvSaline))
      }
  }
  df_out_top5 <- result_out %>% arrange(desc(abs(DBSvSham)), Pair) %>% slice_head(n = 5) %>% mutate(Pair_Clean = str_replace_all(Pair, 
      "_", "-"))
  sourcedata_c <- rbind(df_in_top5, df_out_top5) %>% pivot_longer(cols = c(ShamvSaline, DBSvSham), names_to = "Comparison", 
      values_to = "Delta_Probability") %>% mutate(Comparison = ifelse(Comparison == "ShamvSaline", "Sham vs Saline", 
      "DBS vs Sham"), Pair = Pair_Clean) %>% dplyr::select(Direction, Pair, Comparison, Delta_Probability)
  write.csv(sourcedata_c, file.path(output_dir_data, "3c.csv"), row.names = FALSE)
  cat(sprintf("Exported Source Data: %s (Rows: %d)\n", file.path(output_dir_data, "3c.csv"), nrow(sourcedata_c)))
  get_y_colors <- function(levels_vec) {
      ifelse(levels_vec == "Nrxn1-Nlgn1", "#800000", "black")
  }
  df_in_long <- df_in_top5 %>% mutate(Pair_Clean = factor(Pair_Clean, levels = rev(unique(Pair_Clean)))) %>% 
      pivot_longer(cols = c(DBSvSham, ShamvSaline), names_to = "Condition", values_to = "Score") %>% mutate(Condition = factor(Condition, 
      levels = c("ShamvSaline", "DBSvSham"), labels = c("Sham vs Saline", "DBS vs Sham")))
  in_levels <- levels(df_in_long$Pair_Clean)
  in_colors <- get_y_colors(in_levels)
  p_in <- ggplot(df_in_long, aes(x = Score, y = Pair_Clean, fill = Condition)) + geom_bar(stat = "identity", 
      position = "dodge", color = "black", linewidth = 0.35, alpha = 0.85) + scale_fill_manual(values = c(`Sham vs Saline` = "#ac00cc", 
      `DBS vs Sham` = "#ffd700")) + theme_bw(base_size = 11) + labs(title = "Input to ATN", x = NULL, y = NULL, 
      fill = NULL) + theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold", color = in_colors, 
      size = 10), axis.text.x = element_text(color = "black", size = 9), plot.title = element_text(face = "bold", 
      size = 11, hjust = 0), legend.position = "top", legend.justification = "right")
  df_out_long <- df_out_top5 %>% mutate(Pair_Clean = factor(Pair_Clean, levels = rev(unique(Pair_Clean)))) %>% 
      pivot_longer(cols = c(DBSvSham, ShamvSaline), names_to = "Condition", values_to = "Score") %>% mutate(Condition = factor(Condition, 
      levels = c("ShamvSaline", "DBSvSham"), labels = c("Sham vs Saline", "DBS vs Sham")))
  out_levels <- levels(df_out_long$Pair_Clean)
  out_colors <- get_y_colors(out_levels)
  p_out <- ggplot(df_out_long, aes(x = Score, y = Pair_Clean, fill = Condition)) + geom_bar(stat = "identity", 
      position = "dodge", color = "black", linewidth = 0.35, alpha = 0.85) + scale_fill_manual(values = c(`Sham vs Saline` = "#ac00cc", 
      `DBS vs Sham` = "#ffd700")) + theme_bw(base_size = 11) + labs(title = "Output from ATN", x = "<U+0394> CCC Probability", 
      y = NULL, fill = NULL) + theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold", 
      color = out_colors, size = 10), axis.text.x = element_text(color = "black", size = 9), axis.title.x = element_text(face = "bold", 
      size = 10), plot.title = element_text(face = "bold", size = 11, hjust = 0), legend.position = "none")
  p_combined <- (p_in/p_out) + plot_layout(heights = c(1, 1))
  fig_pdf <- file.path(output_dir_fig, "3c.pdf")
  fig_png <- file.path(output_dir_fig, "3c.png")
  cat(sprintf("Exported Figures: %s and %s\n", fig_pdf, fig_png))
  ggsave(fig_png, p_combined, width = 5, height = 7, dpi = 600)
  invisible(as.list(environment()))
}
