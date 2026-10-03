# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ggplot2)
  metric_csv <- data_path("results/260623_squidff/plots_v7/metric_comparison.csv")
  out_pdf <- output_path("figures/figure 7/7e.pdf")
  out_png <- output_path("figures/figure 7/7e.png")
  out_csv <- output_path("sourcedata/figure 7/7e.csv")
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  region_cols <- c(ATN = "#1460ff", RSPd = "#ff1a71", DG = "#16f2f2", GPe = "#b199ff", CP = "#8d7acc", 
      LSc = "#3283fe", RT = "#ff6600", TRS = "#7609b1", MH = "#faa307", LH = "#c68105", RE = "#1340ff", 
      PF = "#0a4093", CM = "#08306d")
  corrected <- Sys.getenv("DBSEQ_WRS_DIR")
  if (nzchar(corrected)) {
      raw_scores <- supplied_inputs[["raw_scores_2"]][[as.character(file.path(corrected, "regional_sensitivity_scores.csv"))]]
      bootstrap <- supplied_inputs[["bootstrap_3"]][[as.character(file.path(corrected, "donor_bootstrap_scores.csv"))]]
      interval <- bootstrap %>% group_by(region) %>% summarise(Lower = quantile(weighted_reversal_score, 
          0.025), Upper = quantile(weighted_reversal_score, 0.975), .groups = "drop")
      df <- raw_scores %>% transmute(Region = region, wrs_mean = weighted_reversal_score) %>% left_join(interval, 
          by = c(Region = "region"))
      note <- "Encoder-gradient sensitivity; donor-bootstrap 95% interval; calibration cohort"
      write.csv(raw_scores, out_csv, row.names = FALSE)
      write.csv(bootstrap, output_path("statistics/figure 7/7e_bootstrap.csv"), row.names = FALSE)
  } else {
      df <- supplied_inputs[["df_4"]][[as.character(metric_csv)]]
      df$Lower <- df$wrs_mean - df$wrs_sem
      df$Upper <- df$wrs_mean + df$wrs_sem
      note <- "Historical calibration score; error bars are cell-bootstrap uncertainty, not animal SEM"
      write.csv(df %>% select(Region, wrs_mean, wrs_sem), out_csv, row.names = FALSE)
  }
  df$Region <- factor(df$Region, levels = names(region_cols))
  atn_wrs <- df$wrs_mean[df$Region == "ATN"]
  p_wrs <- ggplot(df, aes(x = wrs_mean, y = reorder(Region, wrs_mean), fill = Region)) + geom_col(colour = "black", 
      linewidth = 0.3, width = 0.7) + geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.3, colour = "black", 
      linewidth = 0.4) + scale_fill_manual(values = region_cols, guide = "none") + geom_vline(xintercept = atn_wrs, 
      colour = "#1460ff", linetype = "dashed", linewidth = 0.8, alpha = 0.7) + annotate("text", x = atn_wrs, 
      y = Inf, label = "ATN", vjust = -0.3, hjust = 0.5, colour = "#1460ff", size = 2.8) + labs(x = "Weighted Reversal Score (WRS)", 
      y = NULL, caption = note) + theme_classic(base_size = 10) + coord_cartesian(clip = "off") + theme(axis.text.y = element_text(size = 9), 
      plot.margin = margin(16, 6, 6, 6))
  ggsave(out_png, p_wrs, width = 5, height = 5, dpi = 600)
  cat(">>> Figure 7e complete.\n")
  invisible(as.list(environment()))
}
