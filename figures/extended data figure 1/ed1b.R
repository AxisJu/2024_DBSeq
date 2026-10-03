# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(patchwork)
      library(scales)
  })
  theme_qc <- function() theme_classic(base_size = 5, base_family = "Arial") + theme(plot.title = element_text(size = 5.5, 
      face = "bold", hjust = 0.5, margin = margin(b = 2)), plot.subtitle = element_text(size = 4.5, hjust = 0.5, 
      margin = margin(b = 3)), axis.title = element_text(size = 5, face = "bold"), axis.text = element_text(size = 4.8, 
      color = "black"), axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35), 
      axis.ticks.length = unit(0.8, "mm"), strip.background = element_rect(fill = "grey95", color = "black", 
          linewidth = 0.35), strip.text = element_text(size = 4.8, face = "bold", margin = margin(t = 2, 
          b = 2, unit = "pt")), panel.grid = element_blank(), plot.margin = margin(3, 3, 3, 3, "pt"))
  if (!nzchar(input_root)) stop("Set DBSEQ_INPUT_PACKAGE_ROOT to the source data package.")
  d <- supplied_inputs[["d_2"]][[as.character(file.path(input_root, "sourcedata/extended data figure 1/ed1b.csv"))]]
  stopifnot(all(is.finite(d$doublet_scores)), all(d$doublet_scores >= 0 & d$doublet_scores <= 1))
  export_data(d, "ed1b")
  order <- c("DBS_I_1", "DBS_I_2", "DBS_C_1", "Sham_I_1", "Sham_I_2", "Sham_I_3", "Saline_I_1", "Saline_I_2", 
      "KA_1")
  d$display_sample <- factor(d$display_sample, levels = order)
  thresholds <- d %>% distinct(display_sample, automatic_threshold_from_notebook, prefilter_detected_doublet_percent_from_notebook)
  names(thresholds)[2:3] <- c("Threshold", "Detected_Rate")
  p <- ggplot(d, aes(doublet_scores)) + geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#bdd7e7", 
      color = "#6baed6", linewidth = 0.22) + geom_density(color = "#08519c", linewidth = 0.4) + geom_vline(data = thresholds, 
      aes(xintercept = Threshold), color = "#e31a1c", linetype = "dashed", linewidth = 0.35) + geom_text(data = thresholds, 
      aes(x = ifelse(Threshold > 0.7, Threshold - 0.04, Threshold + 0.04), y = Inf, hjust = ifelse(Threshold > 
          0.7, 1, 0), label = sprintf("T:%.2f\n%.1f%%", Threshold, Detected_Rate)), vjust = 1.25, size = 4.2/.pt, 
      fontface = "bold", color = "#b30000") + facet_wrap(~display_sample, nrow = 1, scales = "free_y") + 
      coord_cartesian(xlim = c(0, 1)) + scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", 
      "1")) + labs(title = "Single-Nucleus Doublet Score Distributions", x = "Scrublet Doublet Score", 
      y = "Density") + theme_qc() + theme(plot.title = element_text(hjust = 0))
  save_plot(p, "ed1b", 172, 52)
  cat("Extended Data Figure 1b rendered with original ggplot methods.\n")
  invisible(as.list(environment()))
}
