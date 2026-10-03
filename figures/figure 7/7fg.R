# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(ggplot2)
  })
  path <- Sys.getenv("DBSEQ_IED_CSV")
  if (!nzchar(path)) path <- input_path("sourcedata/figure 7/7fg_ISD.csv")
  if (!supplied_inputs[["data_2"]]) stop("Missing original animal-by-hour IED counts: ", path)
  d <- supplied_inputs[["d_3"]][[as.character(path)]]
  stopifnot(all(c("Mouse_ID", "ZT", "Condition", "Counts") %in% names(d)), !anyDuplicated(d[c("Mouse_ID", 
      "ZT", "Condition")]), all(is.finite(d$Counts)), all(d$Counts >= 0), setequal(unique(d$Condition), 
      c("Baseline", "14d_TRS_DBS")))
  zt <- sort(unique(d$ZT))
  if (identical(as.numeric(zt), as.numeric(0:23))) d$phase <- d$ZT else if (identical(as.numeric(zt), as.numeric(1:24))) d$phase <- d$ZT - 
      1 else stop("ZT must identify exactly 24 hourly bins, 0:23 or 1:24.")
  stopifnot(all(table(d$Mouse_ID, d$Condition) == 24))
  d$Condition <- factor(d$Condition, levels = c("Baseline", "14d_TRS_DBS"))
  paired <- d %>% select(Mouse_ID, ZT, Condition, Counts) %>% pivot_wider(names_from = Condition, values_from = Counts)
  stopifnot(!anyNA(paired))
  hour_stats <- paired %>% group_by(ZT) %>% summarise(N = n(), P = t.test(Baseline, `14d_TRS_DBS`, paired = TRUE)$p.value, 
      .groups = "drop") %>% mutate(FDR = p.adjust(P, "BH"))
  summary <- d %>% group_by(ZT, Condition) %>% summarise(Mean = mean(Counts), SEM = sd(Counts)/sqrt(n()), 
      N = n(), .groups = "drop")
  write.csv(d, output_path("sourcedata/figure 7/7fg.csv"), row.names = FALSE)
  write.csv(hour_stats, output_path("statistics/figure 7/7f_statistics.csv"), row.names = FALSE)
  colors <- c(Baseline = "#a0a0a0", `14d_TRS_DBS` = "#ec8a64")
  ceiling <- max(summary$Mean + summary$SEM) * 1.15
  hour_stats <- hour_stats %>% mutate(Label = ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**", ifelse(FDR < 
      0.05, "*", ""))))
  pf <- ggplot(summary, aes(ZT, Mean, color = Condition)) + geom_line(linewidth = 0.7) + geom_errorbar(aes(ymin = Mean - 
      SEM, ymax = Mean + SEM), width = 0.4, linewidth = 0.35) + geom_point(shape = 21, fill = "white", 
      size = 2.5) + geom_text(data = hour_stats, aes(ZT, ceiling, label = Label), inherit.aes = FALSE, 
      size = 3) + scale_color_manual(values = colors, labels = c(sprintf("Baseline (N=%s)", n_distinct(d$Mouse_ID)), 
      sprintf("14d TRS-DBS (N=%s)", n_distinct(d$Mouse_ID)))) + annotate("segment", x = 12, xend = 24, 
      y = 0, yend = 0, linewidth = 3) + labs(x = "Zeitgeber time (hour)", y = "IED count per hour", color = NULL, 
      subtitle = "Paired animal tests; stars use BH q across 24 hours") + theme_classic(base_size = 10)
  for (ext in c("png")) ggsave(output_path(paste0("figures/figure 7/7f.", ext)), pf, width = 8, height = 4.5, 
      dpi = 400)
  period <- d %>% group_by(Mouse_ID, Condition) %>% summarise(`24 hr` = mean(Counts), Day = mean(Counts[phase < 
      12]), Night = mean(Counts[phase >= 12]), .groups = "drop") %>% pivot_longer(c(`24 hr`, Day, Night), 
      names_to = "Period", values_to = "Frequency") %>% mutate(Period = factor(Period, levels = c("24 hr", 
      "Day", "Night")))
  period_pair <- period %>% pivot_wider(names_from = Condition, values_from = Frequency)
  stopifnot(!anyNA(period_pair))
  period_stats <- period_pair %>% group_by(Period) %>% summarise(N = n(), P = t.test(Baseline, `14d_TRS_DBS`, 
      paired = TRUE)$p.value, .groups = "drop") %>% mutate(FDR = p.adjust(P, "BH"))
  write.csv(period, output_path("sourcedata/figure 7/7g.csv"), row.names = FALSE)
  write.csv(period_stats, output_path("statistics/figure 7/7g_statistics.csv"), row.names = FALSE)
  pg <- ggplot(period, aes(Condition, Frequency)) + geom_line(aes(group = Mouse_ID), color = "grey65", 
      linewidth = 0.4) + geom_point(aes(color = Condition), size = 2) + scale_color_manual(values = colors, 
      guide = "none") + facet_wrap(~Period, nrow = 1) + geom_text(data = period_stats, aes(x = 1.5, y = Inf, 
      label = sprintf("P=%.4f\nBH q=%.4f", P, FDR)), inherit.aes = FALSE, vjust = 1.2, size = 3) + labs(x = NULL, 
      y = "Frequency (IED/hour)") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + theme_classic(base_size = 10)
  for (ext in c("png")) ggsave(output_path(paste0("figures/figure 7/7g.", ext)), pg, width = 7, height = 4.5, 
      dpi = 400)
  invisible(as.list(environment()))
}
