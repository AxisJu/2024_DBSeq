# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(ggplot2)
      library(dplyr)
      library(tidyr)
      library(ggrepel)
      library(cowplot)
  })
  if (.Platform$OS.type == "windows") grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
  package <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(package)) stop("Set DBSEQ_OUTPUT_ROOT")
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  code_base <- normalizePath(file.path(dirname(normalizePath(entry)), "../.."), winslash = "/")
  out_base <- normalizePath(package, winslash = "/", mustWork = FALSE)
  if (identical(tolower(out_base), tolower(code_base)) || startsWith(tolower(out_base), paste0(tolower(code_base), 
      "/"))) stop("Output must be outside the code package")
  workspace <- package
  layout_dir <- file.path(workspace, "cache/extended_data")
  dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
  cowplot::set_null_device("png")
  grDevices::png(file.path(layout_dir, "ed13-layout-preview.png"), width = 1600, height = 1600, res = 200, 
      type = "cairo")
  input <- supplied_inputs[["input_1"]][[as.character(file.path(package, "sourcedata/extended data figure 13/ed13b.csv"))]]
  metrics <- c("FC", "ED", "SC", "SC2")
  plot_data <- input %>% pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>% mutate(Metric = factor(Metric, 
      levels = metrics)) %>% group_by(Metric) %>% mutate(Color = (Value - min(Value))/(max(Value) - min(Value))) %>% 
      ungroup()
  labels <- lapply(metrics, function(m) {
      test <- cor.test(input$commonness, input[[m]], method = "spearman", exact = FALSE)
      data.frame(Metric = m, Rho = unname(test$estimate), P = test$p.value)
  }) %>% bind_rows() %>% mutate(Metric = factor(Metric, levels = metrics))
  labels$FDR <- p.adjust(labels$P, "BH")
  labels$Label <- sprintf("Spearman <U+03C1> = %.2f\nBH FDR = %.3g", labels$Rho, labels$FDR)
  dir.create(file.path(package, "statistics/extended data figure 13"), recursive = TRUE, showWarnings = FALSE)
  write.csv(labels, file.path(package, "statistics/extended data figure 13/ed13b_tests.csv"), row.names = FALSE)
  figure <- ggplot(plot_data, aes(commonness, Value, color = Color)) + geom_point(alpha = 0.3, size = 0.8) + 
      scale_color_gradientn(colors = c("#7fcdbb", "#2c7fb8", "#253494")) + geom_smooth(method = "loess", 
      span = 1, color = "black", fill = "#cccccc", alpha = 0.6, linewidth = 0.45, linetype = "dashed") + 
      geom_text(data = labels, aes(x = -Inf, y = Inf, label = Label), inherit.aes = FALSE, hjust = -0.06, 
          vjust = 1.1, size = 2.2) + facet_wrap(~Metric, scales = "free_y", ncol = 2, labeller = as_labeller(c(FC = "Functional Connectivity", 
      ED = "Euclidean Distance", SC = "1st Structural Connectivity", SC2 = "2nd Structural Connectivity"))) + 
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) + labs(x = "TRG commonness (neuronal types)", 
      y = NULL) + theme_bw(base_size = 7, base_family = "Arial") + theme(legend.position = "none", panel.grid.minor = element_blank(), 
      strip.background = element_rect(fill = "grey95"))
  ggsave(file.path(package, "figures/extended data figure 13/ed13b.png"), figure, width = 5, height = 5, 
      units = "in", dpi = 400)
  cat("Extended Data Figure 13b: original quadratic LOESS curve and confidence band\n")
  source_dir <- file.path(package, "sourcedata/extended data figure 13")
  output <- file.path(package, "figures/extended data figure 13")
  read_panel <- function(panel) supplied_inputs[["read_panel_2"]][[as.character(file.path(source_dir, paste0("ed13", 
      panel, ".csv")))]]
  save_panel <- function(plot, panel, w, h) ggsave(file.path(output, paste0("ed13", panel, ".png")), plot, 
      width = w, height = h, units = "mm", dpi = 600)
  names_metric <- c(FC = "Functional Connectivity", ED = "Euclidean Distance", SC = "1st Structural Connectivity", 
      SC2 = "2nd Structural Connectivity", SPE = "Shortest Path Efficiency (SPE)", DE = "Diffusion Efficiency (DE)", 
      NE = "Navigation Efficiency (NE)", CMY = "Communicability (CMY)")
  colors <- c(HY = "#e74538", CNU = "#9bd5f4", Isocortex = "#18a799", HPF = "#84c551", TH = "#e0abce")
  regions <- c(RT = "#ff6600", CM = "#08306d", CP = "#8d7acc", RE = "#1340ff", GPe = "#b199ff", RSPd = "#ff1a71", 
      LSc = "#3283fe", TRS = "#7609b1", MH = "#faa307", LH = "#c68105", PF = "#0a4093", DG = "#16f2f2")
  tiers <- c("2-4", "5-7", "8-10", "11-13", "14-16", "17-19", "20+")
  cdata <- supplied_inputs[["cdata_3"]][[as.character("c")]] %>% distinct(bin, Region, count, proportion) %>% 
      mutate(bin = factor(bin, levels = tiers), Region = factor(Region, levels = rev(names(regions))))
  labels_region <- cdata %>% filter(bin == "20+") %>% arrange(match(as.character(Region), names(regions))) %>% 
      mutate(Position = cumsum(proportion) - proportion/2)
  pc <- ggplot(cdata, aes(bin, proportion, fill = Region)) + geom_col(width = 0.7, color = "black", linewidth = 0.5) + 
      scale_fill_manual(values = regions) + geom_text(data = labels_region, aes(x = 7.6, y = Position, 
      label = Region, color = Region), inherit.aes = FALSE, hjust = 0, size = 2.2) + scale_color_manual(values = regions, 
      guide = "none") + scale_x_discrete(expand = expansion(add = c(0.5, 2.1))) + scale_y_continuous(labels = function(x) paste0(x * 
      100, "%"), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) + labs(title = "Proportion of Regions", 
      x = "Commonness of TRG", y = NULL) + theme_classic(base_size = 7, base_family = "Arial") + theme(axis.text.x = element_text(angle = 45, 
      hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
  save_panel(pc, "c", 52, 95)
  plot_region <- function(panel, metric_names, partial) {
      frame <- supplied_inputs[["frame_4"]][[as.character(panel)]]
      test_rows <- lapply(metric_names, function(metric) {
          if (partial) {
              rx <- residuals(lm(rank(frame$commonness) ~ rank(frame$ED)))
              ry <- residuals(lm(rank(frame[[metric]]) ~ rank(frame$ED)))
              rho <- cor(rx, ry)
              p <- 2 * pt(-abs(rho) * sqrt((nrow(frame) - 3)/max(1 - rho * rho, 1e-15)), df = nrow(frame) - 
                  3)
          }
          else {
              test <- cor.test(frame$commonness, frame[[metric]], method = "spearman", exact = FALSE)
              rho <- unname(test$estimate)
              p <- test$p.value
          }
          data.frame(Metric = metric, Rho = rho, P = p, N = nrow(frame))
      }) %>% bind_rows()
      test_rows$FDR <- p.adjust(test_rows$P, "BH")
      write.csv(test_rows, file.path(package, "statistics/extended data figure 13", paste0("ed13", panel, 
          "_tests.csv")), row.names = FALSE)
      plots <- lapply(seq_along(metric_names), function(i) {
          metric <- metric_names[i]
          if (partial) {
              rx <- residuals(lm(rank(frame$commonness) ~ rank(frame$ED)))
              ry <- residuals(lm(rank(frame[[metric]]) ~ rank(frame$ED)))
              rho <- cor(rx, ry)
              p <- 2 * pt(-abs(rho) * sqrt((nrow(frame) - 3)/(1 - rho * rho)), df = nrow(frame) - 3)
          }
          else {
              test <- cor.test(frame$commonness, frame[[metric]], method = "spearman", exact = FALSE)
              rho <- unname(test$estimate)
              p <- test$p.value
          }
          label <- paste0(if (partial) 
              "Partial Spearman's <U+03C1> = "
          else "Spearman's <U+03C1> = ", sprintf("%.2f", rho), paste0("\nBH FDR = ", format(test_rows$FDR[i], 
              digits = 2)))
          d <- frame %>% mutate(Value = .data[[metric]], Metric = names_metric[[metric]])
          ggplot(d, aes(commonness, Value)) + geom_smooth(method = "lm", color = "black", fill = "#cccccc", 
              alpha = 0.25, linewidth = 0.6, linetype = "dashed") + geom_point(aes(color = Division), size = 2, 
              alpha = 0.8) + scale_color_manual(values = colors) + geom_text_repel(aes(label = Region), 
              seed = 260924, size = 2.1, color = "grey20", box.padding = 0.2, max.overlaps = Inf) + annotate("text", 
              x = -Inf, y = Inf, label = label, hjust = -0.05, vjust = 1.1, size = 2.1) + facet_wrap(~Metric, 
              scales = "free_y") + scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) + labs(x = if (i > 
              2) 
              "Weighted Mean Commonness of Region"
          else NULL, y = NULL) + theme_bw(base_size = 7, base_family = "Arial") + theme(panel.grid.minor = element_blank(), 
              panel.grid.major = element_line(color = "#f3f3f3", linewidth = 0.25), strip.background = element_rect(fill = "grey93", 
                  color = "black", linewidth = 0.5), strip.text = element_text(face = "bold", size = 6.3), 
              legend.position = "none", axis.title.x = element_text(size = 6.3, face = "bold"))
      })
      save_panel(plot_grid(plotlist = plots, nrow = 2, ncol = 2, align = "hv"), panel, 100, 100)
  }
  plot_region("d", c("FC", "ED", "SC", "SC2"), FALSE)
  plot_region("e", c("SPE", "DE", "NE", "CMY"), TRUE)
  cat("Extended Data Figure 13 c-e: restored bars, facet strips, labels and original statistics\n")
  invisible(grDevices::dev.off())
  invisible(as.list(environment()))
}
