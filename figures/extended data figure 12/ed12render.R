# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(ggplot2)
      library(dplyr)
      library(tidyr)
      library(ggrepel)
      library(cowplot)
  })
  package <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(package)) stop("Set DBSEQ_OUTPUT_ROOT")
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  code_base <- normalizePath(file.path(dirname(normalizePath(entry)), "../.."), winslash = "/")
  out_base <- normalizePath(package, winslash = "/", mustWork = FALSE)
  if (identical(tolower(out_base), tolower(code_base)) || startsWith(tolower(out_base), paste0(tolower(code_base), 
      "/"))) stop("Output must be outside the code package")
  if (.Platform$OS.type == "windows") grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
  workspace <- package
  layout_dir <- file.path(workspace, "cache/extended_data")
  dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
  cowplot::set_null_device("png")
  grDevices::png(file.path(layout_dir, "ed12-layout-preview.png"), width = 1600, height = 1600, res = 200, 
      type = "cairo")
  source_dir <- file.path(package, "sourcedata/extended data figure 12")
  output <- file.path(package, "figures/extended data figure 12")
  read_panel <- function(panel) supplied_inputs[["read_panel_1"]][[as.character(file.path(source_dir, paste0("ed12", 
      panel, ".csv")))]]
  save_panel <- function(plot, panel, w, h) ggsave(file.path(output, paste0("ed12", panel, ".png")), plot, 
      width = w, height = h, units = "mm", dpi = 600)
  theme_formal <- function() theme_classic(base_size = 5, base_family = "Arial") + theme(axis.line = element_line(linewidth = 0.5), 
      axis.ticks = element_line(linewidth = 0.5), panel.grid.major = element_line(color = "#f8f9fa", linewidth = 0.25), 
      plot.margin = margin(3, 3, 2, 3, "pt"))
  shorten <- function(x) {
      x <- sub(" \\(GO:[0-9]+\\)", "", x)
      x <- sub("Endoplasmic Reticulum To Golgi Vesicle-Mediated Transport", "ER to Golgi Vesicle Transport", 
          x)
      x <- sub("Regulation Of mRNA Splicing, Via Spliceosome", "Reg. of mRNA Splicing", x)
      x <- sub("Protein Insertion Into Mitochondrial Membrane", "Protein Insertion into Mito. Membrane", 
          x)
      x
  }
  correlations <- function(frame) frame %>% group_by(pathway) %>% summarise(rho = cor(TRG_hits, SC2, method = "spearman"), 
      p = cor.test(TRG_hits, SC2, method = "spearman", exact = FALSE)$p.value, .groups = "drop")
  all_data <- supplied_inputs[["all_data_2"]][[as.character("b")]]
  tests <- supplied_inputs[["tests_3"]][[as.character(file.path(package, "statistics/extended data figure 12/ed12_pathway_tests.csv"))]]
  all_cor <- tests[tests$ATN_hits >= 6, ]
  all_cor$p <- all_cor$q_sc2
  pie <- data.frame(Category = factor(c("Positive; BH FDR < 0.05", "Positive Cor.", "Negative Cor."), levels = c("Positive; BH FDR < 0.05", 
      "Positive Cor.", "Negative Cor.")), Count = c(sum(all_cor$rho > 0 & all_cor$p < 0.05), sum(all_cor$rho > 
      0 & all_cor$p >= 0.05), sum(all_cor$rho < 0)))
  pie$Percent <- pie$Count/sum(pie$Count)
  pie$Position <- cumsum(pie$Percent) - pie$Percent/2
  pa <- ggplot(pie, aes(x = 1, y = Percent, fill = Category)) + geom_col(width = 1, color = "white", linewidth = 0.5) + 
      coord_polar("y", start = 0) + scale_fill_manual(values = c("#18a799", "#9bd5f4", "#e74538"), drop = FALSE) + 
      geom_text(data = pie[pie$Count > 0, ], aes(y = Position, label = sprintf("%.2f%%", 100 * Percent)), 
          fontface = "bold", size = 2.8, color = c("white", "#1a365d")) + labs(title = paste0("All Involved Pathways\nN = ", 
      nrow(all_cor)), fill = "Spearman correlation with 2nd SC (BH)") + theme_void(base_size = 5, base_family = "Arial") + 
      theme(plot.title = element_text(hjust = 0.2, face = "bold"), legend.position = "right")
  save_panel(pa, "b", 85, 54)
  bc <- tests %>% filter(pathway %in% supplied_inputs[["bc_4"]][[as.character("a")]]$pathway) %>% mutate(p = q_sc2) %>% 
      arrange(desc(rho), pathway)
  bc$pathway <- factor(bc$pathway, levels = bc$pathway)
  pb <- ggplot(bc, aes(rho, pathway)) + geom_segment(aes(x = 0, xend = rho, yend = pathway), linewidth = 0.5, 
      color = "#555555") + geom_point(aes(size = abs(rho), fill = p < 0.05), shape = 21, color = "#343a40", 
      stroke = 0.5) + scale_size_continuous(range = c(2.5, 4.5), guide = "none") + scale_fill_manual(values = c(`FALSE` = "#dee2e6", 
      `TRUE` = "#18a799"), labels = c(`FALSE` = "Not Significant", `TRUE` = "Significant")) + scale_y_discrete(labels = function(x) vapply(shorten(x), 
      function(v) paste(strwrap(v, width = 42), collapse = "\n"), character(1))) + scale_x_continuous(limits = c(-1, 
      1), expand = expansion(mult = c(0.02, 0.02))) + labs(x = "Spearman correlation with 2nd SC", y = NULL, 
      fill = "BH FDR < 0.05") + theme_formal() + theme(axis.text = element_text(size = 6), legend.position = "top", 
      legend.justification = "right")
  save_panel(pb, "a", 145, 90)
  cdata <- supplied_inputs[["cdata_5"]][[as.character("c")]]
  terms <- c("Aerobic Electron Transport Chain (GO:0019646)", "Spliceosomal Complex Assembly (GO:0000245)", 
      "protein-DNA Complex Organization (GO:0071824)", "Endoplasmic Reticulum To Golgi Vesicle-Mediated Transport (GO:0006888)", 
      "Regulation Of mRNA Splicing, Via Spliceosome (GO:0048024)", "Regulation Of mRNA Stability (GO:0043488)", 
      "Protein Insertion Into Mitochondrial Membrane (GO:0051204)", "Microtubule Depolymerization (GO:0007019)")
  colors <- c(HY = "#e74538", CNU = "#9bd5f4", Isocortex = "#18a799", HPF = "#84c551", TH = "#e0abce")
  plots <- lapply(seq_along(terms), function(i) {
      d <- cdata[cdata$pathway == terms[i], ]
      test <- tests[match(terms[i], tests$pathway), ]
      ggplot(d, aes(SC2, TRG_hits)) + geom_smooth(method = "lm", se = TRUE, color = "#343a40", fill = "#ced4da", 
          alpha = 0.3, linewidth = 0.5) + geom_point(aes(fill = Division), shape = 21, size = 1.9, stroke = 0.35, 
          color = "black") + geom_text_repel(aes(label = Region), seed = 260924, size = 1.45, box.padding = 0.08, 
          point.padding = 0.08, force = 3, max.iter = 3000, max.overlaps = 35, segment.size = 0.2, segment.color = "#ced4da") + 
          scale_fill_manual(values = colors) + scale_x_continuous(breaks = c(0, 1, 2), limits = c(-0.1, 
          2.85), expand = expansion(mult = c(0.04, 0.06))) + scale_y_continuous(expand = expansion(mult = c(0.08, 
          0.16))) + labs(title = shorten(terms[i]), subtitle = sprintf("<U+03C1> = %.2f, BH FDR = %.3f", 
          test$rho, test$q_sc2), x = "2nd Structural Connectivity", y = if (i%%4 == 1) 
          "TRG Hits"
      else NULL) + theme_formal() + theme(plot.title = element_text(size = 4.8, face = "bold", hjust = 0.5), 
          plot.subtitle = element_text(size = 4.4, face = "italic", hjust = 0.5), legend.position = "none")
  })
  save_panel(plot_grid(plotlist = plots, nrow = 2, ncol = 4, align = "hv"), "c", 172, 75)
  ddata <- supplied_inputs[["ddata_6"]][[as.character("d")]]
  metrics <- c("ED", "FC", "SC", "SC2", "SPE", "DE", "NE", "CMY")
  heat <- supplied_inputs[["heat_7"]][[as.character(file.path(package, "statistics/extended data figure 12/ed12_all_metric_tests.csv"))]] %>% 
      filter(pathway %in% ddata$pathway) %>% mutate(p = FDR)
  order <- heat %>% filter(Metric == "SC2") %>% arrange(desc(rho), desc(pathway)) %>% pull(pathway)
  heat <- heat %>% mutate(pathway = factor(pathway, levels = rev(order)), Metric = factor(Metric, levels = metrics), 
      Label = paste0(sprintf("%.2f", rho), ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  pd <- ggplot(heat, aes(Metric, pathway, fill = rho)) + geom_tile(color = "white", linewidth = 0.5) + 
      geom_text(aes(label = Label, color = abs(rho) > 0.6), size = 1.6, fontface = "bold") + scale_color_manual(values = c("black", 
      "white"), guide = "none") + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", 
      midpoint = 0, limits = c(-1, 1), name = "Spearman <U+03C1>") + scale_y_discrete(labels = function(x) sub(" \\(GO:[0-9]+\\)", 
      "", x)) + scale_x_discrete(labels = c("ED", "FC", "1st SC", "2nd SC", "SPE", "DE", "NE", "CMY")) + 
      labs(x = NULL, y = NULL) + theme_minimal(base_size = 5, base_family = "Arial") + theme(panel.grid = element_blank(), 
      axis.text.x = element_text(face = "bold"), axis.text.y = element_text(color = "black"), legend.key.height = grid::unit(5, 
          "mm"))
  save_panel(pd, "d", 172, 75)
  cat("Extended Data Figure 12: original display examples and full publication layout\n")
  invisible(grDevices::dev.off())
  invisible(as.list(environment()))
}
