# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(ComplexHeatmap)
      library(circlize)
      library(grid)
      library(ggrepel)
      library(patchwork)
  })
  if (Sys.getenv("DBSEQ_REBUILD_SOURCE") == "1" || !supplied_inputs[["data_2"]]) {
      trg_file <- Sys.getenv("DBSEQ_TRG_FILE", data_path("results/Table/CosSim_v2/plot_df.qs"))
      x <- if (grepl("\\.qs$", trg_file)) 
          supplied_inputs[["x_3"]][[as.character(trg_file)]]
      else supplied_inputs[["x_4"]][[as.character(trg_file)]]
      x <- unique(x[x$is_TRG %in% TRUE, c("Celltype", "genes")])
      names(x) <- c("celltype", "gene")
      export_data(x, "ed5a")
  }
  x <- supplied_inputs[["x_5"]][[as.character(file.path(source_dir, "ed5a.csv"))]]
  sets <- split(x$gene, x$celltype)
  cts <- names(sets)
  m <- outer(cts, cts, Vectorize(function(a, b) length(intersect(sets[[a]], sets[[b]]))/length(union(sets[[a]], 
      sets[[b]]))))
  dimnames(m) <- list(cts, cts)
  cl <- ifelse(grepl("^ExN", cts), "ExN", ifelse(grepl("^InN", cts), "InN", "Non-neuronal"))
  h <- Heatmap(m, name = "Jaccard index", col = colorRamp2(c(0, 0.05, 0.15, 0.3, max(m[m < 1], 0.4)), c("#f7fcf0", 
      "#ccebc5", "#7bccc4", "#2b8cbe", "#084081")), left_annotation = rowAnnotation(Class = cl, col = list(Class = c(ExN = "#ba1369", 
      InN = "#56ad00", `Non-neuronal` = "#a58946")), width = unit(3, "mm"), annotation_name_gp = gpar(fontsize = 5, 
      fontface = "bold"), annotation_legend_param = list(title_gp = gpar(fontsize = 5, fontface = "bold"), 
      labels_gp = gpar(fontsize = 5), grid_height = unit(2, "mm"), nrow = 1)), cluster_rows = TRUE, cluster_columns = TRUE, 
      show_column_dend = FALSE, column_names_rot = 90, row_dend_width = unit(5, "mm"), row_dend_gp = gpar(lwd = 0.5), 
      border = TRUE, border_gp = gpar(col = "black", lwd = 0.5), row_names_gp = gpar(fontsize = 5, fontfamily = "Arial"), 
      column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"), heatmap_legend_param = list(title_gp = gpar(fontsize = 5), 
          labels_gp = gpar(fontsize = 5), direction = "horizontal", legend_width = unit(35, "mm")))
  ragg::agg_png(file.path(figure_dir, "ed5a.png"), width = 90, height = 100, units = "mm", res = 600, background = "white")
  draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()
  if (Sys.getenv("DBSEQ_REBUILD_SOURCE") == "1" || !all(supplied_inputs[["data_6"]])) {
      folders <- c(MEMENTO = "results/Table/MEMENTO_v2/Level3", EdgeR = "results/Table/EdgeR_v2/Level3")
      specs <- data.frame(method = c("MEMENTO", "MEMENTO", "EdgeR", "EdgeR"), comparison = c("SS", "DS", 
          "SS", "DS"), prefix = c("MEMENTO_1d_ShamvSaline_", "MEMENTO_1d_DBSvSham_", "EdgeR_ShamIvSalineI_", 
          "EdgeR_DBSIvShamI_"), cutoff = c(0.01, 0.01, 0.05, 0.05))
      available <- lapply(seq_len(nrow(specs)), function(i) {
          files <- supplied_inputs[["files_7"]]
          sub("\\.csv$", "", sub(specs$prefix[i], "", files, fixed = TRUE))
      })
      cts <- sort(Reduce(intersect, available))
      if (!length(cts)) 
          stop("No cell types have all four differential-expression inputs.")
      counts <- list()
      intersection_counts <- list()
      for (ct in cts) {
          genes <- list()
          for (i in seq_len(nrow(specs))) {
              f <- file.path(Sys.getenv("DBSEQ_RESULTS_ROOT", file.path(data_root, "results")), sub("^results/", 
                  "", folders[[specs$method[i]]]), paste0(specs$prefix[i], ct, ".csv"))
              z <- supplied_inputs[["z_8"]][[as.character(f)]]
              gene_values <- if (specs$method[i] == "MEMENTO") 
                  z$gene
              else z[[1]]
              if (length(gene_values) != nrow(z)) 
                  stop(paste("Missing gene identifiers:", f))
              g <- unique(as.character(gene_values[!is.na(z$FDR) & z$FDR < specs$cutoff[i]]))
              genes[[paste(specs$method[i], specs$comparison[i])]] <- g
              counts[[length(counts) + 1]] <- data.frame(celltype = ct, method = specs$method[i], comparison = specs$comparison[i], 
                  deg_count = length(g))
          }
          intersection_counts[[ct]] <- data.frame(celltype = ct, SS = length(intersect(genes[["MEMENTO SS"]], 
              genes[["EdgeR SS"]])), DS = length(intersect(genes[["MEMENTO DS"]], genes[["EdgeR DS"]])))
      }
      obj <- supplied_inputs[["obj_9"]][[as.character(data_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]]
      meta <- obj@meta.data
      region4 <- c(IsoCortex = "Isocortex", HPF = "HPF", STRd = "CNU", STRv = "CNU", PALd = "CNU", PALc = "CNU", 
          LSX = "CNU", CLA = "CNU", MBsen = "Midbrain", DORpm = "TH")
      meta$region <- unname(region4[as.character(meta$brainregion_level4)])
      meta$region[is.na(meta$region)] <- ifelse(meta$brainregion_level3[is.na(meta$region)] %in% c("HY", 
          "TH"), as.character(meta$brainregion_level3[is.na(meta$region)]), "Nonspecific")
      meta$region[!grepl("^(ExN|InN)", meta$celltype_level3)] <- "Non-neuronal"
      region <- meta %>% count(celltype_level3, region) %>% group_by(celltype_level3) %>% arrange(desc(n), 
          region) %>% slice_head(n = 1) %>% ungroup() %>% transmute(celltype = as.character(celltype_level3), 
          region)
      sizes <- meta %>% count(celltype_level3, name = "cell_count") %>% rename(celltype = celltype_level3)
      export_data(left_join(bind_rows(counts), region, by = "celltype"), "ed5b")
      export_data(bind_rows(intersection_counts) %>% left_join(sizes, by = "celltype") %>% left_join(region, 
          by = "celltype"), "ed5c")
  }
  scatter <- function(d, x, y, xlab, ylab, method = "pearson", title = NULL) {
      z <- cor.test(d[[x]], d[[y]], method = method)
      p_text <- if (z$p.value < 0.001) 
          sprintf("P = %.2e", z$p.value)
      else sprintf("P = %.3f", z$p.value)
      label <- sprintf("%s's correlation:\n%s = %.2f, %s", tools::toTitleCase(method), ifelse(method == 
          "pearson", "R", "<U+03C1>"), unname(z$estimate), p_text)
      top <- max(d[[y]]) * ifelse(method == "pearson", 1.38, 1.35)
      ggplot(d, aes(.data[[x]], .data[[y]])) + geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, linetype = "dashed", 
          color = "#aaaaaa", fill = "#dddddd", alpha = 0.25) + geom_point(aes(color = region), size = 1.6) + 
          scale_color_manual(values = colors_brainregion, guide = "none") + geom_text_repel(data = d[d$celltype %in% 
          c("ExN_ATN", "InN_HYa", "ExN_RSP_L45IT", "ExN_PF", "InN_Nonspecific", "VCs"), ], aes(label = celltype), 
          size = 5/.pt, fontface = "bold", family = "Arial", seed = 42, max.overlaps = Inf, box.padding = 0.25, 
          point.padding = 0.2, segment.color = "grey50", segment.size = 0.35, min.segment.length = 0) + 
          annotate("text", x = min(d[[x]]), y = top * 0.98, label = label, hjust = 0, vjust = 1, size = 5/.pt, 
              family = "Arial") + scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) + scale_y_continuous(expand = expansion(mult = c(0.06, 
          0.04))) + coord_cartesian(ylim = c(0, top)) + labs(x = xlab, y = ylab, title = title) + theme_classic(base_family = "Arial", 
          base_size = 5) + theme(axis.title = element_text(size = 5, face = "bold", color = "black"), axis.text = element_text(size = 5, 
          color = "black"), axis.line = element_line(linewidth = 0.5), axis.ticks = element_line(linewidth = 0.5), 
          axis.ticks.length = unit(1, "mm"), plot.title = element_text(size = 6, face = "bold", hjust = 0.5), 
          plot.margin = margin(2, 2, 2, 2, "mm"))
  }
  b <- supplied_inputs[["b_10"]][[as.character(file.path(source_dir, "ed5b.csv"))]]
  ps <- lapply(c("MEMENTO", "EdgeR"), function(method) {
      z <- b[b$method == method, ]
      z <- merge(z[z$comparison == "SS", c("celltype", "region", "deg_count")], z[z$comparison == "DS", 
          c("celltype", "deg_count")], by = "celltype", suffixes = c("_SS", "_DS"))
      title <- if (method == "MEMENTO") 
          "Single-cell Level DEG (MEMENTO)"
      else "Pseudo-bulk Level DEG (EdgeR)"
      scatter(z, "deg_count_SS", "deg_count_DS", "DEG Num. (Sham vs Saline)", "DEG Num. (DBS vs Sham)", 
          "spearman", title)
  })
  save_plot(wrap_plots(ps, nrow = 1), "ed5b", 100, 57)
  c <- supplied_inputs[["c_11"]][[as.character(file.path(source_dir, "ed5c.csv"))]]
  c$cell_count_thousands <- c$cell_count/1000
  save_plot(scatter(c, "DS", "cell_count_thousands", "DEG Num. (DBS vs Sham)", "Cell Number (<U+00D7>1,000)") | 
      scatter(c, "SS", "cell_count_thousands", "DEG Num. (Sham vs Saline)", "Cell Number (<U+00D7>1,000)"), 
      "ed5c", 100, 57)
  invisible(as.list(environment()))
}
