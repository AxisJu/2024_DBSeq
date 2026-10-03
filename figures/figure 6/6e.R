# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(GenomicRanges)
      library(rtracklayer)
      library(tidyverse)
      library(grid)
      library(jsonlite)
      library(stringr)
      library(patchwork)
  })
  cat("=== Step 1: Initializing Environments & Output Paths ===\n")
  fig_out_dir <- output_path("figures/figure 6")
  data_out_dir <- output_path("sourcedata/figure 6")
  if (!dir.exists(fig_out_dir)) dir.create(fig_out_dir, recursive = TRUE)
  if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)
  coords_file <- data_path("results/260716_Nr3c1/configs/gene_coords.json")
  exons_file <- data_path("results/260716_Nr3c1/configs/gene_exons.json")
  motifs_file <- data_path("results/260716_Nr3c1/data/target_genes_exact_motifs_mm10_real.csv")
  ccres_file <- data_path("results/260716_Nr3c1/data/target_genes_encode_ccres_mm10.csv")
  peaks_file <- data_path("results/260716_Nr3c1/data/Nr3c1_AP1_genome_wide_cooccupancy.csv")
  cat("=== Step 2: Loading Genomic Coordinates and UCSC mm10 Exon Models ===\n")
  coords_all <- supplied_inputs[["coords_all_2"]][[as.character(coords_file)]]$mm10
  exons_js <- supplied_inputs[["exons_js_3"]][[as.character(exons_file)]]$mm10
  parse_exons <- function(ex_obj) {
      bind_rows(lapply(ex_obj$exons, as.data.frame)) %>% arrange(start)
  }
  cat("=== Step 3: Loading Exact AP-1 Motifs, Peak Overlap, and ENCODE ccREs ===\n")
  motifs_real <- supplied_inputs[["motifs_real_4"]][[as.character(motifs_file)]]
  peaks_df <- supplied_inputs[["peaks_df_5"]][[as.character(peaks_file)]]
  peaks_df$has_Jun_peak <- as.logical(peaks_df$has_Jun_peak)
  peaks_df$has_Fos_peak <- as.logical(peaks_df$has_Fos_peak)
  ccres_real <- supplied_inputs[["ccres_real_6"]][[as.character(ccres_file)]]
  signal_audit <- list()
  signal_replicates <- list()
  extract_signal_fast <- function(rep_files, chrom, start, end, n_bins = 2500, scale_fac = 5) {
      stopifnot(start >= 1, end >= start, n_bins >= 1, n_bins <= end - start + 1)
      breaks <- floor(seq(start, end + 1, length.out = n_bins + 1))
      bin_starts <- head(breaks, -1)
      bin_ends <- tail(breaks, -1) - 1
      gr_query <- GRanges(chrom, IRanges(start, end))
      gr_bins <- GRanges(chrom, IRanges(bin_starts, bin_ends))
      available <- supplied_inputs[["available_7"]]
      if (!any(available)) 
          stop("No available BigWig replicates for ", chrom, ":", start, "-", end)
      rep_mat <- matrix(NA_real_, nrow = length(rep_files), ncol = n_bins)
      baselines <- rep(NA_real_, length(rep_files))
      for (i in which(available)) {
          s_rep <- supplied_inputs[["s_rep_8"]][[as.character(rep_files[[i]])]]
          if (any(!is.finite(s_rep$score))) 
              stop("Nonfinite BigWig score: ", rep_files[[i]])
          if (length(s_rep) && !isDisjoint(s_rep, ignore.strand = TRUE)) 
              stop("Overlapping BigWig intervals: ", rep_files[[i]])
          raw_val <- numeric(n_bins)
          hits <- findOverlaps(gr_bins, s_rep, ignore.strand = TRUE)
          if (length(hits)) {
              q <- queryHits(hits)
              t <- subjectHits(hits)
              widths <- width(pintersect(gr_bins[q], s_rep[t], ignore.strand = TRUE))
              sums <- rowsum(s_rep$score[t] * widths, q, reorder = FALSE)
              ids <- as.integer(rownames(sums))
              raw_val[ids] <- as.numeric(sums)/width(gr_bins)[ids]
          }
          raw_val <- raw_val * scale_fac
          baselines[i] <- if (any(raw_val > 0)) 
              unname(quantile(raw_val[raw_val > 0], 0.15))
          else 0
          rep_mat[i, ] <- pmax(0, raw_val - baselines[i])
      }
      signal_audit[[length(signal_audit) + 1L]] <<- data.frame(File = rep_files, Chrom = chrom, Start = start, 
          End = end, Included = available, Baseline = baselines, Scale_Factor = scale_fac, Bin_Count = n_bins, 
          Status = ifelse(available, "included", "missing_file_excluded"))
      rep_df <- as.data.frame(t(rep_mat))
      names(rep_df) <- paste0("Replicate_", seq_along(rep_files))
      signal_replicates[[length(signal_replicates) + 1L]] <<- data.frame(Track_ID = length(signal_audit), 
          Chrom = chrom, Bin_Start = bin_starts, Bin_End = bin_ends, Available_Replicates = sum(available), 
          rep_df, Mean_Signal = colMeans(rep_mat[available, , drop = FALSE]))
      list(x_pos = (bin_starts + bin_ends)/2, mean_signal = colMeans(rep_mat[available, , drop = FALSE]), 
          rep_matrix = rep_mat, bin_start = bin_starts, bin_end = bin_ends, n_replicates = sum(available))
  }
  cat("=== Step 5: Processing Nrxn1 and Sf3b1 Genomic Intervals ===\n")
  genes_data <- list()
  source_records <- list()
  for (gname in c("Nrxn1", "Sf3b1")) {
      chr <- coords_all[[gname]]$chrom
      g_start <- coords_all[[gname]]$start
      g_end <- coords_all[[gname]]$end
      strand <- coords_all[[gname]]$strand
      reg_start <- ifelse(gname == "Nrxn1", 89980000, 54980000)
      reg_end <- ifelse(gname == "Nrxn1", 91150000, 55030000)
      t1_files <- c(data_path("results/260716_Nr3c1/data/GSE306261/GSM9195899_Vgat_GR_R1.bw"), data_path("results/260716_Nr3c1/data/GSE306261/GSM9195900_Vgat_GR_R2.bw"), 
          data_path("results/260716_Nr3c1/data/GSE306261/GSM9195902_Vgat_GR_R3.bw"))
      s1 <- extract_signal_fast(t1_files, chr, reg_start, reg_end, scale_fac = 5)
      t2_files <- c(data_path("results/260716_Nr3c1/data/GSE306261/GSM9195903_Vglut1_GR_R1.bw"), data_path("results/260716_Nr3c1/data/GSE306261/GSM9195905_Vglut1_GR_R2.bw"), 
          data_path("results/260716_Nr3c1/data/GSE306261/GSM9195906_Vglut1_GR_R3.bw"))
      s2 <- extract_signal_fast(t2_files, chr, reg_start, reg_end, scale_fac = 5)
      m_blocks <- motifs_real %>% filter(gene == gname & pos >= reg_start & pos <= reg_end) %>% mutate(start = pos - 
          15, end = pos + 15, type = ifelse(str_detect(motif_type, "FOS-JUN"), "FOS-Jun", "Jun-only"), 
          color = ifelse(type == "FOS-Jun", "#4575b4", "#d73027"))
      p_blocks <- peaks_df %>% filter(chrom == chr & start <= reg_end & end >= reg_start) %>% mutate(type = case_when(has_Jun_peak & 
          has_Fos_peak ~ "FOS-Jun", has_Jun_peak & !has_Fos_peak ~ "Jun-only", !has_Jun_peak & has_Fos_peak ~ 
          "FOS-Jun", TRUE ~ "Nr3c1-only"), color = ifelse(type == "FOS-Jun", "#4575b4", "#d73027")) %>% 
          filter(type %in% c("FOS-Jun", "Jun-only"))
      c_blocks <- ccres_real %>% filter(gene == gname & start <= reg_end & end >= reg_start) %>% mutate(color = case_when(ccre_type == 
          "PLS" ~ "#d73027", ccre_type == "pELS" ~ "#fe9929", ccre_type == "dELS" ~ "#fec44f", ccre_type == 
          "CTCF" ~ "#41b6c4", TRUE ~ "#999999"))
      sub_exons <- parse_exons(exons_js[[gname]])
      genes_data[[gname]] <- list(gname = gname, chr = chr, start = reg_start, end = reg_end, strand = strand, 
          s1 = s1, s2 = s2, m_blocks = m_blocks, p_blocks = p_blocks, c_blocks = c_blocks, exons = sub_exons)
      df_sig1 <- data.frame(Gene = gname, Chrom = chr, Position = s1$x_pos, Track_Name = "Nr3c1 CUT&RUN (Vgat+ neurons)", 
          Signal_Mean = s1$mean_signal, ElementType = NA, ElementStart = NA, ElementEnd = NA, ElementClass = NA)
      df_sig2 <- data.frame(Gene = gname, Chrom = chr, Position = s2$x_pos, Track_Name = "Nr3c1 CUT&RUN (Vglut1+ neurons)", 
          Signal_Mean = s2$mean_signal, ElementType = NA, ElementStart = NA, ElementEnd = NA, ElementClass = NA)
      df_m <- data.frame(Gene = gname, Chrom = chr, Position = m_blocks$pos, Track_Name = "Annotation: Motif match", 
          Signal_Mean = NA, ElementType = "Motif", ElementStart = m_blocks$start, ElementEnd = m_blocks$end, 
          ElementClass = m_blocks$type)
      df_p <- data.frame(Gene = gname, Chrom = chr, Position = (p_blocks$start + p_blocks$end)/2, Track_Name = "Annotation: Peak overlap", 
          Signal_Mean = NA, ElementType = "Peak", ElementStart = p_blocks$start, ElementEnd = p_blocks$end, 
          ElementClass = p_blocks$type)
      df_c <- data.frame(Gene = gname, Chrom = chr, Position = (c_blocks$start + c_blocks$end)/2, Track_Name = "Annotation: ENCODE ccREs", 
          Signal_Mean = NA, ElementType = "ccRE", ElementStart = c_blocks$start, ElementEnd = c_blocks$end, 
          ElementClass = c_blocks$ccre_type)
      df_ex <- data.frame(Gene = gname, Chrom = chr, Position = (sub_exons$start + sub_exons$end)/2, Track_Name = "Annotation: Exon", 
          Signal_Mean = NA, ElementType = "Exon", ElementStart = sub_exons$start, ElementEnd = sub_exons$end, 
          ElementClass = "CodingExon")
      source_records[[gname]] <- bind_rows(df_sig1, df_sig2, df_m, df_p, df_c, df_ex)
  }
  sourcedata_fig6e <- bind_rows(source_records)
  write.csv(sourcedata_fig6e, file.path(data_out_dir, "6e.csv"), row.names = FALSE)
  cat(sprintf("-> Saved Source Data to: %s\n", file.path(data_out_dir, "6e.csv")))
  cat("=== Step 6: Assembling Publication Quality Graphic Panels ===\n")
  render_fig6e_panel <- function(g_obj) {
      gname <- g_obj$gname
      xmin <- g_obj$start
      xmax <- g_obj$end
      y_max1 <- ifelse(gname == "Nrxn1", 4, 5)
      y_max2 <- ifelse(gname == "Nrxn1", 8, 8)
      plot_signal_track <- function(s_obj, y_max, fill_color = "#c6dbef", line_color = "#2b7bba") {
          df_mean <- data.frame(x = s_obj$x_pos, y = pmin(y_max, s_obj$mean_signal))
          p <- ggplot()
          if (!is.null(s_obj$rep_matrix) && nrow(s_obj$rep_matrix) > 1) {
              for (r in 1:nrow(s_obj$rep_matrix)) {
                  df_rep <- data.frame(x = s_obj$x_pos, y = pmin(y_max, s_obj$rep_matrix[r, ]))
                  p <- p + geom_line(data = df_rep, aes(x = x, y = y), color = line_color, alpha = 0.35, 
                    linewidth = 0.28)
              }
          }
          p <- p + geom_area(data = df_mean, aes(x = x, y = y), fill = fill_color, alpha = 0.65) + geom_line(data = df_mean, 
              aes(x = x, y = y), color = line_color, linewidth = 0.45) + scale_x_continuous(limits = c(xmin, 
              xmax), expand = c(0, 0)) + scale_y_continuous(limits = c(0, y_max), breaks = c(0, y_max), 
              expand = c(0, 0)) + theme_void() + theme(axis.line.y.left = element_line(color = "black", 
              linewidth = 0.55), axis.text.y.left = element_text(size = 7.5, color = "black", margin = margin(r = 3)), 
              axis.ticks.y.left = element_line(color = "black", linewidth = 0.45), axis.ticks.length.y.left = unit(2, 
                  "pt"), axis.line.x = element_line(color = "black", linewidth = 0.4), plot.margin = margin(t = 2, 
                  r = 6, b = 2, l = 4))
          return(p)
      }
      p1 <- plot_signal_track(g_obj$s1, y_max1)
      p2 <- plot_signal_track(g_obj$s2, y_max2)
      plot_discrete_block_row <- function(blocks_df) {
          p <- ggplot()
          if (nrow(blocks_df) > 0) {
              b_df <- blocks_df %>% mutate(w = pmax(xmax - xmin, 1) * 0.0035, xmin_draw = pmax(xmin, start - 
                  w), xmax_draw = pmin(xmax, end + w))
              p <- p + geom_rect(data = b_df, aes(xmin = xmin_draw, xmax = xmax_draw, ymin = 0.15, ymax = 0.85), 
                  fill = b_df$color, color = NA)
          }
          p <- p + scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 
              1), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(t = 1, r = 6, b = 1, l = 4))
          return(p)
      }
      p_motif <- plot_discrete_block_row(g_obj$m_blocks)
      p_peak <- plot_discrete_block_row(g_obj$p_blocks)
      p_ccre <- plot_discrete_block_row(g_obj$c_blocks)
      plot_gene_model <- function(exons_df, gname, strand) {
          p <- ggplot() + geom_segment(aes(x = min(exons_df$start), xend = max(exons_df$end), y = 0.5, 
              yend = 0.5), color = "black", linewidth = 0.55) + geom_segment(aes(x = min(exons_df$start), 
              xend = min(exons_df$start), y = 0.2, yend = 0.8), color = "black", linewidth = 0.55) + geom_segment(aes(x = max(exons_df$end), 
              xend = max(exons_df$end), y = 0.2, yend = 0.8), color = "black", linewidth = 0.55) + geom_rect(data = exons_df, 
              aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8), fill = "black", color = NA) + annotate("text", 
              x = xmax, y = -0.35, label = gname, fontface = "bold.italic", size = 3.6, hjust = 1) + scale_x_continuous(limits = c(xmin, 
              xmax), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.8, 1.2), expand = c(0, 0)) + 
              theme_void() + theme(plot.margin = margin(t = 2, r = 6, b = 4, l = 4))
          return(p)
      }
      p_gene <- plot_gene_model(g_obj$exons, gname, g_obj$strand)
      return(list(t1 = p1, t2 = p2, motif = p_motif, peak = p_peak, ccre = p_ccre, gene = p_gene))
  }
  p_nrxn1 <- render_fig6e_panel(genes_data$Nrxn1)
  p_sf3b1 <- render_fig6e_panel(genes_data$Sf3b1)
  lbl_nr3c1 <- ggplot() + annotate("segment", x = 6, xend = 6, y = 1, yend = 19, color = "black", linewidth = 0.6) + 
      annotate("text", x = 2.5, y = 10, label = "Nr3c1 CUT&RUN", angle = 90, size = 2.7, fontface = "bold", 
          hjust = 0.5) + annotate("text", x = 9.2, y = 15, label = "Vgat+\nneurons", size = 2.4, hjust = 0.5, 
      lineheight = 0.85) + annotate("text", x = 9.2, y = 5, label = "Vglut1+\nneurons", size = 2.4, hjust = 0.5, 
      lineheight = 0.85) + scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 
      20), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(l = 8, r = 0, t = 0, b = 0))
  lbl_motif <- ggplot() + annotate("text", x = 11.5, y = 5, label = "Motif match", fontface = "bold", size = 2.4, 
      hjust = 1) + scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 
      10), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(l = 8, r = 0, t = 0, b = 0))
  lbl_peak <- ggplot() + annotate("text", x = 11.5, y = 5, label = "Peak overlap", fontface = "bold", size = 2.4, 
      hjust = 1) + scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 
      10), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(l = 8, r = 0, t = 0, b = 0))
  lbl_ccre <- ggplot() + annotate("text", x = 11.5, y = 5, label = "ENCODE ccREs", fontface = "bold", size = 2.4, 
      hjust = 1) + scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 
      10), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(l = 8, r = 0, t = 0, b = 0))
  lbl_gene <- ggplot() + theme_void()
  p_legend <- ggplot() + annotate("text", x = 0.04, y = 0.5, label = "AP-1 Elements:", fontface = "bold", 
      size = 2.8, hjust = 0) + annotate("rect", xmin = 0.16, xmax = 0.175, ymin = 0.35, ymax = 0.65, fill = "#4575b4") + 
      annotate("text", x = 0.18, y = 0.5, label = "FOS-Jun", size = 2.7, hjust = 0) + annotate("rect", 
      xmin = 0.26, xmax = 0.275, ymin = 0.35, ymax = 0.65, fill = "#d73027") + annotate("text", x = 0.28, 
      y = 0.5, label = "Jun-only", size = 2.7, hjust = 0) + annotate("text", x = 0.4, y = 0.5, label = "ENCODE ccREs:", 
      fontface = "bold", size = 2.8, hjust = 0) + annotate("rect", xmin = 0.52, xmax = 0.535, ymin = 0.35, 
      ymax = 0.65, fill = "#d73027") + annotate("text", x = 0.54, y = 0.5, label = "PLS: Promoter", size = 2.7, 
      hjust = 0) + annotate("rect", xmin = 0.66, xmax = 0.675, ymin = 0.35, ymax = 0.65, fill = "#fe9929") + 
      annotate("text", x = 0.68, y = 0.5, label = "pELS: Proximal Enhancer", size = 2.7, hjust = 0) + annotate("rect", 
      xmin = 0.83, xmax = 0.845, ymin = 0.35, ymax = 0.65, fill = "#fec44f") + annotate("text", x = 0.85, 
      y = 0.5, label = "dELS: Distal Enhancer", size = 2.7, hjust = 0) + annotate("rect", xmin = 0.97, 
      xmax = 0.985, ymin = 0.35, ymax = 0.65, fill = "#41b6c4") + annotate("text", x = 0.99, y = 0.5, label = "CTCF-only", 
      size = 2.7, hjust = 0) + scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 
      1), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(t = 2, r = 10, b = 2, l = 10))
  r_nr3c1 <- wrap_plots(lbl_nr3c1, (p_nrxn1$t1/p_nrxn1$t2), (p_sf3b1$t1/p_sf3b1$t2), ncol = 3, widths = c(0.22, 
      1, 1))
  r_motif <- wrap_plots(lbl_motif, p_nrxn1$motif, p_sf3b1$motif, ncol = 3, widths = c(0.22, 1, 1))
  r_peak <- wrap_plots(lbl_peak, p_nrxn1$peak, p_sf3b1$peak, ncol = 3, widths = c(0.22, 1, 1))
  r_ccre <- wrap_plots(lbl_ccre, p_nrxn1$ccre, p_sf3b1$ccre, ncol = 3, widths = c(0.22, 1, 1))
  r_gene <- wrap_plots(lbl_gene, p_nrxn1$gene, p_sf3b1$gene, ncol = 3, widths = c(0.22, 1, 1))
  p_fig6e <- wrap_plots(r_nr3c1, r_motif, r_peak, r_ccre, r_gene, p_legend, ncol = 1, heights = c(2, 0.25, 
      0.25, 0.25, 0.45, 0.35))
  cat("=== Step 7: Exporting High-Resolution PNG ===\n")
  pdf(NULL, width = 12, height = 3.6, useDingbats = FALSE)
  print(p_fig6e)
  dev.off()
  ggsave(file.path(fig_out_dir, "6e.png"), plot = p_fig6e, width = 12, height = 3.6, dpi = 300)
  cat(sprintf("\nFigure 6e complete: %s\n", file.path(fig_out_dir, "6e.png")))
  write.csv(bind_rows(signal_audit), output_path("statistics/figure 6/6e_signal_audit.csv"), row.names = FALSE)
  write.csv(bind_rows(signal_replicates), file.path(data_out_dir, "6e_signal_replicates.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
