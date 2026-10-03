# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
      library(tidyr)
      library(stringr)
      library(ggplot2)
      library(patchwork)
      library(scales)
      library(RColorBrewer)
      library(ggrastr)
      library(grid)
      library(gtable)
  })
  args <- commandArgs(trailingOnly = FALSE)
  entry <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
  code_root <- normalizePath(file.path(dirname(entry), "../.."), winslash = "/")
  package_root <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(package_root)) stop("Set DBSEQ_OUTPUT_ROOT outside the code package.")
  output_resolved <- normalizePath(package_root, winslash = "/", mustWork = FALSE)
  if (identical(tolower(output_resolved), tolower(code_root)) || startsWith(tolower(output_resolved), paste0(tolower(code_root), 
      "/"))) stop("Outputs must remain outside the public code package")
  workspace_root <- package_root
  data_root <- Sys.getenv("DBSEQ_DATA_ROOT")
  input_root <- Sys.getenv("DBSEQ_INPUT_PACKAGE_ROOT")
  stopifnot(nzchar(data_root), nzchar(workspace_root))
  panel_dir <- function(kind, n) {
      p <- file.path(package_root, kind, paste("extended data figure", n))
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
      p
  }
  save_data <- function(df, n, panel) {
      stat_columns <- intersect(names(df), c("Pval", "FDR", "LogP", "LogFDR", "padj", "pval", "NES", "Sig_Stars"))
      if (length(stat_columns)) 
          write.csv(df, file.path(panel_dir("statistics", n), paste0("ed", n, panel, ".csv")), row.names = FALSE)
      write.csv(df[, setdiff(names(df), stat_columns), drop = FALSE], file.path(panel_dir("sourcedata", 
          n), paste0("ed", n, panel, ".csv")), row.names = FALSE)
  }
  save_png <- function(plot, n, panel, w, h) {
      png(file.path(panel_dir("figures", n), paste0("ed", n, panel, ".png")), width = w, height = h, units = "mm", 
          res = 600, type = "cairo")
      if (inherits(plot, "grob") || inherits(plot, "gtable")) 
          grid::grid.draw(plot)
      else print(plot)
      dev.off()
  }
  set.seed(42)
  save_multidir_pdf_mm <- function(p, filename_base, width_mm, height_mm) {
      keys <- c(FigS4_A_TRG_ORA_Streamlined_Bubble_172mm.pdf = "7a", FigS4_B_TRG_Commonness_4Subplots_172mm.pdf = "7b", 
          FigS4_C1_TRG_Count_Stacked_40x98mm.pdf = "6a", FigS4_E_Level1_TRG_Direction_ORA_Barplots.pdf = "6b", 
          FigS4_D1_Specific_TRG_Count_40x83mm.pdf = "6c", FigS4_D2_Private_TRGs_ORA_Master_Combined_172mm.pdf = "6d")
      if (!(filename_base %in% names(keys))) 
          return(invisible(NULL))
      key <- keys[[filename_base]]
      save_png(p, as.integer(substr(key, 1, 1)), substr(key, 2, 2), width_mm, height_mm)
  }
  pdf(file = NULL)
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  list2env(supplied_inputs[["helpers_2"]], envir = environment())
  cat(">>> [Step 1] Loading master data assets...\n")
  f_data_path <- file.path(input_root, "sourcedata", "figure 2", "2f.csv")
  trg_file <- Sys.getenv("DBSEQ_TRG_FILE")
  if (nzchar(trg_file)) {
      f_df <- if (grepl("\\.qs$", trg_file)) 
          supplied_inputs[["f_df_3"]][[as.character(trg_file)]]
      else supplied_inputs[["f_df_4"]][[as.character(trg_file)]]
      f_df <- dplyr::rename(f_df, Gene = genes, Is_TRG = is_TRG, CoD_DBSvSham = CoD_DBSvSham_cohend, CoD_ShamvSal = CoD_ShamvSal_cohend)
  } else {
      f_df <- supplied_inputs[["f_df_5"]][[as.character(file.path(input_root, "sourcedata/figure 2/2f.csv"))]]
  }
  gmt_file <- file.path(data_root, "bin/MSigDB_Mm/m5.go.bp.v2023.2.Mm.symbols.gmt")
  lines <- supplied_inputs[["lines_6"]][[as.character(gmt_file)]]
  pathways_raw <- list()
  for (l in lines) {
      parts <- strsplit(l, "\t")[[1]]
      if (length(parts) >= 3) {
          pathways_raw[[parts[1]]] <- parts[3:length(parts)]
      }
  }
  bg_genes <- unique(f_df$Gene)
  N_total <- length(bg_genes)
  pathways_filtered <- lapply(pathways_raw, function(g) intersect(g, bg_genes))
  pathways_filtered <- pathways_filtered[sapply(pathways_filtered, length) >= 15 & sapply(pathways_filtered, 
      length) <= 500]
  cat(sprintf(">>> Background genes: %d, Valid GO-BP pathways: %d\n", N_total, length(pathways_filtered)))
  ordered_celltypes <- c("ExN_RSP_L45IT", "ExN_RSP_L6CT", "ExN_DG", "InN_RSP_MGE", "InN_STRns_MGE", "InN_GPe_MGE", 
      "InN_Nonspecific", "InN_STRns_LGE", "InN_CP_D1", "InN_CP_D2", "InN_LSc", "InN_HYa", "InN_RT", "ExN_TRS", 
      "ExN_MH", "ExN_LH", "ExN_THns", "ExN_ATN", "ExN_RE", "ExN_PF", "ExN_CM", "Astro", "OPC", "Oligo", 
      "VCs", "Micro")
  ordered_celltypes <- c(intersect(ordered_celltypes, unique(f_df$Celltype)), sort(setdiff(unique(f_df$Celltype), 
      ordered_celltypes)))
  non_neuronal_types <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  neuronal_types <- setdiff(ordered_celltypes, non_neuronal_types)
  cat(">>> [Part 1] Computing and selecting streamlined pathways (Top 5 Both, Top 5 Neuron-Only, Top 5 Non-Neuron-Only)...\n")
  ora_results_list <- list()
  for (ct in ordered_celltypes) {
      ct_trgs <- f_df %>% filter(Celltype == ct & Is_TRG == TRUE) %>% pull(Gene) %>% unique()
      n_trgs <- length(ct_trgs)
      if (n_trgs < 3) 
          next
      pvals <- numeric(length(pathways_filtered))
      names(pvals) <- names(pathways_filtered)
      overlaps <- integer(length(pathways_filtered))
      names(overlaps) <- names(pathways_filtered)
      for (pname in names(pathways_filtered)) {
          pgenes <- pathways_filtered[[pname]]
          k <- length(intersect(ct_trgs, pgenes))
          overlaps[pname] <- k
          if (k >= 2) {
              M <- length(pgenes)
              pvals[pname] <- phyper(k - 1, M, N_total - M, n_trgs, lower.tail = FALSE)
          }
          else {
              pvals[pname] <- 1
          }
      }
      fdr_vals <- p.adjust(pvals, method = "BH")
      df_ct <- data.frame(Celltype = ct, Pathway = names(pathways_filtered), Overlap = overlaps, Pval = pvals, 
          FDR = fdr_vals, stringsAsFactors = FALSE)
      ora_results_list[[ct]] <- df_ct
  }
  ora_all <- bind_rows(ora_results_list)
  ora_all$Clean_Term <- gsub("^GOBP_", "", ora_all$Pathway) %>% gsub("_", " ", .) %>% str_to_title()
  sig_records <- ora_all %>% filter(FDR < 0.05)
  pw_summary <- sig_records %>% group_by(Pathway, Clean_Term) %>% summarise(n_total_sig = n_distinct(Celltype), 
      n_neuron_sig = sum(Celltype %in% neuronal_types), n_nonneuron_sig = sum(Celltype %in% non_neuronal_types), 
      min_fdr = min(FDR), mean_fdr = mean(FDR), sig_celltypes = list(unique(Celltype)), .groups = "drop")
  non_n_stats <- ora_all %>% filter(Celltype %in% non_neuronal_types) %>% group_by(Pathway, Clean_Term) %>% 
      summarise(min_nonn_fdr = min(FDR), min_nonn_pval = min(Pval), .groups = "drop")
  neuron_stats <- ora_all %>% filter(Celltype %in% neuronal_types) %>% group_by(Pathway, Clean_Term) %>% 
      summarise(min_neuron_fdr = min(FDR), min_neuron_pval = min(Pval), .groups = "drop")
  tier1_both <- pw_summary %>% filter(n_neuron_sig >= 1 & n_nonneuron_sig >= 1) %>% arrange(desc(n_total_sig), 
      mean_fdr) %>% head(5) %>% mutate(Category = "Both Neurons & Non-Neurons")
  tier2_neuron <- pw_summary %>% inner_join(non_n_stats, by = c("Pathway", "Clean_Term")) %>% filter(n_neuron_sig >= 
      1 & n_nonneuron_sig == 0 & min_nonn_fdr >= 0.1) %>% arrange(desc(n_neuron_sig), mean_fdr) %>% head(5) %>% 
      mutate(Category = "Neuron-Exclusive")
  tier3_nonneuron <- pw_summary %>% inner_join(neuron_stats, by = c("Pathway", "Clean_Term")) %>% filter(n_nonneuron_sig >= 
      1 & n_neuron_sig == 0 & min_neuron_fdr >= 0.1) %>% arrange(desc(n_nonneuron_sig), mean_fdr) %>% head(5) %>% 
      mutate(Category = "Non-Neuron-Exclusive")
  cat(sprintf(">>> Selected 15 pathways: Both=%d, Neuron-only=%d, Non-neuron-only=%d\n", nrow(tier1_both), 
      nrow(tier2_neuron), nrow(tier3_nonneuron)))
  selected_pathways_df <- bind_rows(tier1_both, tier2_neuron, tier3_nonneuron)
  streamlined_terms_ordered <- selected_pathways_df$Clean_Term
  df_bubble_streamlined <- ora_all %>% filter(Clean_Term %in% streamlined_terms_ordered) %>% mutate(Clean_Term = factor(Clean_Term, 
      levels = rev(streamlined_terms_ordered)), Celltype = factor(Celltype, levels = ordered_celltypes), 
      LogFDR = -log10(pmax(FDR, 1e-12)), Plot_Size = ifelse(Overlap >= 2 & FDR < 0.05, Overlap, NA)) %>% 
      filter(!is.na(Plot_Size))
  y_div_nonn <- nrow(tier3_nonneuron) + 0.5
  y_div_neur <- y_div_nonn + nrow(tier2_neuron)
  p1_bubble <- ggplot(df_bubble_streamlined, aes(x = Celltype, y = Clean_Term)) + geom_hline(yintercept = c(y_div_nonn, 
      y_div_neur), linetype = "dashed", color = "grey60", linewidth = 0.5/2.8346) + geom_point(aes(size = Plot_Size, 
      fill = LogFDR), shape = 21, color = "black", stroke = 0.4) + scale_fill_gradientn(colors = c("#f7fcf0", 
      "#7fcdbb", "#2c7fb8", "#084081"), name = expression(-log[10](italic(FDR))), limits = c(1.3, 12), 
      oob = squish) + scale_size_continuous(range = c(1, 3.8), name = "TRG Count") + labs(title = "Cell Type Common TRG Enriched Pathways", 
      subtitle = "Shared Core (Top 5), Neuron-Exclusive (Top 5), and Non-Neuron-Exclusive (Top 5)", x = "Cell Types (Level 3)", 
      y = "GO Biological Process Pathways") + theme_s4_v2(base_size = 5) + theme(axis.text.x = element_text(angle = 45, 
      hjust = 1, vjust = 1, face = "bold", size = 5), axis.text.y = element_text(size = 4.8, color = "black", 
      hjust = 1), legend.position = "right", legend.key.height = unit(2.8, "mm"), legend.key.width = unit(2.8, 
      "mm"), panel.grid.major = element_line(color = "grey92", linewidth = 0.25/2.8346, linetype = "dotted"), 
      plot.margin = margin(3, 4, 3, 6, "pt"))
  save_multidir_pdf_mm(p1_bubble, "FigS4_A_TRG_ORA_Streamlined_Bubble_172mm.pdf", width_mm = 172, height_mm = 72)
  cat(">>> [Part 2] Generating TRG Commonness 4-in-a-row figure (172 mm width, vertical Y titles, 300dpi rasterized points)...\n")
  trg_df <- f_df %>% filter(Is_TRG == TRUE)
  neuron_common_map <- trg_df %>% filter(Celltype %in% neuronal_types) %>% group_by(Gene) %>% summarise(commonness = n_distinct(Celltype), 
      .groups = "drop")
  nonneuron_common_map <- trg_df %>% filter(Celltype %in% non_neuronal_types) %>% group_by(Gene) %>% summarise(commonness = n_distinct(Celltype), 
      .groups = "drop")
  trg_neuron_df <- trg_df %>% filter(Celltype %in% neuronal_types) %>% left_join(neuron_common_map, by = "Gene") %>% 
      mutate(Abs_Delta_CohenD = abs(Delta_CohenD), Abs_DBSvSham_CohenD = abs(CoD_DBSvSham)) %>% group_by(Celltype) %>% 
      mutate(CohenD_Rank = rank(-Abs_Delta_CohenD, ties.method = "average")) %>% ungroup()
  trg_nonneuron_df <- trg_df %>% filter(Celltype %in% non_neuronal_types) %>% left_join(nonneuron_common_map, 
      by = "Gene") %>% mutate(Abs_Delta_CohenD = abs(Delta_CohenD), Abs_DBSvSham_CohenD = abs(CoD_DBSvSham)) %>% 
      group_by(Celltype) %>% mutate(CohenD_Rank = rank(-Abs_Delta_CohenD, ties.method = "average")) %>% 
      ungroup()
  p2_1 <- ggplot(trg_neuron_df, aes(x = factor(commonness), y = CohenD_Rank)) + geom_boxplot(fill = "#18a799", 
      color = "black", linewidth = 0.4, outlier.size = 0.2, alpha = 0.7, width = 0.65) + scale_x_discrete(breaks = c("1", 
      "5", "10", "15", "20")) + labs(title = "Neurons: Commonness vs Rank", x = "Commonness (Subtypes)", 
      y = "Cohen's d Rank (1 = Top)") + theme_s4_v2(base_size = 5) + theme(axis.title.y = element_text(angle = 90, 
      size = 5, face = "bold", vjust = 1, color = "black"), axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 4.5))
  r_res_n <- cor.test(trg_neuron_df$commonness, trg_neuron_df$Abs_Delta_CohenD, method = "pearson")
  p_label_n <- if (r_res_n$p.value < 1e-300) "P < 1.00e-300" else sprintf("P = %.2e", r_res_n$p.value)
  p2_2 <- ggplot(trg_neuron_df, aes(x = commonness, y = Abs_Delta_CohenD)) + geom_jitter_rast(color = "grey70", 
      size = 0.25, alpha = 0.35, width = 0.2, raster.dpi = 300) + geom_smooth(method = "lm", color = "#d90429", 
      fill = "#ffccd5", linewidth = 0.5, se = TRUE) + annotate("text", x = 1, y = max(trg_neuron_df$Abs_Delta_CohenD) * 
      0.95, label = sprintf("r = %.3f\n%s", r_res_n$estimate, p_label_n), size = 4.5/2.8346, fontface = "bold", 
      hjust = 0) + labs(title = "Neurons: Correlation", x = "Commonness", y = expression(paste("|", Delta, 
      " Cohen's d|"))) + theme_s4_v2(base_size = 5) + theme(axis.title.y = element_text(angle = 90, size = 5, 
      face = "bold", vjust = 1, color = "black"), axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 4.5))
  p2_3 <- ggplot(trg_nonneuron_df, aes(x = factor(commonness), y = CohenD_Rank)) + geom_boxplot(fill = "#a58946", 
      color = "black", linewidth = 0.4, outlier.size = 0.2, alpha = 0.7, width = 0.6) + labs(title = "Non-Neurons: Commonness vs Rank", 
      x = "Commonness (Types)", y = "Cohen's d Rank (1 = Top)") + theme_s4_v2(base_size = 5) + theme(axis.title.y = element_text(angle = 90, 
      size = 5, face = "bold", vjust = 1, color = "black"), axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 4.5))
  r_res_nn <- cor.test(trg_nonneuron_df$commonness, trg_nonneuron_df$Abs_Delta_CohenD, method = "pearson")
  p_label_nn <- if (r_res_nn$p.value < 1e-300) "P < 1.00e-300" else sprintf("P = %.2e", r_res_nn$p.value)
  p2_4 <- ggplot(trg_nonneuron_df, aes(x = commonness, y = Abs_Delta_CohenD)) + geom_jitter_rast(color = "grey70", 
      size = 0.25, alpha = 0.35, width = 0.15, raster.dpi = 300) + geom_smooth(method = "lm", color = "#d90429", 
      fill = "#ffccd5", linewidth = 0.5, se = TRUE) + annotate("text", x = 1, y = max(trg_nonneuron_df$Abs_Delta_CohenD) * 
      0.95, label = sprintf("r = %.3f\n%s", r_res_nn$estimate, p_label_nn), size = 4.5/2.8346, fontface = "bold", 
      hjust = 0) + labs(title = "Non-Neurons: Correlation", x = "Commonness", y = expression(paste("|", 
      Delta, " Cohen's d|"))) + theme_s4_v2(base_size = 5) + theme(axis.title.y = element_text(angle = 90, 
      size = 5, face = "bold", vjust = 1, color = "black"), axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 4.5), 
      plot.margin = margin(2, 5, 2, 2, "pt"))
  p2_combined <- (p2_1 | p2_2 | p2_3 | p2_4) + plot_annotation(title = "TRG Commonness, Differential Expression Rank and Effect Size", 
      theme = theme(plot.title = element_text(size = 5.5, face = "bold", hjust = 0.5)))
  save_multidir_pdf_mm(p2_combined, "FigS4_B_TRG_Commonness_4Subplots_172mm.pdf", width_mm = 172, height_mm = 43)
  cat(">>> [Part 3.1] Generating Stacked Horizontal Barplot of TRG Counts (40x98 mm, right-aligned Y)...\n")
  trg_counts_df <- trg_df %>% mutate(Direction = ifelse(CoD_DBSvSham < 0, "Up in DBS", "Down in DBS")) %>% 
      group_by(Celltype, Direction) %>% summarise(Count = n_distinct(Gene), .groups = "drop")
  total_counts <- trg_counts_df %>% group_by(Celltype) %>% summarise(Total = sum(Count), .groups = "drop") %>% 
      arrange(Total)
  trg_counts_df <- trg_counts_df %>% mutate(Celltype = factor(Celltype, levels = total_counts$Celltype), 
      Direction = factor(Direction, levels = c("Down in DBS", "Up in DBS")))
  p3_stacked_bars <- ggplot(trg_counts_df, aes(x = Count, y = Celltype, fill = Direction)) + geom_col(position = "stack", 
      color = "black", linewidth = 0.4, width = 0.72) + scale_fill_manual(values = c(`Up in DBS` = "#d90429", 
      `Down in DBS` = "#3182bd"), name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.08))) + 
      labs(title = "TRG Counts per Cell Type", x = "Number of TRGs", y = NULL) + theme_s4_v2(base_size = 5) + 
      theme(axis.text.y = element_text(size = 4.8, face = "bold", hjust = 1), axis.text.x = element_text(size = 4.8), 
          legend.position = "top", legend.direction = "horizontal", legend.margin = margin(b = -2), plot.margin = margin(2, 
              4, 2, 2, "pt"))
  save_multidir_pdf_mm(p3_stacked_bars, "FigS4_C1_TRG_Count_Stacked_40x98mm.pdf", width_mm = 40, height_mm = 98)
  cat(">>> [Part 3.2] Computing and plotting Level 1 Directional GOBP ORA...\n")
  lineage_rows <- list()
  lineage_ora <- list()
  lvl1_cts <- c("ExN", "InN", "Astro", "OPC", "Oligo", "Micro", "VCs")
  lvl1_base_colors <- c(ExN = "#ba1369", InN = "#56ad00", Astro = "#a58946", OPC = "#5953ff", Oligo = "#201e5a", 
      Micro = "#a87c5a", VCs = "#858881")
  compute_lvl1_data <- function(ct_name) {
      csv_file <- file.path(Sys.getenv("DBSEQ_RESULTS_ROOT", file.path(data_root, "results")), "Table/CosSim_v2/Level1", 
          paste0(ct_name, "_TRG.csv"))
      if (!supplied_inputs[["compute_lvl1_data_7"]]) 
          stop("Missing lineage data: ", csv_file)
      df_lvl1 <- supplied_inputs[["df_lvl1_8"]][[as.character(csv_file)]]
      lineage_rows[[ct_name]] <<- data.frame(Celltype = ct_name, Gene = df_lvl1$genes, log2_fold_change = df_lvl1$PB_DBSvSham_logFC)
      g_up <- unique(df_lvl1$genes[df_lvl1$PB_DBSvSham_logFC > 0])
      g_dn <- unique(df_lvl1$genes[df_lvl1$PB_DBSvSham_logFC < 0])
      run_hyper <- function(gset) {
          gset <- intersect(gset, bg_genes)
          n_g <- length(gset)
          if (n_g < 2) 
              return(data.frame())
          pvals <- numeric(length(pathways_filtered))
          names(pvals) <- names(pathways_filtered)
          overlaps <- integer(length(pathways_filtered))
          names(overlaps) <- names(pathways_filtered)
          for (pname in names(pathways_filtered)) {
              pgenes <- pathways_filtered[[pname]]
              k <- length(intersect(gset, pgenes))
              overlaps[pname] <- k
              if (k >= 2) {
                  M <- length(pgenes)
                  pvals[pname] <- phyper(k - 1, M, N_total - M, n_g, lower.tail = FALSE)
              }
              else {
                  pvals[pname] <- 1
              }
          }
          df_res <- data.frame(Pathway = names(pathways_filtered), Overlap = overlaps, Pval = pvals, stringsAsFactors = FALSE) %>% 
              mutate(FDR = p.adjust(Pval, method = "BH")) %>% filter(Overlap >= 2 & FDR < 0.05) %>% arrange(FDR) %>% 
              head(5)
          if (nrow(df_res) == 0) 
              return(data.frame())
          raw_c <- gsub("^GOBP_", "", df_res$Pathway) %>% gsub("_", " ", .) %>% str_to_title()
          df_res$Clean_Term <- make.unique(raw_c)
          df_res$LogP <- -log10(pmax(df_res$FDR, 1e-25))
          return(df_res)
      }
      res_up <- run_hyper(g_up)
      res_dn <- run_hyper(g_dn)
      lineage_ora[[ct_name]] <<- bind_rows(list(`Up in DBS` = res_up, `Down in DBS` = res_dn), .id = "Direction")
      return(list(Up = res_up, Down = res_dn))
  }
  lvl1_data_all <- list()
  for (ct in lvl1_cts) {
      d <- compute_lvl1_data(ct)
      if (!is.null(d)) 
          lvl1_data_all[[ct]] <- d
  }
  lvl1_plots_all <- list()
  for (ct in names(lvl1_data_all)) {
      res_up <- lvl1_data_all[[ct]]$Up
      res_dn <- lvl1_data_all[[ct]]$Down
      if (nrow(res_up) == 0 && nrow(res_dn) == 0) 
          next
      max_count <- max(c(res_up$Overlap, res_dn$Overlap, 4), na.rm = TRUE)
      max_x <- ceiling(max_count * 1.15)
      all_logp <- c(res_up$LogP, res_dn$LogP)
      min_logp <- min(all_logp, na.rm = TRUE)
      max_logp <- max(all_logp, na.rm = TRUE)
      base_col <- lvl1_base_colors[[ct]]
      tint_col <- colorRampPalette(c("#ffffff", base_col))(10)[3]
      build_lvl1_slot_panel <- function(sub_res, title_text, is_up = TRUE) {
          n_actual <- nrow(sub_res)
          slots_df <- data.frame(Slot = 1:5, Clean_Term = paste0("dummy_", 1:5), Overlap = 0, LogP = min_logp, 
              Is_Real = FALSE, stringsAsFactors = FALSE)
          if (n_actual > 0) {
              sub_res <- sub_res %>% arrange(Overlap)
              for (i in 1:n_actual) {
                  idx <- 5 - n_actual + i
                  slots_df$Clean_Term[idx] <- sub_res$Clean_Term[i]
                  slots_df$Overlap[idx] <- sub_res$Overlap[i]
                  slots_df$LogP[idx] <- sub_res$LogP[i]
                  slots_df$Is_Real[idx] <- TRUE
              }
          }
          slots_df$Term_Factor <- factor(slots_df$Clean_Term, levels = slots_df$Clean_Term)
          label_map <- setNames(ifelse(slots_df$Is_Real, slots_df$Clean_Term, ""), slots_df$Clean_Term)
          p <- ggplot(slots_df, aes(x = Overlap, y = Term_Factor, fill = LogP)) + geom_col(data = filter(slots_df, 
              Is_Real == TRUE), color = "black", linewidth = 0.4, width = 0.72) + scale_fill_gradient(low = tint_col, 
              high = base_col, limits = c(min_logp, max_logp), name = expression(-log[10](FDR))) + scale_y_discrete(drop = FALSE, 
              labels = label_map) + scale_x_continuous(limits = c(0, max_x), expand = expansion(mult = c(0, 
              0.05))) + labs(title = title_text, x = if (is_up) 
              NULL
          else "TRG Counts", y = NULL) + theme_s4_v2(base_size = 5) + theme(plot.title = element_text(size = 4.8, 
              face = "bold", hjust = 0), axis.text.y = element_text(size = 4.2, color = "black", hjust = 1), 
              axis.text.x = if (is_up) 
                  element_blank()
              else element_text(size = 4.5), axis.ticks.x = if (is_up) 
                  element_blank()
              else element_line(color = "black", linewidth = 0.5), axis.title.x = if (is_up) 
                  element_blank()
              else element_text(size = 4.8, face = "bold"), legend.position = "none", plot.margin = margin(t = if (is_up) 
                  2
              else 0, r = 4, b = if (is_up) 
                  0
              else 2, l = 2, "pt"))
          return(p)
      }
      p_up <- build_lvl1_slot_panel(res_up, sprintf("%s: Up in DBS", ct), is_up = TRUE)
      p_dn <- build_lvl1_slot_panel(res_dn, sprintf("%s: Down in DBS", ct), is_up = FALSE)
      p_stacked <- (p_up/p_dn) + plot_layout(heights = c(1, 1))
      g_up_info <- set_panel_size_custom(p_up, width_mm = 18, height_mm = 18)
      g_dn_info <- set_panel_size_custom(p_dn, width_mm = 18, height_mm = 18)
      g_combined <- rbind(g_up_info$grob, g_dn_info$grob, size = "max")
      total_w_mm <- as.numeric(convertWidth(sum(g_combined$widths), "mm"))
      total_h_mm <- as.numeric(convertHeight(sum(g_combined$heights), "mm"))
      out_name <- sprintf("FigS4_C2_Level1_Direction_ORA_%s.pdf", ct)
      save_multidir_pdf_mm(g_combined, out_name, width_mm = total_w_mm, height_mm = total_h_mm)
      lvl1_plots_all[[ct]] <- p_stacked
  }
  if (length(lvl1_plots_all) > 0) {
      p_lvl1_master <- build_lineage_master(lineage_ora, lvl1_base_colors) + plot_annotation(title = "ORA of TRGs per Cell Type", 
          theme = theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5)))
      save_multidir_pdf_mm(p_lvl1_master, "FigS4_E_Level1_TRG_Direction_ORA_Barplots.pdf", width_mm = 172, 
          height_mm = 140)
  }
  cat(">>> [Part 4.1] Generating Cell-Type Specific TRG Count Barplot (Proportional Height to C1: 40x83 mm)...\n")
  gene_recurrence_all <- trg_df %>% group_by(Gene) %>% summarise(n_celltypes = n_distinct(Celltype), .groups = "drop")
  global_specific_genes <- gene_recurrence_all %>% filter(n_celltypes == 1) %>% pull(Gene)
  cat(sprintf(">>> Total global cell-type specific TRGs: %d\n", length(global_specific_genes)))
  spec_trgs_detail <- trg_df %>% filter(Gene %in% global_specific_genes)
  spec_counts_per_ct <- spec_trgs_detail %>% group_by(Celltype) %>% summarise(Specific_Count = n_distinct(Gene), 
      .groups = "drop")
  df_spec_nonzero <- spec_counts_per_ct %>% filter(Specific_Count > 0) %>% arrange(Specific_Count)
  df_spec_nonzero$Celltype_Factor <- factor(df_spec_nonzero$Celltype, levels = df_spec_nonzero$Celltype)
  p4_specific_counts <- ggplot(df_spec_nonzero, aes(x = Specific_Count, y = Celltype_Factor, fill = Celltype)) + 
      geom_col(color = "black", linewidth = 0.4, width = 0.72, show.legend = FALSE) + geom_text(aes(label = as.character(Specific_Count)), 
      hjust = -0.15, size = 4.2/2.8346, fontface = "bold", color = "black") + scale_fill_manual(values = colors_celltype_level3) + 
      scale_x_continuous(expand = expansion(mult = c(0, 0.22))) + labs(title = "Private TRGs per Cell Type", 
      x = "Specific TRG Count", y = NULL) + theme_s4_v2(base_size = 5) + theme(axis.text.y = element_text(size = 4.8, 
      face = "bold", hjust = 1), axis.text.x = element_text(size = 4.8), plot.margin = margin(2, 6, 2, 
      2, "pt"))
  height_d1_mm <- 98 * (nrow(df_spec_nonzero)/26)
  cat(sprintf(">>> FigS4_D1 exact proportional height: %.2f mm\n", height_d1_mm))
  save_multidir_pdf_mm(p4_specific_counts, "FigS4_D1_Specific_TRG_Count_40x83mm.pdf", width_mm = 40, height_mm = height_d1_mm)
  save_multidir_pdf_mm(p4_specific_counts, "FigS4_D1_Specific_TRG_Count_40x63mm.pdf", width_mm = 40, height_mm = height_d1_mm)
  cat(">>> [Part 4.2] Generating Master Combined Private TRGs ORA Figure (Top-to-Bottom, Left-to-Right, Count >= 5)...\n")
  cts_ordered_d1 <- rev(levels(df_spec_nonzero$Celltype_Factor))
  ora_top5_count5_list <- list()
  for (ct in cts_ordered_d1) {
      ct_private_genes <- spec_trgs_detail %>% filter(Celltype == ct) %>% pull(Gene) %>% unique()
      n_g <- length(ct_private_genes)
      if (n_g < 5) 
          next
      pvals <- numeric(length(pathways_filtered))
      names(pvals) <- names(pathways_filtered)
      overlaps <- integer(length(pathways_filtered))
      names(overlaps) <- names(pathways_filtered)
      for (pname in names(pathways_filtered)) {
          pgenes <- pathways_filtered[[pname]]
          k <- length(intersect(ct_private_genes, pgenes))
          overlaps[pname] <- k
          if (k >= 5) {
              M <- length(pgenes)
              pvals[pname] <- phyper(k - 1, M, N_total - M, n_g, lower.tail = FALSE)
          }
          else {
              pvals[pname] <- 1
          }
      }
      df_private_ora <- data.frame(Pathway = names(pathways_filtered), Overlap = overlaps, Pval = pvals, 
          stringsAsFactors = FALSE) %>% mutate(FDR = p.adjust(Pval, method = "BH")) %>% filter(Overlap >= 
          5 & FDR < 0.05) %>% arrange(FDR)
      if (nrow(df_private_ora) > 0) {
          raw_c <- gsub("^GOBP_", "", df_private_ora$Pathway) %>% gsub("_", " ", .) %>% str_to_title()
          df_private_ora$Clean_Term <- make.unique(raw_c)
          df_private_ora$LogP <- -log10(pmax(df_private_ora$FDR, 1e-25))
          ora_top5_count5_list[[ct]] <- head(df_private_ora, 5)
      }
  }
  cat(sprintf(">>> Total Cell Types with high-confidence pathways (Count >= 5): %d\n", length(ora_top5_count5_list)))
  d2_subpanels_list <- list()
  display_private_types <- c("ExN_ATN", "Micro", "Oligo", "InN_HYa", "InN_CP_D1", "InN_RT", "ExN_RSP_L45IT", 
      "InN_CP_D2", "Astro")
  for (ct in display_private_types) {
      sub_res <- ora_top5_count5_list[[ct]]
      if (is.null(sub_res) || nrow(sub_res) == 0) {
          d2_subpanels_list[[ct]] <- ggplot() + theme_void() + annotate("text", x = 0, y = 0, label = paste(ct, 
              "No pathway at BH FDR < 0.05", sep = "\n"), size = 2)
          next
      }
      n_actual <- nrow(sub_res)
      n_priv_trg <- df_spec_nonzero$Specific_Count[df_spec_nonzero$Celltype == ct]
      ct_col <- if (ct %in% names(colors_celltype_level3)) 
          colors_celltype_level3[[ct]]
      else "#18a799"
      tint_col <- colorRampPalette(c("#ffffff", ct_col))(10)[3]
      min_logp <- min(sub_res$LogP)
      max_logp <- max(sub_res$LogP)
      if (min_logp == max_logp) 
          min_logp <- max_logp * 0.5
      max_cnt <- max(sub_res$Overlap, 5)
      max_x <- ceiling(max_cnt * 1.15)
      slots_df <- data.frame(Slot = 1:5, Clean_Term = paste0("dummy_", 1:5), Overlap = 0, LogP = min_logp, 
          Is_Real = FALSE, stringsAsFactors = FALSE)
      sub_res <- sub_res %>% arrange(Overlap)
      for (i in 1:n_actual) {
          idx <- 5 - n_actual + i
          slots_df$Clean_Term[idx] <- sub_res$Clean_Term[i]
          slots_df$Overlap[idx] <- sub_res$Overlap[i]
          slots_df$LogP[idx] <- sub_res$LogP[i]
          slots_df$Is_Real[idx] <- TRUE
      }
      slots_df$Term_Factor <- factor(slots_df$Clean_Term, levels = slots_df$Clean_Term)
      label_map <- setNames(ifelse(slots_df$Is_Real, slots_df$Clean_Term, ""), slots_df$Clean_Term)
      p_sub <- ggplot(slots_df, aes(x = Overlap, y = Term_Factor, fill = LogP)) + geom_col(data = filter(slots_df, 
          Is_Real == TRUE), color = "black", linewidth = 0.4, width = 0.72) + scale_fill_gradient(low = tint_col, 
          high = ct_col) + scale_y_discrete(drop = FALSE, labels = label_map) + scale_x_continuous(limits = c(0, 
          max_x), expand = expansion(mult = c(0, 0.05))) + labs(title = sprintf("%s (%d)", ct, n_priv_trg), 
          x = "TRG Counts") + theme_s4_v2(base_size = 5) + theme(plot.title = element_text(size = 4.8, 
          face = "bold", hjust = 0), axis.title.y = element_blank(), axis.text.y = element_text(size = 4.2, 
          color = "black", hjust = 1), axis.text.x = element_text(size = 4.2), axis.title.x = element_text(size = 4.5, 
          face = "bold"), legend.position = "none", plot.margin = margin(2, 4, 2, 2, "pt"))
      d2_subpanels_list[[ct]] <- p_sub
      panel_h_mm <- n_actual * 3.6
      p_indiv_info <- set_panel_size_custom(p_sub, width_mm = 18, height_mm = panel_h_mm)
      out_single <- sprintf("FigS4_D2_Private_TRGs_ORA_%s.pdf", ct)
      save_multidir_pdf_mm(p_indiv_info$grob, out_single, width_mm = p_indiv_info$width_mm, height_mm = p_indiv_info$height_mm)
  }
  p_master_combined_180 <- wrap_plots(d2_subpanels_list, ncol = 3, byrow = FALSE) + plot_annotation(title = "ORA of Private TRGs per Cell Type", 
      subtitle = "Top 5 GOBP pathways for cell types with high-confidence private therapeutic programs (Overlap >= 5, FDR < 0.05)", 
      theme = theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 4.8, 
          hjust = 0.5)))
  save_multidir_pdf_mm(p_master_combined_180, "FigS4_D2_Private_TRGs_ORA_Master_Combined.pdf", width_mm = 180, 
      height_mm = 125)
  save_multidir_pdf_mm(p_master_combined_180, "FigS4_D2_Private_TRGs_ORA_Master_Combined_172mm.pdf", width_mm = 172, 
      height_mm = 125)
  cat(">>> =====================================================================\n")
  cat(">>> [SUCCESS] All Figure S4 V2 Panels Reconstructed and Exported Successfully!\n")
  cat(">>> =====================================================================\n")
  save_data(trg_df[, c("Gene", "Celltype", "CoD_DBSvSham")], 6, "a")
  save_data(bind_rows(lineage_rows), 6, "b_genes")
  save_data(bind_rows(lineage_ora, .id = "Celltype"), 6, "b")
  save_data(spec_trgs_detail[, c("Gene", "Celltype")], 6, "c")
  private_plot <- bind_rows(ora_top5_count5_list[display_private_types], .id = "Celltype")
  private_plot$Gene_IDs <- vapply(seq_len(nrow(private_plot)), function(i) paste(intersect(spec_trgs_detail$Gene[spec_trgs_detail$Celltype == 
      private_plot$Celltype[i]], pathways_filtered[[private_plot$Pathway[i]]]), collapse = ";"), character(1))
  save_data(private_plot, 6, "d")
  bubble_source <- df_bubble_streamlined[, c("Celltype", "Pathway", "Overlap", "FDR")]
  bubble_source$Gene_IDs <- vapply(seq_len(nrow(bubble_source)), function(i) paste(intersect(trg_df$Gene[trg_df$Celltype == 
      bubble_source$Celltype[i]], pathways_filtered[[bubble_source$Pathway[i]]]), collapse = ";"), character(1))
  save_data(bubble_source, 7, "a")
  rank_source <- bind_rows(list(Neuron = trg_neuron_df, NonNeuron = trg_nonneuron_df), .id = "Lineage")
  save_data(rank_source[, c("Gene", "Celltype", "Lineage", "commonness", "Abs_Delta_CohenD", "CohenD_Rank")], 
      7, "b")
  invisible(as.list(environment()))
}
