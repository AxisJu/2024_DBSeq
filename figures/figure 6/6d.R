# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(tidyverse)
      library(tidygraph)
      library(ggraph)
      library(ggnewscale)
      library(qs)
  })
  dir.create(output_path("figures/figure 6"), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_path("sourcedata/figure 6"), recursive = TRUE, showWarnings = FALSE)
  CACHE_RDS <- data_path("results/260616 pySCENIC/auc_with_meta.rds")
  TF_DB_PATH <- data_path("bin/Mouse_TFs_Kinases_webpage-3-30-2017.csv")
  TRG_TABLE <- data_path("results/Table/CosSim_v2/plot_df.qs")
  REG_CSV <- data_path("data/snRNAseq_mouse/processed/intermediate/pySCENIC/reg.csv")
  CRITERION4 <- data_path("results/260616 pySCENIC/TRG_TF_screening/criterion4_targets.csv")
  colVars_fast <- function(x) {
      n <- nrow(x)
      if (n <= 1) 
          return(rep(0, ncol(x)))
      (colSums(x^2) - (colSums(x)^2)/n)/(n - 1)
  }
  clean_tf_name <- function(x) trimws(gsub("\\s*\\(.*?\\)", "", x))
  cat("Loading cache...\n")
  list2env(supplied_inputs[["helpers_2"]], envir = environment())
  data_obj <- aggregate_donor_auc(supplied_inputs[["data_obj_3"]][[as.character(CACHE_RDS)]])
  global_auc_mat <- data_obj$auc_mat
  global_meta <- data_obj$meta
  tf_db <- supplied_inputs[["tf_db_4"]][[as.character(TF_DB_PATH)]]
  colnames(tf_db)[1] <- "Symbol"
  pure_tfs <- unique(tf_db$Symbol[tf_db$Class == "TF"])
  current_regulons <- colnames(global_auc_mat)
  cleaned_symbols <- clean_tf_name(current_regulons)
  valid_tf_idx <- cleaned_symbols %in% pure_tfs
  global_auc_mat <- global_auc_mat[, valid_tf_idx, drop = FALSE]
  atn_cell_idx <- which(global_meta$final_region == "ATN")
  atn_auc_mat <- global_auc_mat[atn_cell_idx, , drop = FALSE]
  atn_meta <- global_meta[atn_cell_idx, ]
  n_saline <- sum(atn_meta$Group_L2 == "Saline_I")
  n_sham <- sum(atn_meta$Group_L2 == "Sham_I")
  n_dbs <- sum(atn_meta$Group_L2 == "DBS_I")
  m_saline <- colMeans(atn_auc_mat[atn_meta$Group_L2 == "Saline_I", , drop = FALSE])
  m_sham <- colMeans(atn_auc_mat[atn_meta$Group_L2 == "Sham_I", , drop = FALSE])
  m_dbs <- colMeans(atn_auc_mat[atn_meta$Group_L2 == "DBS_I", , drop = FALSE])
  v_saline <- colVars_fast(atn_auc_mat[atn_meta$Group_L2 == "Saline_I", , drop = FALSE])
  v_sham <- colVars_fast(atn_auc_mat[atn_meta$Group_L2 == "Sham_I", , drop = FALSE])
  v_dbs <- colVars_fast(atn_auc_mat[atn_meta$Group_L2 == "DBS_I", , drop = FALSE])
  eps <- 1e-10
  pooled_sd_dbs_sham <- sqrt(((n_dbs - 1) * v_dbs + (n_sham - 1) * v_sham)/(n_dbs + n_sham - 2))
  pooled_sd_dbs_sham[pooled_sd_dbs_sham == 0] <- eps
  pooled_sd_sal_sham <- sqrt(((n_saline - 1) * v_saline + (n_sham - 1) * v_sham)/(n_saline + n_sham - 2))
  pooled_sd_sal_sham[pooled_sd_sal_sham == 0] <- eps
  p_vals_dbs_sham <- numeric(ncol(atn_auc_mat))
  names(p_vals_dbs_sham) <- colnames(atn_auc_mat)
  idx_dbs <- which(atn_meta$Group_L2 == "DBS_I")
  idx_sham <- which(atn_meta$Group_L2 == "Sham_I")
  for (reg in colnames(atn_auc_mat)) {
      wt <- suppressWarnings(wilcox.test(atn_auc_mat[idx_dbs, reg], atn_auc_mat[idx_sham, reg], exact = FALSE))
      p_vals_dbs_sham[reg] <- wt$p.value
  }
  filtered_tfs_df <- tibble(Regulon = colnames(atn_auc_mat), TF = clean_tf_name(Regulon), DBSvSham_cohend = (m_dbs - 
      m_sham)/pooled_sd_dbs_sham, SalinevSham_cohend = (m_saline - m_sham)/pooled_sd_sal_sham, DBSvSham_p_donor = p_vals_dbs_sham) %>% 
      mutate(DBSvSham_FDR = p.adjust(DBSvSham_p_donor, method = "BH")) %>% filter(!is.na(DBSvSham_cohend) & 
      !is.na(SalinevSham_cohend) & !is.infinite(DBSvSham_cohend) & !is.infinite(SalinevSham_cohend)) %>% 
      filter(sign(DBSvSham_cohend) == sign(SalinevSham_cohend) & DBSvSham_FDR < 0.05) %>% arrange(desc(DBSvSham_cohend))
  write.csv(filtered_tfs_df, output_path("statistics/figure 6/6d_selected_TFs.csv"), row.names = FALSE)
  if (nrow(filtered_tfs_df) == 0) {
      blank <- ggplot() + annotate("text", x = 0, y = 0, label = "No TF passes donor-level BH FDR < 0.05", 
          size = 4) + theme_void()
      ggsave(output_path("figures/figure 6/6d.png"), blank, width = 7, height = 4, dpi = 300, bg = "white")
      write.csv(filtered_tfs_df[, c("TF", "DBSvSham_cohend", "SalinevSham_cohend")], output_path("sourcedata/figure 6/6d.csv"), 
          row.names = FALSE)
      return(invisible(NULL))
  }
  all_valid_tfs_inherited <- unique(filtered_tfs_df$TF)
  cat("Part 2...\n")
  trg_data <- supplied_inputs[["trg_data_5"]][[as.character(TRG_TABLE)]]
  stopifnot(all(c("genes", "is_TRG") %in% names(trg_data)))
  global_trgs <- unique(trg_data$genes[!is.na(trg_data$is_TRG) & trg_data$is_TRG])
  splicing_genes <- if (supplied_inputs[["splicing_genes_6"]]) {
      supplied_inputs[["splicing_genes_7"]][[as.character(CRITERION4)]]$trg_names_splicing %>% str_split(";") %>% 
          unlist() %>% unique() %>% setdiff("")
  } else c()
  cat("Part 3...\n")
  raw_csv <- supplied_inputs[["raw_csv_8"]][[as.character(REG_CSV)]]
  header_row_idx <- which(raw_csv$X1 == "TF" | raw_csv$X2 == "TF")[1]
  data_rows <- raw_csv[(header_row_idx + 1):nrow(raw_csv), ]
  target_col_idx <- which(apply(data_rows, 2, function(col) any(str_detect(col, "^\\s*\\[\\("), na.rm = TRUE)))[1]
  if (is.na(target_col_idx)) target_col_idx <- ncol(data_rows)
  reg_targets <- tibble(TF = data_rows[[1]], TargetGenes = data_rows[[target_col_idx]]) %>% mutate(TF = trimws(gsub("\\s*\\(.*?\\)", 
      "", TF))) %>% filter(TF %in% all_valid_tfs_inherited)
  network_edges_raw <- map_dfr(1:nrow(reg_targets), function(r) {
      tf_name <- reg_targets$TF[r]
      matches <- str_match_all(reg_targets$TargetGenes[r], "'([A-Za-z0-9.-]+)',\\s*([0-9.-]+)")[[1]]
      if (nrow(matches) > 0) {
          tibble(from = tf_name, to = matches[, 2], weight = as.numeric(matches[, 3]))
      }
      else NULL
  }) %>% distinct()
  CONFIDENCE_QUANTILE <- 0.75
  edge_audit <- network_edges_raw %>% group_by(from) %>% mutate(Threshold = quantile(weight, CONFIDENCE_QUANTILE, 
      na.rm = TRUE), Included = is.finite(weight) & weight >= Threshold) %>% ungroup()
  network_edges_filtered <- edge_audit %>% filter(Included) %>% select(from, to, weight)
  write_csv(edge_audit, output_path("statistics/figure 6/6d_edge_audit.csv"))
  write_csv(filtered_tfs_df %>% mutate(Selection = "ATN DBS-Sham donor-level BH FDR<0.05 and concordant DBS-Sham/Saline-Sham effect directions"), 
      output_path("statistics/figure 6/6d_TF_selection.csv"))
  edges_main <- network_edges_filtered
  tfs_main <- all_valid_tfs_inherited
  cat("Part 4...\n")
  tf_atn_d_vec <- filtered_tfs_df$DBSvSham_cohend
  names(tf_atn_d_vec) <- filtered_tfs_df$TF
  functional_target_pool <- c("Nlgn1", "Nrxn1", splicing_genes, global_trgs)
  build_tidy_graph <- function(edge_df, valid_tfs) {
      if (nrow(edge_df) == 0) 
          return(NULL)
      nodes_vec <- unique(c(edge_df$from, edge_df$to))
      target_nes_df <- edge_df %>% dplyr::group_by(gene = to) %>% dplyr::summarise(NES = max(weight, na.rm = TRUE), 
          .groups = "drop")
      meta_df <- tibble(gene = nodes_vec) %>% dplyr::left_join(target_nes_df, by = "gene") %>% dplyr::mutate(is_tf = gene %in% 
          valid_tfs, CohenD = dplyr::case_when(is_tf ~ coalesce(tf_atn_d_vec[gene], 0), TRUE ~ 0), Abs_CohenD = abs(CohenD), 
          NES = coalesce(NES, 0), Classification = dplyr::case_when(gene %in% c("Nlgn1", "Nrxn1") ~ "Nlgn1/Nrxn1", 
              gene %in% splicing_genes ~ "RNA splicing", gene %in% global_trgs ~ "Other TRGs", TRUE ~ "Non-TRGs"))
      graph_obj <- as_tbl_graph(edge_df) %>% activate(nodes) %>% dplyr::inner_join(meta_df, by = c(name = "gene")) %>% 
          activate(edges) %>% dplyr::mutate(to_node = .N()$name[to])
      return(graph_obj)
  }
  grn_graph_main <- build_tidy_graph(edges_main, tfs_main)
  stopifnot(!is.null(grn_graph_main))
  all_nodes_combined <- grn_graph_main %>% activate(nodes) %>% as_tibble()
  trg_nodes_all <- all_nodes_combined %>% filter(!is_tf)
  global_trg_nes_limits <- if (nrow(trg_nodes_all) > 0) range(trg_nodes_all$NES, na.rm = TRUE) else c(0, 
      1)
  if (global_trg_nes_limits[1] == global_trg_nes_limits[2]) global_trg_nes_limits <- global_trg_nes_limits + 
      c(-0.1, 0.1)
  tf_nodes_all <- all_nodes_combined %>% filter(is_tf)
  global_tf_cohend_limits <- if (nrow(tf_nodes_all) > 0) range(tf_nodes_all$CohenD, na.rm = TRUE) else c(-1.5, 
      1.5)
  if (global_tf_cohend_limits[1] == global_tf_cohend_limits[2]) global_tf_cohend_limits <- global_tf_cohend_limits + 
      c(-0.1, 0.1)
  global_tf_abs_limits <- if (nrow(tf_nodes_all) > 0) range(tf_nodes_all$Abs_CohenD, na.rm = TRUE) else c(0, 
      1.5)
  if (global_tf_abs_limits[1] == global_tf_abs_limits[2]) global_tf_abs_limits <- global_tf_abs_limits + 
      c(-0.1, 0.1)
  draw_publication_grn <- function(graph_object, title_suffix, trg_nes_lim = global_trg_nes_limits, tf_cohend_lim = global_tf_cohend_limits, 
      tf_abs_lim = global_tf_abs_limits) {
      if (is.null(graph_object)) 
          return(NULL)
      coords <- igraph::layout_components(graph_object, layout = igraph::layout_with_fr, weights = NA)
      p <- ggraph(graph_object, layout = "manual", x = coords[, 1], y = coords[, 2]) + geom_edge_diagonal0(aes(filter = !(to_node %in% 
          functional_target_pool)), color = "#eeeeee", alpha = 0.5, linewidth = 0.25, arrow = arrow(length = unit(0.04, 
          "cm"), type = "closed"), end_cap = circle(0.1, "cm")) + geom_edge_diagonal0(aes(filter = to_node %in% 
          functional_target_pool), color = "#eeeeee", alpha = 0.5, linewidth = 0.25, arrow = arrow(length = unit(0.04, 
          "cm"), type = "closed"), end_cap = circle(0.1, "cm")) + geom_node_point(aes(fill = NES, size = NES, 
          filter = !is_tf & Classification == "Non-TRGs"), shape = 21, color = "#ffffff", stroke = 0.1, 
          alpha = 0.75) + geom_node_point(aes(fill = NES, size = NES, filter = !is_tf & Classification == 
          "Other TRGs"), shape = 21, color = "#000000", stroke = 0.15, alpha = 0.95) + geom_node_point(aes(fill = NES, 
          size = NES, filter = !is_tf & Classification == "RNA splicing"), shape = 21, color = "#ffd700", 
          stroke = 1, alpha = 0.95) + geom_node_point(aes(fill = NES, size = NES, filter = !is_tf & Classification == 
          "Nlgn1/Nrxn1"), shape = 21, color = "#18a799", stroke = 1.2, alpha = 0.95) + scale_fill_gradientn(colors = c("#f4f9f4", 
          "#7fcdbb", "#2c7fb8", "#253494"), limits = trg_nes_lim, name = "Target\nNES") + scale_size_continuous(range = c(2, 
          5), limits = trg_nes_lim, name = "Target\nNES") + new_scale_color() + new_scale_fill() + new_scale("size") + 
          geom_node_point(aes(fill = CohenD, size = Abs_CohenD, filter = is_tf), shape = 21, stroke = 1.2, 
              color = "#111111") + scale_fill_gradient2(low = "#4393c3", mid = "white", high = "#d6604d", 
          midpoint = 0, limits = tf_cohend_lim, name = "TF Specific\nCohen's d") + scale_size_continuous(range = c(5, 
          9), limits = tf_abs_lim, name = "TF\n|Cohen's d|") + geom_node_text(aes(label = name, filter = is_tf), 
          repel = TRUE, color = "#000000", fontface = "bold.italic", size = 3.8, point.padding = unit(0.4, 
              "lines"), segment.color = "#222222", segment.size = 0.38, max.overlaps = Inf) + geom_node_text(aes(label = name, 
          filter = !is_tf & (Classification %in% c("Nlgn1/Nrxn1", "RNA splicing"))), repel = TRUE, color = "#2c3e50", 
          fontface = "bold", size = 3, point.padding = unit(0.3, "lines"), segment.color = "#7f8c8d", segment.size = 0.35, 
          max.overlaps = Inf) + labs(title = sprintf("ATN Transcriptional Regulatory Network - %s", title_suffix), 
          subtitle = "Globally synchronized colormaps & sizes: Target mapped to NES, TF mapped to Cohen's d.") + 
          theme_void(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 11, margin = margin(b = 2)), 
          plot.subtitle = element_text(color = "grey30", size = 7.5, margin = margin(b = 12)), legend.position = "right", 
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5))
      return(p)
  }
  set.seed(42)
  p1 <- draw_publication_grn(grn_graph_main, "All qualifying TFs")
  ggsave(output_path("figures/figure 6/6d.png"), p1, width = 12, height = 12, dpi = 400, bg = "white")
  write_csv(network_edges_filtered, output_path("sourcedata/figure 6/6d.csv"))
  write_csv(all_nodes_combined, output_path("sourcedata/figure 6/6d_nodes.csv"))
  cat("Finished.\n")
  invisible(as.list(environment()))
}
