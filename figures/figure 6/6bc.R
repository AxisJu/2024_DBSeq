# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(tidyverse)
      library(ComplexHeatmap)
      library(circlize)
  })
  CACHE_RDS <- data_path("results/260616 pySCENIC/auc_with_meta.rds")
  TF_DB_PATH <- data_path("bin/Mouse_TFs_Kinases_webpage-3-30-2017.csv")
  MOTIF_DIR <- data_path("results/260616 pySCENIC/motifs_aerts")
  out_pdf <- output_path("figures/figure 6/6bc.pdf")
  out_png <- output_path("figures/figure 6/6bc.png")
  out_csv <- output_path("sourcedata/figure 6/6bc.csv")
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  GROUPS <- c("Saline_I", "Sham_I", "DBS_I")
  colVars_fast <- function(x) {
      n <- nrow(x)
      if (n <= 1) 
          return(rep(0, ncol(x)))
      (colSums(x^2) - (colSums(x)^2)/n)/(n - 1)
  }
  clean_tf_name <- function(x) trimws(gsub("\\s*\\(.*?\\)", "", x))
  list2env(supplied_inputs[["helpers_2"]], envir = environment())
  data_obj <- aggregate_donor_auc(supplied_inputs[["data_obj_3"]][[as.character(CACHE_RDS)]])
  global_auc_mat <- data_obj$auc_mat
  global_meta <- data_obj$meta
  tf_db <- supplied_inputs[["tf_db_4"]][[as.character(TF_DB_PATH)]]
  colnames(tf_db)[1] <- "Symbol"
  pure_tfs <- tf_db %>% filter(Class == "TF") %>% pull(Symbol) %>% unique()
  current_regulons <- colnames(global_auc_mat)
  cleaned_symbols <- clean_tf_name(current_regulons)
  valid_tf_idx <- cleaned_symbols %in% pure_tfs
  global_auc_mat <- global_auc_mat[, valid_tf_idx, drop = FALSE]
  eps <- 1e-10
  atn_idx <- which(global_meta$final_region == "ATN")
  atn_auc <- global_auc_mat[atn_idx, , drop = FALSE]
  atn_meta <- global_meta[atn_idx, ]
  m_sal <- colMeans(atn_auc[atn_meta$Group_L2 == "Saline_I", , drop = FALSE])
  m_shm <- colMeans(atn_auc[atn_meta$Group_L2 == "Sham_I", , drop = FALSE])
  m_dbs <- colMeans(atn_auc[atn_meta$Group_L2 == "DBS_I", , drop = FALSE])
  v_sal <- colVars_fast(atn_auc[atn_meta$Group_L2 == "Saline_I", , drop = FALSE])
  v_shm <- colVars_fast(atn_auc[atn_meta$Group_L2 == "Sham_I", , drop = FALSE])
  v_dbs <- colVars_fast(atn_auc[atn_meta$Group_L2 == "DBS_I", , drop = FALSE])
  n_sal <- sum(atn_meta$Group_L2 == "Saline_I")
  n_shm <- sum(atn_meta$Group_L2 == "Sham_I")
  n_dbs <- sum(atn_meta$Group_L2 == "DBS_I")
  sd_ds <- sqrt(((n_dbs - 1) * v_dbs + (n_shm - 1) * v_shm)/(n_dbs + n_shm - 2))
  sd_ds[sd_ds == 0] <- eps
  sd_ss <- sqrt(((n_sal - 1) * v_sal + (n_shm - 1) * v_shm)/(n_sal + n_shm - 2))
  sd_ss[sd_ss == 0] <- eps
  p_vals <- setNames(numeric(ncol(atn_auc)), colnames(atn_auc))
  idx_d <- which(atn_meta$Group_L2 == "DBS_I")
  idx_s <- which(atn_meta$Group_L2 == "Sham_I")
  for (reg in colnames(atn_auc)) {
      wt <- suppressWarnings(wilcox.test(atn_auc[idx_d, reg], atn_auc[idx_s, reg], exact = FALSE))
      p_vals[reg] <- wt$p.value
  }
  filtered_tfs_df <- tibble(Regulon = colnames(atn_auc), TF = clean_tf_name(Regulon), DBSvSham_cohend = (m_dbs - 
      m_shm)/sd_ds, SalinevSham_cohend = (m_sal - m_shm)/sd_ss, DBSvSham_p_donor = p_vals) %>% mutate(DBSvSham_FDR = p.adjust(DBSvSham_p_donor, 
      method = "BH")) %>% filter(!is.na(DBSvSham_cohend) & !is.na(SalinevSham_cohend) & !is.infinite(DBSvSham_cohend) & 
      !is.infinite(SalinevSham_cohend)) %>% filter(sign(DBSvSham_cohend) == sign(SalinevSham_cohend) & 
      DBSvSham_FDR < 0.05) %>% arrange(desc(DBSvSham_cohend)) %>% mutate(Direction = if_else(DBSvSham_cohend > 
      0, "Upregulated", "Downregulated"))
  write.csv(filtered_tfs_df, output_path("statistics/figure 6/6bc_selected_TFs.csv"), row.names = FALSE)
  if (nrow(filtered_tfs_df) == 0) {
      blank <- ggplot() + annotate("text", x = 0, y = 0, label = "No TF passes donor-level BH FDR < 0.05", 
          size = 4) + theme_void()
      ggsave(output_path("figures/figure 6/6bc.png"), blank, width = 7, height = 4, dpi = 300, bg = "white")
      write.csv(filtered_tfs_df[, c("TF", "DBSvSham_cohend", "SalinevSham_cohend")], output_path("sourcedata/figure 6/6bc.csv"), 
          row.names = FALSE)
      return(invisible(NULL))
  }
  target_regulons <- filtered_tfs_df$Regulon
  region_order <- c("ATN", "CP", "RSPd", "RT", "RE", "GPe", "LH", "PF", "CM", "DG", "LSc", "MH", "TRS")
  region_colors <- c(RSPd = "#ff1a71", DG = "#16f2f2", GPe = "#b199ff", CP = "#8d7acc", LSc = "#3283fe", 
      RT = "#ff6600", TRS = "#7609b1", MH = "#faa307", LH = "#c68105", ATN = "#1460ff", RE = "#1340ff", 
      PF = "#0a4093", CM = "#08306d")
  comp_colors <- c(`DBS vs Sham` = "#ffd700", `Saline vs Sham` = "#18a799")
  n_regions <- length(region_order)
  cohend_mat <- matrix(NA_real_, nrow = length(target_regulons), ncol = 2 * n_regions)
  rownames(cohend_mat) <- target_regulons
  meta_columns <- tibble(Region = rep(region_order, each = 2), Comparison = rep(c("DBS vs Sham", "Saline vs Sham"), 
      times = n_regions))
  for (i in seq_along(region_order)) {
      r <- region_order[i]
      ridx <- which(global_meta$final_region == r)
      if (length(ridx) < 5) 
          next
      rauc <- global_auc_mat[ridx, target_regulons, drop = FALSE]
      rmeta <- global_meta[ridx, ]
      s_sal <- rauc[rmeta$Group_L2 == "Saline_I", , drop = FALSE]
      s_shm <- rauc[rmeta$Group_L2 == "Sham_I", , drop = FALSE]
      s_dbs <- rauc[rmeta$Group_L2 == "DBS_I", , drop = FALSE]
      n1 <- nrow(s_sal)
      m1 <- colMeans(s_sal)
      v1 <- colVars_fast(s_sal)
      n2 <- nrow(s_shm)
      m2 <- colMeans(s_shm)
      v2 <- colVars_fast(s_shm)
      n3 <- nrow(s_dbs)
      m3 <- colMeans(s_dbs)
      v3 <- colVars_fast(s_dbs)
      sd1 <- sqrt(((n3 - 1) * v3 + (n2 - 1) * v2)/(n3 + n2 - 2))
      sd1[sd1 == 0] <- eps
      sd2 <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2)/(n1 + n2 - 2))
      sd2[sd2 == 0] <- eps
      cohend_mat[, (2 * i - 1)] <- (m3 - m2)/sd1
      cohend_mat[, (2 * i)] <- (m1 - m2)/sd2
  }
  colnames(cohend_mat) <- paste(meta_columns$Region, meta_columns$Comparison, sep = "_")
  cd_long <- as_tibble(cohend_mat, rownames = "Regulon") %>% pivot_longer(-Regulon, names_to = "Region_Comparison", 
      values_to = "CohenD") %>% separate(Region_Comparison, into = c("Region", "Comparison"), sep = "_", 
      extra = "merge") %>% mutate(TF = clean_tf_name(Regulon))
  write_csv(cd_long %>% select(TF, Region, Comparison, CohenD), out_csv)
  generate_text_badge <- function(tf_name, output_path) {
      png(output_path, width = 150, height = 50, bg = "transparent")
      par(mar = c(0, 0, 0, 0))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      rect(0.05, 0.1, 0.95, 0.9, col = "#f0f0f2", border = "#dddddf", lwd = 1)
      text(x = 0.5, y = 0.5, labels = tf_name, cex = 1.1, col = "#555555", font = 3)
      dev.off()
  }
  tf_names <- filtered_tfs_df$TF
  motif_files <- character(length(tf_names))
  names(motif_files) <- tf_names
  for (tf in tf_names) {
      img_path <- file.path(MOTIF_DIR, sprintf("%s.png", tf))
      if (supplied_inputs[["data_5"]]) {
          motif_files[tf] <- img_path
      }
      else {
          badge_path <- output_path("cache", sprintf("%s_badge.png", tf))
          if (!supplied_inputs[["data_6"]]) 
              generate_text_badge(tf, badge_path)
          motif_files[tf] <- badge_path
      }
  }
  top_anno <- HeatmapAnnotation(Region = factor(meta_columns$Region, levels = region_order), Comparison = meta_columns$Comparison, 
      col = list(Region = region_colors, Comparison = comp_colors), show_legend = TRUE, annotation_name_side = "left", 
      simple_anno_size = unit(0.35, "cm"))
  right_motif_anno <- rowAnnotation(Motif = anno_image(image = motif_files, width = unit(1.8, "cm"), height = unit(0.5, 
      "cm"), border = FALSE, gp = gpar(fill = "transparent")), show_annotation_name = FALSE)
  row_order_vec <- match(target_regulons, rownames(cohend_mat))
  row_split_factor <- factor(filtered_tfs_df$Direction, levels = c("Upregulated", "Downregulated"))
  max_abs_d <- max(abs(cohend_mat), na.rm = TRUE)
  rdbu_color_fun <- colorRamp2(c(-max_abs_d, 0, max_abs_d), c("#2b63a3", "#ffffff", "#a2163e"))
  ht <- Heatmap(matrix = cohend_mat, name = "Cohen's d", col = rdbu_color_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
      row_order = row_order_vec, row_split = row_split_factor, column_split = factor(meta_columns$Region, 
          levels = region_order), column_title_gp = gpar(fontsize = 8, fontface = "bold"), column_gap = unit(1.5, 
          "mm"), top_annotation = top_anno, right_annotation = right_motif_anno, show_column_names = FALSE, 
      row_names_gp = gpar(fontsize = 7, fontface = "italic"), row_names_side = "left", row_title_gp = gpar(fontsize = 9, 
          fontface = "bold"), use_raster = FALSE)
  pdf(NULL, width = 9, height = 5)
  draw(ht)
  dev.off()
  png(out_png, width = 9, height = 5, units = "in", res = 600)
  draw(ht)
  dev.off()
  cat(">>> Figure 6b/c complete.\n")
  invisible(as.list(environment()))
}
