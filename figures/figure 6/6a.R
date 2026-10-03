# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(tidyverse)
      library(ggrepel)
  })
  CACHE_RDS <- data_path("results/260616 pySCENIC/auc_with_meta.rds")
  TF_DB_PATH <- data_path("bin/Mouse_TFs_Kinases_webpage-3-30-2017.csv")
  out_pdf <- output_path("figures/figure 6/6a.pdf")
  out_png <- output_path("figures/figure 6/6a.png")
  out_csv <- output_path("sourcedata/figure 6/6a.csv")
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
  subset_idx <- global_meta$final_region == "ATN"
  sub_auc <- global_auc_mat[subset_idx, , drop = FALSE]
  sub_meta <- global_meta[subset_idx, ]
  sub_group <- sub_meta$Group_L2
  mat_saline <- sub_auc[sub_group == "Saline_I", , drop = FALSE]
  mat_sham <- sub_auc[sub_group == "Sham_I", , drop = FALSE]
  mat_dbs <- sub_auc[sub_group == "DBS_I", , drop = FALSE]
  n_sal <- nrow(mat_saline)
  m_sal <- colMeans(mat_saline)
  v_sal <- colVars_fast(mat_saline)
  n_shm <- nrow(mat_sham)
  m_shm <- colMeans(mat_sham)
  v_shm <- colVars_fast(mat_sham)
  n_dbs <- nrow(mat_dbs)
  m_dbs <- colMeans(mat_dbs)
  v_dbs <- colVars_fast(mat_dbs)
  eps <- 1e-10
  pooled_sd_dbs <- sqrt(((n_dbs - 1) * v_dbs + (n_shm - 1) * v_shm)/(n_dbs + n_shm - 2))
  pooled_sd_dbs[pooled_sd_dbs == 0] <- eps
  pooled_sd_sal <- sqrt(((n_sal - 1) * v_sal + (n_shm - 1) * v_shm)/(n_sal + n_shm - 2))
  pooled_sd_sal[pooled_sd_sal == 0] <- eps
  idx_dbs <- which(sub_group == "DBS_I")
  idx_sham <- which(sub_group == "Sham_I")
  p_vals <- setNames(numeric(ncol(sub_auc)), colnames(sub_auc))
  for (reg in colnames(sub_auc)) {
      wt <- suppressWarnings(wilcox.test(sub_auc[idx_dbs, reg], sub_auc[idx_sham, reg], exact = FALSE))
      p_vals[reg] <- wt$p.value
  }
  tf_metrics <- tibble(Regulon = colnames(sub_auc), TF = clean_tf_name(Regulon), DBSvSham_cohend = (m_dbs - 
      m_shm)/pooled_sd_dbs, SalinevSham_cohend = (m_sal - m_shm)/pooled_sd_sal, DBSvSham_p_donor = p_vals) %>% 
      mutate(DBSvSham_FDR = p.adjust(DBSvSham_p_donor, method = "BH")) %>% filter(!is.na(DBSvSham_cohend) & 
      !is.na(SalinevSham_cohend) & !is.infinite(DBSvSham_cohend) & !is.infinite(SalinevSham_cohend)) %>% 
      mutate(dir_not_same = sign(DBSvSham_cohend) != sign(SalinevSham_cohend), sc_not_sig = DBSvSham_FDR >= 
          0.05, is_gray = dir_not_same | sc_not_sig)
  plot_df <- tf_metrics %>% arrange(desc(DBSvSham_cohend)) %>% mutate(Regulon = factor(Regulon, levels = Regulon))
  df_gray <- plot_df %>% filter(is_gray == TRUE)
  df_colored <- plot_df %>% filter(is_gray == FALSE)
  write_csv(plot_df %>% select(TF, DBSvSham_cohend, SalinevSham_cohend, is_gray), out_csv)
  write_csv(plot_df, output_path("statistics/figure 6/6a_statistics.csv"))
  PALETTE <- c(`DBS vs Sham` = "#ffd700", `Saline vs Sham` = "#18a799")
  p <- ggplot(plot_df, aes(y = Regulon)) + geom_segment(data = df_gray, aes(yend = Regulon, x = SalinevSham_cohend, 
      xend = DBSvSham_cohend), color = "#f3f3f3", linewidth = 0.3) + geom_segment(data = df_colored, aes(yend = Regulon, 
      x = SalinevSham_cohend, xend = DBSvSham_cohend), color = "#bbbbbb", linewidth = 0.5) + geom_point(data = df_gray, 
      aes(x = SalinevSham_cohend, size = abs(SalinevSham_cohend)), color = "#e5e5e5", alpha = 0.4) + geom_point(data = df_gray, 
      aes(x = DBSvSham_cohend, size = abs(DBSvSham_cohend)), color = "#e5e5e5", alpha = 0.4) + geom_point(data = df_colored, 
      aes(x = SalinevSham_cohend, color = "Saline vs Sham", size = abs(SalinevSham_cohend)), alpha = 0.95) + 
      geom_point(data = df_colored, aes(x = DBSvSham_cohend, color = "DBS vs Sham", size = abs(DBSvSham_cohend)), 
          alpha = 0.95) + geom_vline(xintercept = 0, linetype = "dashed", color = "#777777", linewidth = 0.4) + 
      geom_text_repel(data = df_colored, aes(x = DBSvSham_cohend, label = TF), color = "#111111", fontface = "italic", 
          size = 2.5, box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"), segment.color = "#444444", 
          segment.size = 0.3, min.segment.length = 0, max.overlaps = Inf, force = 3.5, direction = "both") + 
      scale_y_discrete(limits = rev(levels(plot_df$Regulon)), drop = FALSE) + scale_size_continuous(range = c(3, 
      3), name = "|Cohen's d|") + scale_color_manual(values = PALETTE) + labs(y = NULL, x = "Standardized Effect Size (Cohen's d)", 
      color = NULL) + theme_bw(base_size = 9) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
      panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "top", 
      legend.box = "horizontal", legend.box.spacing = unit(0.05, "cm"))
  n_tfs <- nrow(plot_df)
  calc_height <- max(5, n_tfs * 0.08 + 2)
  ggsave(out_png, p, width = 7, height = calc_height, dpi = 600, limitsize = FALSE)
  cat(">>> Figure 6a complete.\n")
  invisible(as.list(environment()))
}
