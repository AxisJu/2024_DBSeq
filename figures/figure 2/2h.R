# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(qs)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(cowplot)
  proj_dir <- data_root
  out_fig_dir <- output_path("figures/figure 2")
  out_data_dir <- output_path("sourcedata/figure 2")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  gene_groups <- list(Common = c("Rnpc3", "Crebzf", "Ndfip1", "Ckb", "Pura"), Neuronal = c("Ube3a", "Bsg", 
      "Matk", "Timm9", "Myl12b"), `Non-neuronal` = c("Mertk", "Ccnd3", "Rab31", "Csrp1", "Sox10"))
  all_target_genes <- unlist(gene_groups)
  ordered_ct <- c("Glu.N", "GABA.N", "Astro.", "OPCs", "Oligo.", "Endo.", "Microglia")
  col_split <- factor(ifelse(ordered_ct %in% c("Glu.N", "GABA.N"), "Neuron", "Non-Neuron"), levels = c("Neuron", 
      "Non-Neuron"))
  ct_colors <- c(Glu.N = "#ba1369", GABA.N = "#56ad00", Astro. = "#a58946", OPCs = "#5953ff", Oligo. = "#201e5a", 
      Endo. = "#858881", Microglia = "#a87c5a", Not_TRG = "#f4f9f4")
  mouse_df <- supplied_inputs[["mouse_df_2"]]
  mouse_non_neuron <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  mouse_neuron <- setdiff(unique(mouse_df$Celltype), mouse_non_neuron)
  mouse_mapped <- mouse_df %>% mutate(Human_Map = case_when(Celltype %in% mouse_neuron & str_detect(Celltype, 
      "ExN") ~ "Glu.N", Celltype %in% mouse_neuron & str_detect(Celltype, "InN") ~ "GABA.N", Celltype == 
      "Astro" ~ "Astro.", Celltype == "OPC" ~ "OPCs", Celltype == "Oligo" ~ "Oligo.", Celltype == "Micro" ~ 
      "Microglia", Celltype == "VCs" ~ "Endo.", TRUE ~ "Other")) %>% filter(Human_Map != "Other")
  mouse_sub <- mouse_mapped %>% filter(genes %in% all_target_genes) %>% group_by(genes, Human_Map) %>% 
      summarise(Mouse_DBSvSham = -median(CoD_DBSvSham_cohend, na.rm = TRUE), Mouse_ShamvSaline = -median(CoD_ShamvSal_cohend, 
          na.rm = TRUE), is_TRG = any(is_TRG == TRUE, na.rm = TRUE), .groups = "drop")
  trg_df <- mouse_sub %>% select(genes, Human_Map, is_TRG) %>% pivot_wider(names_from = Human_Map, values_from = is_TRG, 
      values_fill = FALSE)
  mouse_cohend <- mouse_sub %>% select(-is_TRG) %>% pivot_longer(cols = c(Mouse_DBSvSham, Mouse_ShamvSaline), 
      names_to = "Condition", values_to = "CohenD")
  ortho_source <- Sys.getenv("DBSEQ_ORTHOLOG_FILE", file.path(data_root, "results/Table/human_cohort_de/conserved_trg.csv"))
  orthologs_1to1 <- supplied_inputs[["orthologs_1to1_3"]][[as.character(ortho_source)]] %>% filter(Has_Human_Homolog, 
      Mouse_Gene %in% all_target_genes) %>% select(Mouse_Gene, Human_Gene) %>% distinct()
  stopifnot(setequal(orthologs_1to1$Mouse_Gene, all_target_genes), !anyDuplicated(orthologs_1to1$Mouse_Gene), 
      !anyDuplicated(orthologs_1to1$Human_Gene))
  write.csv(orthologs_1to1 %>% mutate(Source = ortho_source, Ensembl_Release = 115), file.path(out_data_dir, 
      "2h_orthologs.csv"), row.names = FALSE)
  mouse_support <- mouse_sub %>% transmute(Mouse_Gene = genes, Human_Celltype = Human_Map, Mouse_available = TRUE)
  cohort_names <- paste0("cohort", 1:4)
  all_human_res <- list()
  audit_rows <- list()
  for (ct in ordered_ct) {
      ct_mouse_ref <- mouse_support %>% filter(Human_Celltype == ct) %>% inner_join(orthologs_1to1, by = "Mouse_Gene")
      for (cohort_name in cohort_names) {
          sc_path <- file.path(data_root, "results/Table/human_cohort_de/sc", paste0("SC_DE_", cohort_name, 
              "_", ct, ".csv"))
          pb_path <- results_path("Table/human_cohort_de/pb", paste0("PB_DE_", cohort_name, "_", ct, ".csv"))
          sc <- data.frame(Human_Gene = character(), SC_FDR = numeric(), SC_CoD = numeric(), SC_available = logical())
          pb <- data.frame(Human_Gene = character(), PB_P = numeric(), PB_FDR = numeric(), PB_donor_validated = logical(), 
              PB_CoD = numeric(), PB_available = logical())
          if (supplied_inputs[["data_4"]]) 
              sc <- supplied_inputs[["sc_5"]][[as.character(sc_path)]] %>% transmute(Human_Gene = Gene, 
                  SC_FDR = p_val_adj, SC_CoD = CohenD, SC_available = TRUE)
          if (supplied_inputs[["data_6"]]) 
              pb <- supplied_inputs[["pb_7"]][[as.character(pb_path)]] %>% mutate(Patient_validated = if ("biological_unit" %in% 
                  names(.)) 
                  biological_unit == "Patient"
              else FALSE, Adjusted = if ("FDR" %in% names(.)) 
                  FDR
              else p_val_adj, Effect = if ("effect_direction" %in% names(.)) 
                  -CohenD
              else CohenD) %>% transmute(Human_Gene = Gene, PB_P = PValue, PB_FDR = Adjusted, PB_donor_validated = Patient_validated, 
                  PB_CoD = Effect, PB_available = TRUE)
          stopifnot(!anyDuplicated(sc$Human_Gene), !anyDuplicated(pb$Human_Gene))
          all_human_res[[paste(ct, cohort_name)]] <- full_join(sc, pb, by = "Human_Gene") %>% inner_join(ct_mouse_ref, 
              by = "Human_Gene") %>% mutate(Cohort = cohort_name)
          audit_rows[[paste(ct, cohort_name)]] <- orthologs_1to1 %>% left_join(sc, by = "Human_Gene") %>% 
              left_join(pb, by = "Human_Gene") %>% mutate(Human_Celltype = ct, Cohort = cohort_name) %>% 
              left_join(mouse_support, by = c("Mouse_Gene", "Human_Celltype")) %>% mutate(across(c(SC_available, 
              PB_available, Mouse_available), ~coalesce(.x, FALSE)), Displayed = Mouse_available & (is.finite(SC_CoD) | 
              is.finite(PB_CoD)), Raw_CohenD = -ifelse(PB_donor_validated & is.finite(PB_CoD), PB_CoD, 
              SC_CoD), Display_Reason = case_when(!Mouse_available ~ "no_matching_mouse_reference", PB_donor_validated & 
              is.finite(PB_CoD) ~ "donor_level_effect", !is.finite(SC_CoD) ~ "nonfinite_effect", TRUE ~ 
              "historical_cell_effect_descriptive"))
      }
  }
  human_sub <- bind_rows(all_human_res) %>% mutate(Condition = paste0("Human_", Cohort), CohenD = -ifelse(coalesce(PB_donor_validated, 
      FALSE) & is.finite(PB_CoD), PB_CoD, SC_CoD), sig_sc = !is.na(SC_FDR) & SC_FDR < 0.01, sig_pb = Cohort != 
      "cohort2" & PB_donor_validated & !is.na(PB_FDR) & PB_FDR < 0.05, Sig_Label = ifelse(sig_pb, "*", 
      ""))
  stopifnot(!anyDuplicated(human_sub[c("Mouse_Gene", "Human_Celltype", "Condition")]))
  write.csv(bind_rows(audit_rows), output_path("statistics/figure 2/2h_human_display_audit.csv"), row.names = FALSE)
  heat_wide <- bind_rows(mouse_cohend, human_sub %>% select(genes = Mouse_Gene, Human_Map = Human_Celltype, 
      Condition, CohenD)) %>% mutate(ColID = paste(Human_Map, genes, sep = "___")) %>% select(Condition, 
      ColID, CohenD) %>% pivot_wider(names_from = ColID, values_from = CohenD, values_fill = NA_real_) %>% 
      as.data.frame()
  rownames(heat_wide) <- heat_wide$Condition
  heat_wide$Condition <- NULL
  sig_wide <- human_sub %>% mutate(ColID = paste(Human_Celltype, Mouse_Gene, sep = "___")) %>% select(Condition, 
      ColID, Sig_Label) %>% pivot_wider(names_from = ColID, values_from = Sig_Label, values_fill = "", 
      values_fn = function(x) x[1]) %>% as.data.frame()
  rownames(sig_wide) <- sig_wide$Condition
  sig_wide$Condition <- NULL
  row_order <- c("Mouse_DBSvSham", "Mouse_ShamvSaline", paste0("Human_cohort", 1:4))
  for (rn in row_order) {
      if (!rn %in% rownames(heat_wide)) 
          heat_wide[rn, ] <- NA_real_
  }
  for (rn in paste0("Human_cohort", 1:4)) {
      if (!rn %in% rownames(sig_wide)) 
          sig_wide[rn, ] <- ""
  }
  heat_wide <- heat_wide[row_order, ]
  all_expected_cols <- as.vector(outer(ordered_ct, all_target_genes, function(x, y) paste(x, y, sep = "___")))
  missing_cols <- setdiff(all_expected_cols, colnames(heat_wide))
  if (length(missing_cols) > 0) heat_wide[missing_cols] <- NA_real_
  missing_sig_cols <- setdiff(all_expected_cols, colnames(sig_wide))
  if (length(missing_sig_cols) > 0) sig_wide[missing_sig_cols] <- ""
  source_rows <- list()
  for (cat_name in names(gene_groups)) {
      for (g in gene_groups[[cat_name]]) {
          g_trg_status <- if (g %in% trg_df$genes) 
              trg_df[trg_df$genes == g, ]
          else data.frame()
          trg_vec <- sapply(ordered_ct, function(ct) {
              if (ct %in% colnames(g_trg_status) && isTRUE(g_trg_status[[ct]])) 
                  1
              else 0
          })
          r_trg <- as.data.frame(t(trg_vec))
          colnames(r_trg) <- ordered_ct
          r_trg$Gene <- g
          r_trg$Category <- cat_name
          r_trg$Condition <- "is_TRG"
          r_trg$Significance <- ""
          source_rows[[length(source_rows) + 1]] <- r_trg
          for (cond in row_order) {
              cols <- paste(ordered_ct, g, sep = "___")
              vals <- as.numeric(heat_wide[cond, cols])
              sigs <- if (cond %in% rownames(sig_wide)) 
                  as.character(sig_wide[cond, cols])
              else rep("", length(ordered_ct))
              sigs[is.na(sigs)] <- ""
              r_cond <- as.data.frame(t(vals))
              colnames(r_cond) <- ordered_ct
              r_cond$Gene <- g
              r_cond$Category <- cat_name
              r_cond$Condition <- cond
              r_cond$Significance <- paste(sigs, collapse = ",")
              source_rows[[length(source_rows) + 1]] <- r_cond
          }
      }
  }
  source_h <- bind_rows(source_rows) %>% select(Category, Gene, Condition, all_of(ordered_ct), Significance)
  write.csv(source_h, file.path(out_data_dir, "2h.csv"), row.names = FALSE)
  cat("Saved source data: sourcedata/figure 2/2h.csv\n")
  get_breaks_cols <- function(mat, exp_sign, end_cols) {
      mx <- max(mat, na.rm = TRUE)
      mn <- min(mat, na.rm = TRUE)
      if (is.infinite(mx) || is.na(mx)) 
          mx <- 0.01
      if (is.infinite(mn) || is.na(mn)) 
          mn <- -0.01
      if (mx <= 1e-05) 
          mx <- 0.01
      if (mn >= -1e-05) 
          mn <- -0.01
      base_col <- "#f4f9f4"
      opp_col <- "#dddddd"
      len <- length(end_cols)
      if (exp_sign == 1) {
          pos_brks <- seq(0, mx, length.out = len + 1)[-1]
          brks <- c(mn, -1e-06, 0, pos_brks)
          cls <- c(opp_col, opp_col, base_col, end_cols)
      }
      else {
          neg_brks <- seq(mn, 0, length.out = len + 1)[-(len + 1)]
          brks <- c(neg_brks, 0, 1e-06, mx)
          cls <- c(rev(end_cols), base_col, opp_col, opp_col)
      }
      return(colorRamp2(brks, cls))
  }
  draw_2h <- function() {
      grid.newpage()
      grid.rect(gp = gpar(fill = "white", col = NA))
      display_ct <- c("ExN", "InN", "Astro", "OPC", "Oligo", "VC", "Micro")
      cell_x <- seq(0.045, by = 0.088, length.out = 7) + c(0, 0, rep(0.022, 5))
      row_y <- c(0.84, 0.71, 0.58, 0.43, 0.34, 0.25, 0.16)
      legend_bar <- function(col_fun, limits, y, height) {
          vals <- seq(limits[1], limits[2], length.out = 81)
          yy <- y - height/2 + (seq_len(80) - 0.5) * height/80
          grid.rect(x = 0.745, y = yy, width = 0.032, height = height/80 + 5e-04, gp = gpar(fill = col_fun((vals[-1] + 
              vals[-81])/2), col = NA))
          grid.text(format(signif(limits[c(2, 1)], 2), trim = TRUE), x = 0.779, y = c(y + height/2, y - 
              height/2), just = "left", gp = gpar(fontsize = 7.5))
      }
      for (r in seq_along(gene_groups)) {
          cy <- 0.04 + 0.96 * (1 - (r - 0.5)/3)
          pushViewport(viewport(x = 0.071, y = cy, width = 0.138, height = 0.96/3))
          grid.text(names(gene_groups)[r], x = 0.03, y = 0.96, just = "left", gp = gpar(fontsize = 12, 
              fontface = "bold", col = c("#d8af00", "#159b95", "black")[r]))
          grid.text(c("isTRG", "DBS vs Sham", "Sham vs Saline", paste("DRE Cohort", 1:4)), x = 0.98, y = row_y, 
              just = "right", gp = gpar(fontsize = 10))
          popViewport()
          for (j in seq_along(gene_groups[[r]])) {
              g <- gene_groups[[r]][j]
              pushViewport(viewport(x = 0.145 + (j - 0.5) * 0.17, y = cy, width = 0.17, height = 0.96/3))
              cols <- paste(ordered_ct, g, sep = "___")
              mat <- as.matrix(heat_wide[, cols, drop = FALSE])
              trg <- sapply(ordered_ct, function(ct) {
                  z <- mouse_sub %>% filter(genes == g, Human_Map == ct)
                  nrow(z) > 0 && isTRUE(z$is_TRG[1])
              })
              main_dir <- sign(median(mat[2, trg], na.rm = TRUE))
              if (!is.finite(main_dir) || main_dir == 0) 
                  main_dir <- 1
              funs <- list(get_breaks_cols(mat[1, , drop = FALSE], -main_dir, "#f2d016"), get_breaks_cols(mat[2, 
                  , drop = FALSE], main_dir, "#7e418f"), get_breaks_cols(mat[3:6, , drop = FALSE], main_dir, 
                  c("#7fcdbb", "#2c7fb8", "#253494")))
              grid.text(g, x = 0.005, y = 0.96, just = "left", gp = gpar(fontsize = 12, fontface = "italic"))
              fills <- ifelse(trg, unname(ct_colors[ordered_ct]), "#f4f9f4")
              grid.rect(x = cell_x, y = row_y[1], width = 0.088, height = 0.09, gp = gpar(fill = fills, 
                  col = "white", lwd = 0.5))
              for (k in 1:6) {
                  fun <- funs[[min(k, 3)]]
                  fills <- rep("#f4f9f4", 7)
                  ok <- is.finite(mat[k, ])
                  fills[ok] <- fun(mat[k, ok])
                  grid.rect(x = cell_x, y = row_y[k + 1], width = 0.088, height = 0.09, gp = gpar(fill = fills, 
                    col = "white", lwd = 0.5))
                  if (k >= 3) {
                    sigs <- as.character(sig_wide[row_order[k], cols])
                    sigs[is.na(sigs)] <- ""
                    grid.text(sigs, x = cell_x, y = row_y[k + 1] - 0.006, gp = gpar(fontsize = 12))
                  }
              }
              for (block in list(c(row_y[1], 0.09), c(row_y[2], 0.09), c(row_y[3], 0.09), c(mean(row_y[4:7]), 
                  0.36))) {
                  grid.rect(x = mean(cell_x[1:2]), y = block[1], width = 2 * 0.088, height = block[2], 
                    gp = gpar(fill = NA, col = "grey60", lwd = 0.8))
                  grid.rect(x = mean(cell_x[3:7]), y = block[1], width = 5 * 0.088, height = block[2], 
                    gp = gpar(fill = NA, col = "grey60", lwd = 0.8))
              }
              for (k in 1:3) {
                  z <- if (k < 3) 
                    mat[k, ]
                  else mat[3:6, ]
                  direction <- if (k == 1) 
                    -main_dir
                  else main_dir
                  lim <- if (direction > 0) 
                    c(0, max(c(z, 0.01), na.rm = TRUE))
                  else c(min(c(z, -0.01), na.rm = TRUE), 0)
                  legend_bar(funs[[k]], lim, c(0.71, 0.58, 0.325)[k], c(0.09, 0.09, 0.23)[k])
              }
              grid.rect(x = 0.745, y = 0.165, width = 0.032, height = 0.027, gp = gpar(fill = "#dddddd", 
                  col = NA))
              grid.text(ifelse(main_dir < 0, ">0", "<0"), x = 0.779, y = 0.165, just = "left", gp = gpar(fontsize = 7.5))
              if (r == 3) 
                  for (k in 1:7) grid.text(display_ct[k], x = cell_x[k], y = 0.091, rot = 48, just = "right", 
                    gp = gpar(fontsize = 10))
              popViewport()
          }
      }
  }
  pdf(NULL, width = 16, height = 8.4, useDingbats = FALSE)
  draw_2h()
  dev.off()
  png(file.path(out_fig_dir, "2h.png"), width = 4800, height = 2520, res = 300, bg = "white")
  draw_2h()
  dev.off()
  cat("Panel 2h exported: original SC/PB/mouse-reference joins applied using Ensembl 115 mappings.\n")
  invisible(as.list(environment()))
}
