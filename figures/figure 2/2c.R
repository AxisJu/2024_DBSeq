# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ggplot2)
  library(eulerr)
  library(patchwork)
  proj_dir <- Sys.getenv("DBSEQ_DATA_ROOT", unset = "")
  if (!nzchar(proj_dir)) stop("Set DBSEQ_DATA_ROOT to the original DBSeq data root.")
  package_root <- output_root
  if (!nzchar(package_root)) {
      script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      if (!length(script_arg)) 
          stop("Set DBSEQ_PACKAGE_ROOT when sourcing this script.")
      script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
      package_root <- dirname(dirname(dirname(script_path)))
  }
  out_fig_dir <- file.path(package_root, "figures", "figure 2")
  out_data_dir <- file.path(package_root, "sourcedata", "figure 2")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  load_de_data <- function(lvl, mode, comp, ct) {
      ct_safe <- gsub("/", "-", ct)
      file_path <- results_path("Table/EdgeR_v2", lvl, paste0("EdgeR_", comp, "_", ct_safe, ".csv"))
      if (!supplied_inputs[["load_de_data_2"]]) 
          return(NULL)
      df <- supplied_inputs[["df_3"]][[as.character(file_path)]]
      gene_col <- if ("X" %in% colnames(df)) 
          "X"
      else colnames(df)[1]
      data.frame(gene = df[[gene_col]], PValue = df$PValue, FDR = df$FDR, logFC = df$logFC)
  }
  get_celltype_venn <- function(ct) {
      de1 <- supplied_inputs[["de1_4"]][[as.character("Level1")]] %>% filter(!is.na(PValue), !is.na(logFC)) %>% 
          distinct(gene, .keep_all = TRUE)
      de2 <- supplied_inputs[["de2_5"]][[as.character("Level1")]] %>% filter(!is.na(PValue), !is.na(logFC)) %>% 
          distinct(gene, .keep_all = TRUE)
      common_genes <- intersect(de1$gene, de2$gene)
      de1 <- de1 %>% filter(gene %in% common_genes)
      de2 <- de2 %>% filter(gene %in% common_genes) %>% mutate(logFC = logFC * -1)
      up1 <- de1$gene[de1$FDR < 0.05 & de1$logFC > 0]
      down1 <- de1$gene[de1$FDR < 0.05 & de1$logFC < 0]
      up2 <- de2$gene[de2$FDR < 0.05 & de2$logFC > 0]
      down2 <- de2$gene[de2$FDR < 0.05 & de2$logFC < 0]
      list(up_dbs = up1, up_sham = up2, down_dbs = down1, down_sham = down2)
  }
  inn_sets <- get_celltype_venn("InN")
  exn_sets <- get_celltype_venn("ExN")
  add_to_source <- function(ct, direction, dbs_set, sham_set) {
      all_genes <- union(dbs_set, sham_set)
      data.frame(Celltype = ct, Direction = direction, Gene = all_genes, In_DBS_vs_Sham = all_genes %in% 
          dbs_set, In_Sham_vs_Saline = all_genes %in% sham_set, Category = case_when(all_genes %in% dbs_set & 
          all_genes %in% sham_set ~ "Overlap", all_genes %in% dbs_set ~ "DBS_vs_Sham_only", TRUE ~ "Sham_vs_Saline_only"))
  }
  source_df <- bind_rows(add_to_source("InN", "Up", inn_sets$up_dbs, inn_sets$up_sham), add_to_source("InN", 
      "Down", inn_sets$down_dbs, inn_sets$down_sham), add_to_source("ExN", "Up", exn_sets$up_dbs, exn_sets$up_sham), 
      add_to_source("ExN", "Down", exn_sets$down_dbs, exn_sets$down_sham))
  write.csv(source_df, file.path(out_data_dir, "2c.csv"), row.names = FALSE)
  cat("Saved source data: sourcedata/figure 2/2c.csv\n")
  plot_single_venn <- function(dbs_set, sham_set, main_title = NULL, show_labels = FALSE) {
      fit <- euler(list(`DBS vs Sham` = dbs_set, `Sham vs Saline` = sham_set))
      p <- plot(fit, fills = list(fill = c("#ffd700", "#18a799"), alpha = 0.5), edges = list(col = "black", 
          lwd = 1), labels = if (show_labels) 
          list(font = 1, cex = 0.8)
      else FALSE, quantities = list(font = 1, cex = 0.85), main = if (!is.null(main_title)) 
          list(label = main_title, font = 3, cex = 0.9)
      else NULL)
      wrap_elements(p)
  }
  p_inn_up <- plot_single_venn(inn_sets$up_dbs, inn_sets$up_sham, main_title = "Up", show_labels = FALSE)
  p_inn_dn <- plot_single_venn(inn_sets$down_dbs, inn_sets$down_sham, main_title = "Down", show_labels = FALSE)
  p_exn_up <- plot_single_venn(exn_sets$up_dbs, exn_sets$up_sham, main_title = "Up", show_labels = FALSE)
  p_exn_dn <- plot_single_venn(exn_sets$down_dbs, exn_sets$down_sham, main_title = "Down", show_labels = FALSE)
  final_plot <- (p_inn_up/p_inn_dn/p_exn_up/p_exn_dn) + plot_layout(heights = c(1, 1, 1, 1))
  pdf_path <- file.path(out_fig_dir, "2c.pdf")
  png_path <- file.path(out_fig_dir, "2c.png")
  ggsave(png_path, final_plot, width = 2.8, height = 7.5, dpi = 300)
  cat("Saved figures: figures/figure 2/2c.pdf and 2c.png\n")
  invisible(as.list(environment()))
}
