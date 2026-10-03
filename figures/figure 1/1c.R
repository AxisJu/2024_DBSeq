# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(Seurat)
      library(ggplot2)
      library(dplyr)
      library(patchwork)
      library(RColorBrewer)
  })
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
  script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/")
  package_root <- output_root
  data_root <- Sys.getenv("DBSEQ_DATA_ROOT")
  if (!nzchar(data_root)) stop("Set DBSEQ_DATA_ROOT to the original analysis root.")
  fig_dir <- file.path(package_root, "figures", "figure 1")
  data_dir <- file.path(package_root, "sourcedata", "figure 1")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  seurat_path <- file.path(data_root, "data/snRNAseq_mouse/processed/matrix/running_all_250704.qs")
  obj <- supplied_inputs[["obj_2"]][[as.character(seurat_path)]]
  meta_df <- obj@meta.data
  meta_df$cell_id <- rownames(meta_df)
  stopifnot("Group_L2" %in% names(meta_df))
  meta_df$Group_display <- factor(case_when(as.character(meta_df$Group_L2) %in% c("DBS", "DBS_I", "DBS_C") ~ 
      "DBS", as.character(meta_df$Group_L2) %in% c("Sham", "Sham_I") ~ "Sham", as.character(meta_df$Group_L2) %in% 
      c("Saline", "Saline_I") ~ "Saline", TRUE ~ "Others"), levels = c("DBS", "Sham", "Saline", "Others"))
  level4_dict <- c(IsoCortex = "Isocortex", HPF = "HPF", STRd = "CNU", STRv = "CNU", PALd = "CNU", PALc = "CNU", 
      LSX = "CNU", CLA = "CNU", MBsen = "Midbrain", DORpm = "TH")
  meta_df$brainregion_abc <- unname(case_when(as.character(meta_df$brainregion_level4) %in% names(level4_dict) ~ 
      unname(level4_dict[as.character(meta_df$brainregion_level4)]), as.character(meta_df$brainregion_level3) == 
      "HY" ~ "HY", as.character(meta_df$brainregion_level3) == "TH" ~ "TH", as.character(meta_df$brainregion_level3) == 
      "STR" ~ "CNU", TRUE ~ "Nonspecific"))
  colors_celltype_level3 <- c(Astro = "#a58946", EPCs = "#594a26", OPC = "#5953ff", Oligo = "#201e5a", 
      VCs = "#858881", CHPCs = "#cfd4c9", Micro = "#a87c5a", ExN_RSP_L23IT = "#ff1a71", ExN_RSP_L45IT = "#ed1986", 
      ExN_CLA = "#ba1369", ExN_CA3 = "#860d4c", ExN_CA1 = "#53082f", ExN_RSP_L6CT = "#61e2a4", ExN_HPF_CajalRetzius = "#d00000", 
      ExN_DG = "#16f2f2", InN_Immature = "#1b4332", InN_RSP_MGE = "#f954ee", InN_STRns_MGE = "#70e000", 
      InN_GPe_MGE = "#56ad00", InN_GPi_MGE = "#2f6000", InN_Nonspecific = "#f0a0ff", InN_STRns_LGE = "#d899ff", 
      InN_GPe_LGE = "#b199ff", InN_CP_D1 = "#8d7acc", InN_CP_D2 = "#695b99", InN_OT = "#4e4372", InN_LSc = "#3283fe", 
      InN_HYa = "#450099", InN_RT = "#ff6600", ExN_HYa = "#aa0dfe", ExN_TRS = "#7609b1", ExN_BST = "#420564", 
      ExN_MH = "#faa307", ExN_LH = "#c68105", ExN_THns = "#0d47a1", ExN_ATN = "#1460ff", ExN_RE = "#1340ff", 
      ExN_PF = "#0a4093", ExN_CM = "#08306d", InN_SCsg = "#9ef01a")
  colors_mapmycells <- c(`01 IT-ET Glut` = "#fa0087", `02 NP-CT-L6b Glut` = "#61e2a4", `03 OB-CR Glut` = "#a31515", 
      `04 DG-IMN Glut` = "#16f2f2", `05 OB-IMN GABA` = "#1b4332", `06 CTX-CGE GABA` = "#ccff33", `07 CTX-MGE GABA` = "#f954ee", 
      `08 CNU-MGE GABA` = "#70e000", `09 CNU-LGE GABA` = "#b199ff", `10 LSX GABA` = "#3283fe", `11 CNU-HYa GABA` = "#450099", 
      `12 HY GABA` = "#ff6600", `17 MH-LH Glut` = "#faa307", `18 TH Glut` = "#0d47a1", `19 MB Glut` = "#007200", 
      `20 MB GABA` = "#9ef01a", `30 Astro-Ependymal` = "#594a26", `31 OPC-Oligo` = "#03045e", `33 Vascular` = "#858881", 
      `34 Immune` = "#825f45", `Low_confidence (r<=0.5)` = "#f0f0f0", Others = "#dddddd")
  colors_brainregion <- c(HY = "#e74538", CNU = "#9bd5f4", Isocortex = "#18a799", HPF = "#84c551", TH = "#e0abce", 
      Midbrain = "#fcc396", Nonspecific = "#cccccc")
  colors_group <- c(DBS = "#ffd700", Sham = "#ac00cc", Saline = "#18a799", Others = "#dddddd")
  level3_order <- names(colors_celltype_level3)
  cluster_dominant <- meta_df %>% group_by(leiden_res_2.0) %>% count(celltype_level3) %>% slice_max(n, 
      n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(celltype_level3 = factor(celltype_level3, levels = level3_order), 
      leiden_num = as.numeric(as.character(leiden_res_2.0))) %>% arrange(celltype_level3, leiden_num)
  meta_df$leiden_res_2.0 <- factor(meta_df$leiden_res_2.0, levels = as.character(cluster_dominant$leiden_res_2.0))
  map_csv_path <- file.path(data_root, "data/snRNAseq_mouse/processed/intermediate/MapMyCells/merged_preprocessed_10xWholeMouseBrain(CCN20230722)_CorrelationMapping_UTC_1740584377620/merged_preprocessed_10xWholeMouseBrain(CCN20230722)_CorrelationMapping_UTC_1740584377620.csv")
  map_df <- supplied_inputs[["map_df_3"]][[as.character(map_csv_path)]]
  meta_df <- left_join(meta_df, map_df[, c("cell_id", "class_name", "class_correlation_coefficient")], 
      by = "cell_id")
  valid_prefixes <- c(sprintf("%02d", 1:12), sprintf("%02d", 17:20), sprintf("%02d", 30:34))
  meta_df <- meta_df %>% mutate(prefix = substr(class_name, 1, 2), step_a_class = case_when(is.na(class_correlation_coefficient) | 
      class_correlation_coefficient <= 0.5 ~ "Low_confidence (r<=0.5)", !(prefix %in% valid_prefixes) ~ 
      "Others", TRUE ~ class_name))
  valid_cluster_props <- meta_df %>% filter(!step_a_class %in% c("Low_confidence (r<=0.5)", "Others")) %>% 
      group_by(leiden_res_2.0, step_a_class) %>% tally() %>% group_by(leiden_res_2.0) %>% mutate(prop = n/sum(n)) %>% 
      filter(prop > 0.1) %>% select(leiden_res_2.0, step_a_class) %>% mutate(keep = TRUE)
  meta_df <- left_join(meta_df, valid_cluster_props, by = c("leiden_res_2.0", "step_a_class")) %>% mutate(final_map_class = case_when(step_a_class == 
      "Low_confidence (r<=0.5)" ~ "Low_confidence (r<=0.5)", step_a_class == "Others" ~ "Others", is.na(keep) ~ 
      "Others", TRUE ~ step_a_class))
  df_source <- meta_df %>% select(cell_id, leiden_res_2.0, celltype_level3, final_map_class, brainregion_abc, 
      Group_L1, Group_L2, Group_display)
  write.csv(df_source, file.path(data_dir, "1c.csv"), row.names = FALSE)
  cat(sprintf("Source data exported: %d cells\n", nrow(df_source)))
  theme_nature <- theme_classic() + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 9, 
          face = "bold"), legend.position = "right", legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 7), 
      plot.margin = margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"))
  p1 <- ggplot(meta_df, aes(x = leiden_res_2.0, fill = celltype_level3)) + geom_bar(position = "stack", 
      width = 0.9) + scale_fill_manual(values = colors_celltype_level3) + scale_y_continuous(expand = expansion(mult = c(0, 
      0.05))) + labs(y = "Number of\nCells") + theme_nature
  p2 <- ggplot(meta_df, aes(x = leiden_res_2.0, fill = final_map_class)) + geom_bar(position = "fill", 
      width = 0.9) + scale_fill_manual(values = colors_mapmycells) + scale_y_continuous(expand = c(0, 0)) + 
      labs(y = "ABC Atlas Class\nFraction") + theme_nature
  p3 <- ggplot(meta_df, aes(x = leiden_res_2.0, fill = brainregion_abc)) + geom_bar(position = "fill", 
      width = 0.9) + scale_fill_manual(values = colors_brainregion) + scale_y_continuous(expand = c(0, 
      0)) + labs(y = "Major Brain Region\nFraction") + theme_nature
  p4 <- ggplot(meta_df, aes(x = leiden_res_2.0, fill = Group_display)) + geom_bar(position = "fill", width = 0.9) + 
      scale_fill_manual(values = colors_group, drop = FALSE, name = "Group") + scale_y_continuous(expand = c(0, 
      0)) + labs(y = "Group\nFraction", x = "Cluster (leiden_res_2.0)") + theme_nature
  final_plot <- p1/p2/p3/p4 + plot_layout(heights = c(2, 2, 2, 0.7), guides = "collect")
  pdf_out <- file.path(fig_dir, "1c.pdf")
  png_out <- file.path(fig_dir, "1c.png")
  pdf(NULL, width = 12, height = 7)
  print(final_plot)
  dev.off()
  png(png_out, width = 3600, height = 2100, res = 300)
  print(final_plot)
  dev.off()
  cat("Panel 1c successfully rendered.\n")
  invisible(as.list(environment()))
}
