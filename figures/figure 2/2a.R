# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
      library(ggplot2)
      library(ggrepel)
  })
  fig_dir <- output_path("figures/figure 2")
  data_dir <- output_path("sourcedata/figure 2")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  seurat_path <- data_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs")
  obj <- supplied_inputs[["obj_2"]][[as.character(seurat_path)]]
  non_neuronal_list <- c("Oligo", "Astro", "Micro", "OPC", "VCs", "EPCs", "CHPCs")
  level4_dict <- c(IsoCortex = "Isocortex", HPF = "HPF", STRd = "CNU", STRv = "CNU", PALd = "CNU", PALc = "CNU", 
      LSX = "CNU", CLA = "CNU", MBsen = "Midbrain", DORpm = "TH")
  obj$brainregion_abc <- unname(case_when(as.character(obj$celltype_level3) %in% non_neuronal_list ~ "Non-neuronal", 
      as.character(obj$brainregion_level4) %in% names(level4_dict) ~ unname(level4_dict[as.character(obj$brainregion_level4)]), 
      as.character(obj$brainregion_level3) == "HY" ~ "HY", as.character(obj$brainregion_level3) == "TH" ~ 
          "TH", as.character(obj$brainregion_level3) == "STR" ~ "CNU", TRUE ~ "Nonspecific"))
  region_map <- obj@meta.data %>% group_by(celltype_level3, brainregion_abc) %>% tally() %>% slice_max(n, 
      n = 1, with_ties = FALSE) %>% select(Celltype = celltype_level3, brainregion_abc)
  colors_brainregion <- c(HY = "#e74538", CNU = "#9bd5f4", Isocortex = "#18a799", HPF = "#84c551", TH = "#e0abce", 
      Midbrain = "#fcc396", `Non-neuronal` = "#000000", Nonspecific = "#cccccc")
  celltype_list_L3 <- c("ExN_THns", "ExN_PF", "ExN_MH", "Oligo", "Astro", "InN_Nonspecific", "ExN_LH", 
      "ExN_CM", "VCs", "ExN_ATN", "InN_HYa", "InN_CP_D1", "ExN_RSP_L6CT", "OPC", "ExN_RE", "InN_GPe_MGE", 
      "Micro", "InN_RT", "InN_CP_D2", "InN_STRns_MGE", "InN_GPe_LGE", "ExN_TRS", "ExN_RSP_L45IT", "InN_RSP_MGE", 
      "InN_STRns_LGE", "InN_LSc", "ExN_DG", "InN_Immature", "ExN_CA3", "ExN_CA1")
  de_res_SS_L3 <- NULL
  for (ct in celltype_list_L3) {
      f1 <- paste0(data_path("results/Table/MEMENTO_v2/Level3/MEMENTO_1d_ShamvSaline_"), gsub("/", "-", 
          ct), ".csv")
      f2 <- paste0(data_path("results/Table/EdgeR_v2/Level3/EdgeR_ShamIvSalineI_"), gsub("/", "-", ct), 
          ".csv")
      if (supplied_inputs[["data_3"]] && supplied_inputs[["data_4"]]) {
          de1 <- supplied_inputs[["de1_5"]][[as.character(f1)]] %>% filter(FDR < 0.01)
          de2 <- supplied_inputs[["de2_6"]][[as.character(f2)]] %>% filter(FDR < 0.05)
          de <- intersect(de1$gene, de2$X)
          de_res_SS_L3 <- rbind(de_res_SS_L3, data.frame(Celltype = ct, DE_Num_SS = length(de)))
      }
  }
  celltype_list_DS <- c("ExN_THns", "ExN_PF", "ExN_MH", "Oligo", "Astro", "InN_Nonspecific", "ExN_LH", 
      "ExN_CM", "VCs", "ExN_ATN", "InN_HYa", "InN_CP_D1", "ExN_RSP_L6CT", "OPC", "ExN_RE", "InN_GPe_MGE", 
      "Micro", "InN_RT", "InN_CP_D2", "InN_STRns_MGE", "ExN_TRS", "ExN_RSP_L45IT", "InN_RSP_MGE", "ExN_CLA", 
      "InN_STRns_LGE", "InN_LSc", "ExN_DG", "ExN_BST")
  de_res_DS_L3 <- NULL
  for (ct in celltype_list_DS) {
      f1 <- paste0(data_path("results/Table/MEMENTO_v2/Level3/MEMENTO_1d_DBSvSham_"), gsub("/", "-", ct), 
          ".csv")
      f2 <- paste0(data_path("results/Table/EdgeR_v2/Level3/EdgeR_DBSIvShamI_"), gsub("/", "-", ct), ".csv")
      if (supplied_inputs[["data_7"]] && supplied_inputs[["data_8"]]) {
          de1 <- supplied_inputs[["de1_9"]][[as.character(f1)]] %>% filter(FDR < 0.01)
          de2 <- supplied_inputs[["de2_10"]][[as.character(f2)]] %>% filter(FDR < 0.05)
          de <- intersect(de1$gene, de2$X)
          de_res_DS_L3 <- rbind(de_res_DS_L3, data.frame(Celltype = ct, DE_Num_DS = length(de)))
      }
  }
  df_merged <- inner_join(de_res_SS_L3, de_res_DS_L3, by = "Celltype") %>% left_join(region_map, by = "Celltype")
  write.csv(df_merged, file.path(data_dir, "2a.csv"), row.names = FALSE)
  cat(sprintf("Source data exported: %d cell types\n", nrow(df_merged)))
  cor_res <- cor.test(df_merged$DE_Num_SS, df_merged$DE_Num_DS, method = "spearman")
  stat_label <- sprintf("Spearman's correlation:\nrho = %.2f, P = %.2e", cor_res$estimate, cor_res$p.value)
  key_labels <- c("ExN_ATN", "ExN_RSP_L45IT", "InN_HYa", "ExN_PF", "InN_Nonspecific", "VCs")
  df_label <- df_merged %>% filter(Celltype %in% key_labels)
  p <- ggplot(df_merged, aes(x = DE_Num_SS, y = DE_Num_DS)) + geom_point(aes(color = brainregion_abc), 
      size = 3.5, alpha = 0.95) + geom_smooth(method = "lm", se = TRUE, color = "#aaaaaa", fill = "#dddddd", 
      alpha = 0.2, linewidth = 0.8, linetype = "dashed") + scale_color_manual(values = colors_brainregion) + 
      annotate("text", x = 100, y = 3100, label = stat_label, hjust = 0, vjust = 1, size = 3.8, fontface = "plain") + 
      geom_text_repel(data = df_label, aes(label = Celltype), size = 3.5, box.padding = 0.5, point.padding = 0.3, 
          segment.color = "grey60", segment.size = 0.4, min.segment.length = 0) + labs(x = "DEG Num. (Sham vs Saline)", 
      y = "DEG Num. (DBS vs Sham)", color = "Brain Region") + theme_classic() + theme(axis.title = element_text(size = 11, 
      face = "bold", color = "black"), axis.text = element_text(size = 10, color = "black"), axis.line = element_line(linewidth = 0.5), 
      legend.position = "right")
  pdf_out <- file.path(fig_dir, "2a.pdf")
  png_out <- file.path(fig_dir, "2a.png")
  pdf(NULL, width = 6.5, height = 5.5)
  print(p)
  dev.off()
  png(png_out, width = 1950, height = 1650, res = 300)
  print(p)
  dev.off()
  cat("Panel 2a successfully rendered.\n")
  invisible(as.list(environment()))
}
