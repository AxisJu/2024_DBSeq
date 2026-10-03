# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(Seurat)
      library(Matrix)
      library(patchwork)
      library(scales)
      library(cowplot)
      library(scattermore)
      library(ggh4x)
  })
  raw_file <- data_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs")
  if (Sys.getenv("DBSEQ_REBUILD_SOURCE") == "1" || !supplied_inputs[["data_2"]]) {
      obj <- supplied_inputs[["obj_3"]][[as.character(raw_file)]]
      obj <- subset(obj, cells = colnames(obj)[!obj$donor_id %in% c("unassigned", "doublet")])
      target_subtypes <- c("Astro", "EPCs", "OPC", "Oligo", "VCs", "CHPCs", "Micro", "ExN_RSP_L23IT", "ExN_RSP_L45IT", 
          "ExN_CLA", "ExN_CA3", "ExN_CA1", "ExN_RSP_L6CT", "ExN_HPF_CajalRetzius", "ExN_DG", "InN_Immature", 
          "InN_RSP_MGE", "InN_STRns_MGE", "InN_GPe_MGE", "InN_GPi_MGE", "InN_Nonspecific", "InN_STRns_LGE", 
          "InN_GPe_LGE", "InN_CP_D1", "InN_CP_D2", "InN_OT", "InN_LSc", "InN_HYa", "InN_RT", "ExN_HYa", 
          "ExN_TRS", "ExN_BST", "ExN_MH", "ExN_LH", "ExN_THns", "ExN_ATN", "ExN_RE", "ExN_PF", "ExN_CM", 
          "InN_SCsg")
      target_subtypes <- intersect(target_subtypes, unique(as.character(obj$celltype_level3)))
      marker_map <- c(Astro = "Aqp4", EPCs = "Foxj1", OPC = "Pdgfra", Oligo = "Mbp", VCs = "Cldn5", CHPCs = "Ttr", 
          Micro = "Cx3cr1", ExN_RSP_L23IT = "Tshz2", ExN_RSP_L45IT = "Phactr1", ExN_CLA = "Ptprt", ExN_CA3 = "Zbtb20", 
          ExN_CA1 = "Epha6", ExN_RSP_L6CT = "Hs3st4", ExN_HPF_CajalRetzius = "Reln", ExN_DG = "Ppp3ca", 
          InN_Immature = "Sox2ot", InN_RSP_MGE = "Zeb2", InN_STRns_MGE = "Zfp536", InN_GPe_MGE = "Cntnap5a", 
          InN_GPi_MGE = "Elavl2", InN_Nonspecific = "Gm42418", InN_STRns_LGE = "Foxp2", InN_GPe_LGE = "Slc44a5", 
          InN_CP_D1 = "Kcnq5", InN_CP_D2 = "Foxp1", InN_OT = "Sgcd", InN_LSc = "Slc8a1", InN_HYa = "Dscam", 
          InN_RT = "Fign", ExN_HYa = "Shisa6", ExN_TRS = "Asic2", ExN_BST = "Ebf2", ExN_MH = "Nwd2", ExN_LH = "Nlgn1", 
          ExN_THns = "Pcp4", ExN_ATN = "Kcnc2", ExN_RE = "Nell1", ExN_PF = "Lypd6b", ExN_CM = "Necab1", 
          InN_SCsg = "Otx2os1")
      selected_features <- unname(marker_map[target_subtypes])
      stopifnot(length(selected_features) == 40, all(selected_features %in% rownames(obj)))
      Idents(obj) <- "celltype_level3"
      d <- DotPlot(obj, features = selected_features, idents = target_subtypes)$data
      d$marker_order <- match(d$features.plot, selected_features)
      d$marker_for <- target_subtypes[d$marker_order]
      d$cell_count <- as.integer(table(obj$celltype_level3)[as.character(d$id)])
      d$displayed_bubble <- d$pct.exp >= 3
      names(d)[names(d) == "features.plot"] <- "gene"
      names(d)[names(d) == "id"] <- "celltype"
      names(d)[names(d) == "avg.exp"] <- "mean_linear_normalized_expression"
      names(d)[names(d) == "pct.exp"] <- "percent_expressed"
      names(d)[names(d) == "avg.exp.scaled"] <- "scaled_mean_log_expression"
      export_data(d, "ed2b")
  }
  d <- supplied_inputs[["d_4"]][[as.character(file.path(source_dir, "ed2b.csv"))]]
  selected_features <- unique(d$gene[order(d$marker_order)])
  target_subtypes <- unique(d$marker_for[order(d$marker_order)])
  df_plot <- d
  df_plot$features.plot <- df_plot$gene
  df_plot$id <- df_plot$celltype
  df_plot$pct_display <- ifelse(df_plot$percent_expressed < 3, NA, df_plot$percent_expressed)
  df_plot$avg.exp.scaled <- df_plot$scaled_mean_log_expression
  p_bubble <- ggplot(df_plot, aes(x = features.plot, y = id)) + geom_point(data = subset(df_plot, is.na(pct_display)), 
      aes(x = features.plot, y = id), color = "grey93", size = 0.25, shape = 16) + geom_point(data = subset(df_plot, 
      !is.na(pct_display)), aes(size = pct_display, fill = avg.exp.scaled), shape = 21, color = "black", 
      stroke = 0.35) + scale_x_discrete(limits = selected_features) + scale_y_discrete(limits = rev(target_subtypes), 
      labels = function(x) sub("ExN_HPF_CajalRetzius", "ExN_HPF_CR", x, fixed = TRUE)) + scale_fill_gradientn(colors = c("#f4f9f4", 
      "#7fcdbb", "#2c7fb8", "#253494"), name = "Mean Scaled\nExpression", limits = c(-2, 2.5), oob = squish) + 
      scale_size_continuous(range = c(0.3, 2.2), limits = c(0, 100), breaks = c(25, 50, 75, 100), name = "Percent\nExpressed (%)") + 
      labs(x = "Subtype Marker Genes", y = "Cell Subtype (Level 3)") + theme_classic(base_family = "Arial") + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, face = "italic", size = 5, color = "black"), 
          axis.text.y = element_text(face = "bold", size = 5, color = "black"), axis.title = element_text(size = 5, 
              face = "bold", color = "black"), axis.line = element_line(linewidth = 0.5, color = "black"), 
          axis.ticks = element_line(linewidth = 0.5, color = "black"), axis.ticks.length = unit(1, "mm"), 
          legend.title = element_text(size = 5, face = "bold", color = "black"), legend.text = element_text(size = 5, 
              color = "black"), legend.key.size = unit(2.8, "mm"), legend.spacing.y = unit(1, "mm"), panel.grid.major = element_line(color = "grey94", 
              linewidth = 0.25, linetype = "dashed"), plot.margin = margin(t = 3, r = 4, b = 3, l = 4, 
              unit = "mm"))
  save_plot(p_bubble, "ed2b", 170, 110)
  invisible(as.list(environment()))
}
