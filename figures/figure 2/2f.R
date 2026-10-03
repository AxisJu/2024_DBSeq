# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggrastr)
  proj_dir <- data_root
  out_fig_dir <- output_path("figures/figure 2")
  out_data_dir <- output_path("sourcedata/figure 2")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  colors_celltype_level3 <- c(Astro = "#a58946", EPCs = "#594a26", OPC = "#5953ff", Oligo = "#201e5a", 
      VCs = "#858881", CHPCs = "#cfd4c9", Micro = "#a87c5a", ExN_RSP_L23IT = "#ff1a71", ExN_RSP_L45IT = "#ed1986", 
      ExN_CLA = "#ba1369", ExN_CA3 = "#860d4c", ExN_CA1 = "#53082f", ExN_RSP_L6CT = "#61e2a4", ExN_HPF_CajalRetzius = "#d00000", 
      ExN_DG = "#16f2f2", InN_Immature = "#1b4332", InN_RSP_MGE = "#f954ee", InN_STRns_MGE = "#70e000", 
      InN_GPe_MGE = "#56ad00", InN_GPi_MGE = "#2f6000", InN_Nonspecific = "#f0a0ff", InN_STRns_LGE = "#d899ff", 
      InN_GPe_LGE = "#b199ff", InN_CP_D1 = "#8d7acc", InN_CP_D2 = "#695b99", InN_OT = "#4e4372", InN_LSc = "#3283fe", 
      InN_HYa = "#450099", InN_RT = "#ff6600", ExN_HYa = "#aa0dfe", ExN_TRS = "#7609b1", ExN_BST = "#420564", 
      ExN_MH = "#faa307", ExN_LH = "#c68105", ExN_THns = "#0d47a1", ExN_ATN = "#1460ff", ExN_RE = "#1340ff", 
      ExN_PF = "#0a4093", ExN_CM = "#08306d", InN_SCsg = "#9ef01a")
  plot_df <- supplied_inputs[["plot_df_2"]]
  stopifnot(sum(grepl("^median_cos", names(plot_df))) == 1)
  names(plot_df)[grepl("^median_cos", names(plot_df))] <- "median_costheta"
  if (!"x_pos" %in% names(plot_df)) {
      types <- unique(as.character(plot_df$Celltype))
      plot_df$Celltype <- factor(plot_df$Celltype, levels = types)
      plot_df <- plot_df %>% group_by(Celltype) %>% arrange(genes, .by_group = TRUE) %>% mutate(x_pos = as.numeric(Celltype) + 
          (row_number() - (n() + 1)/2)/max(n(), 1) * 0.7) %>% ungroup()
  }
  source_f <- plot_df %>% select(Gene = genes, Celltype, CoD_DBSvSham = CoD_DBSvSham_cohend, CoD_ShamvSaline = CoD_ShamvSal_cohend, 
      Delta_CohenD = delta_CohenD, Median_CosTheta = median_costheta, Is_TRG = is_TRG)
  write.csv(source_f, file.path(out_data_dir, "2f.csv"), row.names = FALSE)
  cat("Saved source data: sourcedata/figure 2/2f.csv\n")
  top5_trgs <- plot_df %>% filter(is_TRG == TRUE) %>% group_by(Celltype) %>% slice_max(order_by = median_costheta, 
      n = 5, with_ties = FALSE) %>% ungroup()
  valid_celltype_order <- levels(plot_df$Celltype)
  axis_breaks <- 1:length(valid_celltype_order)
  p <- ggplot(plot_df, aes(x = x_pos, y = delta_CohenD)) + geom_hline(yintercept = 0, color = "grey60", 
      linetype = "solid", linewidth = 0.5) + rasterise(geom_point(data = filter(plot_df, !is_TRG), color = "grey80", 
      size = 0.01, alpha = 0.1), dpi = 300) + geom_point(data = filter(plot_df, is_TRG), aes(color = Celltype), 
      size = 0.1, alpha = 0.9) + geom_text_repel(data = top5_trgs, aes(label = genes, color = Celltype), 
      size = 3, box.padding = 1, point.padding = 0.1, segment.color = "grey50", segment.size = 0.3, max.overlaps = 50, 
      min.segment.length = 0, show.legend = FALSE) + scale_color_manual(values = colors_celltype_level3) + 
      scale_x_continuous(breaks = axis_breaks, labels = valid_celltype_order, expand = expansion(add = c(0.6, 
          0.6))) + labs(x = "Cell Types", y = expression(paste(Delta, " Cohen's d (Sham vs Saline - DBS vs Sham)")), 
      title = "Transcriptional Shifts and TRGs across Cell Types") + theme_classic(base_size = 14) + theme(axis.text.x = element_text(angle = 45, 
      hjust = 1, vjust = 1, size = 10, color = "black"), axis.text.y = element_text(color = "black"), axis.ticks.x = element_blank(), 
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
      legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.margin = margin(t = 10, 
          r = 20, b = 10, l = 10))
  pdf_path <- file.path(out_fig_dir, "2f.pdf")
  png_path <- file.path(out_fig_dir, "2f.png")
  ggsave(png_path, p, width = 21.3, height = 6.7, dpi = 300)
  cat("Saved figures: figures/figure 2/2f.pdf and 2f.png\n")
  invisible(as.list(environment()))
}
