# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(dplyr)
      library(ggplot2)
  })
  fig_dir <- output_path("figures/figure 2")
  data_dir <- output_path("sourcedata/figure 2")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  get_deg_lists <- function(file_path, flip = FALSE, cutoff = 0.05) {
      if (!supplied_inputs[["get_deg_lists_2"]]) 
          stop("Missing DEG input: ", file_path)
      df <- supplied_inputs[["df_3"]][[as.character(file_path)]]
      if (!"FDR" %in% colnames(df)) 
          stop("Missing FDR column: ", file_path)
      gene_col <- if ("X" %in% colnames(df)) 
          "X"
      else if ("gene" %in% colnames(df)) 
          "gene"
      else colnames(df)[1]
      fc_col <- if ("logFC" %in% colnames(df)) 
          "logFC"
      else if ("de_coef" %in% colnames(df)) 
          "de_coef"
      else NULL
      if (is.null(fc_col)) 
          stop("Missing effect column: ", file_path)
      if (flip) 
          df[[fc_col]] <- df[[fc_col]] * -1
      up_genes <- df[[gene_col]][!is.na(df$FDR) & df$FDR < cutoff & df[[fc_col]] > 0]
      dn_genes <- df[[gene_col]][!is.na(df$FDR) & df$FDR < cutoff & df[[fc_col]] < 0]
      return(list(up = up_genes, dn = dn_genes))
  }
  ct_list <- c("ExN", "InN", "Astro", "OPC", "Oligo", "VCs", "Micro")
  overlap_stats <- data.frame()
  for (ct in ct_list) {
      ct_safe <- gsub("/", "-", ct)
      file_pb_1 <- paste0(data_path("results/Table/EdgeR_v2/Level1/EdgeR_DBSIvShamI_"), ct_safe, ".csv")
      file_pb_2 <- paste0(data_path("results/Table/EdgeR_v2/Level1/EdgeR_ShamIvSalineI_"), ct_safe, ".csv")
      file_sc_1 <- paste0(data_path("results/Table/MEMENTO_v2/Level1/MEMENTO_1d_DBSvSham_"), ct_safe, ".csv")
      file_sc_2 <- paste0(data_path("results/Table/MEMENTO_v2/Level1/MEMENTO_1d_ShamvSaline_"), ct_safe, 
          ".csv")
      pb1 <- get_deg_lists(file_pb_1, flip = FALSE)
      pb2 <- get_deg_lists(file_pb_2, flip = TRUE)
      sc1 <- get_deg_lists(file_sc_1, flip = FALSE, cutoff = 0.01)
      sc2 <- get_deg_lists(file_sc_2, flip = TRUE, cutoff = 0.01)
      calc_overlap <- function(setA, setB, modality, direction) {
          intersect_n <- length(intersect(setA, setB))
          union_n <- length(union(setA, setB))
          jaccard <- ifelse(union_n == 0, 0, intersect_n/union_n)
          data.frame(Celltype = ct, Modality = modality, Direction = direction, Category = paste(modality, 
              direction), Overlap_Count = intersect_n, Jaccard_Index = jaccard)
      }
      overlap_stats <- bind_rows(overlap_stats, calc_overlap(sc1$up, sc2$up, "SC Level", "Up"), calc_overlap(sc1$dn, 
          sc2$dn, "SC Level", "Down"), calc_overlap(pb1$up, pb2$up, "PB Level", "Up"), calc_overlap(pb1$dn, 
          pb2$dn, "PB Level", "Down"))
  }
  overlap_stats$Celltype <- factor(overlap_stats$Celltype, levels = rev(ct_list))
  overlap_stats$Category <- factor(overlap_stats$Category, levels = c("SC Level Up", "SC Level Down", "PB Level Up", 
      "PB Level Down"))
  write.csv(overlap_stats, file.path(data_dir, "2b.csv"), row.names = FALSE)
  cat(sprintf("Source data exported: %d records\n", nrow(overlap_stats)))
  overlap_stats$Direction <- factor(overlap_stats$Direction, levels = c("Up", "Down"))
  overlap_stats$Modality <- factor(overlap_stats$Modality, levels = c("SC Level", "PB Level"))
  p <- ggplot(overlap_stats, aes(x = Direction, y = Celltype)) + geom_hline(yintercept = 1:7, colour = "grey85", 
      linetype = "dashed", linewidth = 0.3) + geom_point(aes(size = Overlap_Count, fill = Jaccard_Index), 
      shape = 21, color = "black", stroke = 0.5) + scale_fill_gradientn(colors = c("#f4f9f4", "#7fcdbb", 
      "#2c7fb8", "#253494"), name = "Jaccard\nSimilarity", limits = c(0, max(overlap_stats$Jaccard_Index)), 
      breaks = c(0, 0.4, 0.8)) + scale_size_continuous(range = c(2, 9), name = "Overlapped\nGene Num.", 
      limits = c(0, max(overlap_stats$Overlap_Count)), breaks = c(0, 1000, 2000)) + facet_grid(. ~ Modality, 
      scales = "free_x", space = "free_x", labeller = as_labeller(c(`SC Level` = "SC Level\nMEMENTO\nFDR < 0.05", 
          `PB Level` = "PB Level\nEdgeR\nFDR < 0.05"))) + guides(fill = guide_colourbar(direction = "horizontal", 
      title.position = "top", barwidth = unit(18, "mm")), size = guide_legend(nrow = 1, title.position = "top")) + 
      theme_classic() + labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 9, angle = 45, 
      hjust = 1, color = "black"), axis.text.y = element_text(size = 10, face = "bold", color = "black"), 
      axis.line = element_blank(), axis.ticks = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"), 
      legend.position = "bottom", legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8))
  pdf_out <- file.path(fig_dir, "2b.pdf")
  png_out <- file.path(fig_dir, "2b.png")
  pdf(NULL, width = 3.8, height = 6)
  print(p)
  dev.off()
  png(png_out, width = 1140, height = 1800, res = 300)
  print(p)
  dev.off()
  cat("Panel 2b successfully rendered.\n")
  invisible(as.list(environment()))
}
