# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(RColorBrewer)
  library(scales)
  cache_file <- cache_path("enrichr_cache.qs")
  if (!supplied_inputs[["data_2"]] || Sys.getenv("DBSEQ_REBUILD_SOURCE") == "1") {
      list2env(supplied_inputs[["helpers_3"]], envir = environment())
      prepare_4_ora()
  }
  enrich_cache <- supplied_inputs[["enrich_cache_4"]][[as.character(cache_file)]]
  res_atn_go <- enrich_cache$atn
  res_com_go <- enrich_cache$common
  top1_atn <- "Cell-Cell Junction Assembly (GO:0007043)"
  top5_com <- c("RNA Processing (GO:0006396)", "RNA Splicing (GO:0008380)", "mRNA Splicing, Via Spliceosome (GO:0000398)", 
      "Regulation Of mRNA Splicing, Via Spliceosome (GO:0048024)", "RNA Splicing, Via Transesterification Reactions With Bulged Adenosine As Nucleophile (GO:0000377)")
  top5_shared <- c("Golgi Vesicle Transport (GO:0048193)", "Endoplasmic Reticulum To Golgi Vesicle-Mediated Transport (GO:0006888)", 
      "Intracellular Protein Transport (GO:0006886)", "intra-Golgi Vesicle-Mediated Transport (GO:0006891)", 
      "Golgi Organization (GO:0007030)")
  target_terms <- c(top1_atn, top5_com, top5_shared)
  clean_label_dict <- c(`Cell-Cell Junction Assembly (GO:0007043)` = "Cell-Cell Junction Assembly", `RNA Processing (GO:0006396)` = "RNA Processing", 
      `RNA Splicing (GO:0008380)` = "RNA Splicing", `mRNA Splicing, Via Spliceosome (GO:0000398)` = "mRNA Splicing,\nVia Spliceosome", 
      `Regulation Of mRNA Splicing, Via Spliceosome (GO:0048024)` = "Regulation Of mRNA Splicing,\nVia Spliceosome", 
      `RNA Splicing, Via Transesterification Reactions With Bulged Adenosine As Nucleophile (GO:0000377)` = "RNA Splicing,\nVia Transesterification Reactions", 
      `Golgi Vesicle Transport (GO:0048193)` = "Golgi Vesicle Transport", `Endoplasmic Reticulum To Golgi Vesicle-Mediated Transport (GO:0006888)` = "Endoplasmic Reticulum to Golgi\nVesicle-Mediated Transport", 
      `Intracellular Protein Transport (GO:0006886)` = "Intracellular Protein Transport", `intra-Golgi Vesicle-Mediated Transport (GO:0006891)` = "intra-Golgi\nVesicle-Mediated Transport", 
      `Golgi Organization (GO:0007030)` = "Golgi Organization")
  df_bubble <- bind_rows(res_atn_go %>% filter(Term %in% target_terms) %>% mutate(TRG_Group = "ATN\nSpecific"), 
      res_com_go %>% filter(Term %in% target_terms) %>% mutate(TRG_Group = "Neuronal\nCommon")) %>% mutate(Count = as.numeric(sub("/.*", 
      "", Overlap)), Clean_Term = unname(clean_label_dict[Term]), Pathway_Type = case_when(Term %in% top1_atn ~ 
      "ATN", Term %in% top5_com ~ "Neuronal-Common", Term %in% top5_shared ~ "Shared"), Sig_Stars = case_when(Adjusted.P.value < 
      0.001 ~ "***", Adjusted.P.value < 0.01 ~ "**", Adjusted.P.value < 0.05 ~ "*", TRUE ~ ""))
  ordered_terms <- rev(unname(clean_label_dict[target_terms]))
  df_bubble <- df_bubble %>% mutate(Clean_Term = factor(Clean_Term, levels = ordered_terms), TRG_Group = factor(TRG_Group, 
      levels = c("ATN\nSpecific", "Neuronal\nCommon")), Pathway_Type = factor(Pathway_Type, levels = c("ATN", 
      "Neuronal-Common", "Shared")))
  ylgnbu_cols <- rev(brewer.pal(9, "YlGnBu")[3:7])
  p_bubble <- ggplot(df_bubble, aes(x = TRG_Group, y = Clean_Term)) + geom_point(aes(size = Count, fill = Adjusted.P.value), 
      shape = 21, color = "black", stroke = 0.8) + geom_text(aes(label = Sig_Stars), color = "white", fontface = "bold", 
      size = 3.3, vjust = 0.75) + labs(caption = "Stars: BH FDR < 0.05 (*), 0.01 (**), 0.001 (***)\nColour: FDR") + 
      facet_grid(Pathway_Type ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_fill_gradientn(colors = ylgnbu_cols, 
      limits = c(0, 0.05), breaks = c(0, 0.05), labels = c("0", ">=0.05"), oob = scales::squish, na.value = rev(brewer.pal(9, 
          "YlGnBu")[3:7])[5], name = "FDR", guide = guide_colorbar(title.position = "top", title.theme = element_text(face = "bold", 
          size = 10), barwidth = unit(2.8, "cm"), barheight = unit(0.35, "cm"))) + scale_size_continuous(range = c(4.5, 
      12.5), limits = c(0, 95), breaks = c(30, 60, 90), name = "Counts", guide = guide_legend(title.position = "top", 
      title.theme = element_text(face = "bold", size = 10), override.aes = list(fill = "white", color = "black", 
          stroke = 0.8))) + scale_x_discrete(expand = expansion(mult = c(0.4, 0.4))) + theme_bw(base_size = 11) + 
      theme(strip.placement = "outside", strip.background = element_rect(fill = "white", color = "black", 
          linewidth = 0.8), strip.text.y.left = element_text(angle = 90, face = "bold", color = "black", 
          size = 10), axis.text.y = element_text(face = "bold", color = "black", size = 9.5, hjust = 1), 
          axis.text.x = element_text(face = "bold", color = "black", size = 10), axis.title = element_blank(), 
          axis.ticks = element_line(color = "black", linewidth = 0.5), panel.border = element_rect(color = "black", 
              fill = NA, linewidth = 0.8), panel.grid.major = element_line(color = "grey92", linewidth = 0.4), 
          panel.grid.minor = element_blank(), legend.position = "bottom", legend.box = "vertical", plot.caption = element_text(size = 7, 
              hjust = 0.5), legend.spacing.x = unit(0.6, "cm"), plot.margin = margin(t = 10, r = 10, b = 10, 
              l = 10))
  out_fig_dir <- output_path("figures/figure 4")
  out_data_dir <- output_path("sourcedata/figure 4")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_fig_dir, "4b.png"), plot = p_bubble, width = 4.8, height = 7.2, dpi = 600)
  df_source <- df_bubble %>% select(Pathway_Type, Term, Clean_Term, TRG_Group, Count, Overlap, P.value, 
      Adjusted.P.value, Sig_Stars, Genes) %>% rename(Pathway_Category = Pathway_Type, Pathway_Term = Term, 
      Pathway_Name = Clean_Term, TRG_Subset = TRG_Group, Enriched_Gene_Count = Count, Pathway_Overlap = Overlap, 
      P_Value = P.value, FDR = Adjusted.P.value, Significance = Sig_Stars, Enriched_Genes = Genes)
  df_source$Significance_Basis <- "Full-library BH FDR within each TRG subset"
  write.csv(df_source, output_path("statistics/figure 4/4b_enrichment.csv"), row.names = FALSE)
  write.csv(df_source %>% select(-P_Value, -FDR, -Significance, -Significance_Basis), file.path(out_data_dir, 
      "4b.csv"), row.names = FALSE)
  cat("Figure 4b generated successfully.\n")
  invisible(as.list(environment()))
}
