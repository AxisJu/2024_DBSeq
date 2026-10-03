# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(stringr)
  library(ggtext)
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
  df_atn_sub <- res_atn_go %>% filter(Term %in% target_terms) %>% mutate(Count_ATN = as.numeric(sub("/.*", 
      "", Overlap)), Genes_ATN = Genes) %>% select(Term, Count_ATN, Genes_ATN)
  df_com_sub <- res_com_go %>% filter(Term %in% target_terms) %>% mutate(Count_Common = as.numeric(sub("/.*", 
      "", Overlap)), Genes_Common = Genes) %>% select(Term, Count_Common, Genes_Common)
  df_merged <- tibble(Term = target_terms) %>% left_join(df_atn_sub, by = "Term") %>% left_join(df_com_sub, 
      by = "Term") %>% mutate(Count_ATN = ifelse(is.na(Count_ATN), 0, Count_ATN), Count_Common = ifelse(is.na(Count_Common), 
      0, Count_Common), Total_Count = Count_ATN + Count_Common, Clean_Term = unname(clean_label_dict[Term]), 
      Pathway_Type = case_when(Term %in% top1_atn ~ "ATN", Term %in% top5_com ~ "Neuronal-Common", Term %in% 
          top5_shared ~ "Shared"))
  format_atn_genes <- function(genes_str, width = 65) {
      if (is.na(genes_str) || genes_str == "") 
          return("")
      genes <- unlist(strsplit(genes_str, ";"))
      genes_title <- genes
      genes_formatted <- sapply(genes_title, function(g) {
          if (g == "Zbtb7a") {
              "<span style='color:#E64B35;'><b><i>Zbtb7a</i></b></span>"
          }
          else {
              paste0("<i>", g, "</i>")
          }
      })
      raw_text <- paste(genes_formatted, collapse = ", ")
      if (nchar(paste(genes_title, collapse = ", ")) > width) {
          mid <- floor(length(genes_formatted)/2)
          line1 <- paste(genes_formatted[1:mid], collapse = ", ")
          line2 <- paste(genes_formatted[(mid + 1):length(genes_formatted)], collapse = ", ")
          return(paste0(line1, ",<br>", line2))
      }
      return(raw_text)
  }
  df_merged <- df_merged %>% mutate(Gene_Text_Rich = sapply(Genes_ATN, format_atn_genes))
  df_bar_long <- df_merged %>% tidyr::pivot_longer(cols = c(Count_Common, Count_ATN), names_to = "TRG_Group", 
      values_to = "Count") %>% mutate(TRG_Group = factor(recode(TRG_Group, Count_Common = "Neuronal-Common", 
      Count_ATN = "ATN-Specific"), levels = c("ATN-Specific", "Neuronal-Common")))
  ordered_terms <- rev(unname(clean_label_dict[target_terms]))
  df_bar_long <- df_bar_long %>% mutate(Clean_Term = factor(Clean_Term, levels = ordered_terms), Pathway_Type = factor(Pathway_Type, 
      levels = c("ATN", "Neuronal-Common", "Shared")))
  df_labels <- df_merged %>% mutate(Clean_Term = factor(Clean_Term, levels = ordered_terms), Pathway_Type = factor(Pathway_Type, 
      levels = c("ATN", "Neuronal-Common", "Shared")))
  p_bar <- ggplot() + geom_bar(data = df_bar_long, aes(x = Clean_Term, y = Count, fill = TRG_Group), stat = "identity", 
      color = "black", linewidth = 0.5, width = 0.65) + geom_richtext(data = df_labels, aes(x = Clean_Term, 
      y = Total_Count + 2, label = Gene_Text_Rich), hjust = 0, size = 2.65, lineheight = 1.05, fill = NA, 
      label.color = NA, label.padding = unit(0, "lines")) + coord_flip() + facet_grid(Pathway_Type ~ ., 
      scales = "free_y", space = "free_y") + scale_fill_manual(values = c(`Neuronal-Common` = "#56B4E9", 
      `ATN-Specific` = "#E64B35"), breaks = c("Neuronal-Common", "ATN-Specific")) + scale_y_continuous(limits = c(0, 
      115), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0.05))) + theme_bw(base_size = 11) + 
      theme(strip.background = element_blank(), strip.text = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(face = "bold", 
              color = "black", size = 10), axis.title.x = element_text(face = "bold", color = "black", 
              size = 11, margin = margin(t = 8)), panel.border = element_rect(color = "black", fill = NA, 
              linewidth = 0.8), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4), legend.position = "bottom", 
          legend.title = element_blank(), legend.text = element_text(face = "bold", size = 10), legend.key.size = unit(0.4, 
              "cm"), plot.margin = margin(t = 10, r = 25, b = 10, l = 0)) + labs(y = "Num. of enriched TRGs")
  out_fig_dir <- output_path("figures/figure 4")
  out_data_dir <- output_path("sourcedata/figure 4")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_fig_dir, "4c.png"), plot = p_bar, width = 6.8, height = 7.2, dpi = 600)
  df_source_c <- df_merged %>% select(Pathway_Type, Term, Clean_Term, Count_Common, Count_ATN, Total_Count, 
      Genes_ATN, Genes_Common) %>% rename(Pathway_Category = Pathway_Type, Pathway_Term = Term, Pathway_Name = Clean_Term, 
      Neuronal_Common_Count = Count_Common, ATN_Specific_Count = Count_ATN, Total_Enriched_Count = Total_Count, 
      ATN_Specific_Genes = Genes_ATN, Neuronal_Common_Genes = Genes_Common)
  write.csv(df_source_c, file.path(out_data_dir, "4c.csv"), row.names = FALSE)
  cat("Figure 4c generated successfully.\n")
  invisible(as.list(environment()))
}
