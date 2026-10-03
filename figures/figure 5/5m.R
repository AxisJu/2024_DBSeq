# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  dir.create(output_path("figures/figure 5"), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_path("sourcedata/figure 5"), recursive = TRUE, showWarnings = FALSE)
  plot_df <- supplied_inputs[["plot_df_2"]]
  mouse_non_neuron <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  mouse_neuron <- setdiff(unique(plot_df$Celltype), mouse_non_neuron)
  trg_class <- plot_df %>% filter(is_TRG == TRUE, Celltype %in% mouse_neuron) %>% group_by(genes) %>% summarise(n_neurons = n_distinct(Celltype), 
      is_in_ATN = "ExN_ATN" %in% Celltype, .groups = "drop") %>% filter(is_in_ATN == TRUE) %>% mutate(TRG_Group = ifelse(n_neurons == 
      1, "ATN-Specific", "Neuron-Common"))
  gnomad_file <- data_path("bin/gnomad.v4.1.1.constraint_metrics.tsv")
  df_gnomad_raw <- supplied_inputs[["df_gnomad_raw_3"]][[as.character(gnomad_file)]]
  ortholog_source <- Sys.getenv("DBSEQ_ORTHOLOG_FILE", file.path(data_root, "results/Table/human_cohort_de/conserved_trg.csv"))
  orthologs <- supplied_inputs[["orthologs_4"]][[as.character(ortholog_source)]] %>% filter(Has_Human_Homolog, 
      !is.na(Human_Gene), nzchar(Human_Gene)) %>% distinct(Mouse_Gene, Human_Gene)
  orthologs <- orthologs %>% group_by(Mouse_Gene) %>% filter(n() == 1) %>% ungroup() %>% group_by(Human_Gene) %>% 
      filter(n() == 1) %>% ungroup()
  df_gnomad <- df_gnomad_raw %>% select(Human_Gene = gene, gene_id, transcript, canonical, mane_select, 
      LOEUF = lof.oe_ci.upper) %>% mutate(LOEUF = as.numeric(LOEUF), is_mane = tolower(as.character(mane_select)) == 
      "true", is_canonical = tolower(as.character(canonical)) == "true", is_ensembl = grepl("^ENST", transcript)) %>% 
      filter(is_mane | is_canonical) %>% arrange(Human_Gene, desc(is_mane), desc(is_canonical), desc(is_ensembl), 
      transcript) %>% distinct(Human_Gene, .keep_all = TRUE)
  mapping_audit <- trg_class %>% left_join(orthologs, by = c(genes = "Mouse_Gene")) %>% left_join(df_gnomad, 
      by = "Human_Gene") %>% mutate(Included = !is.na(Human_Gene) & is.finite(LOEUF) & LOEUF > 0, Reason = case_when(is.na(Human_Gene) ~ 
      "No unique versioned ortholog", is.na(transcript) ~ "No MANE or canonical transcript", !is.finite(LOEUF) | 
      LOEUF <= 0 ~ "Unavailable positive LOEUF", TRUE ~ "Included"))
  write.csv(mapping_audit, output_path("statistics/figure 5/5m_gene_mapping.csv"), row.names = FALSE)
  df_loeuf <- mapping_audit %>% filter(Included) %>% mutate(log10_LOEUF = log10(LOEUF))
  stopifnot(!anyDuplicated(df_loeuf$Human_Gene), !anyDuplicated(df_loeuf$genes))
  csv_out <- df_loeuf %>% select(Gene = genes, Human_Gene, Category = TRG_Group, gene_id, transcript, canonical, 
      mane_select, LOEUF, log10_LOEUF)
  write.csv(csv_out, output_path("sourcedata/figure 5/5m.csv"), row.names = FALSE)
  p_value <- wilcox.test(log10_LOEUF ~ TRG_Group, df_loeuf, exact = FALSE)$p.value
  write.csv(data.frame(unit = "gene with one MANE/canonical transcript", n_records = nrow(df_loeuf), n_genes = n_distinct(df_loeuf$genes), 
      test = "Wilcoxon rank-sum", p_value = p_value, reproduced_reference_stars = FALSE), output_path("statistics/figure 5/5m_statistics.csv"), 
      row.names = FALSE)
  df_loeuf$Category_plot <- factor(df_loeuf$TRG_Group, levels = c("ATN-Specific", "Neuron-Common"), labels = c("ATN\nSpecific", 
      "Neuronal\nCommon"))
  p_loeuf <- ggplot(df_loeuf, aes(x = Category_plot, y = log10_LOEUF, fill = Category_plot)) + geom_boxplot(color = "black", 
      alpha = 0.75, outlier.shape = 21, outlier.alpha = 0, width = 0.5) + scale_fill_manual(values = c(`ATN\nSpecific` = "#E64B35", 
      `Neuronal\nCommon` = "#4DBBD5")) + annotate("segment", x = c(1, 1, 2), xend = c(1, 2, 2), y = c(0.73, 
      0.8, 0.8), yend = c(0.8, 0.8, 0.73)) + annotate("text", x = 1.5, y = 0.86, label = as.character(symnum(p_value, 
      corr = FALSE, na = FALSE, cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", 
          "**", "*", "ns")))) + scale_y_continuous(breaks = c(-1, 0, 1), expand = expansion(mult = c(0, 
      0))) + coord_cartesian(ylim = c(-1.8, 1.2)) + theme_bw(base_size = 14) + theme(legend.position = "none", 
      plot.title = element_text(size = 11), axis.title.x = element_blank(), axis.text = element_text(face = "bold", 
          color = "black"), panel.grid.minor = element_blank()) + labs(y = expression(log[10](LOEUF)), 
      title = "Tolerance to\nLoss-of-Function Mutations")
  ggsave(output_path("figures/figure 5/5m.png"), p_loeuf, width = 3.5, height = 4.5, dpi = 600)
  cat(">>> Figure 5m complete.\n")
  invisible(as.list(environment()))
}
