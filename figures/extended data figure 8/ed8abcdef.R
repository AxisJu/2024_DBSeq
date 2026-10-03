# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
      library(tidyr)
      library(stringr)
      library(ggplot2)
      library(patchwork)
      library(scales)
      library(fgsea)
      library(grid)
  })
  args <- commandArgs(trailingOnly = FALSE)
  entry <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
  code_root <- normalizePath(file.path(dirname(entry), "../.."), winslash = "/")
  package_root <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(package_root)) stop("Set DBSEQ_OUTPUT_ROOT outside the code package.")
  output_resolved <- normalizePath(package_root, winslash = "/", mustWork = FALSE)
  if (identical(tolower(output_resolved), tolower(code_root)) || startsWith(tolower(output_resolved), paste0(tolower(code_root), 
      "/"))) stop("Outputs must remain outside the public code package")
  workspace_root <- package_root
  data_root <- Sys.getenv("DBSEQ_DATA_ROOT")
  input_root <- Sys.getenv("DBSEQ_INPUT_PACKAGE_ROOT")
  stopifnot(nzchar(data_root), nzchar(workspace_root))
  panel_dir <- function(kind, n) {
      p <- file.path(package_root, kind, paste("extended data figure", n))
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
      p
  }
  save_data <- function(df, n, panel) {
      stat_columns <- intersect(names(df), c("Pval", "FDR", "LogP", "LogFDR", "padj", "pval", "NES", "Sig_Stars"))
      if (length(stat_columns)) 
          write.csv(df, file.path(panel_dir("statistics", n), paste0("ed", n, panel, ".csv")), row.names = FALSE)
      write.csv(df[, setdiff(names(df), stat_columns), drop = FALSE], file.path(panel_dir("sourcedata", 
          n), paste0("ed", n, panel, ".csv")), row.names = FALSE)
  }
  save_png <- function(plot, n, panel, w, h) {
      png(file.path(panel_dir("figures", n), paste0("ed", n, panel, ".png")), width = w, height = h, units = "mm", 
          res = 600, type = "cairo")
      if (inherits(plot, "grob") || inherits(plot, "gtable")) 
          grid::grid.draw(plot)
      else print(plot)
      dev.off()
  }
  set.seed(42)
  pdf(file = NULL)
  trg_file <- Sys.getenv("DBSEQ_TRG_FILE")
  if (nzchar(trg_file)) {
      f_df <- if (grepl("\\.qs$", trg_file)) 
          supplied_inputs[["f_df_1"]][[as.character(trg_file)]]
      else supplied_inputs[["f_df_2"]][[as.character(trg_file)]]
      f_df <- dplyr::rename(f_df, Gene = genes, Is_TRG = is_TRG, CoD_DBSvSham = CoD_DBSvSham_cohend, CoD_ShamvSal = CoD_ShamvSal_cohend)
  } else {
      f_df <- supplied_inputs[["f_df_3"]][[as.character(file.path(input_root, "sourcedata/figure 2/2f.csv"))]]
  }
  ordered_celltypes <- c("ExN_RSP_L45IT", "ExN_RSP_L6CT", "ExN_DG", "InN_RSP_MGE", "InN_STRns_MGE", "InN_GPe_MGE", 
      "InN_Nonspecific", "InN_STRns_LGE", "InN_CP_D1", "InN_CP_D2", "InN_LSc", "InN_HYa", "InN_RT", "ExN_TRS", 
      "ExN_MH", "ExN_LH", "ExN_THns", "ExN_ATN", "ExN_RE", "ExN_PF", "ExN_CM", "Astro", "OPC", "Oligo", 
      "VCs", "Micro")
  ordered_celltypes <- c(intersect(ordered_celltypes, unique(f_df$Celltype)), sort(setdiff(unique(f_df$Celltype), 
      ordered_celltypes)))
  non_neuronal_types <- c("Astro", "OPC", "Oligo", "VCs", "Micro")
  neuronal_types <- setdiff(ordered_celltypes, non_neuronal_types)
  trg_df <- f_df %>% filter(Is_TRG == TRUE)
  gene_stats <- trg_df %>% group_by(Gene) %>% summarise(n_neuron = n_distinct(Celltype[Celltype %in% neuronal_types]), 
      n_nonneuron = n_distinct(Celltype[Celltype %in% non_neuronal_types]), .groups = "drop") %>% mutate(Group = case_when(n_neuron >= 
      1 & n_nonneuron >= 1 ~ "Neuron & Non-Neuron Common", n_neuron >= 2 & n_nonneuron == 0 ~ "Neuron Common", 
      n_nonneuron >= 2 & n_neuron == 0 ~ "Non-Neuron Common", TRUE ~ "Specific / Other"), ratio_neuron = (n_neuron + 
      0.1)/(n_nonneuron + 0.1), ratio_nonneuron = (n_nonneuron + 0.1)/(n_neuron + 0.1))
  common_genes_df <- gene_stats %>% filter(Group %in% c("Neuron & Non-Neuron Common", "Neuron Common", 
      "Non-Neuron Common"))
  cat(">>> Gene count per Common Group:\n")
  print(table(common_genes_df$Group))
  cat(">>> Generating Panel A (Height = 70 mm)...\n")
  tf_file <- file.path(data_root, "bin/Mouse_TFs_Kinases_webpage-3-30-2017.csv")
  tf_df <- supplied_inputs[["tf_df_4"]][[as.character(tf_file)]]
  colnames(tf_df)[1] <- "Symbol"
  known_tfs <- tf_df %>% filter(Class == "TF") %>% pull(Symbol) %>% tolower() %>% unique()
  common_genes_df$Is_TF <- ifelse(tolower(common_genes_df$Gene) %in% known_tfs, "TF", "Non-TF")
  common_genes_df$Is_TF <- factor(common_genes_df$Is_TF, levels = c("TF", "Non-TF"))
  cc_file <- file.path(data_root, "bin/MSigDB_Mm/m5.go.cc.v2023.2.Mm.symbols.gmt")
  pathways_cc <- gmtPathways(cc_file)
  all_cc_genes <- tolower(unique(unlist(pathways_cc)))
  common_genes_df$Coding_Type <- ifelse((grepl("Rik$", common_genes_df$Gene) | grepl("^Gm\\d+", common_genes_df$Gene)) & 
      !(tolower(common_genes_df$Gene) %in% all_cc_genes), "Non-Coding", "Protein-Coding")
  common_genes_df$Coding_Type <- factor(common_genes_df$Coding_Type, levels = c("Protein-Coding", "Non-Coding"))
  axon_syn_genes <- tolower(unique(unlist(pathways_cc[grepl("SYNAPSE|AXON|DENDRIT|POSTSYNAP|PRESYNAP", 
      names(pathways_cc))])))
  nuc_genes <- tolower(unique(unlist(pathways_cc[grepl("NUCLEUS|CHROMATIN|NUCLEOPLASM|NUCLEAR", names(pathways_cc))])))
  cyto_genes <- tolower(unique(unlist(pathways_cc[grepl("CYTOPLASM|CYTOSOL|RIBOSOM|GOLGI|ENDOPLASMIC", 
      names(pathways_cc))])))
  mito_genes <- tolower(unique(unlist(pathways_cc[grepl("MITOCHONDR", names(pathways_cc))])))
  pm_genes <- tolower(unique(unlist(pathways_cc[grepl("PLASMA_MEMBRANE|CELL_SURFACE", names(pathways_cc))])))
  assign_cc <- function(g) {
      gl <- tolower(g)
      if (gl %in% axon_syn_genes) 
          return("Synapse & Axon")
      if (gl %in% nuc_genes) 
          return("Nucleus")
      if (gl %in% mito_genes) 
          return("Mitochondria")
      if (gl %in% cyto_genes) 
          return("Cytoplasm & Organelle")
      if (gl %in% pm_genes) 
          return("Plasma Membrane")
      return("Other / Unspecified")
  }
  common_genes_df$CC_Compartment <- sapply(common_genes_df$Gene, assign_cc)
  common_genes_df$CC_Compartment <- factor(common_genes_df$CC_Compartment, levels = c("Synapse & Axon", 
      "Nucleus", "Cytoplasm & Organelle", "Mitochondria", "Plasma Membrane", "Other / Unspecified"))
  group_counts <- table(common_genes_df$Group)
  group_levels <- c("Neuron & Non-Neuron Common", "Neuron Common", "Non-Neuron Common")
  group_labels <- paste0(c("Common", "Neuronal", "Non-Neuronal"), "\n(n=", format(as.integer(group_counts[group_levels]), 
      big.mark = ",", trim = TRUE), ")")
  common_genes_df$Group_Plot <- factor(common_genes_df$Group, levels = group_levels, labels = group_labels)
  theme_panel_a <- theme_classic(base_size = 5, base_family = "Arial") + theme(text = element_text(family = "Arial", 
      size = 5, color = "black"), axis.line = element_line(linewidth = 0.4, color = "black"), axis.ticks = element_line(linewidth = 0.4, 
      color = "black"), axis.ticks.length = unit(0.8, "mm"), axis.text.x = element_text(size = 4.8, color = "black", 
      angle = 0, hjust = 0.5, vjust = 1, lineheight = 1), axis.text.y = element_text(size = 5, color = "black"), 
      axis.title.y = element_text(size = 5.5, face = "bold", color = "black"), axis.title.x = element_blank(), 
      plot.title = element_text(size = 5.5, face = "bold", hjust = 0.5, margin = margin(b = 3)), legend.position = "bottom", 
      legend.justification = "center", legend.box.margin = margin(0, 0, -1, 0), legend.margin = margin(0, 
          0, 0, 0), legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 4.2, family = "Arial"), 
      legend.title = element_text(size = 4.8, face = "bold", family = "Arial"), plot.margin = margin(4, 
          4, 4, 4))
  p_a1 <- ggplot(common_genes_df, aes(x = Group_Plot, fill = Coding_Type)) + geom_bar(position = "fill", 
      width = 0.6, color = "black", linewidth = 0.3) + scale_y_continuous(labels = label_percent(), expand = c(0, 
      0), limits = c(0, 1.02)) + scale_fill_manual(values = c(`Protein-Coding` = "#2c7fb8", `Non-Coding` = "#bdbdbd"), 
      name = "Biotype:") + labs(title = "Gene Biotype", y = "Proportion") + theme_panel_a + theme(plot.margin = margin(2, 
      3, 2, 6))
  p_a2 <- ggplot(common_genes_df, aes(x = Group_Plot, fill = CC_Compartment)) + geom_bar(position = "fill", 
      width = 0.6, color = "black", linewidth = 0.3) + scale_y_continuous(labels = label_percent(), expand = c(0, 
      0), limits = c(0, 1.02)) + scale_fill_manual(values = c(`Synapse & Axon` = "#ff1a71", Nucleus = "#8d7acc", 
      `Cytoplasm & Organelle` = "#ffd700", Mitochondria = "#e64b35", `Plasma Membrane` = "#18a799", `Other / Unspecified` = "#d9d9d9"), 
      name = "CC Compartment:", guide = guide_legend(ncol = 2, byrow = TRUE, title.position = "top")) + 
      labs(title = "Subcellular Localization (CC)", y = "") + theme_panel_a
  p_a3 <- ggplot(common_genes_df, aes(x = Group_Plot, fill = Is_TF)) + geom_bar(position = "fill", width = 0.6, 
      color = "black", linewidth = 0.3) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), 
      limits = c(0, 1.02)) + scale_fill_manual(values = c(TF = "#e74c3c", `Non-TF` = "#95a5a6"), name = "TF Status:") + 
      labs(title = "Transcription Factor Ratio", y = "") + theme_panel_a
  panel_a_combined <- (p_a1 | p_a2 | p_a3) + plot_layout(widths = c(1, 1.35, 1))
  save_png(p_a1, 8, "a", 55, 70)
  save_png(p_a2, 8, "b", 65, 70)
  save_png(p_a3, 8, "c", 55, 70)
  save_data(common_genes_df[, c("Gene", "Group", "Coding_Type")], 8, "a")
  save_data(common_genes_df[, c("Gene", "Group", "CC_Compartment")], 8, "b")
  save_data(common_genes_df[, c("Gene", "Group", "Is_TF")], 8, "c")
  cat(">>> Generating Panel B (Top 10 ORA, X=Count, Fill=-log10P)...\n")
  gmt_file <- file.path(data_root, "bin/MSigDB_Mm/m5.go.bp.v2023.2.Mm.symbols.gmt")
  pathways_bp <- gmtPathways(gmt_file)
  bg_genes <- unique(f_df$Gene)
  N_bg <- length(bg_genes)
  pathways_bp <- lapply(pathways_bp, function(g) intersect(g, bg_genes))
  pathways_bp <- pathways_bp[sapply(pathways_bp, length) >= 15 & sapply(pathways_bp, length) <= 500]
  ora_top10_list <- list()
  groups_order <- c("Neuron & Non-Neuron Common", "Neuron Common", "Non-Neuron Common")
  for (grp in groups_order) {
      gset <- common_genes_df %>% filter(Group == grp) %>% pull(Gene)
      gset <- intersect(gset, bg_genes)
      n_g <- length(gset)
      pvals <- numeric(length(pathways_bp))
      names(pvals) <- names(pathways_bp)
      overlaps <- integer(length(pathways_bp))
      names(overlaps) <- names(pathways_bp)
      for (pname in names(pathways_bp)) {
          pgenes <- pathways_bp[[pname]]
          k <- length(intersect(gset, pgenes))
          overlaps[pname] <- k
          if (k >= 3) {
              M <- length(pgenes)
              pvals[pname] <- phyper(k - 1, M, N_bg - M, n_g, lower.tail = FALSE)
          }
          else {
              pvals[pname] <- 1
          }
      }
      df_ora <- data.frame(Group = grp, Pathway = names(pathways_bp), Overlap = overlaps, Pval = pvals, 
          FDR = p.adjust(pvals, method = "BH"), stringsAsFactors = FALSE) %>% mutate(FDR = p.adjust(Pval, 
          method = "BH")) %>% filter(Overlap >= 3 & FDR < 0.05) %>% arrange(FDR) %>% head(10)
      df_ora$Clean_Term <- gsub("^GOBP_", "", df_ora$Pathway) %>% gsub("_", " ", .) %>% str_to_title()
      df_ora$Clean_Term <- gsub("Mrna", "mRNA", df_ora$Clean_Term)
      df_ora$Clean_Term <- gsub("Rna", "RNA", df_ora$Clean_Term)
      df_ora$Clean_Term <- gsub("Atp", "ATP", df_ora$Clean_Term)
      df_ora$Clean_Term <- gsub("Dna", "DNA", df_ora$Clean_Term)
      df_ora$LogP <- -log10(pmax(df_ora$FDR, 1e-25))
      df_ora$Wrapped_Term <- str_wrap(df_ora$Clean_Term, width = 24)
      ora_top10_list[[grp]] <- df_ora
  }
  ora_colors <- c("#f4f9f4", "#7fcdbb", "#2c7fb8", "#253494")
  all_logp <- unlist(lapply(ora_top10_list, function(d) d$LogP))
  min_lp <- 3
  max_lp <- 7.5
  create_top10_ora_bar <- function(df, title_str, x_limit = NULL) {
      df$Wrapped_Term <- factor(df$Wrapped_Term, levels = rev(df$Wrapped_Term))
      p <- ggplot(df, aes(x = Overlap, y = Wrapped_Term, fill = LogP)) + geom_col(width = 0.62, color = "black", 
          linewidth = 0.28) + scale_fill_gradientn(colors = ora_colors, limits = c(min_lp, max_lp), breaks = c(3, 
          4, 5, 6, 7), name = expression(-log[10](FDR))) + scale_x_continuous(expand = expansion(mult = c(0, 
          0.1)), limits = if (is.null(x_limit)) 
          NULL
      else c(0, x_limit), breaks = pretty_breaks(n = 3)) + labs(title = title_str, x = "Gene Count", y = "") + 
          theme_classic(base_size = 5, base_family = "Arial") + theme(text = element_text(family = "Arial", 
          size = 5, color = "black"), axis.line = element_line(linewidth = 0.4, color = "black"), axis.ticks = element_line(linewidth = 0.4, 
          color = "black"), axis.ticks.length = unit(0.7, "mm"), axis.text.y = element_text(size = 4.4, 
          color = "black", hjust = 1, lineheight = 0.82), axis.text.x = element_text(size = 4.8, color = "black"), 
          axis.title.x = element_text(size = 5.2, face = "bold", color = "black"), plot.title = element_text(size = 5.5, 
              face = "bold", hjust = 0.5, margin = margin(b = 2)), legend.position = "none", plot.margin = margin(2, 
              3, 2, 2))
      return(p)
  }
  p_b1 <- create_top10_ora_bar(ora_top10_list[["Neuron & Non-Neuron Common"]], "Neuron & Non-Neuron Common")
  p_b2 <- create_top10_ora_bar(ora_top10_list[["Neuron Common"]], "Neuron Common")
  p_b3 <- create_top10_ora_bar(ora_top10_list[["Non-Neuron Common"]], "Non-Neuron Common") + theme(legend.position = "right", 
      legend.key.height = unit(13, "mm"), legend.key.width = unit(2.2, "mm"), legend.title = element_text(size = 4.8, 
          face = "bold", angle = 90, hjust = 0.5), legend.text = element_text(size = 4.2), legend.box.margin = margin(0, 
          0, 0, -2))
  p_b3 <- p_b3 + theme(legend.title = element_text(size = 4.8, angle = 0)) + guides(fill = guide_colorbar(barheight = unit(32, 
      "mm"), barwidth = unit(2.2, "mm")))
  panel_b_combined <- (p_b1 | p_b2 | p_b3) + plot_layout(widths = c(1, 1, 1.22))
  save_png(panel_b_combined, 8, "d", 172, 70)
  ora_source <- bind_rows(ora_top10_list, .id = "Group_ID")
  ora_source$Gene_IDs <- vapply(seq_len(nrow(ora_source)), function(i) paste(intersect(common_genes_df$Gene[common_genes_df$Group == 
      ora_source$Group[i]], pathways_bp[[ora_source$Pathway[i]]]), collapse = ";"), character(1))
  save_data(ora_source, 8, "d")
  cat(">>> Generating Panel C (Scatter + Neuron Selectivity GSEA + Non-Neuron Selectivity GSEA)...\n")
  p_c1 <- ggplot(gene_stats, aes(x = n_nonneuron, y = n_neuron, color = ratio_neuron)) + geom_jitter(width = 0.22, 
      height = 0.22, size = 0.75, alpha = 0.75, stroke = 0.2, shape = 16) + scale_color_gradientn(colors = ora_colors, 
      trans = "log10", breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"), name = "Neuron Selectivity\nRatio") + 
      scale_x_continuous(breaks = 0:5, limits = c(-0.5, 5.5), expand = c(0, 0)) + scale_y_continuous(breaks = seq(0, 
      21, by = 5), limits = c(-0.5, 22), expand = c(0, 0)) + labs(title = "TRG Commonness & Selectivity", 
      x = "Non-Neuronal Commonness", y = "Neuronal Commonness") + theme_classic(base_size = 5, base_family = "Arial") + 
      theme(text = element_text(family = "Arial", size = 5, color = "black"), axis.line = element_line(linewidth = 0.4, 
          color = "black"), axis.ticks = element_line(linewidth = 0.4, color = "black"), axis.ticks.length = unit(0.8, 
          "mm"), axis.text = element_text(size = 5, color = "black"), axis.title = element_text(size = 5.5, 
          face = "bold", color = "black"), plot.title = element_text(size = 5.5, face = "bold", hjust = 0.5, 
          margin = margin(b = 2)), legend.position = "right", legend.key.height = unit(13, "mm"), legend.key.width = unit(2.2, 
          "mm"), legend.title = element_text(size = 4.5, face = "bold", lineheight = 0.85, margin = margin(b = 2)), 
          legend.text = element_text(size = 4.2), legend.box.margin = margin(2, 0, 0, -2), plot.margin = margin(5, 
              4, 1, 2))
  all_genes_df <- data.frame(Gene = unique(f_df$Gene), stringsAsFactors = FALSE) %>% left_join(gene_stats, 
      by = "Gene") %>% mutate(n_neuron = ifelse(is.na(n_neuron), 0, n_neuron), n_nonneuron = ifelse(is.na(n_nonneuron), 
      0, n_nonneuron))
  set.seed(42)
  j_v <- runif(nrow(all_genes_df), -1e-06, 1e-06)
  ranks_neuron_ratio <- sort(setNames((all_genes_df$n_neuron + 0.1)/(all_genes_df$n_nonneuron + 0.1) + 
      j_v, all_genes_df$Gene), decreasing = TRUE)
  ranks_nonneuron_ratio <- sort(setNames((all_genes_df$n_nonneuron + 0.1)/(all_genes_df$n_neuron + 0.1) + 
      j_v, all_genes_df$Gene), decreasing = TRUE)
  calc_gsea_profile <- function(pathway_genes, stats, gseaParam = 1) {
      gene_names <- names(stats)
      N <- length(stats)
      hits <- gene_names %in% pathway_genes
      hit_indices <- which(hits)
      N_H <- length(hit_indices)
      if (N_H == 0) 
          return(list(df_curve = data.frame(x = 1:N, y = 0), hits = integer(0), max_es = 0))
      hit_weights <- abs(stats[hits])^gseaParam
      NR <- sum(hit_weights)
      if (NR == 0) 
          NR <- 1
      step_hit <- numeric(N)
      step_hit[hits] <- hit_weights/NR
      step_miss <- numeric(N)
      step_miss[!hits] <- 1/(N - N_H)
      running_score <- cumsum(step_hit - step_miss)
      return(list(df_curve = data.frame(x = 1:N, y = running_score), hits = hit_indices, max_es = max(running_score)))
  }
  c2_terms <- c("GOBP_GLUCOSE_METABOLIC_PROCESS", "GOBP_MONOSACCHARIDE_METABOLIC_PROCESS", "GOBP_MONOSACCHARIDE_BIOSYNTHETIC_PROCESS", 
      "GOBP_MONOSACCHARIDE_CATABOLIC_PROCESS", "GOBP_CARBOHYDRATE_CATABOLIC_PROCESS")
  c3_terms <- c("GOBP_IMMUNE_EFFECTOR_PROCESS", "GOBP_ENSHEATHMENT_OF_NEURONS", "GOBP_CELL_ACTIVATION", 
      "GOBP_INFLAMMATORY_RESPONSE", "GOBP_ENDOCYTOSIS")
  gsea_source <- list()
  line_palette <- c("#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a")
  clean_term_name <- function(term) {
      name <- gsub("^GOBP_", "", term) %>% gsub("_", " ", .) %>% str_to_title()
      name <- gsub("Mrna", "mRNA", name)
      name <- gsub("Rna", "RNA", name)
      name <- gsub("Atp", "ATP", name)
      name <- gsub("Dna", "DNA", name)
      return(name)
  }
  build_gsea_subplot <- function(terms, ranks, col_title, is_middle = TRUE) {
      N_genes <- length(ranks)
      curve_list <- list()
      hit_list <- list()
      for (i in seq_along(terms)) {
          pname <- terms[i]
          pgenes <- pathways_bp[[pname]]
          prof <- calc_gsea_profile(pgenes, ranks)
          gsea_source[[paste(col_title, pname)]] <<- data.frame(Lineage = col_title, Pathway = pname, Gene = names(ranks), 
              Rank = seq_along(ranks), Selectivity_ratio = as.numeric(ranks), In_pathway = names(ranks) %in% 
                  pgenes, Enrichment_score = prof$df_curve$y)
          clean_lbl <- clean_term_name(pname)
          df_c <- prof$df_curve
          df_c$Pathway <- factor(clean_lbl, levels = clean_term_name(terms))
          curve_list[[i]] <- df_c
          if (length(prof$hits) > 0) {
              df_h <- data.frame(x = prof$hits, Pathway_Idx = i, Pathway = factor(clean_lbl, levels = clean_term_name(terms)))
              hit_list[[i]] <- df_h
          }
      }
      all_curves <- bind_rows(curve_list)
      all_hits <- bind_rows(hit_list)
      named_palette <- setNames(line_palette, clean_term_name(terms))
      p_curve <- ggplot(all_curves, aes(x = x, y = y, color = Pathway)) + geom_hline(yintercept = 0, color = "grey75", 
          linewidth = 0.35, linetype = "dashed") + geom_line(linewidth = 0.5) + scale_color_manual(values = named_palette) + 
          scale_x_continuous(expand = c(0, 0), limits = c(1, N_genes)) + scale_y_continuous(expand = expansion(mult = c(0.02, 
          0.08)), breaks = pretty_breaks(n = 4)) + labs(title = col_title, x = "", y = if (is_middle) 
          "Enrichment score"
      else "") + theme_classic(base_size = 5, base_family = "Arial") + theme(text = element_text(family = "Arial", 
          size = 5, color = "black"), axis.line = element_line(linewidth = 0.4, color = "black"), axis.ticks = element_line(linewidth = 0.4, 
          color = "black"), axis.ticks.length = unit(0.8, "mm"), axis.text.y = element_text(size = 5, color = "black"), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size = 5.5, 
              face = "bold", color = "black"), plot.title = element_text(size = 5.5, face = "bold", hjust = 0.5, 
              margin = margin(b = 2)), legend.position = "none", plot.margin = margin(2, 4, 1, if (is_middle) 
              2
          else 4))
      barcode_plots <- list()
      for (i in seq_along(terms)) {
          c_name <- clean_term_name(terms[i])
          sub_hits <- all_hits %>% filter(Pathway_Idx == i)
          c_color <- line_palette[i]
          p_bar <- ggplot() + geom_segment(data = sub_hits, aes(x = x, xend = x, y = 0, yend = 1), color = c_color, 
              linewidth = 0.35) + scale_x_continuous(expand = c(0, 0), limits = c(1, N_genes)) + scale_y_continuous(expand = c(0, 
              0), limits = c(0, 1)) + labs(title = c_name) + theme_void(base_size = 5, base_family = "Arial") + 
              theme(plot.title = element_text(family = "Arial", size = 4.8, color = "black", hjust = 0, 
                  margin = margin(t = 2.5, b = 1)), panel.border = element_blank(), plot.margin = margin(0, 
                  4, 1, if (is_middle) 
                    2
                  else 4))
          barcode_plots[[i]] <- p_bar
      }
      col_combined <- (p_curve/barcode_plots[[1]]/barcode_plots[[2]]/barcode_plots[[3]]/barcode_plots[[4]]/barcode_plots[[5]]) + 
          plot_layout(heights = c(2.6, 0.38, 0.38, 0.38, 0.38, 0.38))
      return(col_combined)
  }
  p_c2 <- build_gsea_subplot(c2_terms, ranks_neuron_ratio, "Neuron selectivity: descriptive ES profiles", 
      is_middle = TRUE)
  p_c3 <- build_gsea_subplot(c3_terms, ranks_nonneuron_ratio, "Non-neuron selectivity: descriptive ES profiles", 
      is_middle = FALSE)
  panel_c_combined <- (p_c1 | p_c2 | p_c3) + plot_layout(widths = c(1.15, 1.05, 1))
  save_png(p_c1, 8, "e", 65, 78)
  save_png(p_c2 | p_c3, 8, "f", 110, 78)
  save_data(gene_stats[, c("Gene", "n_neuron", "n_nonneuron", "ratio_neuron")], 8, "e")
  save_data(bind_rows(gsea_source), 8, "f")
  invisible(as.list(environment()))
}
