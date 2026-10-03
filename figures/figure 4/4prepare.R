# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
  })
  read_trg <- function() {
      path <- Sys.getenv("DBSEQ_TRG_FILE", data_path("results/Table/CosSim_v2/plot_df.qs"))
      supplied_inputs[["read_trg_2"]]
  }
  prepare_4a <- function() {
      suppressPackageStartupMessages({
          library(Seurat)
          library(Matrix)
      })
      object <- supplied_inputs[["object_3"]][[as.character(data_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]]
      meta <- object@meta.data
      keep <- meta$celltype_level3 == "ExN_ATN" & meta$Group_L2 %in% c("DBS_I", "Sham_I", "Saline_I") & 
          !meta$donor_id %in% c("unassigned", "doublet")
      keep[is.na(keep)] <- FALSE
      meta <- meta[keep, , drop = FALSE]
      expression <- GetAssayData(object, assay = "RNA", layer = "data")[, rownames(meta), drop = FALSE]
      means <- sapply(c("DBS_I", "Sham_I", "Saline_I"), function(group) Matrix::rowMeans(expression[, meta$Group_L2 == 
          group, drop = FALSE]))
      means <- as.data.frame(means)
      means$gene <- rownames(means)
      trg <- supplied_inputs[["trg_4"]] %>% filter(Celltype == "ExN_ATN")
      for (prefix in c("median_cos", "mean_cos")) {
          index <- grep(paste0("^", prefix), names(trg))
          if (length(index) > 1) 
              stop("Ambiguous cosine columns: ", prefix)
          if (length(index) == 1) 
              names(trg)[index] <- paste0(prefix, "ine")
      }
      if (!"mean_cosine" %in% names(trg)) 
          trg$mean_cosine <- NA_real_
      stopifnot(!anyDuplicated(trg$genes))
      trg <- trg %>% select(genes, is_TRG, CoD_DBSvSham_cohend, CoD_ShamvSal_cohend, median_cosine, mean_cosine)
      combined <- left_join(means, trg, by = c(gene = "genes")) %>% mutate(is_TRG = coalesce(is_TRG, FALSE), 
          delta_CohenD = CoD_DBSvSham_cohend - CoD_ShamvSal_cohend)
      total <- rowSums(combined[, c("DBS_I", "Sham_I", "Saline_I")])
      write.csv(data.frame(gene = combined$gene, total_expression = total, included = total > 0), output_path("statistics/figure 4/4a_expression_eligibility.csv"), 
          row.names = FALSE)
      combined <- combined[is.finite(total) & total > 0, ]
      qsave(combined, cache_path("df_combined_tern.qs"))
  }
  prepare_4_ora <- function() {
      mouse <- supplied_inputs[["mouse_5"]]
      library <- supplied_inputs[["library_6"]][[as.character(input_path("sourcedata/figure 5/inputs/go_bp_2023.qs"))]]
      background <- sort(unique(mouse$genes[!is.na(mouse$genes)]))
      library <- lapply(library, function(genes) intersect(genes, background))
      library <- library[lengths(library) > 0 & lengths(library) < length(background)]
      neuron <- grepl("^(ExN|InN)", mouse$Celltype)
      classified <- mouse[neuron & mouse$is_TRG, ] %>% group_by(genes) %>% summarise(n_neurons = n_distinct(Celltype), 
          is_in_ATN = any(Celltype == "ExN_ATN"), .groups = "drop") %>% filter(is_in_ATN)
      specific <- intersect(classified$genes[classified$n_neurons == 1], background)
      common <- intersect(classified$genes[classified$n_neurons > 1], background)
      all_atn <- union(specific, common)
      ora <- function(query) {
          query <- intersect(query, background)
          k <- vapply(library, function(genes) length(intersect(query, genes)), integer(1))
          sizes <- lengths(library)
          p <- phyper(k - 1, sizes, length(background) - sizes, length(query), lower.tail = FALSE)
          data.frame(Term = names(library), P.value = p, Adjusted.P.value = p.adjust(p, "BH"), Overlap = paste0(k, 
              "/", sizes), Gene_Count = k, Genes = vapply(library, function(genes) paste(sort(intersect(query, 
              genes)), collapse = ";"), character(1)), stringsAsFactors = FALSE)
      }
      atn_res <- ora(specific)
      common_res <- ora(common)
      all_res <- ora(all_atn)
      qsave(list(atn = atn_res, common = common_res), cache_path("enrichr_cache.qs"))
      for (name in c("specific", "common", "all")) {
          result <- switch(name, specific = atn_res, common = common_res, all = all_res)
          write.csv(result, output_path("statistics/figure 4", paste0("4_ORA_", name, ".csv")), row.names = FALSE)
      }
      selected <- all_res$Term[all_res$Adjusted.P.value < 0.05 & all_res$Gene_Count > 5]
      if (!length(selected)) 
          stop("No pathway passes full-library BH FDR < 0.05 and overlap > 5.")
      pathways <- library[selected]
      stats <- data.frame(Term = selected, Clean_Term = sub(" \\(GO:[0-9]+\\)", "", selected), Count_ATN = vapply(pathways, 
          function(g) length(intersect(g, specific)), integer(1)), Count_Common = vapply(pathways, function(g) length(intersect(g, 
          common)), integer(1)), Total_Genes = lengths(pathways))
      stats$Total_Enriched <- stats$Count_ATN + stats$Count_Common
      stats$Prop_ATN <- stats$Count_ATN/stats$Total_Genes
      stats$Prop_Common <- stats$Count_Common/stats$Total_Genes
      stats$Category <- ifelse(stats$Count_ATN == 0, "Only Neuron-Common", ifelse(stats$Count_Common == 
          0, "Only ATN-specific", "Both Shared"))
      union_genes <- sort(unique(unlist(pathways)))
      gene_category <- data.frame(Gene = union_genes, Category = ifelse(union_genes %in% specific, "ATN-specific", 
          ifelse(union_genes %in% common, "Neuronal-Common", "Other Genes")))
      counts <- as.data.frame(table(factor(gene_category$Category, levels = c("Other Genes", "Neuronal-Common", 
          "ATN-specific"))))
      names(counts) <- c("Category", "Count")
      weights <- vapply(all_atn, function(g) mean(vapply(pathways, function(p) as.numeric(g %in% p)/length(p), 
          numeric(1))), numeric(1))
      observed <- sum(weights[all_atn %in% specific])
      set.seed(261002)
      null <- replicate(9999, sum(sample(weights, length(specific), replace = FALSE)))
      permutation_p <- (1 + sum(null >= observed - 1e-14))/(1 + length(null))
      master <- list(df_pathway_stats = stats, grn_counts = counts, bg_rate = mean(null), p_value = permutation_p, 
          test = "Conditional ATN-specific/common gene-label permutation; 9999 draws")
      qsave(master, cache_path("fig4_def_master.qs"))
      write.csv(stats, output_path("statistics/figure 4/4_pathway_counts.csv"), row.names = FALSE)
      write.csv(gene_category, output_path("sourcedata/figure 4/4e_gene_membership.csv"), row.names = FALSE)
      write.csv(data.frame(Term = rep(names(pathways), lengths(pathways)), Gene = unlist(pathways, use.names = FALSE)), 
          output_path("sourcedata/figure 4/4d_pathway_membership.csv"), row.names = FALSE)
      write.csv(data.frame(Observed = observed, Null_mean = mean(null), P = permutation_p, B = length(null), 
          Universe = length(background), ATN_specific = length(specific), ATN_common = length(common)), 
          output_path("statistics/figure 4/4f_permutation.csv"), row.names = FALSE)
  }
  if (basename(sub("^--file=", "", entry_arg)) == "4prepare.R") {
      prepare_4_ora()
      prepare_4a()
  }
  invisible(as.list(environment()))
}
