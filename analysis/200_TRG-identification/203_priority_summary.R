# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(dplyr)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  scores <- supplied_inputs[["scores_4"]][[as.character(input_path(opt$scdrs))]]
  if (anyDuplicated(rownames(scores))) stop("Duplicate scDRS cell IDs")
  idx <- match(rownames(obj@meta.data), rownames(scores))
  meta <- obj@meta.data
  meta$scdrs <- scores$norm_score[idx]
  comparisons <- list(SS = c("ShamvSaline", "salinesham"), DS = c("DBSvSham", "dbssham"))
  for (level in strsplit(opt$levels, ",")[[1]]) {
      column <- paste0("celltype_level", level)
      result <- meta %>% group_by(Celltype = .data[[column]]) %>% summarise(scDRS_nScore = mean(scdrs, 
          na.rm = TRUE), .groups = "drop")
      for (name in names(comparisons)) {
          codes <- comparisons[[name]]
          augur <- supplied_inputs[["augur_5"]][[as.character(input_path(paste0("data/snRNAseq_mouse/processed/intermediate/Augur/augur_", 
              codes[2], "_level", level, "_true.qs")))]]$AUC
          names(augur)[1:2] <- c("Celltype", paste0("AUC_", name))
          rows <- lapply(as.character(result$Celltype), function(ct) {
              path <- input_path(paste0("results/Table/MEMENTO_v2/Level", level, "/MEMENTO_1d_", codes[1], 
                  "_", ct, ".csv"))
              if (!supplied_inputs[["rows_6"]]) 
                  return(NULL)
              de <- supplied_inputs[["de_7"]][[as.character(path)]]
              row <- data.frame(Celltype = ct, n = sum(de$FDR < 0.01, na.rm = TRUE))
              names(row)[2] <- paste0("DE_Num_", name)
              row
          })
          result <- result %>% left_join(augur, by = "Celltype") %>% left_join(bind_rows(rows), by = "Celltype")
      }
      write.csv(result, output_path(paste0("results/Table/celltype_priority_Level", level, ".csv")), row.names = FALSE)
      variables <- setdiff(names(result), "Celltype")
      comparisons_out <- combn(variables, 2, simplify = FALSE)
      tests <- lapply(comparisons_out, function(pair) {
          valid <- complete.cases(result[, pair]) & is.finite(result[[pair[1]]]) & is.finite(result[[pair[2]]])
          if (sum(valid) < 3 || sd(result[[pair[1]]][valid]) == 0 || sd(result[[pair[2]]][valid]) == 0) 
              return(NULL)
          test <- cor.test(result[[pair[1]]][valid], result[[pair[2]]][valid], method = "pearson")
          data.frame(first = pair[1], second = pair[2], n_celltypes = sum(valid), correlation = unname(test$estimate), 
              p_value = test$p.value)
      })
      table <- bind_rows(tests)
      if (nrow(table)) {
          table$fdr <- p.adjust(table$p_value, "BH")
          write.csv(table, output_path(paste0("results/Table/celltype_priority_correlations_Level", level, 
              ".csv")), row.names = FALSE)
      }
  }
  invisible(as.list(environment()))
}
