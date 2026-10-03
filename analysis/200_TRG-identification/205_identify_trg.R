# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(dplyr)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  columns <- if (opt$scope == "atn") c(ANT = "celltype") else setNames(paste0("celltype_level", strsplit(opt$levels, 
      ",")[[1]]), paste0("Level", strsplit(opt$levels, ",")[[1]]))
  read_de <- function(relative, prefix) {
      path <- input_path(relative)
      if (!supplied_inputs[["read_de_4"]]) 
          stop("Missing prerequisite: ", relative)
      df <- supplied_inputs[["df_5"]][[as.character(path)]]
      names(df)[1] <- if (names(df)[1] %in% c("", "X", "gene")) 
          "genes"
      else names(df)[1]
      if ("gene" %in% names(df)) 
          names(df)[names(df) == "gene"] <- "genes"
      stopifnot("genes" %in% names(df), !anyDuplicated(df$genes))
      names(df)[names(df) != "genes"] <- paste0(prefix, "_", names(df)[names(df) != "genes"])
      df
  }
  for (level in names(columns)) {
      all_results <- list()
      trgs <- list()
      types <- sort(unique(as.character(obj@meta.data[[columns[[level]]]])))
      for (ct in types) {
          filename <- file.path("results/Table/CosSim", level, paste0("CosSim_", ct, ".csv"))
          if (!supplied_inputs[["data_6"]]) 
              next
          result <- supplied_inputs[["result_7"]][[as.character(input_path(filename))]]
          memroot <- if (opt$scope == "atn") 
              "results/Table/MEMENTO/ANT"
          else file.path("results/Table/MEMENTO_v2", level)
          pbroot <- if (opt$scope == "atn") 
              "results/Table/EdgeR/ANT"
          else file.path("results/Table/EdgeR_v2", level)
          result <- result %>% left_join(supplied_inputs[["result_8"]][[as.character(file.path(memroot, 
              paste0("MEMENTO_1d_DBSvSham_", ct, ".csv")))]], by = "genes") %>% left_join(supplied_inputs[["result_9"]][[as.character(file.path(memroot, 
              paste0("MEMENTO_1d_ShamvSaline_", ct, ".csv")))]], by = "genes") %>% left_join(supplied_inputs[["result_10"]][[as.character(file.path(pbroot, 
              paste0("EdgeR_DBSIvShamI_", ct, ".csv")))]], by = "genes") %>% left_join(supplied_inputs[["result_11"]][[as.character(file.path(pbroot, 
              paste0("EdgeR_ShamIvSalineI_", ct, ".csv")))]], by = "genes")
          if (opt$scope != "atn") {
              result <- result %>% left_join(supplied_inputs[["result_12"]][[as.character(file.path("results/Table/CohenD", 
                  level, paste0("CohenD_DBSIvShamI_", ct, ".csv")))]], by = "genes") %>% left_join(supplied_inputs[["result_13"]][[as.character(file.path("results/Table/CohenD", 
                  level, paste0("CohenD_ShamIvSalineI_", ct, ".csv")))]], by = "genes")
              opposite <- result$CoD_DBSvSham_cohend * result$CoD_ShamvSal_cohend < 0
          }
          else opposite <- result$Mem_DBSvSham_de_coef * result$Mem_ShamvSal_de_coef < 0
          result$is_TRG <- with(result, bootstrap_p < 0.05 & shuffle_p < 0.05 & Mem_DBSvSham_FDR < 0.01 & 
              Mem_ShamvSal_FDR < 0.01 & PB_DBSvSham_FDR < 0.05 & PB_ShamvSal_FDR < 0.05) & opposite
          result$is_TRG[is.na(result$is_TRG)] <- FALSE
          all_results[[ct]] <- result
          trgs[[ct]] <- result[result$is_TRG, ]
          write.csv(result, output_path(file.path("results/Table/CosSim_v2", level, paste0(ct, ".csv"))), 
              row.names = FALSE)
          write.csv(trgs[[ct]], output_path(file.path("results/Table/CosSim_v2", level, paste0(ct, "_TRG.csv"))), 
              row.names = FALSE)
      }
      if (!length(all_results)) 
          stop("No CosSim inputs for ", level)
      base <- "data/snRNAseq_mouse/processed/intermediate/CosSim_v2"
      qsave(all_results, output_path(file.path(base, paste0(level, "_res.qs"))))
      qsave(trgs, output_path(file.path(base, paste0(level, "_res_TRG.qs"))))
      plot_df <- bind_rows(all_results, .id = "Celltype")
      qsave(plot_df, output_path(file.path("results/Table/CosSim_v2", paste0(level, "_plot_df.qs"))))
  }
  if (opt$scope == "all") {
      files <- file.path(opt$output_root, "results/Table/CosSim_v2", c("Level1_plot_df.qs", "Level3_plot_df.qs"))
      if (all(supplied_inputs[["data_14"]])) {
          level1 <- supplied_inputs[["level1_15"]][[as.character(files[1])]]
          level3 <- supplied_inputs[["level3_16"]][[as.character(files[2])]]
          plot_df <- bind_rows(level1[!level1$Celltype %in% c("ExN", "InN", "CHPCs", "EPCs"), ], level3[grepl("^(ExN|InN)_", 
              level3$Celltype), ])
          plot_df$analysis_mode <- "reanalysis"
          plot_df$delta_CohenD <- plot_df$CoD_ShamvSal_cohend - plot_df$CoD_DBSvSham_cohend
          qsave(plot_df, output_path("results/Table/CosSim_v2/plot_df.qs"))
          write.csv(plot_df, output_path("results/Table/CosSim_v2/plot_df.csv"), row.names = FALSE)
      }
  }
  invisible(as.list(environment()))
}
