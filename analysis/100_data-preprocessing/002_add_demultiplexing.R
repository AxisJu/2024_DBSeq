# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  base <- "data/snRNAseq_mouse/processed/metadata"
  annotation <- supplied_inputs[["annotation_3"]][[as.character(input_path(file.path(base, "metadata_obs_101.csv")))]]
  mapping <- supplied_inputs[["mapping_4"]][[as.character(input_path(file.path(base, "metadata_cellnamemapping.csv")))]]
  stopifnot(!anyDuplicated(mapping$unique_cellname))
  idx <- match(rownames(annotation), mapping$unique_cellname)
  stopifnot(!anyNA(idx))
  annotation$raw_cellname <- mapping$raw_cellname[idx]
  libraries <- c(DBS_I_1 = "DBS-1", DBS_I_2 = "DBS-2", DBS_C_1 = "DBSDC-1", PTZ_1 = "PTZ-old", Saline_I_1 = "Saline-1", 
      Saline_I_2 = "Saline-2", Sham_I_1 = "Sham-1", Sham_I_2 = "Sham-2", Sham_I_3 = "Sham-3")
  donors <- list(DBS_I_1 = c("DBS_M1", "DBS_M2", "DBS_M3"), DBS_I_2 = c("DBS_M4", "DBS_M5"), DBS_C_1 = c("DBS_M1", 
      "DBS_M2", "DBS_M3"), Saline_I_1 = c("Saline_M1", "Saline_M2", "Saline_M3"), Saline_I_2 = c("Saline_M3", 
      "Saline_M1", "Saline_M2"), Sham_I_1 = c("Sham_M1", "Sham_M2", "Sham_M3"), Sham_I_2 = c("Sham_M2", 
      "Sham_M1", "Sham_M3"), Sham_I_3 = c("Sham_M3", "Sham_M1", "Sham_M2"))
  results <- lapply(names(libraries), function(sample) {
      selected <- annotation[annotation$Sample == sample, , drop = FALSE]
      if (sample == "PTZ_1") {
          return(data.frame(donor_id = rep("PTZ_M1", nrow(selected)), prob_max = NA_real_, prob_doublet = NA_real_, 
              n_vars = NA_integer_, row.names = rownames(selected)))
      }
      path <- input_path(paste0("data/snRNAseq_mouse/processed/intermediate/vireoSNP/", libraries[[sample]], 
          "/donor_ids.tsv"))
      dat <- supplied_inputs[["dat_5"]][[as.character(path)]]
      stopifnot(!anyDuplicated(rownames(dat)), all(selected$raw_cellname %in% rownames(dat)))
      dat <- dat[selected$raw_cellname, c("donor_id", "prob_max", "prob_doublet", "n_vars"), drop = FALSE]
      lookup <- setNames(donors[[sample]], paste0("donor", seq_along(donors[[sample]]) - 1))
      mapped <- lookup[dat$donor_id]
      dat$donor_id[!is.na(mapped)] <- mapped[!is.na(mapped)]
      rownames(dat) <- rownames(selected)
      dat
  })
  combined <- do.call(rbind, results)
  stopifnot(setequal(rownames(combined), rownames(annotation)))
  write.csv(combined[rownames(annotation), ], output_path(file.path(base, "metadata_demultiplexing.csv")))
  invisible(as.list(environment()))
}
