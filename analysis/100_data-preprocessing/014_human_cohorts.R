# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  values <- supplied_inputs[["values_4"]][[as.character(input_path(opt$log_expression))]]
  stopifnot(setequal(rownames(values), rownames(obj)), setequal(colnames(values), colnames(obj)))
  LayerData(obj, assay = "RNA", layer = "data") <- values[rownames(obj), colnames(obj)]
  meta <- obj@meta.data
  selections <- list(cohort1_MTLE_brainbenigh_temporalcortex = meta$Cohort == "YC", cohort2_FCD2b_selfcontrol_frontalcortex = meta$Patient == 
      "Epi_13", cohort3_MTLE_autopsy_hippocampus = meta$BrainRegion == "Hippocampus" & meta$Cohort != "Tran", 
      cohort4_MTLE_autopsy_amygdala = meta$BrainRegion == "Amygdala" & meta$Cohort != "Tran" & meta$Disease != 
          "Epilepsy, LEAT")
  for (name in names(selections)) {
      cells <- rownames(meta)[which(selections[[name]])]
      if (!length(cells)) 
          stop("Empty cohort: ", name)
      qsave(obj[, cells], output_path(paste0("data/snRNAseq_human_DRE/", name, ".qs")))
  }
  invisible(as.list(environment()))
}
