# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(qs)
      library(Seurat)
      library(dplyr)
  })
  cohort_paths <- list(`Cohort 1 (Temporal)` = file.path(data_root, "data/snRNAseq_human_DRE/cohort1_MTLE_brainbenigh_temporalcortex.qs"), 
      `Cohort 2 (Frontal)` = file.path(data_root, "data/snRNAseq_human_DRE/cohort2_FCD2b_selfcontrol_frontalcortex.qs"), 
      `Cohort 3 (Hippocampus)` = file.path(data_root, "data/snRNAseq_human_DRE/cohort3_MTLE_autopsy_hippocampus.qs"), 
      `Cohort 4 (Amygdala)` = file.path(data_root, "data/snRNAseq_human_DRE/cohort4_MTLE_autopsy_amygdala.qs"))
  target_human_cts <- c("Glu.N", "GABA.N", "Astro.", "OPCs", "Oligo.", "Microglia", "Endo.")
  human_markers <- c("SLC17A7", "CAMK2A", "GAD1", "GAD2", "AQP4", "GFAP", "PDGFRA", "VCAN", "MBP", "MOG", 
      "CX3CR1", "P2RY12", "CLDN5", "PECAM1")
  qc_meta_list <- list()
  dotplot_data_list <- list()
  for (cname in names(cohort_paths)) {
      cat(sprintf(">>> Loading and processing %s...\n", cname))
      cpath <- cohort_paths[[cname]]
      tmp_obj <- supplied_inputs[["tmp_obj_1"]][[as.character(cpath)]]
      valid_cells <- Cells(tmp_obj)[tmp_obj$celltype_coarse %in% target_human_cts]
      tmp_obj <- tmp_obj[, valid_cells]
      tmp_obj$celltype_coarse <- factor(tmp_obj$celltype_coarse, levels = target_human_cts)
      tmp_obj$Group <- factor(tmp_obj$Group, levels = c("Control", "Epilepsy"))
      qc_df <- data.frame(Cell_ID = Cells(tmp_obj), Sample = tmp_obj$Sample_snSeq.processed, nFeature_RNA = tmp_obj$nFeature_RNA, 
          nCount_RNA = tmp_obj$nCount_RNA, percent.mt = tmp_obj$percent.mt, percent.rb = tmp_obj$percent.rb, 
          Group = tmp_obj$Group, celltype = tmp_obj$celltype_coarse, Cohort = cname, stringsAsFactors = FALSE)
      qc_meta_list[[cname]] <- qc_df
      Idents(tmp_obj) <- "celltype_coarse"
      present_markers <- intersect(human_markers, rownames(tmp_obj))
      dp_tmp <- DotPlot(tmp_obj, features = present_markers, idents = target_human_cts)
      dp_df <- dp_tmp$data
      dp_df$Cohort <- cname
      dotplot_data_list[[cname]] <- dp_df
      rm(tmp_obj)
      gc()
  }
  s7_cache <- list(qc = bind_rows(qc_meta_list), dotplot = bind_rows(dotplot_data_list), target_cts = target_human_cts, 
      human_markers = human_markers)
  s7_data <- s7_cache
  invisible(as.list(environment()))
}
