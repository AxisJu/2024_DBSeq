# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  aggregate_donor_auc <- function(data_obj) {
      meta <- data_obj$meta
      auc <- data_obj$auc_mat
      stopifnot(nrow(meta) == nrow(auc), !anyDuplicated(meta$cell_id))
      donor <- supplied_inputs[["donor_1"]][[as.character(data_path("data/snRNAseq_mouse/processed/metadata/metadata_demultiplexing.csv"))]]
      ids <- donor[[1]]
      stopifnot(!anyDuplicated(ids), all(meta$cell_id %in% ids))
      meta$donor_id <- donor$donor_id[match(meta$cell_id, ids)]
      keep <- !is.na(meta$donor_id) & !meta$donor_id %in% c("unassigned", "doublet") & meta$Group_L2 %in% 
          c("DBS_I", "Sham_I", "Saline_I")
      meta <- meta[keep, , drop = FALSE]
      auc <- auc[keep, , drop = FALSE]
      consistency <- tapply(meta$Group_L2, meta$donor_id, function(x) length(unique(x)))
      stopifnot(all(consistency == 1))
      key <- paste(meta$final_region, meta$Group_L2, meta$donor_id, sep = "|")
      keys <- unique(key)
      avg <- t(vapply(keys, function(k) colMeans(auc[key == k, , drop = FALSE]), numeric(ncol(auc))))
      colnames(avg) <- colnames(auc)
      rownames(avg) <- keys
      out_meta <- meta[match(keys, key), c("final_region", "Group_L2", "donor_id")]
      out_meta$n_cells <- as.integer(table(factor(key, levels = keys)))
      out_meta$cell_id <- keys
      write.csv(cbind(out_meta, as.data.frame(avg, check.names = FALSE)), output_path("sourcedata/figure 6/donor_mean_auc.csv"), 
          row.names = FALSE)
      list(auc_mat = avg, meta = out_meta)
  }
  invisible(as.list(environment()))
}
