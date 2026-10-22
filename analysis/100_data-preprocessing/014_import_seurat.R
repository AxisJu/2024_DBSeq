# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  library(reticulate)
  np <- import("numpy")
  scipy <- import("scipy.sparse")
  base <- "data/snRNAseq_mouse/processed"
  obs <- supplied_inputs[["obs_3"]][[as.character(input_path(file.path(base, "metadata/metadata_obs_109.csv")))]]
  var <- supplied_inputs[["var_4"]][[as.character(input_path(file.path(base, "metadata/metadata_var.csv")))]]
  donors <- supplied_inputs[["donors_5"]][[as.character(input_path(file.path(base, "metadata/metadata_demultiplexing.csv")))]]
  stopifnot(!anyDuplicated(rownames(obs)), !anyDuplicated(rownames(donors)), all(rownames(obs) %in% rownames(donors)))
  obs <- cbind(obs, donors[rownames(obs), , drop = FALSE])
  counts <- supplied_inputs[["counts_6"]][[as.character(input_path(file.path(base, "matrix/matrix_counts.npz")))]]
  stopifnot(identical(as.integer(dim(counts)), as.integer(c(nrow(var), nrow(obs)))))
  dimnames(counts) <- list(rownames(var), rownames(obs))
  obj <- CreateSeuratObject(counts = counts, meta.data = obs, assay = "RNA", project = "DBSeq", min.cells = 0)
  logdata <- supplied_inputs[["logdata_7"]][[as.character(input_path(file.path(base, "matrix/matrix_logCPM.npz")))]]
  dimnames(logdata) <- dimnames(counts)
  LayerData(obj, assay = "RNA", layer = "data") <- logdata
  qsave(logdata, output_path(file.path(base, "matrix/matrix_logCPM.qs")))
  for (key in c("PCA", "scVI", "scANVI", "UMAP_scANVI")) {
      values <- supplied_inputs[["values_8"]][[as.character(input_path(file.path(base, paste0("metadata/reduction_", 
          key, ".csv"))))]]
      stopifnot(all(rownames(obs) %in% rownames(values)))
      values <- as.matrix(values[rownames(obs), , drop = FALSE])
      name <- c(PCA = "pca", scVI = "scvi", scANVI = "scanvi", UMAP_scANVI = "umap")[[key]]
      colnames(values) <- paste0(name, "_", seq_len(ncol(values)))
      obj[[name]] <- CreateDimReducObject(embeddings = values, key = paste0(name, "_"), assay = "RNA")
  }
  qsave(obj, output_path(file.path(base, "matrix/running_all_250704.qs")))
  magic_file <- input_path(file.path(base, "matrix/matrix_MAGIC.npz"))
  if (supplied_inputs[["data_9"]]) {
      values <- t(supplied_inputs[["values_10"]][[as.character(magic_file)]]["arr_0"])
      dimnames(values) <- dimnames(counts)
      qsave(values, output_path(file.path(base, "matrix/matrix_MAGIC.qs")))
  }
  invisible(as.list(environment()))
}
