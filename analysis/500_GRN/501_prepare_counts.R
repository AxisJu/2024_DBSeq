# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(Seurat)
      library(qs)
      library(Matrix)
  })
  obj <- supplied_inputs[["obj_2"]][[as.character(resolve_input(opt("object", "data/snRNAseq_mouse/processed/matrix/running_all_250704.qs")))]]
  m <- obj[[]]
  required <- c("Group_L2", "donor_id", "Sample", "celltype_level1", "celltype_level3", "brainregion_projection")
  stopifnot(all(required %in% names(m)), !anyDuplicated(rownames(m)))
  regions <- c("RSPd", "DG", "GPe", "CP", "LSc", "RT", "TRS", "MH", "LH", "RE", "PF", "CM")
  keep <- m$Group_L2 %in% c("DBS_I", "Sham_I", "Saline_I") & !is.na(m$donor_id) & !m$donor_id %in% c("unassigned", 
      "doublet", "") & m$celltype_level1 %in% c("ExN", "InN") & (m$celltype_level3 == "ExN_ATN" | m$brainregion_projection %in% 
      regions)
  keep[is.na(keep)] <- FALSE
  m <- m[keep, , drop = FALSE]
  m$Cell_ID <- rownames(m)
  m$final_region <- ifelse(m$celltype_level3 == "ExN_ATN", "ATN", as.character(m$brainregion_projection))
  m$clean_barcode <- sub(".*([ACGT]{16}-[0-9]+)$", "\\1", m$Cell_ID)
  x <- GetAssayData(obj, assay = "RNA", slot = "counts")[, m$Cell_ID, drop = FALSE]
  stopifnot(identical(colnames(x), m$Cell_ID), all(x@x >= 0), all(x@x == round(x@x)))
  out <- file.path(output_root, "500_GRN", "counts")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(x, file.path(out, "counts.mtx"))
  write.csv(m, file.path(out, "metadata.csv"), row.names = FALSE)
  write.csv(data.frame(gene = rownames(x)), file.path(out, "genes.csv"), row.names = FALSE)
  qs::qsave(x, file.path(out, "counts.qs"))
  invisible(as.list(environment()))
}
