# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  obj_all <- supplied_inputs[["obj_all_3"]][[as.character(input_path(opt$input))]]
  comparisons <- list(dbssham = c("DBS_I", "Sham_I"), salinesham = c("Saline_I", "Sham_I"), dbsic = c("DBS_I", 
      "DBS_C"))
  for (name in names(comparisons)) {
      obj <- obj_all[, colnames(obj_all)[obj_all$Group_L2 %in% comparisons[[name]]]]
      obj$Group_L2 <- as.character(obj$Group_L2)
      DefaultAssay(obj) <- "RNA"
      for (level in strsplit(opt$levels, ",")[[1]]) {
          column <- paste0("celltype_level", level)
          if (opt$method == "Augur") {
              result <- Augur::calculate_auc(obj, cell_type_col = column, label_col = "Group_L2", n_subsamples = 50, 
                  subsample_size = 20, folds = 3, rf_params = list(mtry = 2, trees = 50, importance = "accuracy"), 
                  n_threads = as.integer(opt$threads))
              filename <- paste0("Augur/augur_", name, "_level", level, "_true.qs")
          }
          else if (opt$method == "scDist") {
              valid <- !is.na(obj$donor_id) & !obj$donor_id %in% c("unassigned", "doublet")
              obj <- obj[, which(valid)]
              result <- scDist::scDist(normalized_counts = GetAssayData(obj, assay = "RNA", layer = "data"), 
                  meta.data = obj@meta.data, fixed.effects = "Group_L2", random.effects = "donor_id", clusters = column, 
                  d = 30, truncate = TRUE, min.count.per.cell = 50)
              filename <- paste0("scDist/scDist_", name, "_level", level, ".qs")
          }
          else stop("method must be Augur or scDist")
          qsave(result, output_path(file.path("data/snRNAseq_mouse/processed/intermediate", filename)))
      }
  }
  invisible(as.list(environment()))
}
