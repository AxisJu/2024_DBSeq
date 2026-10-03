# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  library(Matrix)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  meta <- obj@meta.data
  keep <- !is.na(meta$donor_id) & !meta$donor_id %in% c("unassigned", "doublet") & !meta$celltype_level1 %in% 
      c("CHPCs", "EPCs")
  obj <- obj[, rownames(meta)[keep]]
  meta <- obj@meta.data
  x <- GetAssayData(obj, assay = "RNA", layer = "data")
  comparisons <- list(DBSIvShamI = c("DBS_I", "Sham_I"), ShamIvSalineI = c("Sham_I", "Saline_I"), DBSIvC = c("DBS_I", 
      "DBS_C"), ShamIvPTZ = c("Sham_I", "PTZ"))
  moments <- function(z) {
      n <- ncol(z)
      m <- rowMeans(z)
      list(n = n, mean = m, variance = pmax((rowSums(z^2) - n * m^2)/(n - 1), 0))
  }
  for (level in strsplit(opt$levels, ",")[[1]]) {
      column <- paste0("celltype_level", level)
      for (ct in sort(unique(as.character(meta[[column]])))) for (name in names(comparisons)) {
          pair <- comparisons[[name]]
          idx1 <- which(meta[[column]] == ct & meta$Group_L2 == pair[1])
          idx0 <- which(meta[[column]] == ct & meta$Group_L2 == pair[2])
          if (min(length(idx1), length(idx0)) < 2) 
              next
          a <- moments(x[, idx0, drop = FALSE])
          b <- moments(x[, idx1, drop = FALSE])
          pooled <- sqrt(((a$n - 1) * a$variance + (b$n - 1) * b$variance)/(a$n + b$n - 2))
          result <- data.frame(gene = rownames(x), cohend = (b$mean - a$mean)/pooled, effect_direction = "comparison_first_minus_second", 
              effect_unit = "cell_distribution_descriptive", n_reference_cells = a$n, n_treatment_cells = b$n)
          result <- result[order(result$cohend, decreasing = TRUE), ]
          write.csv(result, output_path(paste0("results/Table/CohenD/Level", level, "/CohenD_", name, "_", 
              ct, ".csv")), row.names = FALSE)
      }
  }
  invisible(as.list(environment()))
}
