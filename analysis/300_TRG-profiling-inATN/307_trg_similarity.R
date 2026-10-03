# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(dplyr)
  library(tidyr)
  res <- c(supplied_inputs[["res_3"]][[as.character(input_path(opt$results))]], supplied_inputs[["res_4"]][[as.character(input_path(opt$atn_results))]])
  res <- res[!names(res) %in% c("ExN_THns", "InN_Immature", "InN_Nonspecific", "ExN_ATN")]
  if (anyDuplicated(names(res))) stop("Duplicate cell type names")
  all_genes <- Reduce(union, lapply(res, function(x) x$genes))
  matrix <- vapply(res, function(x) {
      stopifnot(!anyDuplicated(x$genes))
      values <- x$observed
      if ("shuffle_mean" %in% names(x)) 
          values <- values - x$shuffle_mean
      setNames(values, x$genes)[all_genes]
  }, numeric(length(all_genes)))
  rownames(matrix) <- all_genes
  trgs <- lapply(res, function(x) unique(x$genes[!is.na(x$is_TRG) & x$is_TRG]))
  if (all(c("AV", "AM") %in% colnames(matrix))) {
      matrix <- cbind(matrix, AV_AM = rowMeans(matrix[, c("AV", "AM"), drop = FALSE]))
      matrix <- matrix[, !colnames(matrix) %in% c("AV", "AM"), drop = FALSE]
      trgs$AV_AM <- union(trgs$AV, trgs$AM)
      trgs <- trgs[!names(trgs) %in% c("AV", "AM")]
  }
  genes <- intersect(Reduce(union, trgs), rownames(matrix))
  correlation <- cor(matrix[genes, , drop = FALSE], method = "spearman", use = "pairwise.complete.obs")
  overlap <- outer(names(trgs), names(trgs), Vectorize(function(a, b) {
      denominator <- min(length(trgs[[a]]), length(trgs[[b]]))
      if (denominator == 0) 
          NA_real_
      else length(intersect(trgs[[a]], trgs[[b]]))/denominator
  }))
  dimnames(overlap) <- list(names(trgs), names(trgs))
  write.csv(correlation, output_path("results/Table/TRG_similarity_celltype.csv"))
  write.csv(overlap, output_path("results/Table/TRG_overlap_coefficient.csv"))
  obj <- supplied_inputs[["obj_5"]][[as.character(input_path(opt$input))]]
  mapping <- unique(obj@meta.data[c("celltype_level3", "brainregion_projection")])
  mapping <- mapping[mapping$brainregion_projection != "/" & mapping$celltype_level3 %in% colnames(matrix), 
      ]
  mapping <- rbind(mapping, data.frame(celltype_level3 = c("AD", "AV_AM"), brainregion_projection = c("AD", 
      "AV_AM")))
  long <- as.data.frame(matrix) %>% tibble::rownames_to_column("gene") %>% pivot_longer(-gene, names_to = "celltype_level3", 
      values_to = "score")
  projected <- inner_join(long, mapping, by = "celltype_level3", relationship = "many-to-many") %>% group_by(gene, 
      brainregion_projection) %>% summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = brainregion_projection, 
      values_from = score)
  projected <- as.data.frame(projected)
  rownames(projected) <- projected$gene
  projected$gene <- NULL
  proj_cor <- cor(as.matrix(projected[intersect(genes, rownames(projected)), , drop = FALSE]), method = "spearman", 
      use = "pairwise.complete.obs")
  write.csv(proj_cor, output_path("results/Table/TRG_similarity_projection.csv"))
  write.csv(projected[intersect(genes, rownames(projected)), , drop = FALSE], output_path("results/Table/TRG_projection_scores.csv"))
  invisible(as.list(environment()))
}
