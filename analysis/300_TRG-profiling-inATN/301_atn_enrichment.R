# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(fgsea)
  all_results <- supplied_inputs[["all_results_3"]][[as.character(input_path(opt$results))]]
  res <- all_results[[opt$celltype]]
  if (is.null(res)) stop("Cell type missing from TRG results")
  stopifnot(!anyDuplicated(res$genes))
  tested <- !is.na(res$Mem_DBSvSham_de_coef) & !is.na(res$Mem_ShamvSal_de_coef)
  rank <- setNames(res[[opt$rank_statistic]][tested], res$genes[tested])
  rank <- sort(rank[is.finite(rank)], decreasing = TRUE)
  pathways <- gmtPathways(input_path(opt$pathways))
  enriched <- fgsea(pathways, rank, minSize = 20, maxSize = 100, nproc = as.integer(opt$threads))
  enriched$leadingEdge <- vapply(enriched$leadingEdge, paste, collapse = ",", character(1))
  write.csv(enriched, output_path(paste0("results/Table/gsea_", opt$celltype, "_CosSim.csv")), row.names = FALSE)
  invisible(as.list(environment()))
}
