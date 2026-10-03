# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(neuroestimator)
  })
  python <- Sys.getenv("DBSEQ_NEUROESTIMATOR_PYTHON")
  if (nzchar(python)) reticulate::use_python(python, required = TRUE)
  counts_path <- opt("counts")
  if (is.null(counts_path)) counts_path <- file.path(output_root, "500_GRN", "counts", "counts.qs")
  x <- supplied_inputs[["x_2"]][[as.character(resolve_input(counts_path))]]
  limit <- as.numeric(opt("max-dense-gb", "16"))
  if (prod(dim(x)) * 8/1024^3 > limit) stop("Dense matrix exceeds --max-dense-gb")
  score <- neuroestimator::neuroestimator(as.matrix(x), species = "mmusculus")
  out <- file.path(output_root, "500_GRN", "neuroestimator")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(score, file.path(out, "scores.csv"), row.names = TRUE)
  invisible(as.list(environment()))
}
