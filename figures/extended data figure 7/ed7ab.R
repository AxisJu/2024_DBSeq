# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args <- commandArgs(trailingOnly = FALSE)
  entry <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
  code_root <- normalizePath(file.path(dirname(entry), "../.."), winslash = "/")
  package_root <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(package_root)) stop("Set DBSEQ_OUTPUT_ROOT outside the code package.")
  workspace_root <- package_root
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  invisible(as.list(environment()))
}
