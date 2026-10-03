# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  entry_dir <- dirname(normalizePath(sub("^--file=", "", entry_arg), winslash = "/"))
  required_root <- function(key) {
      value <- Sys.getenv(key)
      if (!nzchar(value)) 
          stop("Set ", key, " before running this script.")
      normalizePath(value, winslash = "/", mustWork = FALSE)
  }
  data_root <- required_root("DBSEQ_DATA_ROOT")
  output_root <- required_root("DBSEQ_OUTPUT_ROOT")
  code_root <- normalizePath(file.path(entry_dir, "../.."), winslash = "/")
  if (identical(tolower(output_root), tolower(code_root)) || startsWith(tolower(output_root), paste0(tolower(code_root), 
      "/"))) stop("DBSEQ_OUTPUT_ROOT must be outside the public code directory.")
  input_package_root <- Sys.getenv("DBSEQ_INPUT_PACKAGE_ROOT")
  workspace_root <- Sys.getenv("DBSEQ_WORKSPACE_ROOT")
  results_root <- Sys.getenv("DBSEQ_RESULTS_ROOT", file.path(data_root, "results"))
  results_path <- function(...) file.path(results_root, ...)
  data_path <- function(...) {
      relative <- file.path(...)
      if (startsWith(relative, "results/")) 
          results_path(sub("^results/", "", relative))
      else file.path(data_root, relative)
  }
  read_trg_plot <- function() {
      path <- Sys.getenv("DBSEQ_TRG_FILE", results_path("Table/CosSim_v2/plot_df.qs"))
      result <- if (grepl("\\.qs$", path)) 
          supplied_inputs[["result_1"]][[as.character(path)]]
      else supplied_inputs[["result_2"]][[as.character(path)]]
      if ("analysis_mode" %in% names(result)) {
          result$CoD_DBSvSham_cohend <- -result$CoD_DBSvSham_cohend
          result$CoD_ShamvSal_cohend <- -result$CoD_ShamvSal_cohend
      }
      result
  }
  input_path <- function(...) {
      if (!nzchar(input_package_root)) 
          stop("Set DBSEQ_INPUT_PACKAGE_ROOT for frozen inputs.")
      file.path(input_package_root, ...)
  }
  output_path <- function(...) {
      target <- file.path(output_root, ...)
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      target
  }
  cache_path <- function(name) output_path("cache", name)
  code_path <- function(name) file.path(entry_dir, name)
  external_path <- function(key, ...) file.path(required_root(key), ...)
  set.seed(260924)
  if (.Platform$OS.type == "windows") grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
  invisible(as.list(environment()))
}
