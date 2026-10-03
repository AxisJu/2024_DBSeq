# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
      cat("Set DBSEQ_DATA_ROOT and DBSEQ_OUTPUT_ROOT. Optional arguments: --name value.\n")
      return(invisible(NULL))
  }
  opt <- function(name, default = NULL) {
      i <- match(paste0("--", name), args)
      if (is.na(i)) 
          return(default)
      if (i == length(args)) 
          stop("Missing value: ", name)
      args[[i + 1L]]
  }
  data_root <- Sys.getenv("DBSEQ_DATA_ROOT")
  output_root <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(data_root) || !nzchar(output_root)) stop("Set DBSEQ_DATA_ROOT and DBSEQ_OUTPUT_ROOT")
  entry <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), winslash = "/")
  package_root <- dirname(entry)
  while (basename(package_root) != "analysis" && dirname(package_root) != package_root) package_root <- dirname(package_root)
  package_root <- dirname(package_root)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  output_root <- normalizePath(output_root, winslash = "/")
  if (startsWith(tolower(paste0(output_root, "/")), tolower(paste0(package_root, "/")))) stop("Output is inside code package")
  resolve_input <- function(value) {
      if (grepl("^(/|[A-Za-z]:|~)", value)) 
          path.expand(value)
      else file.path(data_root, value)
  }
  invisible(as.list(environment()))
}
