# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  parse_options <- function(defaults = list()) {
      defaults <- modifyList(list(data_root = Sys.getenv("DBSEQ_DATA_ROOT", getwd()), output_root = Sys.getenv("DBSEQ_OUTPUT_ROOT", 
          file.path(path.expand("~"), "DBSeq-results", "publication-analysis")), external_root = Sys.getenv("DBSEQ_EXTERNAL_ROOT", 
          ""), seed = "777", threads = "4", force = FALSE), defaults)
      argv <- commandArgs(trailingOnly = TRUE)
      if ("--help" %in% argv || "-h" %in% argv) {
          cat("Options (value follows each option):\n")
          for (key in names(defaults)) cat("  --", gsub("_", "-", key), "  ", as.character(defaults[[key]]), 
              "\n", sep = "")
          return(invisible(NULL))
      }
      i <- 1L
      while (i <= length(argv)) {
          key <- gsub("-", "_", sub("^--", "", argv[i]))
          if (!key %in% names(defaults)) 
              stop("Unknown option: ", argv[i])
          if (is.logical(defaults[[key]])) {
              defaults[[key]] <- TRUE
              i <- i + 1L
          }
          else {
              if (i == length(argv)) 
                  stop("Missing option value")
              defaults[[key]] <- argv[i + 1L]
              i <- i + 2L
          }
      }
      defaults$data_root <- normalizePath(defaults$data_root, winslash = "/", mustWork = FALSE)
      defaults$output_root <- normalizePath(defaults$output_root, winslash = "/", mustWork = FALSE)
      package <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = FALSE)
      if (startsWith(tolower(paste0(defaults$output_root, "/")), tolower(paste0(package, "/")))) 
          stop("Output root must be outside the code package")
      if (!nzchar(defaults$external_root)) 
          defaults$external_root <- file.path(defaults$data_root, "external")
      set.seed(as.integer(defaults$seed))
      Sys.setenv(OMP_NUM_THREADS = defaults$threads, OPENBLAS_NUM_THREADS = defaults$threads, MKL_NUM_THREADS = defaults$threads)
      defaults
  }
  input_path <- function(relative) {
      if (grepl("^(/|[[:alpha:]]:)", relative)) 
          return(relative)
      if (startsWith(relative, "external/")) 
          return(file.path(opt$external_root, sub("^external/", "", relative)))
      generated <- file.path(opt$output_root, relative)
      if (supplied_inputs[["input_path_1"]]) 
          generated
      else file.path(opt$data_root, relative)
  }
  output_path <- function(relative) {
      if (grepl("^(/|[[:alpha:]]:)", relative) || ".." %in% strsplit(gsub("\\\\", "/", relative), "/", 
          fixed = TRUE)[[1]]) 
          stop("Output paths must remain inside output root")
      path <- file.path(opt$output_root, relative)
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      path
  }
  invisible(as.list(environment()))
}
