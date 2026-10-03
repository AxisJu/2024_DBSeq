# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(RRHO2)
  library(qs)
  if (!nzchar(opt$first) || !nzchar(opt$second)) stop("Provide --first and --second differential-expression CSVs")
  load_de <- function(path) {
      df <- supplied_inputs[["df_3"]][[as.character(input_path(path))]]
      gene <- if ("gene" %in% names(df)) 
          df$gene
      else df[[1]]
      p <- if ("PValue" %in% names(df)) 
          df$PValue
      else df$de_pval
      effect <- if ("logFC" %in% names(df)) 
          df$logFC
      else df$de_coef
      result <- data.frame(gene, p, effect)
      result <- result[is.finite(result$p) & is.finite(result$effect), ]
      if (anyDuplicated(result$gene)) 
          stop("Duplicate genes in DE input")
      result
  }
  a <- supplied_inputs[["a_4"]][[as.character(opt$first)]]
  b <- supplied_inputs[["b_5"]][[as.character(opt$second)]]
  common <- intersect(a$gene, b$gene)
  if (length(common) < as.integer(opt$stepsize) * 2) stop("Insufficient common genes")
  rank <- function(df, flip = FALSE) {
      df <- df[df$gene %in% common, ]
      values <- -log10(df$p + 1e-300) * sign(df$effect) * ifelse(flip, -1, 1)
      out <- data.frame(gene = df$gene, metric = values)
      out[order(out$metric, decreasing = TRUE), ]
  }
  one <- rank(a)
  two <- rank(b, opt$flip_second)
  results <- lapply(c("hyper", "fisher"), function(method) RRHO2_initialize(one, two, labels = c("first", 
      "second"), log10.ind = TRUE, boundary = 0, method = method, stepsize = as.integer(opt$stepsize)))
  names(results) <- c("hyper", "fisher")
  qsave(results, output_path("results/Table/RRHO2/overlap.qs"))
  write.csv(results$hyper$hypermat, output_path("results/Table/RRHO2/hypergeometric.csv"))
  write.csv(results$fisher$hypermat, output_path("results/Table/RRHO2/fisher_odds.csv"))
  invisible(as.list(environment()))
}
