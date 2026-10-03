# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || length(args) != 1L) {
      cat("Usage: Rscript 10_feature_expression.R <pseudobulk directory>; set DBSEQ_OUTPUT_ROOT\n")
      return(invisible(NULL))
  }
  suppressPackageStartupMessages({
      library(edgeR)
      library(Matrix)
  })
  input <- normalizePath(args[[1]], winslash = "/")
  root <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(root)) stop("Set DBSEQ_OUTPUT_ROOT")
  entry <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), winslash = "/")
  public <- dirname(dirname(dirname(dirname(entry))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, winslash = "/")
  if (startsWith(tolower(paste0(root, "/")), tolower(paste0(public, "/")))) stop("Output is inside code package")
  out <- file.path(root, "600_Others", "260909-DOLPHIN", paste0(basename(input), "_expression"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  x <- supplied_inputs[["x_1"]][[as.character(file.path(input, "counts.mtx"))]]
  m <- supplied_inputs[["m_2"]][[as.character(file.path(input, "samples.csv"))]]
  features <- supplied_inputs[["features_3"]][[as.character(file.path(input, "features.csv"))]]
  stopifnot(ncol(x) == nrow(m), nrow(x) == nrow(features), !anyDuplicated(features$feature_id))
  rownames(x) <- features$feature_id
  results <- list()
  for (region in sort(unique(m$final_region))) {
      idx <- which(m$final_region == region & m$Group_L2 %in% c("Sham_I", "DBS_I"))
      mm <- m[idx, , drop = FALSE]
      if (anyDuplicated(mm$donor_id)) 
          stop("Donor duplicated within a region; review pairing")
      group <- factor(mm$Group_L2, levels = c("Sham_I", "DBS_I"))
      if (any(table(group) < 2L)) 
          next
      design <- model.matrix(~group)
      y <- DGEList(as.matrix(x[, idx, drop = FALSE]))
      keep <- filterByExpr(y, design)
      y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
      y <- estimateDisp(y, design, robust = TRUE)
      fit <- glmQLFit(y, design, robust = TRUE)
      table <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
      table$feature_id <- rownames(table)
      table$region <- region
      results[[region]] <- merge(table, features, by = "feature_id", all.x = TRUE)
  }
  if (!length(results)) stop("No region has two donors per condition")
  result <- do.call(rbind, results)
  result$FDR_all_regions_features <- p.adjust(result$PValue, "BH")
  write.csv(result, file.path(out, "feature_expression_tests.csv"), row.names = FALSE)
  if ("gene_id" %in% names(result)) {
      groups <- split(result, interaction(result$region, result$gene_id, drop = TRUE))
      gene <- do.call(rbind, lapply(groups, function(d) {
          data.frame(region = d$region[1], gene_id = d$gene_id[1], n_features = nrow(d), P = min(1, min(d$PValue) * 
              nrow(d)))
      }))
      gene$FDR_all_regions_genes <- p.adjust(gene$P, "BH")
      write.csv(gene, file.path(out, "gene_any_feature_expression_tests.csv"), row.names = FALSE)
  }
  invisible(as.list(environment()))
}
