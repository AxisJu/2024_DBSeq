# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  if (!nzchar(opt$ed)) stop("Provide --ed with the labeled Allen distance matrix")
  scores <- supplied_inputs[["scores_3"]][[as.character(input_path(opt$scores))]]
  targets <- strsplit(opt$targets, ",")[[1]]
  target_metric <- function(path) {
      table <- supplied_inputs[["table_4"]][[as.character(input_path(path))]]
      rownames(table) <- sub("_ipsi$", "", rownames(table))
      names(table) <- sub("_ipsi$", "", names(table))
      if (anyDuplicated(rownames(table)) || anyDuplicated(names(table))) 
          stop("Nonunique metric region labels")
      if (!all(targets %in% names(table)) || !all(targets %in% rownames(table))) 
          stop("Missing target regions")
      values <- rowMeans(table[, targets, drop = FALSE])
      values["AV_AM"] <- mean(values[targets])
      values[!names(values) %in% targets]
  }
  ne <- target_metric(opt$ne)
  ed <- target_metric(opt$ed)
  regions <- Reduce(intersect, list(colnames(scores), names(ne), names(ed)))
  if (length(regions) < 5) stop("Fewer than five common regions")
  fit <- function(gene) {
      df <- data.frame(y = as.numeric(scores[gene, regions]), NE = ne[regions], ED = ed[regions])
      df <- df[complete.cases(df) & apply(df, 1, function(z) all(is.finite(z))), ]
      if (nrow(df) < 5) 
          return(NULL)
      model <- lm(y ~ NE + ED, df)
      if (model$rank != 3) 
          return(NULL)
      co <- summary(model)$coefficients
      t <- co["NE", "t value"]
      data.frame(gene, beta = co["NE", "Estimate"], se = co["NE", "Std. Error"], t = t, p = co["NE", "Pr(>|t|)"], 
          r2_partial = t^2/(t^2 + df.residual(model)), n_regions = nrow(df))
  }
  result <- do.call(rbind, lapply(rownames(scores), fit))
  if (is.null(result)) stop("No estimable gene models")
  result$q <- p.adjust(result$p, method = "BH")
  write.csv(result[order(result$q), ], output_path("results/Table/TRG_projection_gene_models.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
