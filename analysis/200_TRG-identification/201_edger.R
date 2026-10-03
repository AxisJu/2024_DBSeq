# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  library(edgeR)
  library(Matrix)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  meta <- obj@meta.data
  keep <- !is.na(meta$donor_id) & !meta$donor_id %in% c("unassigned", "doublet")
  if ("celltype_level1" %in% names(meta)) keep <- keep & !meta$celltype_level1 %in% c("CHPCs", "EPCs")
  obj <- obj[, rownames(meta)[keep]]
  meta <- obj@meta.data
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  comparisons <- list(DBSIvShamI = c("DBS_I", "Sham_I"), ShamIvSalineI = c("Sham_I", "Saline_I"), DBSIvC = c("DBS_I", 
      "DBS_C"), ShamIvPTZ = c("Sham_I", "PTZ"))
  if (opt$scope == "atn") comparisons <- comparisons[1:2]
  columns <- if (opt$scope == "atn") c(ANT = "celltype") else setNames(paste0("celltype_level", strsplit(opt$levels, 
      ",")[[1]]), paste0("Level", strsplit(opt$levels, ",")[[1]]))
  audit <- list()
  for (level in names(columns)) {
      column <- columns[[level]]
      for (ct in sort(unique(as.character(meta[[column]])))) {
          for (comparison in names(comparisons)) {
              groups <- comparisons[[comparison]]
              idx <- which(meta[[column]] == ct & meta$Group_L2 %in% groups)
              sub <- meta[idx, , drop = FALSE]
              if (!all(groups %in% sub$Group_L2)) 
                  next
              x <- counts[, idx, drop = FALSE]
              key <- paste(sub$Group_L2, sub$donor_id, sep = "::")
              units <- unique(key)
              incidence <- sparseMatrix(i = seq_along(key), j = match(key, units), x = 1)
              exp <- x %*% incidence
              colnames(exp) <- units
              n_cells <- tabulate(match(key, units), nbins = length(units))
              metadata <- sub[match(units, key), c("Group_L2", "donor_id"), drop = FALSE]
              rownames(metadata) <- units
              retained <- n_cells >= as.integer(opt$min_cells)
              exp <- exp[, retained, drop = FALSE]
              metadata <- metadata[retained, , drop = FALSE]
              sizes <- table(factor(metadata$Group_L2, levels = groups))
              if (any(sizes < 2)) {
                  audit[[length(audit) + 1]] <- data.frame(level, ct, comparison, status = "insufficient_donors", 
                    donors0 = sizes[2], donors1 = sizes[1])
                  next
              }
              metadata$group <- factor(metadata$Group_L2, levels = rev(groups))
              y <- DGEList(exp, samples = metadata)
              if (opt$scope == "atn") {
                  genes <- filterByExpr(y, group = metadata$group)
              }
              else {
                  fraction <- as.numeric(opt$fraction)
                  genes <- rowMeans(x[, sub$Group_L2 == groups[1], drop = FALSE] > 0) > fraction & rowMeans(x[, 
                    sub$Group_L2 == groups[2], drop = FALSE] > 0) > fraction
              }
              y <- y[genes, , keep.lib.sizes = FALSE]
              if (!nrow(y)) 
                  next
              repeated <- any(table(metadata$donor_id) > 1)
              paired <- opt$paired == "yes" || (opt$paired == "auto" && repeated)
              design <- if (paired) 
                  model.matrix(~factor(donor_id) + group, metadata)
              else model.matrix(~group, metadata)
              if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) 
                  stop("Non-estimable design: ", ct, " ", comparison)
              y <- calcNormFactors(y)
              y <- estimateDisp(y, design)
              fit <- glmQLFit(y, design, robust = TRUE)
              result <- topTags(glmQLFTest(fit, coef = ncol(design)), n = Inf, adjust.method = "BH")$table
              result <- data.frame(gene = rownames(result), result, row.names = NULL)
              root <- if (opt$scope == "atn") 
                  "results/Table/EdgeR/ANT"
              else file.path("results/Table/EdgeR_v2", level)
              write.csv(result, output_path(file.path(root, paste0("EdgeR_", comparison, "_", ct, ".csv"))), 
                  row.names = FALSE)
              audit[[length(audit) + 1]] <- data.frame(level, ct, comparison, status = ifelse(paired, "paired", 
                  "independent"), donors0 = sizes[2], donors1 = sizes[1])
          }
      }
  }
  if (length(audit)) write.csv(do.call(rbind, audit), output_path(paste0("results/Table/EdgeR_", opt$scope, 
      "_units.csv")), row.names = FALSE)
  invisible(as.list(environment()))
}
