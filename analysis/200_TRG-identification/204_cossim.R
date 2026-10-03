# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  library(Matrix)
  groups <- c("DBS_I", "Sham_I", "Saline_I")
  if (!opt$statistic %in% c("mean", "median")) stop("statistic must be mean or median")
  if (!opt$sampling %in% c("donor-bootstrap", "historical-cells")) stop("Unknown sampling mode")
  if (!opt$permutation %in% c("donor", "exploratory-cells")) stop("Unknown permutation mode")
  summarize_cos <- function(x) {
      if (opt$statistic == "mean") 
          rowMeans(x)
      else apply(x, 1, median)
  }
  cosine_values <- function(x, metadata, cells_per_group) {
      totals <- lapply(groups, function(group) {
          rows <- which(metadata$grp == group)
          donors <- unique(metadata$donor[rows])
          if (length(donors) < 2 && opt$sampling == "donor-bootstrap") 
              stop("Fewer than two donors in a group")
          if (opt$sampling == "donor-bootstrap") {
              draws <- sample(donors, length(donors), replace = TRUE)
              n_each <- ceiling(cells_per_group/length(draws))
              values <- lapply(draws, function(donor) {
                  available <- rows[metadata$donor[rows] == donor]
                  take <- available[sample.int(length(available), n_each, replace = length(available) < 
                    n_each)]
                  rowMeans(x[, take, drop = FALSE])
              })
              Reduce(`+`, values)/length(values) * cells_per_group
          }
          else {
              counts <- table(metadata$donor[rows])
              allocation <- floor(counts/sum(counts) * cells_per_group)
              allocation[which.max(counts)] <- allocation[which.max(counts)] + cells_per_group - sum(allocation)
              take <- unlist(lapply(names(allocation), function(donor) {
                  available <- rows[metadata$donor[rows] == donor]
                  available[sample.int(length(available), min(allocation[[donor]], length(available)))]
              }))
              rowSums(x[, take, drop = FALSE])
          }
      })
      dbs <- totals[[1]] - totals[[2]]
      saline <- totals[[3]] - totals[[2]]
      (-1 + dbs * saline)/sqrt((1 + dbs^2) * (1 + saline^2))
  }
  permute_metadata <- function(metadata) {
      if (opt$permutation == "exploratory-cells") {
          metadata$grp <- sample(metadata$grp)
      }
      else {
          donor_groups <- unique(metadata[c("donor", "grp")])
          if (anyDuplicated(donor_groups$donor)) 
              stop("Each donor must belong to one of DBS_I, Sham_I, Saline_I")
          donor_groups$grp <- sample(donor_groups$grp)
          metadata$grp <- donor_groups$grp[match(metadata$donor, donor_groups$donor)]
      }
      metadata
  }
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  meta <- obj@meta.data
  keep <- meta$Group_L2 %in% groups & !is.na(meta$donor_id) & !meta$donor_id %in% c("unassigned", "doublet")
  if ("celltype_level1" %in% names(meta)) keep <- keep & !meta$celltype_level1 %in% c("CHPCs", "EPCs")
  obj <- obj[, rownames(meta)[keep]]
  meta <- obj@meta.data
  exp <- GetAssayData(obj, assay = "RNA", layer = "data")
  columns <- if (opt$scope == "atn") c(ANT = "celltype") else setNames(paste0("celltype_level", strsplit(opt$levels, 
      ",")[[1]]), paste0("Level", strsplit(opt$levels, ",")[[1]]))
  for (level in names(columns)) {
      column <- columns[[level]]
      celltypes <- sort(unique(as.character(meta[[column]])))
      celltypes <- setdiff(celltypes, strsplit(opt$exclude_celltypes, ",", fixed = TRUE)[[1]])
      for (ct in celltypes) {
          idx <- which(meta[[column]] == ct)
          md <- data.frame(grp = as.character(meta$Group_L2[idx]), donor = as.character(meta$donor_id[idx]))
          if (!all(groups %in% md$grp)) 
              next
          if (opt$sampling == "donor-bootstrap" && any(vapply(groups, function(g) length(unique(md$donor[md$grp == 
              g])), integer(1)) < 2)) 
              next
          x <- exp[, idx, drop = FALSE]
          minimum <- min(table(md$grp))
          n_cells <- if (minimum > 500) 
              100
          else if (minimum > 250) 
              50
          else if (minimum > 100 || level == "Level1") 
              25
          else 10
          n_boots <- as.integer(opt$boots)
          n_shuffles <- as.integer(opt$shuffles)
          observed <- replicate(n_boots, cosine_values(x, md, n_cells))
          observed_stat <- summarize_cos(observed)
          null <- vapply(seq_len(n_shuffles), function(i) {
              shuffled <- permute_metadata(md)
              summarize_cos(replicate(n_boots, cosine_values(x, shuffled, n_cells)))
          }, numeric(nrow(x)))
          stats <- data.frame(genes = rownames(x), statistic = opt$statistic, observed = observed_stat, 
              shuffle_mean = rowMeans(null), mean_cosine = rowMeans(observed), median_cosine = apply(observed, 
                  1, median), ci_lower = apply(observed, 1, quantile, probs = 0.025), ci_upper = apply(observed, 
                  1, quantile, probs = 0.975), bootstrap_p = rowMeans(observed <= 0), shuffle_p = (1 + 
                  rowSums(null >= observed_stat))/(n_shuffles + 1), sampling_unit = opt$sampling, permutation_unit = opt$permutation)
          root <- file.path("data/snRNAseq_mouse/processed/intermediate/CosSim", level)
          qsave(list(genes = rownames(x), observed = observed, null = null, metadata = md, options = opt), 
              output_path(file.path(root, paste0("CosSim_", ct, ".qs"))))
          write.csv(stats, output_path(file.path("results/Table/CosSim", level, paste0("CosSim_", ct, ".csv"))), 
              row.names = FALSE)
      }
  }
  invisible(as.list(environment()))
}
