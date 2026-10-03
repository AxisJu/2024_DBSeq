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
  library(dplyr)
  library(tidyr)
  human_types <- c("Glu.N", "GABA.N", "OPCs", "Oligo.", "Microglia", "Astro.", "Endo.")
  cohorts <- c(cohort1 = "cohort1_MTLE_brainbenigh_temporalcortex", cohort2 = "cohort2_FCD2b_selfcontrol_frontalcortex", 
      cohort3 = "cohort3_MTLE_autopsy_hippocampus", cohort4 = "cohort4_MTLE_autopsy_amygdala")
  effects <- function(x, groups) {
      i <- which(groups == "Control")
      j <- which(groups == "Epilepsy")
      n1 <- length(i)
      n2 <- length(j)
      mean1 <- rowMeans(x[, i, drop = FALSE])
      mean2 <- rowMeans(x[, j, drop = FALSE])
      if (min(n1, n2) < 2) 
          return(rep(NA_real_, nrow(x)))
      var1 <- pmax((rowSums(x[, i, drop = FALSE]^2) - n1 * mean1^2)/(n1 - 1), 0)
      var2 <- pmax((rowSums(x[, j, drop = FALSE]^2) - n2 * mean2^2)/(n2 - 1), 0)
      (mean2 - mean1)/sqrt(((n1 - 1) * var1 + (n2 - 1) * var2)/(n1 + n2 - 2))
  }
  combined <- list()
  units <- list()
  for (cohort in names(cohorts)) {
      obj <- supplied_inputs[["obj_3"]][[as.character(input_path(paste0("data/snRNAseq_human_DRE/", cohorts[[cohort]], 
          ".qs")))]]
      md <- obj@meta.data
      if (!opt$donor_column %in% names(md)) 
          stop("Missing biological donor column: ", opt$donor_column)
      for (ct in human_types) {
          cells <- which(md$celltype_coarse == ct & md$Group %in% c("Control", "Epilepsy") & !is.na(md[[opt$donor_column]]))
          if (length(cells) < 50) 
              next
          meta <- md[cells, , drop = FALSE]
          x <- GetAssayData(obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
          pair <- c("Control", "Epilepsy")
          if (!all(pair %in% meta$Group)) 
              next
          genes <- rowMeans(x[, meta$Group == pair[1], drop = FALSE] > 0) > 0.1 & rowMeans(x[, meta$Group == 
              pair[2], drop = FALSE] > 0) > 0.1
          key <- paste(meta$Group, meta[[opt$donor_column]], sep = "::")
          ids <- unique(key)
          incidence <- sparseMatrix(i = seq_along(key), j = match(key, ids), x = 1)
          pb <- x %*% incidence
          colnames(pb) <- ids
          info <- meta[match(ids, key), , drop = FALSE]
          rownames(info) <- ids
          info$donor <- as.character(info[[opt$donor_column]])
          for (column in intersect(c("Group", "Age", "Gender"), names(info))) {
              if (any(vapply(split(meta[[column]], key), function(z) length(unique(na.omit(z))) > 1, logical(1)))) 
                  stop("Conflicting donor metadata: ", column)
          }
          retained <- tabulate(match(key, ids), nbins = length(ids)) >= 10
          pb <- pb[genes, retained, drop = FALSE]
          info <- info[retained, , drop = FALSE]
          n <- table(factor(info$Group, levels = pair))
          if (!all(n > 0) || !nrow(pb)) 
              next
          descriptive <- data.frame(Gene = rownames(pb), Control_mean_count = rowMeans(pb[, info$Group == 
              "Control", drop = FALSE]), Epilepsy_mean_count = rowMeans(pb[, info$Group == "Epilepsy", 
              drop = FALSE]), Control_donors = n[1], Epilepsy_donors = n[2])
          write.csv(descriptive, output_path(paste0("results/Table/human_cohort_de/descriptive/", cohort, 
              "_", ct, ".csv")), row.names = FALSE)
          status <- if (any(n < 2)) 
              "descriptive_only_insufficient_donors"
          else "eligible"
          units[[length(units) + 1]] <- data.frame(cohort, celltype = ct, n_control = n[1], n_epilepsy = n[2], 
              status)
          if (any(n < 2)) 
              next
          info$Group <- factor(info$Group, levels = pair)
          covariates <- character()
          repeated <- anyDuplicated(info$donor) > 0
          if (repeated) 
              covariates <- "factor(donor)"
          for (name in intersect(c("Age", "Gender"), names(info))) {
              if (anyNA(info[[name]]) || length(unique(info[[name]])) < 2) 
                  next
              if (name == "Age") 
                  info[[name]] <- as.numeric(as.character(info[[name]]))
              trial <- model.matrix(as.formula(paste("~", paste(c(covariates, name, "Group"), collapse = "+"))), 
                  info)
              if (qr(trial)$rank == ncol(trial) && nrow(trial) > ncol(trial)) 
                  covariates <- c(covariates, name)
          }
          design <- model.matrix(as.formula(paste("~", paste(c(covariates, "Group"), collapse = "+"))), 
              info)
          if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) 
              next
          y <- calcNormFactors(DGEList(pb, samples = info))
          y <- estimateDisp(y, design)
          fit <- glmQLFit(y, design, robust = TRUE)
          result <- topTags(glmQLFTest(fit, coef = ncol(design)), n = Inf, adjust.method = "BH")$table
          result$Gene <- rownames(result)
          d <- setNames(effects(cpm(y, log = TRUE), as.character(info$Group)), rownames(y))
          result$CohenD <- d[result$Gene]
          result$effect_direction <- "Epilepsy_minus_Control"
          result$biological_unit <- "Patient"
          result$Cohort <- cohort
          result$Human_Celltype <- ct
          write.csv(result, output_path(paste0("results/Table/human_cohort_de/pb/PB_DE_", cohort, "_", 
              ct, ".csv")), row.names = FALSE)
          combined[[length(combined) + 1]] <- result
          if (opt$exploratory_sc) {
              single <- obj[, rownames(meta)]
              Idents(single) <- "Group"
              exploratory <- FindMarkers(single, ident.1 = "Epilepsy", ident.2 = "Control", test.use = "MAST", 
                  max.cells.per.ident = 2000, logfc.threshold = 0.1, min.pct = 0.1, random.seed = as.integer(opt$seed))
              write.csv(exploratory, output_path(paste0("results/Table/human_cohort_de/exploratory_sc/SC_DE_", 
                  cohort, "_", ct, ".csv")))
          }
      }
  }
  if (length(units)) write.csv(bind_rows(units), output_path("results/Table/human_cohort_de/biological_units.csv"), 
      row.names = FALSE)
  if (!length(combined)) stop("No cohort had an estimable donor-level contrast")
  human <- bind_rows(combined)
  orthologs <- supplied_inputs[["orthologs_4"]][[as.character(input_path(opt$orthologs))]]
  needed <- c("external_gene_name", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type")
  stopifnot(all(needed %in% names(orthologs)))
  orthologs <- orthologs %>% filter(mmusculus_homolog_orthology_type == "ortholog_one2one", external_gene_name != 
      "", mmusculus_homolog_associated_gene_name != "") %>% distinct()
  mapped <- inner_join(human, orthologs, by = c(Gene = "external_gene_name"))
  write.csv(mapped, output_path("results/Table/human_cohort_de/conserved_trg_candidates.csv"), row.names = FALSE)
  non_neuronal <- c(Astro = "Astro.", OPC = "OPCs", Oligo = "Oligo.", Micro = "Microglia", VCs = "Endo.")
  mouse_tables <- lapply(c("Level1", "Level3"), function(level) {
      tables <- supplied_inputs[["tables_5"]][[as.character(input_path(paste0("data/snRNAseq_mouse/processed/intermediate/CosSim_v2/", 
          level, "_res.qs")))]]
      selected <- if (level == "Level1") 
          intersect(names(tables), names(non_neuronal))
      else grep("^(ExN|InN)_", names(tables), value = TRUE)
      corrected <- lapply(tables[selected], function(frame) {
          if (!"CoD_ShamvSal_effect_direction" %in% names(frame)) {
              frame$CoD_ShamvSal_cohend <- -frame$CoD_ShamvSal_cohend
          }
          else if (!all(frame$CoD_ShamvSal_effect_direction == "comparison_first_minus_second")) {
              stop("Unrecognized mouse Cohen d direction")
          }
          frame
      })
      bind_rows(corrected, .id = "Celltype")
  })
  mouse <- bind_rows(mouse_tables) %>% mutate(Human_Celltype = case_when(grepl("^ExN_", Celltype) ~ "Glu.N", 
      grepl("^InN_", Celltype) ~ "GABA.N", TRUE ~ unname(non_neuronal[Celltype])))
  mouse_summary <- mouse %>% group_by(genes, Human_Celltype) %>% summarise(Mouse_ShamVSal_CohenD = median(CoD_ShamvSal_cohend, 
      na.rm = TRUE), N_Mouse_TRG_Celltypes = sum(is_TRG, na.rm = TRUE), N_Mouse_Celltypes = n_distinct(Celltype), 
      .groups = "drop")
  conservation <- inner_join(mapped, mouse_summary, by = c(mmusculus_homolog_associated_gene_name = "genes", 
      "Human_Celltype"))
  conservation$Direction_Concordant <- sign(conservation$CohenD) == sign(conservation$Mouse_ShamVSal_CohenD)
  write.csv(conservation, output_path("results/Table/human_cohort_de/conserved_trg.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
