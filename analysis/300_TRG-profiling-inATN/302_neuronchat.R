# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(Seurat)
  library(NeuronChat)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path(opt$input))]]
  meta <- obj@meta.data
  excluded <- c("CHPCs", "EPCs", "ExN_CLA", "ExN_CA1", "ExN_CA3", "ExN_HPF_CajalRetzius", "ExN_RSP_L23IT", 
      "InN_OT")
  groups <- c("Saline_I", "Sham_I", "DBS_I")
  valid <- meta$Group_L2 %in% groups & !is.na(meta$donor_id) & !meta$donor_id %in% c("unassigned", "doublet") & 
      !meta$celltype_level3 %in% excluded
  obj <- obj[, rownames(meta)[valid]]
  networks <- list()
  output <- list()
  for (group in groups) {
      selected <- obj[, which(obj$Group_L2 == group)]
      model <- createNeuronChat(GetAssayData(selected, assay = "RNA", layer = "data"), DB = "mouse", group.by = selected$celltype_level3)
      model <- run_NeuronChat(model, M = as.integer(opt$permutations))
      qsave(model, output_path(paste0("data/snRNAseq_mouse/processed/intermediate/NeuronChat/", group, 
          ".qs")))
      total <- sum(vapply(model@net, sum, numeric(1)))
      if (!is.finite(total) || total <= 0) 
          stop("Zero communication total for ", group)
      networks[[group]] <- lapply(model@net, function(x) x * 500/total)
      for (pathway in names(networks[[group]])) {
          values <- networks[[group]][[pathway]]
          if (!opt$target %in% rownames(values) || !opt$target %in% colnames(values)) 
              next
          output[[length(output) + 1]] <- data.frame(Group = group, Pathway = pathway, Celltype = colnames(values), 
              Output = values[opt$target, ], Input = values[, opt$target], Inference_unit = "pooled_cells")
      }
  }
  qsave(networks, output_path("data/snRNAseq_mouse/processed/intermediate/NeuronChat/normalized_networks.qs"))
  if (length(output)) write.csv(do.call(rbind, output), output_path("results/Table/NeuronChat/ATN_pathway_connections.csv"), 
      row.names = FALSE)
  invisible(as.list(environment()))
}
