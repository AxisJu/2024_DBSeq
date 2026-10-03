# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  library(qs)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  data_root <- Sys.getenv("DBSEQ_DATA_ROOT")
  workspace <- Sys.getenv("DBSEQ_OUTPUT_ROOT")
  if (!nzchar(data_root) || !nzchar(workspace)) stop("Set DBSEQ_DATA_ROOT and DBSEQ_OUTPUT_ROOT")
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  code_base <- normalizePath(file.path(dirname(normalizePath(entry)), "../.."), winslash = "/")
  out_base <- normalizePath(workspace, winslash = "/", mustWork = FALSE)
  if (identical(tolower(out_base), tolower(code_base)) || startsWith(tolower(out_base), paste0(tolower(code_base), 
      "/"))) stop("Output must be outside the code package")
  output <- file.path(workspace, "cache/extended_data")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  trg <- Sys.getenv("DBSEQ_TRG_FILE", file.path(data_root, "results/Table/CosSim_v2/plot_df.qs"))
  mouse <- if (grepl("\\.qs$", trg)) supplied_inputs[["mouse_1"]][[as.character(trg)]] else supplied_inputs[["mouse_2"]][[as.character(trg)]]
  neuron <- unique(as.character(mouse$Celltype[mouse$is_TRG & grepl("^(ExN|InN)_", mouse$Celltype)]))
  object <- supplied_inputs[["object_3"]][[as.character(file.path(data_root, "data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]]
  meta <- object@meta.data
  rm(object)
  gc()
  meta <- meta %>% filter(Group_L2 %in% c("DBS_I", "Sham_I", "Saline_I"), !donor_id %in% c("unassigned", 
      "doublet"), celltype_level1 %in% c("ExN", "InN")) %>% semi_join(data.frame(celltype_level3 = neuron), 
      by = "celltype_level3")
  mapping <- meta %>% distinct(celltype_level3, brainregion_projection) %>% filter(brainregion_projection != 
      "/", celltype_level3 != "ExN_ATN")
  mapped <- mouse %>% filter(is_TRG, Celltype %in% neuron) %>% left_join(mapping, by = c(Celltype = "celltype_level3")) %>% 
      mutate(Region = ifelse(is.na(brainregion_projection) | brainregion_projection == "/", "Unknown", 
          brainregion_projection)) %>% group_by(genes) %>% mutate(n_neurons = n_distinct(Celltype), is_in_ATN = "ExN_ATN" %in% 
      Celltype) %>% ungroup() %>% filter(is_in_ATN)
  records <- mapped %>% filter(Region != "Unknown", Celltype != "ExN_ATN") %>% select(genes, Celltype, 
      Region, n_neurons) %>% arrange(genes, Celltype, Region)
  write.csv(records, file.path(output, "ed13_membership.csv"), row.names = FALSE)
  cross <- supplied_inputs[["cross_4"]][[as.character(file.path(Sys.getenv("DBSEQ_ANALYSIS_OUTPUT_ROOT", 
      data_root), "data/snRNAseq_mouse/processed/intermediate/CosSim_v2/Level3_res.qs"))]]
  effect <- bind_rows(cross, .id = "Celltype")
  mean_column <- grep("^mean_cos", names(effect), value = TRUE)
  stopifnot(length(mean_column) == 1)
  effect$expr <- effect[[mean_column]] - effect$shuffle_mean
  effect <- effect %>% select(Celltype, genes, expr)
  frontier <- records %>% left_join(effect, by = c("genes", "Celltype")) %>% filter(expr > 0)
  write.csv(frontier, file.path(output, "ed13_frontier_membership.csv"), row.names = FALSE)
  cat("Prepared Extended Data Figure 13 gene membership and propagation weights\n")
  invisible(as.list(environment()))
}
