# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(qs)
      library(Seurat)
      library(dplyr)
  })
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
  trg_path <- Sys.getenv("DBSEQ_TRG_FILE", file.path(data_root, "results/Table/CosSim_v2/plot_df.qs"))
  mouse <- if (grepl("\\.qs$", trg_path)) supplied_inputs[["mouse_1"]][[as.character(trg_path)]] else supplied_inputs[["mouse_2"]][[as.character(trg_path)]]
  object <- supplied_inputs[["object_3"]][[as.character(file.path(data_root, "data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]]
  mapping <- object@meta.data %>% distinct(celltype_level3, brainregion_projection) %>% filter(brainregion_projection != 
      "/", celltype_level3 != "ExN_ATN")
  rm(object)
  gc()
  neuron <- setdiff(unique(mouse$Celltype), c("Astro", "OPC", "Oligo", "VCs", "Micro"))
  mapped <- mouse %>% filter(is_TRG, Celltype %in% neuron) %>% left_join(mapping, by = c(Celltype = "celltype_level3")) %>% 
      mutate(Region = ifelse(is.na(brainregion_projection) | brainregion_projection == "/", "Unknown", 
          brainregion_projection)) %>% group_by(genes) %>% mutate(is_in_ATN = "ExN_ATN" %in% Celltype) %>% 
      ungroup() %>% filter(is_in_ATN) %>% mutate(Region = ifelse(Celltype == "ExN_ATN", "ATN", Region)) %>% 
      filter(Region != "Unknown") %>% select(Mouse_Gene = genes, Region) %>% distinct()
  write.csv(mapped, file.path(output, "ed12_regional_genes.csv"), row.names = FALSE)
  writeLines(unique(mouse$genes), file.path(output, "ed12_background_genes.txt"))
  cat("Prepared regional TRG membership and measured-gene background\n")
  go <- supplied_inputs[["go_4"]][[as.character(file.path(Sys.getenv("DBSEQ_INPUT_PACKAGE_ROOT"), "sourcedata/figure 5/inputs/go_bp_2023.qs"))]]
  write.csv(data.frame(pathway = names(go), genes = vapply(go, paste, collapse = ";", FUN.VALUE = character(1))), 
      file.path(output, "ed12_all_gene_sets.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
