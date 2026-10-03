# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(fgsea)
      library(glue)
      library(qs)
      library(dplyr)
      library(stringr)
      library(aPEAR)
      library(ggplot2)
  })
  base_dir <- data_root
  output_dir_fig <- output_path("figures/figure 3")
  output_dir_data <- output_path("sourcedata/figure 3")
  dir.create(output_dir_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  gmt_path <- file.path(base_dir, "bin/MSigDB_Mm/msigdb.v2023.2.Mm.symbols.gmt")
  memento_dbs_path <- results_path("Table/MEMENTO_v2/Level3/MEMENTO_1d_DBSvSham_ExN_ATN.csv")
  memento_sham_path <- results_path("Table/MEMENTO_v2/Level3/MEMENTO_1d_ShamvSaline_ExN_ATN.csv")
  cossim_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/CosSim_v2/level3_res.qs")
  pathways <- gmtPathways(gmt_path)
  memento_res1 <- supplied_inputs[["memento_res1_2"]][[as.character(memento_dbs_path)]]
  memento_res2 <- supplied_inputs[["memento_res2_3"]][[as.character(memento_sham_path)]]
  expressed_genes <- intersect(rownames(memento_res1), rownames(memento_res2))
  gene_rank <- supplied_inputs[["gene_rank_4"]][[as.character(cossim_path)]][["ExN_ATN"]] %>% rename_with(~ifelse(grepl("^median_cos", 
      .x), "median_costheta", .x)) %>% filter(genes %in% expressed_genes) %>% arrange(desc(median_costheta)) %>% 
      dplyr::select(genes, median_costheta)
  rank <- setNames(gene_rank$median_costheta, gene_rank$genes)
  set.seed(42)
  fgsea_res <- fgsea(pathways, rank, minSize = 20, maxSize = 100, nproc = 1)
  res_tp <- fgsea_res %>% dplyr::filter(str_detect(pathway, "GOBP_")) %>% dplyr::filter(padj < 0.05) %>% 
      dplyr::rowwise() %>% dplyr::mutate(Description = str_remove(pathway, "GOBP_") %>% str_to_title(), 
      pathwayGenes = paste(unlist(leadingEdge), collapse = "/"), colorBy = NES, nodeSize = size) %>% dplyr::ungroup() %>% 
      dplyr::select(pathway, Description, NES, pval, padj, size, pathwayGenes, colorBy, nodeSize) %>% as.data.frame()
  sourcedata_a <- res_tp %>% dplyr::select(Pathway = pathway, Description, NES, PValue = pval, FDR = padj, 
      Size = size, LeadingEdge_Genes = pathwayGenes)
  write.csv(sourcedata_a, file.path(output_dir_data, "3a.csv"), row.names = FALSE)
  cat(sprintf("Exported Source Data: %s (Rows: %d)\n", file.path(output_dir_data, "3a.csv"), nrow(sourcedata_a)))
  set.seed(123)
  p <- enrichmentNetwork(res_tp %>% dplyr::select(Description, pathwayGenes, colorBy, nodeSize), repelLabels = TRUE, 
      fontSize = 3, minClusterSize = 3, colorBy = "colorBy", nodeSize = "nodeSize", drawEllipses = TRUE)
  fig_pdf <- file.path(output_dir_fig, "3a.pdf")
  fig_png <- file.path(output_dir_fig, "3a.png")
  ggsave(fig_png, p, width = 661/72, height = 409/72, dpi = 300)
  invisible(as.list(environment()))
}
