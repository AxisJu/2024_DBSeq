# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(TFBSTools)
      library(motifmatchr)
      library(JASPAR2024)
      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      library(org.Mm.eg.db)
      library(BSgenome.Mmusculus.UCSC.mm10)
      library(GenomicFeatures)
  })
  gene_file <- opt("genes")
  if (is.null(gene_file)) stop("--genes CSV with gene column is required")
  symbols <- unique(supplied_inputs[["symbols_2"]][[as.character(resolve_input(gene_file))]]$gene)
  mapping <- AnnotationDbi::select(org.Mm.eg.db, symbols, "ENTREZID", "SYMBOL")
  mapping <- unique(mapping[complete.cases(mapping), ])
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  tx <- transcripts(txdb, columns = c("tx_id", "gene_id"))
  entrez <- as.character(unlist(tx$gene_id))
  valid <- lengths(tx$gene_id) == 1L
  tx <- tx[valid]
  tx$entrez <- vapply(as.list(tx$gene_id), as.character, character(1))
  tx <- tx[tx$entrez %in% mapping$ENTREZID]
  prom <- promoters(tx, upstream = as.integer(opt("upstream", "2000")), downstream = as.integer(opt("downstream", 
      "500")))
  prom <- trim(prom)
  pfms <- getMatrixSet(JASPAR2024, list(species = 10090, collection = "CORE", all_versions = FALSE))
  hit <- motifmatchr::matchMotifs(pfms, prom, genome = BSgenome.Mmusculus.UCSC.mm10, out = "scores", p.cutoff = as.numeric(opt("p-cutoff", 
      "5e-5")))
  score <- as.matrix(motifmatchr::motifScores(hit))
  idx <- which(score > 0, arr.ind = TRUE)
  result <- data.frame(ENTREZID = tx$entrez[idx[, 1]], transcript_id = tx$tx_id[idx[, 1]], motif_id = names(pfms)[idx[, 
      2]], TF = vapply(pfms, TFBSTools::name, character(1))[idx[, 2]], score = score[idx])
  result <- merge(mapping, result, by = "ENTREZID", all.y = TRUE)
  out <- file.path(output_root, "500_GRN")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(mapping, file.path(out, "promoter_gene_mapping.csv"), row.names = FALSE)
  write.csv(result, file.path(out, "promoter_motif_hits.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
