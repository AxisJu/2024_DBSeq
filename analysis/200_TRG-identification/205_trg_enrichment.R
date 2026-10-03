# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(fgsea)
  tables <- supplied_inputs[["tables_3"]][[as.character(input_path(paste0("data/snRNAseq_mouse/processed/intermediate/CosSim_v2/", 
      opt$level, "_res.qs")))]]
  pathways <- gmtPathways(input_path(opt$pathways))
  for (celltype in names(tables)) {
      table <- tables[[celltype]]
      if (anyDuplicated(table$genes)) 
          stop("Duplicate gene identifiers")
      tested <- !is.na(table$Mem_DBSvSham_de_coef) & !is.na(table$Mem_ShamvSal_de_coef)
      rank <- setNames(table[[opt$rank_statistic]][tested], table$genes[tested])
      rank <- sort(rank[is.finite(rank)], decreasing = TRUE)
      if (length(rank) < 20) 
          next
      result <- fgsea(pathways, rank, minSize = 20, maxSize = 100, nproc = as.integer(opt$threads))
      result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ",", character(1))
      write.csv(result, output_path(paste0("results/Table/fgsea/", opt$level, "/fgsea_cossim_", celltype, 
          ".csv")), row.names = FALSE)
      if (opt$ora) {
          for (direction in c("Up", "Down")) {
              selected <- table$is_TRG & if (direction == "Up") 
                  table$CoD_DBSvSham_cohend < 0
              else table$CoD_DBSvSham_cohend > 0
              genes <- table$genes[!is.na(selected) & selected]
              if (!length(genes)) 
                  next
              enrichment <- enrichR::enrichr(genes, databases = c("GO_Biological_Process_2023", "KEGG_2019_Mouse"))
              for (database in names(enrichment)) write.csv(enrichment[[database]], output_path(paste0("results/Table/ora/", 
                  opt$level, "/", celltype, "_", direction, "_", database, ".csv")), row.names = FALSE)
          }
      }
  }
  invisible(as.list(environment()))
}
