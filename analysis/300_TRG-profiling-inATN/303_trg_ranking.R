# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  library(qs)
  library(dplyr)
  trgs <- supplied_inputs[["trgs_3"]][[as.character(input_path(opt$results))]]
  long <- bind_rows(trgs, .id = "Celltype")
  long <- long[!duplicated(long[c("Celltype", "genes")]), ]
  ranking <- long %>% mutate(Regulation = ifelse(CoD_DBSvSham_cohend > 0, "Down", "Up")) %>% count(Celltype, 
      Regulation)
  counts <- long %>% count(genes, name = "N_celltypes")
  write.csv(ranking, output_path("results/Table/TRG_ranking.csv"), row.names = FALSE)
  write.csv(counts, output_path("results/Table/TRG_celltype_counts.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
