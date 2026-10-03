# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  get_publication_universe <- function() {
      inputs <- input_path("sourcedata/figure 5/inputs")
      go <- supplied_inputs[["go_2"]][[as.character(file.path(inputs, "go_bp_2023.qs"))]]
      enrichment <- supplied_inputs[["enrichment_3"]][[as.character(file.path(inputs, "enrichr_atn_all.qs"))]]
      terms <- enrichment$Term[!is.na(enrichment$Adjusted.P.value) & enrichment$Adjusted.P.value < 0.05]
      universe <- unique(toupper(unlist(go[intersect(terms, names(go))], use.names = FALSE)))
      stopifnot(length(universe) == 2664L)
      universe
  }
  invisible(as.list(environment()))
}
