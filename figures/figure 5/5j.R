# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(igraph)
      library(ggraph)
      library(dplyr)
      library(ggplot2)
      library(patchwork)
  })
  root <- output_root
  data_dir <- file.path(root, "sourcedata/figure 5")
  fig_dir <- file.path(root, "figures/figure 5")
  g <- supplied_inputs[["g_2"]][[as.character(input_path("sourcedata/figure 5/inputs/PPI_TopHubs_50.graphml"))]]
  stopifnot(vcount(g) == 124L)
  membership <- components(g)$membership
  fill_cols <- c(`ATN-Specific TRG` = "#E64B35", `Neuron-Common TRG` = "#4DBBD5", Uninvolved = "grey85")
  line_cols <- c(`ATN-Specific TRG` = "#801a0c", `Neuron-Common TRG` = "#1e687a", Uninvolved = "grey60")
  global_limits <- range(V(g)$internal_degree)
  node_rows <- list()
  edge_rows <- list()
  plots <- list()
  for (i in seq_len(max(membership))) {
      sub <- induced_subgraph(g, which(membership == i))
      nodes <- igraph::as_data_frame(sub, what = "vertices")
      labels <- nodes %>% filter(Category != "Uninvolved") %>% arrange(desc(internal_degree)) %>% slice_head(n = 10) %>% 
          pull(name)
      V(sub)$Label <- ifelse(V(sub)$name %in% labels, paste0(substr(V(sub)$name, 1, 1), tolower(substring(V(sub)$name, 
          2))), "")
      set.seed(42)
      lay <- if (vcount(sub) > 50) 
          create_layout(sub, layout = "graphopt", charge = 0.5)
      else create_layout(sub, layout = "fr")
      node_rows[[i]] <- data.frame(Component = i, Gene = lay$name, Category = lay$Category, Degree_in_full_GRN = lay$internal_degree, 
          Label = lay$Label, x = lay$x, y = lay$y)
      edges <- igraph::as_data_frame(sub, what = "edges")
      edge_rows[[i]] <- data.frame(Component = i, edges)
      plots[[i]] <- ggraph(lay) + geom_edge_link(edge_colour = "grey80", edge_alpha = 0.3, edge_width = 0.25) + 
          geom_node_point(aes(size = internal_degree, fill = Category, color = Category), shape = 21, stroke = 0.5) + 
          geom_node_text(aes(label = Label), repel = TRUE, size = 3.3, fontface = "bold", bg.color = "white", 
              bg.r = 0.12, max.overlaps = Inf, seed = 42) + scale_fill_manual(values = fill_cols, drop = FALSE, 
          limits = names(fill_cols), name = NULL, labels = c("ATN-specific TRGs", "Neuronal-common TRGs", 
              "Other genes")) + scale_color_manual(values = line_cols, guide = "none") + scale_size_continuous(limits = global_limits, 
          range = c(2.5, 11), guide = "none") + theme_void() + theme(plot.margin = margin(6, 8, 6, 8), 
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5), legend.position = "bottom")
  }
  small_components <- wrap_plots(plots[c(8, 9, 3, 5, 7, 6)], ncol = 2)
  right <- (plots[[4]] | plots[[1]])/small_components + plot_layout(heights = c(1, 1.2))
  combined_main <- (plots[[2]] | right) + plot_layout(widths = c(1.1, 1)) & theme(legend.position = "none")
  legend <- grid::grobTree(grid::pointsGrob(x = c(0.14, 0.41, 0.73), y = rep(0.5, 3), pch = 21, size = grid::unit(3, 
      "mm"), gp = grid::gpar(fill = unname(fill_cols), col = unname(line_cols))), grid::textGrob(c("ATN-specific TRGs", 
      "Neuronal-common TRGs", "Other genes"), x = c(0.155, 0.425, 0.745), y = 0.5, just = "left", gp = grid::gpar(fontsize = 11)))
  combined <- combined_main/wrap_elements(full = legend) + plot_layout(heights = c(1, 0.055))
  write.csv(bind_rows(node_rows), file.path(data_dir, "5j.csv"), row.names = FALSE)
  write.csv(bind_rows(edge_rows), file.path(data_dir, "5j_edges.csv"), row.names = FALSE)
  ggsave(file.path(fig_dir, "5j.png"), combined, width = 13, height = 8, dpi = 400, bg = "white")
  cat("5j: ", vcount(g), " nodes in ", max(membership), " components exported.\n")
  invisible(as.list(environment()))
}
