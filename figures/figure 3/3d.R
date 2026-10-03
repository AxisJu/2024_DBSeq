# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  suppressPackageStartupMessages({
      library(qs)
      library(dplyr)
      library(stringr)
      library(NeuronChat)
      library(CellChat)
      library(circlize)
      library(colorRamp2)
      library(RColorBrewer)
  })
  base_dir <- data_root
  output_dir_fig <- output_path("figures/figure 3")
  output_dir_data <- output_path("sourcedata/figure 3")
  dir.create(output_dir_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  dbs_qs_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/NeuronChat/DBS_I.qs")
  sham_qs_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/NeuronChat/Sham_I.qs")
  saline_qs_path <- file.path(base_dir, "data/snRNAseq_mouse/processed/intermediate/NeuronChat/Saline_I.qs")
  colors_celltype_level3 <- c(Astro = "#a58946", EPCs = "#594a26", OPC = "#5953ff", Oligo = "#201e5a", 
      VCs = "#858881", CHPCs = "#cfd4c9", Micro = "#a87c5a", ExN_RSP_L23IT = "#ff1a71", ExN_RSP_L45IT = "#ed1986", 
      ExN_CLA = "#ba1369", ExN_CA3 = "#860d4c", ExN_CA1 = "#53082f", ExN_RSP_L6CT = "#61e2a4", ExN_HPF_CajalRetzius = "#d00000", 
      ExN_DG = "#16f2f2", InN_Immature = "#1b4332", InN_RSP_MGE = "#f954ee", InN_STRns_MGE = "#70e000", 
      InN_GPe_MGE = "#56ad00", InN_GPi_MGE = "#2f6000", InN_Nonspecific = "#f0a0ff", InN_STRns_LGE = "#d899ff", 
      InN_GPe_LGE = "#b199ff", InN_CP_D1 = "#8d7acc", InN_CP_D2 = "#695b99", InN_OT = "#4e4372", InN_LSc = "#3283fe", 
      InN_HYa = "#450099", InN_RT = "#ff6600", ExN_HYa = "#aa0dfe", ExN_TRS = "#7609b1", ExN_BST = "#420564", 
      ExN_MH = "#faa307", ExN_LH = "#c68105", ExN_THns = "#0d47a1", ExN_ATN = "#1460ff", ExN_RE = "#1340ff", 
      ExN_PF = "#0a4093", ExN_CM = "#08306d", InN_SCsg = "#9ef01a")
  x_DBS_I <- supplied_inputs[["x_DBS_I_2"]][[as.character(dbs_qs_path)]]
  x_Sham_I <- supplied_inputs[["x_Sham_I_3"]][[as.character(sham_qs_path)]]
  x_Saline_I <- supplied_inputs[["x_Saline_I_4"]][[as.character(saline_qs_path)]]
  all_samples <- list(x_Saline_I@net, x_Sham_I@net, x_DBS_I@net)
  target_total <- 500
  all_samples_norm <- lapply(all_samples, function(sample) {
      current_total <- sum(sapply(sample, sum))
      scale_factor <- target_total/current_total
      lapply(sample, function(mat) mat * scale_factor)
  })
  net_saline <- all_samples_norm[[1]][["Nrxn1_Nlgn1"]]
  net_sham <- all_samples_norm[[2]][["Nrxn1_Nlgn1"]]
  net_dbs <- all_samples_norm[[3]][["Nrxn1_Nlgn1"]]
  net_dbssham <- net_dbs - net_sham
  net_salinesham <- net_sham - net_saline
  celltypes <- rownames(net_dbssham)
  group <- ifelse(str_detect(celltypes, "ExN"), "Excitatory", ifelse(str_detect(celltypes, "InN"), "Inhibitory", 
      "Non-neuronal"))
  names(group) <- celltypes
  extract_df <- function(mat, dir_name, comp_name) {
      if (dir_name == "Input to ATN") {
          data.frame(Direction = dir_name, Comparison = comp_name, Source_CellType = rownames(mat), Target_CellType = "ExN_ATN", 
              Differential_Probability = mat[, "ExN_ATN"], stringsAsFactors = FALSE)
      }
      else {
          data.frame(Direction = dir_name, Comparison = comp_name, Source_CellType = "ExN_ATN", Target_CellType = colnames(mat), 
              Differential_Probability = mat["ExN_ATN", ], stringsAsFactors = FALSE)
      }
  }
  sourcedata_d <- rbind(extract_df(net_dbssham, "Input to ATN", "DBS vs Sham"), extract_df(net_salinesham, 
      "Input to ATN", "Sham vs Saline"), extract_df(net_dbssham, "Output from ATN", "DBS vs Sham"), extract_df(net_salinesham, 
      "Output from ATN", "Sham vs Saline")) %>% mutate(Is_Self = Source_CellType == Target_CellType, Panel_Row = ifelse(Direction == 
      "Input to ATN", 1L, 2L), Panel_Column = ifelse(Comparison == "DBS vs Sham", 1L, 2L)) %>% arrange(Direction, 
      Comparison, desc(abs(Differential_Probability)))
  write.csv(sourcedata_d, file.path(output_dir_data, "3d.csv"), row.names = FALSE)
  cat(sprintf("Exported Source Data: %s (Rows: %d)\n", file.path(output_dir_data, "3d.csv"), nrow(sourcedata_d)))
  axis_chord <- function(net, edge.transparency = FALSE, color.use = NULL, group = NULL, cell.order = NULL, 
      sources.use = NULL, targets.use = NULL, lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03), 
      remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
      reduce = -1, transparency = 0, link.border = NA, title.cex = 1.2, title.name = NULL, show.legend = FALSE, 
      legend.pos.x = 20, legend.pos.y = 20, ...) {
      if (inherits(x = net, what = c("matrix", "Matrix"))) {
          cell.levels <- union(rownames(net), colnames(net))
          net <- reshape2::melt(net, value.name = "prob")
          colnames(net)[1:2] <- c("source", "target")
      }
      else if (is.data.frame(net)) {
          cell.levels <- as.character(union(net$source, net$target))
      }
      if (!is.null(cell.order)) 
          cell.levels <- cell.order
      net$source <- as.character(net$source)
      net$target <- as.character(net$target)
      if (!is.null(sources.use)) 
          net <- subset(net, source %in% sources.use)
      if (!is.null(targets.use)) 
          net <- subset(net, target %in% targets.use)
      net <- subset(net, prob != 0)
      df <- net
      cells.use <- union(df$source, df$target)
      order.sector <- cell.levels[cell.levels %in% cells.use]
      if (is.null(color.use)) {
          color.use <- scPalette(length(cell.levels))
          names(color.use) <- cell.levels
      }
      grid.col <- color.use[order.sector]
      names(grid.col) <- order.sector
      if (!is.null(group)) 
          group <- group[names(group) %in% order.sector]
      range_value <- 0.2
      col_fun <- colorRamp2(breaks = c(-range_value, 0, range_value), colors = rev(brewer.pal(11, "RdBu")[c(2, 
          6, 10)]))
      df$col <- col_fun(df$prob)
      edge.color <- df$col
      link.arr.type <- ifelse(directional %in% c(0, 2), "triangle", "big.arrow")
      circos.clear()
      chordDiagram(df, order = order.sector, col = edge.color, grid.col = grid.col, transparency = transparency, 
          link.border = link.border, directional = directional, direction.type = c("diffHeight", "arrows"), 
          link.arr.type = link.arr.type, annotationTrack = "grid", annotationTrackHeight = annotationTrackHeight, 
          preAllocateTracks = list(track.height = mm_h(6)), small.gap = small.gap, big.gap = big.gap, scale = TRUE, 
          group = group, link.target.prop = link.target.prop, diffHeight = mm_h(4), target.prop.height = mm_h(2), 
          reduce = reduce, ...)
      circos.track(track.index = 1, panel.fun = function(x, y) {
      }, bg.border = NA)
      if (!is.null(group)) {
          drawn_sectors <- get.all.sector.index()
          unique_groups <- unique(group[names(group) %in% drawn_sectors])
          outer_colors <- c(Excitatory = "#ba1369", Inhibitory = "#56ad00", `Non-neuronal` = "#201e5a")
          for (g in unique_groups) {
              sectors_in_g <- names(group)[group == g & names(group) %in% drawn_sectors]
              if (length(sectors_in_g) > 0) {
                  base_col <- ifelse(g %in% names(outer_colors), outer_colors[g], "#999999")
                  highlight.sector(sector.index = sectors_in_g, track.index = 1, col = adjustcolor(base_col, 
                    alpha.f = 0.25), border = base_col, lwd = 1.5, text = g, cex = lab.cex * 1.5, text.vjust = -0.5, 
                    niceFacing = TRUE)
              }
          }
      }
      circos.clear()
  }
  draw_panel <- function() {
      layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.18))
      par(mar = c(0.5, 3.3, 3.5, 0.5))
      axis_chord(net_dbssham, group = group, targets.use = "ExN_ATN", color.use = colors_celltype_level3)
      title(main = "DBS vs Sham", line = 1)
      mtext("Nrxn1-Nlgn1: Input to ATN", side = 2, line = 1, font = 2)
      axis_chord(net_salinesham, group = group, targets.use = "ExN_ATN", color.use = colors_celltype_level3)
      title(main = "Sham vs Saline", line = 1)
      axis_chord(net_dbssham, group = group, sources.use = "ExN_ATN", color.use = colors_celltype_level3)
      mtext("Nrxn1-Nlgn1: Output from ATN", side = 2, line = 1, font = 2)
      axis_chord(net_salinesham, group = group, sources.use = "ExN_ATN", color.use = colors_celltype_level3)
      par(mar = c(2.5, 6, 1, 6))
      values <- seq(-0.2, 0.2, length.out = 201)
      cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")[c(2, 6, 10)]))(200)
      plot.new()
      plot.window(xlim = c(-0.2, 0.2), ylim = c(0, 1), xaxs = "i")
      rect(values[-201], 0.4, values[-1], 0.8, col = cols, border = NA)
      axis(1, at = c(-0.2, 0, 0.2), pos = 0.4)
      mtext("Differential CCC Probability", side = 1, line = 1.1)
  }
  pdf(NULL, width = 12, height = 12, useDingbats = FALSE)
  draw_panel()
  dev.off()
  png(file.path(output_dir_fig, "3d.png"), width = 3600, height = 3600, res = 300)
  draw_panel()
  dev.off()
  cat("Panel 3d exported with explicit direction/comparison labels.\n")
  invisible(as.list(environment()))
}
