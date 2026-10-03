# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  theme_s4_v2 <- function(base_size = 5, base_family = "Arial") {
      theme_classic(base_size = base_size, base_family = base_family) %+replace% theme(plot.title = element_text(size = 5.2, 
          face = "bold", hjust = 0.5, margin = margin(b = 2)), plot.subtitle = element_text(size = 4.5, 
          hjust = 0.5, margin = margin(b = 2)), axis.title = element_text(size = 5, face = "bold", color = "black"), 
          axis.title.x = element_text(size = 5, face = "bold", color = "black", margin = margin(t = 2)), 
          axis.title.y = element_text(size = 5, face = "bold", color = "black", angle = 90, vjust = 1, 
              margin = margin(r = 2)), axis.text = element_text(size = 5, color = "black"), axis.text.x = element_text(size = 4.8, 
              color = "black"), axis.text.y = element_text(size = 4.8, color = "black"), axis.line = element_line(color = "black", 
              linewidth = 0.5), axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.length = unit(0.8, 
              "mm"), legend.title = element_text(size = 4.8, face = "bold"), legend.text = element_text(size = 4.5), 
          legend.key.size = unit(2.5, "mm"), legend.background = element_blank(), strip.background = element_rect(fill = "grey95", 
              color = "black", linewidth = 0.5), strip.text = element_text(size = 5, face = "bold", color = "black", 
              margin = margin(1.5, 1.5, 1.5, 1.5)), panel.grid = element_blank(), plot.margin = margin(2, 
              3, 2, 3, "pt"))
  }
  build_lineage_master <- function(tables, palette) {
      display_order <- intersect(c("ExN", "InN", "OPC", "Oligo", "Astro", "Micro"), names(tables))
      tables <- tables[display_order]
      plots <- lapply(names(tables), function(ct) {
          d <- tables[[ct]]
          if (!nrow(d)) 
              return(NULL)
          d$Direction <- factor(d$Direction, levels = c("Up in DBS", "Down in DBS"))
          d <- d[order(d$Direction, -d$Overlap, d$Pathway), ]
          d$Term <- stringr::str_wrap(d$Clean_Term, 31)
          d$Term <- factor(d$Term, levels = rev(unique(d$Term)))
          ggplot(d, aes(Overlap, Term, fill = LogP)) + geom_col(width = 0.7, color = "black", linewidth = 0.25) + 
              facet_wrap(~Direction, ncol = 1, scales = "free_y") + scale_fill_gradient(low = "#f6f4f2", 
              high = palette[[ct]]) + scale_x_continuous(expand = expansion(mult = c(0, 0.12))) + labs(title = ct, 
              x = "TRG count", y = NULL) + theme_s4_v2(5) + theme(legend.position = "none", axis.text.y = element_text(size = 4.3, 
              lineheight = 0.9), strip.background = element_blank(), strip.text = element_text(size = 5), 
              plot.title = element_text(size = 5.5, face = "bold", hjust = 0.5))
      })
      patchwork::wrap_plots(Filter(Negate(is.null), plots), ncol = 3)
  }
  set_panel_size_custom <- function(p, width_mm, height_mm) {
      g <- ggplotGrob(p)
      panel_idx_w <- g$layout$l[g$layout$name == "panel"]
      panel_idx_h <- g$layout$t[g$layout$name == "panel"]
      g$widths[panel_idx_w] <- unit(width_mm, "mm")
      g$heights[panel_idx_h] <- unit(height_mm, "mm")
      total_w_mm <- as.numeric(convertWidth(sum(g$widths), "mm"))
      total_h_mm <- as.numeric(convertHeight(sum(g$heights), "mm"))
      return(list(grob = g, width_mm = total_w_mm, height_mm = total_h_mm))
  }
  invisible(as.list(environment()))
}
