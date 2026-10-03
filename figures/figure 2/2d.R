# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(RRHO2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  proj_dir <- Sys.getenv("DBSEQ_DATA_ROOT", unset = "")
  if (!nzchar(proj_dir)) stop("Set DBSEQ_DATA_ROOT to the original DBSeq data root.")
  package_root <- output_root
  if (!nzchar(package_root)) {
      script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      if (!length(script_arg)) 
          stop("Set DBSEQ_PACKAGE_ROOT when sourcing this script.")
      script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
      package_root <- dirname(dirname(dirname(script_path)))
  }
  out_fig_dir <- file.path(package_root, "figures", "figure 2")
  out_data_dir <- file.path(package_root, "sourcedata", "figure 2")
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  load_de_data <- function(lvl, mode, comp, ct) {
      ct_safe <- gsub("/", "-", ct)
      file_path <- results_path("Table/EdgeR_v2", lvl, paste0("EdgeR_", comp, "_", ct_safe, ".csv"))
      df <- supplied_inputs[["df_2"]][[as.character(file_path)]]
      gene_col <- if ("X" %in% colnames(df)) 
          "X"
      else colnames(df)[1]
      data.frame(X = df[[gene_col]], PValue = df$PValue, logFC = df$logFC)
  }
  compute_rrho <- function(de1, de2, name1, name2, step_size = 50, flip_de2 = FALSE) {
      de1 <- de1 %>% filter(!is.na(PValue), !is.na(logFC)) %>% distinct(X, .keep_all = TRUE)
      de2 <- de2 %>% filter(!is.na(PValue), !is.na(logFC)) %>% distinct(X, .keep_all = TRUE)
      common <- intersect(de1$X, de2$X)
      if (length(common) < step_size * 2) 
          return(NULL)
      df1 <- de1 %>% filter(X %in% common) %>% mutate(metric = -log10(PValue + 1e-300) * sign(logFC)) %>% 
          arrange(desc(metric)) %>% select(X, metric) %>% as.data.frame()
      if (flip_de2) {
          df2 <- de2 %>% filter(X %in% common) %>% mutate(metric = -log10(PValue + 1e-300) * sign(logFC) * 
              -1) %>% arrange(desc(metric)) %>% select(X, metric) %>% as.data.frame()
      }
      else {
          df2 <- de2 %>% filter(X %in% common) %>% mutate(metric = -log10(PValue + 1e-300) * sign(logFC)) %>% 
              arrange(desc(metric)) %>% select(X, metric) %>% as.data.frame()
      }
      rrho_hyper <- RRHO2_initialize(df1, df2, labels = c(name1, name2), log10.ind = TRUE, boundary = 0, 
          method = "hyper", stepsize = step_size)
      n_up1 <- sum(df1$metric > 0, na.rm = TRUE)
      n_up2 <- sum(df2$metric > 0, na.rm = TRUE)
      cross1 <- ceiling(n_up1/step_size) + 0.5
      cross2 <- ceiling(n_up2/step_size) + 0.5
      list(mat_pval = rrho_hyper$hypermat, cross1 = cross1, cross2 = cross2, name1 = name1, name2 = name2)
  }
  cat("Computing RRHO2 matrices for InN and ExN...\n")
  res_inn_ss <- compute_rrho(supplied_inputs[["res_inn_ss_3"]][[as.character("Level1")]], supplied_inputs[["res_inn_ss_4"]][[as.character("Level1")]], 
      "DBS vs Sham", "Sham vs Saline", 50, TRUE)
  res_inn_ic <- compute_rrho(supplied_inputs[["res_inn_ic_5"]][[as.character("Level1")]], supplied_inputs[["res_inn_ic_6"]][[as.character("Level1")]], 
      "DBS vs Sham", "DBS Ipsi. vs Contra.", 50, FALSE)
  res_exn_ss <- compute_rrho(supplied_inputs[["res_exn_ss_7"]][[as.character("Level1")]], supplied_inputs[["res_exn_ss_8"]][[as.character("Level1")]], 
      "DBS vs Sham", "Sham vs Saline", 50, TRUE)
  res_exn_ic <- compute_rrho(supplied_inputs[["res_exn_ic_9"]][[as.character("Level1")]], supplied_inputs[["res_exn_ic_10"]][[as.character("Level1")]], 
      "DBS vs Sham", "DBS Ipsi. vs Contra.", 50, FALSE)
  build_source_df <- function(rrho_res, celltype, comparison) {
      df <- reshape2::melt(rrho_res$mat_pval)
      colnames(df) <- c("Bin_X", "Bin_Y", "MinusLog10_PValue")
      df$Celltype <- celltype
      df$Comparison <- comparison
      df[, c("Celltype", "Comparison", "Bin_X", "Bin_Y", "MinusLog10_PValue")]
  }
  source_d <- bind_rows(build_source_df(res_inn_ss, "InN", "DBSvSham_vs_ShamvSaline"), build_source_df(res_inn_ic, 
      "InN", "DBSvSham_vs_DBSIvC"), build_source_df(res_exn_ss, "ExN", "DBSvSham_vs_ShamvSaline"), build_source_df(res_exn_ic, 
      "ExN", "DBSvSham_vs_DBSIvC"))
  write.csv(source_d, file.path(out_data_dir, "2d.csv"), row.names = FALSE)
  cat("Saved source data: sourcedata/figure 2/2d.csv\n")
  plot_single_rrho <- function(rrho_obj, x_label, y_label, show_x = FALSE, show_y = FALSE, show_legend = FALSE) {
      plot_df <- reshape2::melt(rrho_obj$mat_pval)
      p <- ggplot(plot_df, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + scale_fill_gradientn(colors = c("#f4f9f4", 
          "#7fcdbb", "#2c7fb8", "#253494"), name = expression(-log[10] ~ italic(P)), limits = c(0, 600), 
          oob = scales::squish) + geom_vline(xintercept = rrho_obj$cross1, color = "white", linewidth = 1.2) + 
          geom_hline(yintercept = rrho_obj$cross2, color = "white", linewidth = 1.2) + theme_classic() + 
          theme(axis.title.x = if (show_x) 
              element_text(face = "bold", size = 9, color = "#ffd700")
          else element_blank(), axis.title.y = if (show_y) 
              element_text(face = "bold", size = 9)
          else element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), legend.position = if (show_legend) 
                  "right"
              else "none", legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8), 
              plot.margin = margin(2, 2, 2, 2)) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
          0)) + labs(x = x_label, y = y_label)
      return(p)
  }
  p1 <- plot_single_rrho(res_inn_ss, "DBS vs Sham", "Sham vs Saline", show_x = FALSE, show_y = TRUE) + 
      theme(axis.title.y = element_text(face = "bold", size = 9, color = "#18a799"))
  p2 <- plot_single_rrho(res_inn_ic, "DBS vs Sham", "DBS Ipsi. vs Contra.", show_x = FALSE, show_y = TRUE) + 
      theme(axis.title.y = element_text(face = "bold", size = 9, color = "#ffd700"))
  p3 <- plot_single_rrho(res_exn_ss, "DBS vs Sham", "Sham vs Saline", show_x = TRUE, show_y = TRUE) + theme(axis.title.y = element_text(face = "bold", 
      size = 9, color = "#18a799"))
  p4 <- plot_single_rrho(res_exn_ic, "DBS vs Sham", "DBS Ipsi. vs Contra.", show_x = TRUE, show_y = TRUE, 
      show_legend = TRUE) + theme(axis.title.y = element_text(face = "bold", size = 9, color = "#ffd700"))
  final_rrho <- (p1 | p2)/(p3 | p4) + plot_layout(guides = "collect")
  pdf_path <- file.path(out_fig_dir, "2d.pdf")
  png_path <- file.path(out_fig_dir, "2d.png")
  ggsave(png_path, final_rrho, width = 6.2, height = 5.5, dpi = 300)
  cat("Saved figures: figures/figure 2/2d.pdf and 2d.png\n")
  invisible(as.list(environment()))
}
