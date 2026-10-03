# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  entry_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  library(dplyr)
  library(ppcor)
  library(qs)
  library(dplyr)
  library(tibble)
  dir_derivatives <- data_path("data/derivatives/dbseq")
  dir_allen_bin <- required_root("DBSEQ_ALLEN_ROOT")
  load_and_clean_mat <- function(filepath, clean_ipsi = TRUE) {
      mat <- supplied_inputs[["mat_2"]][[as.character(filepath)]]
      first_col <- colnames(mat)[1]
      if (clean_ipsi) {
          mat[[first_col]] <- make.names(gsub("_ipsi", "", mat[[first_col]]), unique = TRUE)
          colnames(mat) <- gsub("_ipsi", "", colnames(mat))
      }
      else {
          mat[[first_col]] <- make.names(mat[[first_col]], unique = TRUE)
      }
      mat <- mat %>% tibble::column_to_rownames(var = first_col)
      return(mat)
  }
  fc_region_mapping <- c(CP = "STRd", LSc = "LSX", TRS = "LSX", CM = "ILM", PF = "ILM", RE = "MTN", GPe = "PALd", 
      MH = "EPI", LH = "LZ", RT = "RT", RSPd = "RSPd", DG = "DG")
  extract_seed_vector <- function(mat, target_regions, is_fc = FALSE) {
      vec <- rep(NA, length(target_regions))
      names(vec) <- target_regions
      if (is_fc) {
          avail_seed <- intersect("ATN", rownames(mat))
          if (length(avail_seed) == 0) 
              return(vec)
          for (i in seq_along(target_regions)) {
              orig_name <- target_regions[i]
              mapped_name <- ifelse(orig_name %in% names(fc_region_mapping), fc_region_mapping[[orig_name]], 
                  orig_name)
              if (mapped_name %in% colnames(mat)) {
                  vec[i] <- mean(mat[avail_seed, mapped_name], na.rm = TRUE)
              }
          }
      }
      else {
          avail_seed <- intersect(c("AV", "AD", "AMd", "AMv"), rownames(mat))
          avail_target <- intersect(target_regions, colnames(mat))
          if (length(avail_seed) > 0 && length(avail_target) > 0) {
              raw_means <- colMeans(mat[avail_seed, avail_target, drop = FALSE], na.rm = TRUE)
              vec[names(raw_means)] <- raw_means
          }
      }
      return(vec)
  }
  cat("Loading and processing matrices...\n")
  sc_w_path <- file.path(dir_allen_bin, "41586_2014_BFnature13186_MOESM71_ESM_W_ipsi.csv")
  sc_p_path <- file.path(dir_allen_bin, "41586_2014_BFnature13186_MOESM71_ESM_PValue_ipsi.csv")
  mat_sc_w <- supplied_inputs[["mat_sc_w_3"]][[as.character(sc_w_path)]]
  mat_sc_p <- supplied_inputs[["mat_sc_p_4"]][[as.character(sc_p_path)]]
  mat_sc_w[[1]] <- make.names(gsub("_ipsi", "", mat_sc_w[[1]]), unique = TRUE)
  mat_sc_p[[1]] <- make.names(gsub("_ipsi", "", mat_sc_p[[1]]), unique = TRUE)
  mat_sc_w <- mat_sc_w %>% tibble::column_to_rownames(var = colnames(mat_sc_w)[1])
  mat_sc_p <- mat_sc_p %>% tibble::column_to_rownames(var = colnames(mat_sc_p)[1])
  colnames(mat_sc_w) <- gsub("_ipsi", "", colnames(mat_sc_w))
  colnames(mat_sc_p) <- gsub("_ipsi", "", colnames(mat_sc_p))
  mat_sc_raw <- mat_sc_w
  mat_sc_raw[mat_sc_p > 0.05] <- 0
  m_sc_raw <- as.matrix(mat_sc_raw)
  mat_sc2 <- m_sc_raw %*% m_sc_raw
  mat_sc <- mat_sc_raw
  mask_nonzero <- mat_sc > 0
  if (any(mask_nonzero)) {
      log_vals <- log10(mat_sc[mask_nonzero])
      shift_val <- abs(min(log_vals)) + 1
      mat_sc[mask_nonzero] <- log_vals + shift_val
  }
  mat_ed <- supplied_inputs[["mat_ed_5"]][[as.character(file.path(dir_allen_bin, "41586_2014_BFnature13186_MOESM72_ESM.csv"))]]
  mat_fc <- supplied_inputs[["mat_fc_6"]][[as.character(file.path(dir_derivatives, "mouse_BOLD_fc.csv"))]]
  mat_ne <- supplied_inputs[["mat_ne_7"]][[as.character(file.path(dir_derivatives, "allen_mousebrainconnectome_NE.csv"))]]
  mat_cmy <- supplied_inputs[["mat_cmy_8"]][[as.character(file.path(dir_derivatives, "allen_mousebrainconnectome_CMY.csv"))]]
  mat_de <- supplied_inputs[["mat_de_9"]][[as.character(file.path(dir_derivatives, "allen_mousebrainconnectome_DE.csv"))]]
  mat_spe <- supplied_inputs[["mat_spe_10"]][[as.character(file.path(dir_derivatives, "allen_mousebrainconnectome_SPE.csv"))]]
  target_12_regions <- c("CP", "GPe", "LSc", "TRS", "CM", "PF", "RE", "RT", "MH", "LH", "RSPd", "DG")
  df_imaging_features_12 <- data.frame(Region = target_12_regions) %>% dplyr::mutate(FC = extract_seed_vector(mat_fc, 
      Region, is_fc = TRUE), SC = extract_seed_vector(mat_sc, Region, is_fc = FALSE), SC2 = extract_seed_vector(mat_sc2, 
      Region, is_fc = FALSE), ED = extract_seed_vector(mat_ed, Region, is_fc = FALSE), NE = extract_seed_vector(mat_ne, 
      Region, is_fc = FALSE), CMY = extract_seed_vector(mat_cmy, Region, is_fc = FALSE), DE = extract_seed_vector(mat_de, 
      Region, is_fc = FALSE), SPE = extract_seed_vector(mat_spe, Region, is_fc = FALSE))
  cat("Extraction complete. Displaying results:\n")
  print(df_imaging_features_12)
  list2env(supplied_inputs[["helpers_11"]], envir = environment())
  library(Seurat)
  mouse <- supplied_inputs[["mouse_12"]]
  object <- supplied_inputs[["object_13"]][[as.character(data_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]]
  mapping <- object@meta.data %>% dplyr::filter(Group_L2 %in% c("DBS_I", "Sham_I", "Saline_I"), !donor_id %in% 
      c("unassigned", "doublet"), brainregion_projection != "/") %>% dplyr::distinct(celltype_level3, brainregion_projection)
  rm(object)
  atn_genes <- mouse$genes[mouse$Celltype == "ExN_ATN" & mouse$is_TRG]
  membership <- mouse %>% dplyr::filter(is_TRG, grepl("^(ExN|InN)", Celltype), genes %in% atn_genes) %>% 
      dplyr::inner_join(mapping, by = c(Celltype = "celltype_level3")) %>% dplyr::transmute(Gene = genes, 
      Region = brainregion_projection) %>% dplyr::distinct()
  go <- supplied_inputs[["go_14"]][[as.character(input_path("sourcedata/figure 5/inputs/go_bp_2023.qs"))]]
  splicing_genes <- go[["RNA Splicing (GO:0008380)"]]
  if (!length(splicing_genes)) stop("RNA splicing gene set is missing.")
  if (!all(target_12_regions %in% membership$Region)) stop("Regional TRG membership is incomplete.")
  df_splicing <- data.frame(Region = target_12_regions, Enriched_Count = vapply(target_12_regions, function(region) length(intersect(membership$Gene[membership$Region == 
      region], splicing_genes)), integer(1)))
  write.csv(membership %>% dplyr::filter(Gene %in% splicing_genes), output_path("sourcedata/figure 4/4jk_splicing_gene_membership.csv"), 
      row.names = FALSE)
  df_jk <- df_imaging_features_12 %>% inner_join(df_splicing, by = "Region")
  print(df_jk)
  region_audit <- data.frame(Region = c("ATN", target_12_regions)) %>% mutate(Has_RNA_splicing_record = Region %in% 
      df_splicing$Region, Is_Imaging_Target = Region %in% target_12_regions, Included = Region %in% df_jk$Region, 
      Reason = case_when(Region == "ATN" ~ "connectivity_seed", !Has_RNA_splicing_record ~ "no_RNA_splicing_record", 
          TRUE ~ "included"))
  data_dir <- output_path("sourcedata/figure 4")
  write.csv(df_jk, output_path("cache/4jk_input.csv"), row.names = FALSE)
  write.csv(region_audit, output_path("statistics/figure 4/4jk_region_audit.csv"), row.names = FALSE)
  invisible(as.list(environment()))
}
