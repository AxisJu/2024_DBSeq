# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
  list2env(supplied_inputs[["helpers_1"]], envir = environment())
  opt <- supplied_inputs[["opt_2"]]
  set.seed(as.integer(opt$seed))
  library(qs)
  library(dplyr)
  library(Seurat)
  library(magrittr)
  obj <- supplied_inputs[["obj_3"]][[as.character(input_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]] %>% 
      subset(Group_L2 %in% c("DBS_I", "Sham_I", "Saline_I")) %>% subset(donor_id %in% c("unassigned", "doublet"), 
      invert = TRUE) %>% subset(celltype_level3 %in% c("ExN_ATN"))
  public_atn_mat_raw <- supplied_inputs[["public_atn_mat_raw_4"]][[as.character(input_path("data/snRNAseq_Mm_ATN_Cembrowski/raw/GSE227627_merged_counts.txt"))]]
  public_atn_mat_anno <- public_atn_mat_raw[, c(1:2)]
  unique_gene_idx <- !duplicated(public_atn_mat_anno[, 2])
  idx <- which(unique_gene_idx)
  public_atn_mat <- public_atn_mat_raw[idx, c(3:ncol(public_atn_mat_raw))]
  rownames(public_atn_mat) <- public_atn_mat_anno[idx, 2]
  public_atn_meta <- supplied_inputs[["public_atn_meta_5"]][[as.character(input_path("data/snRNAseq_Mm_ATN_Cembrowski/raw/GSE227627_metadata.txt"))]]
  seuratobj <- CreateSeuratObject(counts = public_atn_mat, meta.data = public_atn_meta, assay = "RNA", 
      min.cells = 3, min.features = 200, project = "Public_ATN", names.field = 1, names.delim = "_")
  seuratobj[["percent.mt"]] <- PercentageFeatureSet(seuratobj, pattern = "^mt-")
  seuratobj[["percent.rb"]] <- PercentageFeatureSet(seuratobj, pattern = "^Rp[sl]")
  seuratobj[["is.oligo"]] <- (GetAssayData(seuratobj, assay = "RNA", layer = "counts")["Olig1", ] > 0) | 
      (GetAssayData(seuratobj, assay = "RNA", layer = "counts")["Olig2", ] > 0)
  seuratobj[["is.nn"]] <- (GetAssayData(seuratobj, assay = "RNA", layer = "counts")["Snap25", ] < 10)
  seuratobj[["QC"]] <- ifelse((seuratobj$nFeature_RNA > 1500) & (seuratobj$percent.mt < 5) & (!seuratobj$is.oligo) & 
      (!seuratobj$is.nn) & (seuratobj$animalId == "Mark_Anterior"), "Pass", "Failed")
  table(seuratobj$QC)
  seuratobj <- subset(seuratobj, subset = QC == "Pass")
  seuratobj <- seuratobj %>% SCTransform(.) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
      ScaleData(.) %>% RunPCA(.)
  seuratobj <- FindNeighbors(seuratobj, dims = 1:10, reduction = "pca")
  seuratobj <- FindClusters(seuratobj, resolution = 1)
  seuratobj <- RunUMAP(seuratobj, dims = 1:10, reduction = "pca")
  celltype <- seuratobj$SCT_snn_res.1 %>% case_match(c("7") ~ "AD", c("0", "5") ~ "AV", c("1") ~ "AM", 
      .default = "OTH")
  table(celltype)
  seuratobj[["celltype"]] <- factor(celltype, levels = c("AD", "AV", "AM", "OTH"))
  seuratobj <- FindNeighbors(seuratobj, dims = 1:10, reduction = "pca")
  seuratobj <- RunUMAP(seuratobj, dims = 1:10, reduction = "pca")
  qsave(seuratobj, output_path("data/snRNAseq_Mm_ATN_Cembrowski/processed/running_ATN.qs"))
  obj <- supplied_inputs[["obj_6"]][[as.character(input_path("data/snRNAseq_mouse/processed/matrix/running_all_250704.qs"))]] %>% 
      subset(Group_L2 %in% c("DBS_I", "Sham_I", "Saline_I")) %>% subset(donor_id %in% c("unassigned", "doublet"), 
      invert = TRUE) %>% subset(celltype_level3 %in% c("ExN_ATN"))
  obj <- obj %>% SCTransform(.) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
      ScaleData(.) %>% RunPCA(.)
  anchors <- FindTransferAnchors(reference = seuratobj, query = obj, dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = seuratobj$celltype, dims = 1:30)
  obj <- AddMetaData(obj, metadata = predictions)
  table(obj$predicted.id)
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "scanvi", graph.name = "subclust")
  obj <- FindClusters(obj, resolution = 0.1, graph.name = "subclust")
  obj <- RunUMAP(obj, dims = 1:30, min.dist = 1, spread = 1)
  celltype <- obj$subclust_res.0.1 %>% case_match(c("3") ~ "AD", c("2") ~ "AV", c("0", "1") ~ "AM", .default = "OTH")
  table(celltype)
  obj[["celltype"]] <- factor(celltype, levels = c("AD", "AV", "AM", "OTH"))
  obj <- subset(obj, subset = celltype == "OTH", invert = TRUE)
  obj <- RunUMAP(obj, dims = 1:30, min.dist = 1, spread = 1)
  qsave(obj, output_path("data/snRNAseq_mouse/processed/matrix/running_ATN_250704.qs"))
  obj <- supplied_inputs[["obj_7"]][[as.character(input_path("data/snRNAseq_mouse/processed/matrix/running_ATN_250704.qs"))]]
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
  write.csv(as.data.frame(as.matrix(expr_mat)), file = output_path("data/snRNAseq_mouse/processed/matrix/running_ATN_250704_mat.csv"))
  meta_df <- obj@meta.data
  write.csv(meta_df, file = output_path("data/snRNAseq_mouse/processed/metadata/running_ATN_250704_metadata.csv"))
  var_df <- rownames(expr_mat)
  write.csv(var_df, file = output_path("data/snRNAseq_mouse/processed/metadata/running_ATN_250704_var.csv"))
  invisible(as.list(environment()))
}
