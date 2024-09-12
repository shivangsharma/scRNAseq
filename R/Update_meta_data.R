obj <- readRDS("/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/seurat_processed_object/SeuratObject_sqrtNormalized_SCTransform_UMAP_annotated.RData")

obj@meta.data$orig.ident <- NULL

obj@meta.data$Run_name <- "Run5324" %>% as.factor()

obj@meta.data$Slide_name <- str_replace(string = obj@meta.data[["Slide_name"]],
                                        pattern = "_",
                                        replacement = "-") %>% factor()
obj@meta.data$Slide_name <- paste0("Slide", obj@meta.data[["Slide_name"]]) %>% as.factor()
obj@meta.data$Slide_name %>% unique()

fovs <- str_split(string = obj@meta.data[["cell_ID"]], pattern = "_") %>%
  lapply(X = ., FUN = "[", 3) %>%
  unlist() %>%
  as.numeric()
width <- fovs %>%
  log10() %>%
  floor() %>%
  max(.) + 1
fovs <- sprintf(paste0("%0", width, "d"), fovs)
obj@meta.data$fov <- paste0("FOV", fovs) %>% as.factor()
obj@meta.data$FOV <- obj@meta.data$fov
obj@meta.data$fov <- NULL

obj@meta.data$Patient_ID <- NA
obj@meta.data[(obj@meta.data$Slide_name == "23-8") & (obj@meta.data$fov %in% c(1:11)), "Patient_ID"] <- "Patient-8"
obj@meta.data[(obj@meta.data$Slide_name == "23-8") & (obj@meta.data$fov %in% c(12:21)), "Patient_ID"] <- "Patient-23"
obj@meta.data[(obj@meta.data$Slide_name == "30-18") & (obj@meta.data$fov %in% c(1:8)), "Patient_ID"] <- "Patient-30"
obj@meta.data[(obj@meta.data$Slide_name == "30-18") & (obj@meta.data$fov %in% c(9:17)), "Patient_ID"] <- "Patient-18"

obj@meta.data$Tissue_name <- obj@meta.data$Category$Patient_ID

obj@meta.data$Cell_number <- NA
obj@meta.data$Cell_number <- str_split(string = obj@meta.data$cell_ID, pattern = "_") %>%
  lapply(X = ., FUN = "[", 4) %>%
  unlist()
obj@meta.data$Cell_number <- paste0("Cell", obj@meta.data[["Cell_number"]])

obj@meta.data$Organ_name <- factor("Prostate")

obj@meta.data$Cell_ID <- paste(obj@meta.data[["Run_name"]], obj@meta.data[["Slide_name"]], obj@meta.data[["Patient_ID"]], obj@meta.data[["Organ_name"]], obj@meta.data[["Tissue_name"]], obj@meta.data[["FOV"]], obj@meta.data[["Cell_number"]], sep = "_")
colnames(obj) <- obj@meta.data[["Cell_ID"]]

obj@meta.data$tissue <- NULL

obj_alt <- RenameCells(obj, new.names = obj@meta.data$Cell_ID)
obj <- obj_alt

df_annotation <- obj@meta.data[, c("Luminal_Score", "Basal_Score", "Endothelial_Score", "BCell_Score", "TCell_Score", "Myeloid_Score", "Mast_Score", "SMC_Score", "Fibroblast_Score", "CellAnnCoarse", "seurat_clusters")]
obj@meta.data$Annotation <- df_annotation
obj@meta.data[, c("Luminal_Score", "Basal_Score", "Endothelial_Score", "BCell_Score", "TCell_Score", "Myeloid_Score", "Mast_Score", "SMC_Score", "Fibroblast_Score", "CellAnnCoarse", "seurat_clusters")] <- NULL
df_other <- obj@meta.data[, c("ISH.concentration", "Dash", "nb_clus", "leiden_clus")]
obj@meta.data[, c("ISH.concentration", "Dash", "nb_clus", "leiden_clus")] <- NULL
obj@meta.data$Other <- df_other
df_image <- obj@meta.data[, c("Area", "AspectRatio", "Width", "Height", "Mean.CD298", "Max.CD298", "Mean.PanCK", "Max.PanCK", "Mean.CD45", "Max.CD45", "Mean.CD3", "Max.CD3", "Mean.DAPI", "Max.DAPI", "IFcolor")]
obj@meta.data[, c("Area", "AspectRatio", "Width", "Height", "Mean.CD298", "Max.CD298", "Mean.PanCK", "Max.PanCK", "Mean.CD45", "Max.CD45", "Mean.CD3", "Max.CD3", "Mean.DAPI", "Max.DAPI", "IFcolor")] <- NULL
obj@meta.data$Image <- df_image
df_categorical <- obj@meta.data[, c("Run_name", "Slide_name", "Patient_ID", "Organ_name", "Tissue_name", "FOV", "Cell_number")]
obj@meta.data[, c("Run_name", "Slide_name", "Patient_ID", "Organ_name", "Tissue_name", "FOV", "Cell_number")] <- NULL
obj@meta.data$Category <- df_categorical
obj@meta.data$cell_ID <- NULL
obj@meta.data$Cell_ID <- NULL

obj@meta.data$Centroid <- obj@meta.data$locations
obj@meta.data$locations <- NULL

obj@meta.data$Category[["Morphology_name"]] <- NA
meta_data_23_8 <- read.csv(file = "/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/raw_data/Run5324_23_8/Run5324_23_8_annotation_file.csv", sep = ",", header = TRUE)
morpho_dict <- meta_data_23_8$Annotation.group.1
names(morpho_dict) <- meta_data_23_8$FOV..
obj@meta.data$Category$Morphology_name[obj@meta.data$Category$Slide_name %in% "Slide23-8"] <- morpho_dict[obj@meta.data$Category$FOV[obj@meta.data$Category$Slide_name %in% "Slide23-8"]]

meta_data_30_18 <- read.csv(file = "/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/raw_data/Run5324_30_18/Run5324_30_18_annotation_file.csv", sep = ",", header = TRUE)
morpho_dict <- meta_data_30_18$Annotation.group.1
names(morpho_dict) <- meta_data_30_18$FOV..
obj@meta.data$Category$Morphology_name[obj@meta.data$Category$Slide_name %in% "Slide30-18"] <- morpho_dict[obj@meta.data$Category$FOV[obj@meta.data$Category$Slide_name %in% "Slide30-18"]]
obj@meta.data$Category$Morphology_name <- str_replace(obj@meta.data$Category$Morphology_name, pattern = " ", replacement = "_")
obj@meta.data$Category$Morphology_name <- factor(obj@meta.data$Category$Morphology_name, levels = c("Benign_Epithelium", "Tumor", "PIN", "Stroma", "Lymphoid_Aggregates"))
obj@meta.data$Category$Tissue_name <- obj@meta.data$Tissue_name
obj@meta.data$Tissue_name <- NULL
obj@meta.data$CategoryTissue_name <- NULL
obj@meta.data$Category$Patient_ID <- as.factor(obj@meta.data$Category$Patient_ID)
obj@meta.data$Category$Tissue_name <- as.factor(obj@meta.data$Category$Tissue_name)

cell_ID <- paste(obj@meta.data$Category[["Run_name"]], obj@meta.data$Category[["Slide_name"]], obj@meta.data$Category[["Patient_ID"]], obj@meta.data$Category[["Organ_name"]], obj@meta.data$Category[["Tissue_name"]], obj@meta.data$Category[["FOV"]], obj@meta.data$Category[["Morphology_name"]], obj@meta.data$Category[["Cell_number"]], sep = "_")

obj <- RenameCells(obj, new.names = cell_ID)

obj_temp <- ReadNanostring(data.dir = "/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/raw_data/Run5324_23_8",
                           mtx.file = "Run5324_23_8_exprMat_file.csv",
                           metadata.file = "Run5324_23_8_metadata_file.csv",
                           molecules.file = "Run5324_23_8_tx_file.csv",
                           segmentations.file = "Run5324_23_8-polygons.csv",
                           type = "segmentation",
                           mol.type = "pixels",
                           metadata = c("fov", "Area", "Mean.DAPI", "Max.DAPI"),
                           fov = "FOV")

for (name in names(obj@meta.data)) {

  if (is.data.frame(obj@meta.data[[name]]) | is.data.table(obj@meta.data[[name]])) {

    rownames(obj@meta.data[[name]]) <- rownames(obj@meta.data)

  } else {

    names(obj@meta.data[[name]]) <- names(obj@meta.data)

  }


}

saveRDS(obj, file = "/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/seurat_processed_object/SeuratObject_sqrtNormalized_SCTransform_UMAP_annotated_structured.RData")
