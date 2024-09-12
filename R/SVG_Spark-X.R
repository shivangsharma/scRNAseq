library(tidyverse)
library(Seurat)
library(ggplot2)
library(SPARK)
obj <- readRDS("/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/seurat_processed_object/SeuratObject_sqrtNormalized_SCTransform_UMAP_annotated.RData")

info <- obj@meta.data$Centroid
rownames(info) <- rownames(obj@meta.data)
location <- as.matrix(info)
sp_count <- GetAssayData(object = obj, assay = "Nanostring", layer = "counts")

sparkX <- vector(mode = "list", length = obj@meta.data$Category$FOV %>% unique() %>% length())
names(sparkX) <- obj@meta.data$Category$FOV %>% unique()
gene_short_list <- vector(mode = "list", length = obj@meta.data$Category$FOV %>% unique() %>% length())
names(gene_short_list) <- obj@meta.data$Category$FOV %>% unique()
fov_list <- vector(mode = "list", length = obj@meta.data$Category$Slide_name %>% unique() %>% length())
names(fov_list) <- levels(obj@meta.data$Category$Slide_name)

for (slide in levels(obj@meta.data$Category$Slide_name)) {

  fov_list[[slide]] <- obj@meta.data$Category$FOV[obj@meta.data$Category$Slide_name %in% slide] %>% unique()

}

for (slide in names(fov_list)) {

  for (fov in fov_list[[slide]]) {

    cell_in_fov <- colnames(obj)[which((obj@meta.data$Category$Slide_name == slide) & (obj@meta.data$Category$FOV == fov))]
    location_fov <- location[cell_in_fov, ]
    sp_count_fov <- sp_count[, cell_in_fov]
    sparkX[[fov]] <- sparkx(sp_count_fov, location_fov, numCores = 1, option = "mixture")
    cum_gene_expr <- rowSums(sp_count_fov) %>% sort(decreasing = TRUE)
    cum_gene_expr_short_list <- names(cum_gene_expr[cum_gene_expr > mean(cum_gene_expr) + 1*sd(cum_gene_expr)])
    spatial_var_gene_short_list <- rownames(sparkX[[fov]][["res_mtest"]])[sparkX[[fov]]$res_mtest$adjustedPval < 0.001]
    gene_short_list[[fov]] <- intersect(cum_gene_expr_short_list, spatial_var_gene_short_list)

  }

}

gene_vec <- purrr::reduce(.x = gene_short_list, .f = intersect, .dir = "backward")
# > "KLK3"  "MZT2A" "NEAT1" "MALAT1"





cents <- CreateCentroids(coords = obj@meta.data$locations,  nsides = Inf, radius = 15)
segm <-
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = NULL
  )

coords <- CreateFOV(
  coords = segmentations.data,
  type = c("segmentation"),
  molecules = NULL,
  assay = "CosMxSMI"
)


