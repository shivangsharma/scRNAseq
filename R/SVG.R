library(tidyverse)
library(Seurat)
library(ggplot2)
library(SPARK)
obj <- readRDS("/media/singhlab/B684-19A61/CosMxJuly2022Experiment/data/seurat_processed_object/SeuratObject_sqrtNormalized_SCTransform_UMAP_annotated.RData")

info <- obj@meta.data$locations
rownames(info) <- rownames(obj@meta.data)
location <- as.matrix(info)
sp_count <- GetAssayData(object = obj, assay = "Nanostring", layer = "counts")

sparkX <- vector(mode = "list", length = obj@meta.data$fov %>% unique() %>% length())
names(sparkX) <- obj@meta.data$fov %>% unique()
gene_short_list <- vector(mode = "list", length = obj@meta.data$fov %>% unique() %>% length())
names(gene_short_list) <- obj@meta.data$fov %>% unique()
for (fov in unique(obj@meta.data$fov)) {

  cell_in_fov <- colnames(obj)[which(obj@meta.data$fov == fov)]
  location_fov <- location[cell_in_fov, ]
  sp_count_fov <- sp_count[, cell_in_fov]
  sparkX[[fov]] <- sparkx(sp_count_fov, location_fov, numCores = 1, option = "mixture")
  cum_gene_count <- rowSums(sp_count_fov) %>% sort(decreasing = TRUE) %>% log1p()
  cum_gene_count_short_list <- names(cum_gene_count)[cum_gene_count > mean(cum_gene_count) + 2*sd(cum_gene_count)]
  diff_expr_gene_short_list <- rownames(sparkX[[fov]][["res_mtest"]])[-log10(sparkX[[fov]][["res_mtest"]]$adjustedPval) >= 3]
  gene_short_list[[fov]] <- intersect(cum_gene_count_short_list, diff_expr_gene_short_list)

}

gene_list <- reduce(.x = gene_short_list, .f = intersect, .dir = "forward")
# [1] "KLK3"  "MZT2A" "NEAT1" "MSMB"  "RPL21" "HLA-C" "RPL34"



