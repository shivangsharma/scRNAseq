rawNanoStringToStandard <- function(raw_data_dir, exprMat_filename,
                                fov_location_filename, metadata_filename,
                                fov_annotation = FALSE,
                                fov_annotation_input_mode = "auto",
                                annotation_filename, slide_name, output_dir) {
  # raw_data_dir variable takes folder path where the raw data files are stored
  #
  # expr_path variable takes .csv file with expression matrix as input
  #
  # loc_path variable takes offset information of fov. The offset
  # information is provided in fov_position csv file
  #
  # meta_path variable takes the metadata csv file containing relative locations
  # of all the files
  #
  # fov_annotations = TRUE, if you want to add additional annotations to fov data
  #                   FALSE, otherwiswe (default). If TRUE, mode of input needs to be specified
  #
  # fov_annotation_input_mode = "auto", detects number of samples and associated annotations automatically
  #                             "manual", add annotations manually in interactive RShiny app
  #
  # annotation_filename: Takes file containing additional annotations as input.
  # Required only if fov_annotations is set to TRUE and fov_annotation_input_mode is set to "auto"
  #
  # output_dir requires name of the folder where standard raw data will be saved
  #
  # Example
  # rawNanoStringToGiotto(raw_data_dir = "~/Desktop/CosMx/Sample_1",
  #                       exprMat = "Run1008_21_7_exprMat_file.csv",
  #                       fov_pos = "Run1008_21_7_fov_positions_file.csv",
  #                       metadata = "Run1008_21_7_metadata_file.csv",
  #                       fov_annotation = TRUE,
  #                       fov_annotation_input_mode = "auto",
  #                       annotation_filename = "Run1008_21_7_annotation_file.csv",
  #                       slide_name = "Run1008_21_7")

  library(readxl)
  library(Matrix)
  library(tidyr)

  if(!(file.path(raw_data_dir, output_dir) %>% dir.exists())) {
    file.path(raw_data_dir, output_dir) %>% dir.create()
  }

  # 1. Read raw data
  expr_data <- read.csv(file = paste0(raw_data_dir, "/", exprMat_filename),
                        header = TRUE)
  offset_data <- read.csv(file = paste0(raw_data_dir, "/", fov_location_filename),
                          header = TRUE)
  meta_data <- read.csv(file = paste0(raw_data_dir, "/", metadata_filename),
                       header = TRUE)
  annotation_data <- read.csv(file = paste0(raw_data_dir, "/", annotation_filename),
                               header = TRUE)
  colnames(annotation_data)[1] <- "FOV"

  # 2. Remove cell ID 0 from expression data
  expr_data <- expr_data[!(expr_data$cell_ID == 0), ]
  row.names(expr_data) <- 1:nrow(expr_data)

  # 3. Create unique cell IDs
  unique_cell_ids <- paste0("SN-", slide_name, "FOV-", expr_data$fov, "CID-", expr_data$cell_ID)
  expr_data$cell_ID <- unique_cell_ids

  # 4. Generate annotation for each cell
  annotations <- colnames(annotation_data)[2:dim(annotation_data)[2]]
  df_cell_annotation <- data.frame(cell_ID = unique_cell_ids, fov = expr_data$fov)
  rownames(df_cell_annotation) <- df_cell_annotation$cell_ID

  for (x in annotations) {
    df_cell_annotation[x] <- NA
  }

  for (i in 1:dim(df_cell_annotation)[1]) {
    ind <- (annotation_data$FOV %in% df_cell_annotation$fov[i]) %>% which(arr.ind = TRUE)
    df_cell_annotation[i, annotations] <-
      annotation_data[ind, 2:(length(annotations) + 1)] %>%
      as.character() %>% t()
  }

  write.table(x = df_cell_annotation[, 2:dim(df_cell_annotation)[2]],
              file = file.path(raw_data_dir, output_dir, "annotations.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 5. Generate cell names file for each matrix
  write.table(x = df_cell_annotation$cell_ID, file =
                file.path(raw_data_dir, output_dir, "barcodes.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 6. Generate expression matrix without fov and cell ID column
  expr_data <- expr_data[, 3:dim(expr_data)[2]]
  expr_data <- Matrix(data = expr_data %>% as.matrix(),sparse = TRUE)
  writeMM(obj = expr_data, file = file.path(raw_data_dir, output_dir, "matrix.mtx"))

  # 7. Generate list of Cell_type from expression data
  gene_names <- colnames(expr_data)
  write.table(x = gene_names, file = file.path(raw_data_dir, output_dir,
                                               "features.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 8. Generate location data
  write.table(x = meta_data[, c(5, 6)], file = file.path(raw_data_dir, output_dir,
                                              "locations.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 9. Generate offset data
  write.table(x = offset_data, file = file.path(raw_data_dir, output_dir,
                                                         "offsets.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

rawNanoStringToGiotto <- function(my_working_dir, output_dir = "GiottoRawData") {

  library(tidyverse)
  library(stringr)

  exprMat_filename <- list.files(path = my_working_dir, pattern = ".*_exprMat_.*")
  metadata_filename <- list.files(path = my_working_dir, pattern = ".*_metadata_.*")
  annotation_filename <- list.files(path = my_working_dir, pattern = ".*_annotation_.*")
  fov_position_filename <- list.files(path = my_working_dir, pattern = ".*_fov_.*")

  slide_name <- list.files(path = my_working_dir, pattern = ".*_exprMat_.*") %>% str_extract(pattern = "[^exprMat]+")

  expr_data <-read.csv(file = file.path(my_working_dir, exprMat_filename), header = TRUE)
  expr_data <- expr_data[!(expr_data$cell_ID == 0), ]
  rownames(expr_data) <- 1:nrow(expr_data)
  expr_data$cell_ID <- paste0("cell_", rownames(expr_data))
  expr_data$fov <- NULL
  source("~/Desktop/CosMxJuly2022Experiment/Package functions/df_transpose.R")
  expr_data <- df_transpose(expr_data)
  colnames(expr_data) <- expr_data[1, ]
  expr_data <- expr_data[2:dim(expr_data)[1], ]
  colnames(expr_data) <- paste0(slide_name, colnames(expr_data))
  write.table(x = expr_data, file = file.path(my_working_dir, output_dir, paste0(slide_name,"expr_data.csv")), row.names = TRUE, quote = FALSE, sep = " ")

  cell_metadata <- read.csv(file = file.path(my_working_dir, metadata_filename))
  cell_metadata$cell_ID <- colnames(expr_data)
  cell_metadata$fov <- paste0(slide_name, cell_metadata$fov)
  cell_local_pos <- cell_metadata[, c(2,5,6)]
  colnames(cell_local_pos) <- c("ID", "X", "Y")
  write.table(x = cell_local_pos, file = file.path(my_working_dir, output_dir, paste0(slide_name,"cell_local_pos.csv")),  row.names = FALSE, quote = FALSE, sep = "\t")

  cell_ann <- cell_metadata[, -c(5,6)]
  ann_table <- read.csv(file = file.path(my_working_dir, annotation_filename))
  ann_table$FOV.. <- paste0(slide_name,ann_table$FOV..)
  ann_list <- list()
  ann_list[ann_table$FOV..] <- ann_table$Annotation.group.2
  cell_ann[, "PatientID"] <- ann_list[cell_ann$fov] %>% as.character()
  ann_list <- list()
  ann_list[ann_table$FOV..] <- ann_table$Annotation.group.1
  cell_ann[, "TissueSubtype"] <- ann_list[cell_ann$fov] %>% as.character()
  colnames(cell_ann)[2] <- "ID"
  # colnames(cell_ann) <- c("ID", "FOV", colnames(ann_table)[4]) # this step needs to be generalized for future uses
  write.table(x = cell_ann, file = file.path(my_working_dir, output_dir, paste0(slide_name,"cell_ann.csv")),  row.names = FALSE, quote = FALSE, sep = "\t")

  my_offset_data <- read.csv(file = file.path(my_working_dir, fov_position_filename))
  colnames(my_offset_data) <- c("field", "x_offset", "y_offset")
  my_offset_data$field <- paste0(slide_name, my_offset_data$field)
  write.table(x = my_offset_data, file = file.path(my_working_dir, output_dir, paste0(slide_name,"cell_offset_data.csv")),  row.names = FALSE, quote = FALSE, sep = "\t")
}

integrateRawGiottoData <- function(dir1, dir2, output_dir) {

  # dir1 <- "~/Desktop/CosMxJuly2022Experiment/Flat Files/5 Raw data/5 Raw data/Run5324_23_8/GiottoRawData"
  # dir2 <- "~/Desktop/CosMxJuly2022Experiment/Flat Files/5 Raw data/5 Raw data/Run5324_30_18/GiottoRawData"
  # output_dir <- "~/Desktop/CosMxJuly2022Experiment/Flat Files/5 Raw data/5 Raw data/CombinedGiottoData"

  expr_data_1_filename <- list.files(path = dir1, pattern = ".*_expr_data.*")
  cell_local_pos_1_filename <- list.files(path = dir1, pattern = ".*_cell_local_pos.*")
  cell_ann_1_filename <- list.files(path = dir1, pattern = ".*_cell_ann.*")
  my_offset_data_1_filename <- list.files(path = dir1, pattern = ".*_cell_offset_data.*")

  expr_data_2_filename <- list.files(path = dir2, pattern = ".*_expr_data.*")
  cell_local_pos_2_filename <- list.files(path = dir2, pattern = ".*_cell_local_pos.*")
  cell_ann_2_filename <- list.files(path = dir2, pattern = ".*_cell_ann.*")
  my_offset_data_2_filename <- list.files(path = dir2, pattern = ".*_cell_offset_data.*")

  expr_data_1 <- data.table::fread(file = file.path(dir1, expr_data_1_filename))
  rownames_expr_data <- expr_data_1$V1
  rownames(expr_data_1) <- rownames_expr_data
  expr_data_1$V1 <- NULL

  expr_data_2 <- data.table::fread(file = file.path(dir2, expr_data_2_filename))
  rownames_expr_data <- expr_data_2$V1
  rownames(expr_data_2) <- rownames_expr_data
  expr_data_2$V1 <- NULL

  expr_data_combined <- cbind(expr_data_1, expr_data_2)
  rownames(expr_data_combined) <- rownames(expr_data_1)
  write.table(x = expr_data_combined, file = file.path(output_dir, "expr_data.csv"), row.names = TRUE, quote = FALSE, sep = "\t")

  cell_local_pos_1 <- data.table::fread(file = file.path(dir1, cell_local_pos_1_filename))
  cell_local_pos_2 <- data.table::fread(file = file.path(dir2, cell_local_pos_2_filename))
  cell_local_pos_combined <- rbind(cell_local_pos_1, cell_local_pos_2)
  write.table(x = cell_local_pos_combined, file = file.path(output_dir, "cell_local_pos.csv"),  row.names = FALSE, quote = FALSE, sep = "\t")

  cell_ann_1 <- data.table::fread(file = file.path(dir1, cell_ann_1_filename))
  cell_ann_2 <- data.table::fread(file = file.path(dir2, cell_ann_2_filename))
  cell_ann_combined <- rbind(cell_ann_1, cell_ann_2)
  write.table(x = cell_ann_combined, file = file.path(output_dir, "cell_ann.csv"), row.names = FALSE, quote = FALSE, sep = "\t")

  my_offset_data_1 <-  data.table::fread(file = file.path(dir1, my_offset_data_1_filename))
  my_offset_data_2 <-  data.table::fread(file = file.path(dir2, my_offset_data_2_filename))
  my_offset_data_combined <- rbind(my_offset_data_1, my_offset_data_2)
  write.table(x = my_offset_data_combined, file = file.path(output_dir, "offset_data.csv"), row.names = FALSE, quote = FALSE, sep = "\t")
}

# barPlot_Spatial <- function(obj, Cell_type, groups, subtypes) {
#   library(ggplot2)
#   library(Giotto)
#   library(dplyr)
#   library(tidyverse)
#
#   df_combined <- data.frame(Cell_ID = gem@cell_metadata[["rna"]][["cell_ID"]])
#   temp <- strsplit(df_combined$Cell_ID, "_") %>% as.data.frame() %>% df_transpose()
#   df_combined$SlideID <- temp$`2`
#   df_combined$FOV <- temp$`3`
#   df_combined$CellName <- temp$`4`
#   rownames(df_combined) <- df_combined$Cell_ID
#
#   slideID <- df_combined$SlideID
#   fov <- df_combined$FOV
#   patientID <- rep(x = 0, times = length(slideID))
#   tumorType <- rep(x = "NA", times = length(slideID))
#
#   for (i in 1:length(slideID)) {
#     if (slideID[i] == 1) {
#       if (fov[i] %in% c(1,2,3,4,11,12,13,14,15,16)) {
#         tumorType[i] <- "Tumor"
#       }
#       else {
#         tumorType[i] <- "Other"
#       }
#       if (fov[i] %in% c(1:10)) {
#         patientID[i] <- 8
#       }
#       else {
#         patientID[i] <- 23
#       }
#     }
#
#     else {
#       if (fov[i] %in% c(1,2,3,4,5,6,11,12,13,14,15,16)) {
#         tumorType[i] <- "Tumor"
#       }
#       else {
#         tumorType[i] <- "Other"
#       }
#       if (fov[i] %in% c(1:10)) {
#         patientID[i] <- 18
#       }
#       else {
#         patientID[i] <- 30
#       }
#     }
#   }
#
#   df_combined$PatientID <- patientID
#   df_combined$TumorType <- tumorType
#
#   feat_idx <- (rownames(gem@expression[["rna"]][["normalized"]]) %in% c("CD8A", "CD8B", "GZMB")) %>% which()
#   temp <- gem@expression[["rna"]][["normalized"]][feat_idx, ] %>% df_transpose
#   df_combined_CD8A <- df_combined
#   df_combined_CD8A$Expr <- temp$CD8A
#   df_combined_CD8A$Gene <- "CD8A"
#   df_combined_CD8B <- df_combined
#   df_combined_CD8B$Expr <- temp$CD8B
#   df_combined_CD8B$Gene <- "CD8B"
#   df_combined_GZMB <- df_combined
#   df_combined_GZMB$Expr <- temp$GZMB
#   df_combined_GZMB$Gene <- "GZMB"
#   df_combined <- rbind(df_combined_CD8A, df_combined_CD8B, df_combined_GZMB)
#   df_combined$PatientID <- paste0("Patient ", df_combined$PatientID)
#
#   df_tumor <- df_combined[df_combined$TumorType == "Tumor", ]
#   ggplot(df_combined, aes(x = Gene))
#
#   ggplot(df_tumor, aes(fill = PatientID, y = Expr, x = Gene)) + geom_bar(position="dodge", stat="identity")
#   ggplot(df_tumor, aes(fill = PatientID, y = Expr, x = Gene)) + geom_boxplot() + coord_cartesian(ylim = c(0.1, 20))
#
#   df_summary <- data.frame(row.names = c(1:12))
#   df_summary$Gene <- c(rep("CD8A", times = 4), rep("CD8B", times = 4), rep("GZMB", times = 4))
#   df_summary$PatientID <- rep(c("Patient 8", "Patient 23", "Patient 18", "Patient 30"), times = 3)
#   df_summary$Mean <- NA
#   df_summary$CI <- NA
#   for(i in df_summary$Gene) {
#     for(j in df_summary$PatientID) {
#       temp <- df_tumor[((df_tumor$Gene %in% i) & (df_tumor$PatientID %in% j)) %>% which(), ]
#       df_summary[((df_summary$Gene %in% i) & (df_summary$PatientID %in% j)) %>% which(), 3] <- mean(temp$Expr)
#       df_summary[((df_summary$Gene %in% i) & (df_summary$PatientID %in% j)) %>% which(), 4] <- calc_CI(vector = temp$Expr, alpha = 0.05)
#     }
#   }
#
#   ggplot(df_summary, aes(fill = PatientID, y = Mean, x = Gene)) +
#     geom_bar(position=position_dodge(), stat="identity", color = "black") +
#     geom_errorbar(aes(ymin = Mean - CI, ymax = Mean + CI),
#                   width=0.2, colour="black",
#                   alpha=0.9, size=1.5, position = position_dodge(0.9))
# }

barPlot_Spatial <- function(gem, Cell_type) {
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(stringr)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  source("~/Desktop/CosMxJuly2022Experiment/Package functions/functions.R")

  # which columns are 23, 8, 30, 18
  patient_idx <- str_split(string = gem@cell_ID, pattern = "_") %>%
    as.data.frame() %>% df_transpose()
  colnames(patient_idx) <- c("c", "slide_no", "FOV", "cell_name")
  rownames(patient_idx) <- paste(patient_idx$c, patient_idx$slide_no,
                                 patient_idx$FOV, patient_idx$cell_name,
                                 sep = "_")
  patient_idx$data_idx <- c(1:length(gem@cell_ID))
  patient_idx$Patient_ID <- NA
  slide_1_idx <- patient_idx[which(patient_idx$slide_no %in% "1"), ]
  slide_2_idx <- patient_idx[which(patient_idx$slide_no %in% "2"), ]
  patient_8_idx <- slide_1_idx[which(slide_1_idx$FOV %in% c(1:10, 21)), ] %>%
    rownames()
  patient_23_idx <- slide_1_idx[which(slide_1_idx$FOV %in% c(11:20)), ] %>%
    rownames()
  patient_18_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c(1:10)), ] %>%
    rownames()
  patient_30_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c(11:20)), ] %>%
    rownames()
  patient_idx[patient_8_idx, "Patient_ID"] <- "Patient_8"
  patient_idx[patient_23_idx, "Patient_ID"] <- "Patient_23"
  patient_idx[patient_18_idx, "Patient_ID"] <- "Patient_18"
  patient_idx[patient_30_idx, "Patient_ID"] <- "Patient_30"

  patient_idx$Tissue <- NA
  tumor_slide_1_idx <- slide_1_idx[which(slide_1_idx$FOV %in%
                                           c(1:4, 11:16, 21)), ] %>%
    rownames()
  pin_slide_1_idx <- slide_1_idx[which(slide_1_idx$FOV %in% c(5, 10)), ] %>%
    rownames()
  stroma_slide_1_idx <- slide_1_idx[which(slide_1_idx$FOV %in%
                                            c(6, 7, 19, 20)), ] %>% rownames()
  benignEpi_slide_1_idx <- slide_1_idx[which(slide_1_idx$FOV %in%
                                               c(8, 9, 17, 18)), ] %>% rownames()
  tumor_slide_2_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c(1:6, 11:16)), ] %>%
    rownames()
  pin_slide_2_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c()), ] %>%
    rownames()
  stroma_slide_2_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c(9, 10)), ] %>%
    rownames()
  benignEpi_slide_2_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c(7, 8)), ] %>%
    rownames()
  lymphAggr_slide_2_idx <- slide_2_idx[which(slide_2_idx$FOV %in% c(17)), ] %>%
    rownames()
  patient_idx[tumor_slide_1_idx, "Tissue"] <- "Tumor"
  patient_idx[pin_slide_1_idx, "Tissue"] <- "PIN"
  patient_idx[stroma_slide_1_idx, "Tissue"] <- "Stroma"
  patient_idx[benignEpi_slide_1_idx, "Tissue"] <- "BenignEpi"
  patient_idx[tumor_slide_2_idx, "Tissue"] <- "Tumor"
  patient_idx[pin_slide_2_idx, "Tissue"] <- "PIN"
  patient_idx[stroma_slide_2_idx, "Tissue"] <- "Stroma"
  patient_idx[benignEpi_slide_2_idx, "Tissue"] <- "BenignEpi"
  patient_idx[lymphAggr_slide_2_idx, "Tissue"] <- "LymphAggr"
  # extract expression of CD8A, CD8B, GZMB (or any other Cell_type)
  patient_idx[, Cell_type] <- gem@expression$rna$normalized[Cell_type, ] %>%
    df_transpose()
  patient_idx$leiden_clus <- gem@cell_metadata[["rna"]][["leiden_clus"]]
  patient_expr <- patient_idx
  patient_expr[, c("c", "slide_no", "FOV", "cell_name", "data_idx")] <- NULL
  patient_expr <- patient_expr[patient_expr$Tissue %in% "Tumor",]
  patient_expr[, c("Tissue", "leiden_clus")] <- NULL

  # melt the data frame
  patient_expr <- melt(patient_expr,id=c("Patient_ID"), variable.name = "Cell_type",
                       value.name = "Norm_Expr")
  patient_expr$Cell_type <- as.factor(patient_expr$Cell_type) %>%
    fct_relevel(levels = Cell_type)
  patient_expr$Patient_ID <- as.factor(patient_expr$Patient_ID) %>%
    fct_relevel("Patient_8", "Patient_23", "Patient_18", "Patient_30")

  summarized_data <- patient_expr %>%
    group_by(Cell_type, Patient_ID) %>%
    summarize(Mean = mean(Norm_Expr),
              CI = calc_CI(vector = Norm_Expr, alpha = 0.05))

  p <- ggplot(summarized_data, aes(x = Cell_type, y = Mean, fill = Patient_ID)) +
    geom_bar(position = "dodge", stat = "identity", colour = "black") +
    geom_errorbar(aes(x = Cell_type, ymin = Mean - CI, ymax = Mean + CI, ),
                  width=0.2, colour="black", alpha=0.9, size=0.7,
                  position = position_dodge(.9)) +
    theme(
      axis.title = element_text( color="black", size=16, face=1),
      axis.line = element_line(size = 0.5, colour = "black", linetype=1),
      axis.text = element_text(angle = 0, color="black", size=16, face=1),
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(color="black", size=16, face=1),
      legend.title = element_text(color="black", size=14, face=1)
    ) +
    labs(title = "CD8B, GZMB and CD58 Normalized Expression in Tumor", x = "Cell_type",
         y = "Expression", fill = "Patient_ID") +
    scale_y_continuous(expand=c(0,0))


  stat_summarized_data <- patient_expr %>%
    group_by(Cell_type) %>%
    t_test(Norm_Expr ~ Patient_ID) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  stat_summarized_data <- stat_summarized_data %>%
    add_xy_position(x = "Cell_type", dodge = 0.9, fun = "mean_ci")

  p2 <- p + stat_pvalue_manual(stat_summarized_data,
                               label = "p.adj.signif",
                               tip.length = 0.01,
                               bracket.nudge.y = 0, bracket.size = 0.7, size = 6, scales = "fixed",
                               inherit.aes = FALSE) + ylim(0, max(stat_summarized_data$y.position)*1.05)
  p2

  ############################### Alternate code ###############################
  # p3 <- ggbarplot(
  #   patient_expr, x = "Cell_type", y = "Norm_Expr", fill = "Patient_ID",
  #   palette = "npg", add = "mean_ci",
  #   position = position_dodge(0.7)
  # )
  # p4 <- p3 + stat_pvalue_manual(stat_summarized_data,
  #                               label = "p.adj.signif",
  #                               tip.length = 0.00005, scales = "free")
}


barPlot_Spatial <- function(gem, Cell_type) {
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(stringr)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  source("~/Desktop/CosMxJuly2022Experiment/Package functions/data_transformation_functions.R")
  source("~/Desktop/CosMxJuly2022Experiment/Package functions/statistical_functions.R")

  obj <- readRDS("~/Desktop/CosMxJuly2022Experiment/R Seurat Oject/SeuratObject_sqrtNormalized_SCTransform_UMAP.RData")
  patient_expr <- matrix(nrow = dim(obj)[2], ncol = 4) %>% as.data.frame()
  colnames(patient_expr) <- c("CD44", "Patient_ID", "Tissue_Type", "Cell_type")

  patient_expr$CD44 <- gem@expression$rna$sqrtnormalized["CD44", ]
  #patient_expr$CD44 <- obj@assays$RNA@data["CD44", ]
  patient_expr$Patient_ID <- factor(x = obj@meta.data$Patient_ID,
                                    levels = c("Patient_8", "Patient_23",
                                               "Patient_18", "Patient_30"))
  patient_expr$Tissue_Type <- factor(x = obj@meta.data$Tissue_Type,
                                     levels = c("Tumor", "Stroma", "LymphAggr",
                                                "BenignEpi", "PIN"))
  patient_expr$Cell_type <- factor(x = obj@meta.data$CellAnnCoarse,
                                   levels = c("Luminal", "Basal", "Fibroblast",
                                              "Endothelial", "VSMC", "Myeloid",
                                              "TCell", "BCell", "Mast", "SMC"))
  patient_expr <- patient_expr[patient_expr$Tissue_Type %in% "Tumor",]

  patient_expr$Tissue_Type <- NULL
  colnames(patient_expr)[1] <- "Norm_Expr"

  # melt the data frame
  patient_expr <- melt(patient_expr,id=c("Patient_ID"), variable.name = "CD44",
                       value.name = "Norm_Expr")
  patient_expr$Cell_type <- as.factor(patient_expr$Cell_type) %>%
    fct_relevel(levels = Cell_type)
  patient_expr$Patient_ID <- as.factor(patient_expr$Patient_ID) %>%
    fct_relevel("Patient_8", "Patient_23", "Patient_18", "Patient_30")

  summarized_data <- patient_expr %>%
    group_by(Cell_type, Patient_ID) %>%
    summarize(Mean = mean(Norm_Expr),
              CI = calc_CI(vector = Norm_Expr, alpha = 0.05), Counts = n())
  summarized_data <- summarized_data[summarized_data$Counts >= 200, ]

  p <- ggplot(summarized_data, aes(x = Cell_type, y = Mean, fill = Patient_ID)) +
    geom_bar(position = "dodge", stat = "identity", colour = "black") +
    geom_errorbar(aes(x = Cell_type, ymin = Mean - CI, ymax = Mean + CI, ),
                  width=0.2, colour="black", alpha=0.9, size=0.7,
                  position = position_dodge(.9)) +
    theme(
      axis.title = element_text( color="black", size=16, face=1),
      axis.line = element_line(size = 0.5, colour = "black", linetype=1),
      axis.text = element_text(angle = 0, color="black", size=16, face=1),
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(color="black", size=16, face=1),
      legend.title = element_text(color="black", size=14, face=1)
    ) +
    labs(title = "CD44 Normalized Expression in Tumor", x = "Cell_type",
         y = "Expression", fill = "Patient_ID") +
    scale_y_continuous(expand=c(0,0))


  stat_summarized_data <- patient_expr %>%
    group_by(Cell_type) %>%
    t_test(Norm_Expr ~ Patient_ID) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  stat_summarized_data <- stat_summarized_data %>%
    add_xy_position(x = "Cell_type", dodge = 0.9, fun = "mean_ci")

  p2 <- p + stat_pvalue_manual(stat_summarized_data,
                               label = "p.adj.signif",
                               tip.length = 0.01,
                               bracket.nudge.y = 0, bracket.size = 0.7, size = 6, scales = "fixed",
                               inherit.aes = FALSE) + ylim(0, max(stat_summarized_data$y.position)*1.05)
  p2

  ############################### Alternate code ###############################
  # p3 <- ggbarplot(
  #   patient_expr, x = "Cell_type", y = "Norm_Expr", fill = "Patient_ID",
  #   palette = "npg", add = "mean_ci",
  #   position = position_dodge(0.7)
  # )
  # p4 <- p3 + stat_pvalue_manual(stat_summarized_data,
  #                               label = "p.adj.signif",
  #                               tip.length = 0.00005, scales = "free")
}

featurePlot_panel <- function(gem, markers, save_path) {
  plot_data <- as.data.frame(gem@dimension_reduction[["cells"]][["umap"]][["umap"]][["coordinates"]])
  useful_expr <- gem@expression[["rna"]][["normalized"]][markers, ]
  plot_data[, markers] <- useful_expr %>% t()
  plot <- vector(mode = "list")
  for (i in colnames(plot_data)[-c(1,2)]) {
    plot[[i]] <- ggplot(data = plot_data, aes(x = Dim.1, y = Dim.2)) +
      geom_point(aes(color = !!sym(i)), size = 0.3, stroke = 0, shape = 16) +
      scale_colour_gradient(low = "grey80", high = "#00008b") +
      labs(title = i) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0(i, ".tif"),
           plot[[i]],
           path = save_path,
           width = 1600,
           height = 900,
           units = "px",
           device = "tiff"
           )
  }
  #grid.arrange(grobs = plot[1:length(markers)], ncol = 5, nrow = 2)
}

violinPlot <- function(obj, Cell_type, save_path) {



}
