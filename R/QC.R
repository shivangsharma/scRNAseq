plot_lib_sizes <- function(obj) {

  data_obj <- LayerData(obj,
                        assay = "Nanostring",
                        layer = "counts")
  categories <- obj@meta.data$Category %>% colnames()
  df <- data.frame(matrix(data = colSums(data_obj) %>% log1p(),
                          ncol = 1,
                          dimnames = list(colnames(data_obj), c("logCounts"))
                          )
                   )
  df <- cbind(df, obj@meta.data$Category)
  p <- vector(mode = "list", length = length(categories))
  names(p) <- categories
  parent_category <- NA
  for (category in categories) {

    # print(category)
    child_category <- category
    df$y <- NA
    max_y <- 0

    if (is.na(parent_category)) {

      for (child_category_level in levels(df[[child_category]])) {

        df_subset <- df[df[[child_category]] %in% child_category_level, ]
        density_table <- density(df_subset$logCounts)
        max_y <- ifelse(test = max(density_table$y) > max_y,
                        yes = max(density_table$y),
                        no = max_y)

      }
      rel_dist <- max_y * seq(from = 1.05,
                              by = 0.05,
                              length.out = df[[child_category]] %>%
                                levels() %>%
                                length())
      names(rel_dist) <- levels(df[[child_category]])
      df$y <- rel_dist[df[[child_category]]]
      p[[child_category]] <- ggplot(data = df,
                                    mapping = aes(x = logCounts,
                                                  color = !!sym(child_category)
                                    )
      ) +
        geom_density(aes(y = after_stat(density)),
                     show.legend = FALSE) +
        stat_density(geom = "line",
                     position = "identity") +
        geom_boxplot(aes(y = y),
                     width = 0.025 * max_y,
                     show.legend = FALSE,
                     outlier.shape = NA) +
        labs(x = "log2(Counts + 1)",
             y = "Cell count",
             color = str_replace(string = child_category,
                                 pattern = "_",
                                 replacement = " ")) +
        theme_classic()

    } else {

      parent_child_categories <- vector(mode = "list",
                                        length = df[[parent_category]] %>%
                                          levels() %>%
                                          length())
      names(parent_child_categories) <- levels(df[[parent_category]])
      for (parent_category_level in names(parent_child_categories)) {

        df_subset <- df[df[[parent_category]] %in% parent_category_level, ]
        parent_child_categories[[parent_category_level]] <-
          unique(df_subset[[child_category]])

      }

      for (parent_category_level in names(parent_child_categories)) {

        for (child_category_level in parent_child_categories[[parent_category_level]]) {

          density_table <- df[(df[[parent_category]] %in%
                                 parent_category_level) &
                                (df[[child_category]] %in%
                                   child_category_level), ]$logCounts %>%
            density()
          max_y <- ifelse(test = max(density_table$y) > max_y,
                          yes = max(density_table$y),
                          no = max_y)

        }

      }

      rel_dist <- max_y * seq(from = 1.05,
                              by = 0.05,
                              length.out = parent_child_categories %>%
                                unlist() %>%
                                length())
      i <- 1
      for (parent_category_level in names(parent_child_categories)) {

        for (child_category_level in parent_child_categories[[parent_category_level]]) {

          df$y[(df[[parent_category]] %in% parent_category_level) & (df[[child_category]] %in% child_category_level)] <- rel_dist[i]
          i <- i + 1

        }

      }

      p[[child_category]] <- ggplot(data = df,
                                    mapping = aes(x = logCounts,
                                                  color = !!sym(child_category),
                                                  linetype = as.factor(!!sym(parent_category))
                                    )
      ) +
        geom_density(aes(y = after_stat(density)),
                     show.legend = FALSE) +
        stat_density(geom = "line",
                     position = "identity") +
        geom_boxplot(aes(y = y),
                     width = 0.025 * max_y,
                     show.legend = FALSE,
                     outlier.shape = NA) +
        labs(x = "log2(Counts + 1)",
             y = "Cell count",
             color = str_replace(string = child_category,
                                 pattern = "_",
                                 replacement = " "),
             linetype = str_replace(string = parent_category,
                                    pattern = "_",
                                    replacement = " ")
        ) +
        theme_classic()

    }

    parent_category <- child_category

  }

  return(p)

}

plot_gene_abundance <- function(obj, plot_expressed_genes = FALSE) {

  data_obj <- LayerData(obj,
                        assay = "Nanostring",
                        layer = "counts")
  if (plot_expressed_genes == TRUE) {

    data_obj[data_obj > 1] <- 1
    data_obj[data_obj <= 1] <- 0

  } else {

    data_obj <- data_obj

  }
  meta_data <- obj@meta.data$Category
  categories <- colnames(meta_data)
  p <- vector(mode = "list", length = length(categories))
  names(p) <- categories
  parent_category <- NA
  for (category in categories) {

    child_category <- category
    max_y <- 0
    df$y <- NA

    if (is.na(parent_category)) {

      df <- data.frame(matrix(nrow = 0, ncol = 2))
      colnames(df) <- c("logCounts", child_category)
      for (child_category_level in levels(meta_data[[child_category]])) {

        cell_IDs <- rownames(meta_data[meta_data[[child_category]] %in% child_category_level, ])
        data_obj_subset <- data_obj[, cell_IDs]
        gene_abundance <- rowSums(data_obj_subset) %>% log1p()
        df_temp <- data.frame(logCounts = gene_abundance,
                              V1 = child_category_level)
        colnames(df_temp)[2] <- child_category
        df <- rbind(df, df_temp)

      }

      for (child_category_level in levels(meta_data[[child_category]])) {

        df_subset <- df[df[[child_category]] %in% child_category_level, ]
        density_table <- density(df_subset$logCounts)
        max_y <- ifelse(test = max(density_table$y) > max_y,
                        yes = max(density_table$y),
                        no = max_y)

      }

      rel_dist <- max_y * seq(from = 1.05,
                              by = 0.05,
                              length.out = df[[child_category]] %>%
                                levels() %>%
                                length())
      names(rel_dist) <- levels(df[[child_category]])
      df$y <- rel_dist[df[[child_category]]]

      p[[child_category]] <- ggplot(data = df,
                                    mapping = aes(x = logCounts,
                                                  color = !!sym(child_category)
                                                  )
                                    ) +
        geom_density(aes(y = after_stat(density)),
                     show.legend = FALSE) +
        stat_density(geom = "line",
                     position = "identity") +
        geom_boxplot(aes(y = y),
                     width = 0.025 * max_y,
                     show.legend = FALSE,
                     outlier.shape = NA) +
        labs(x = "log2(Counts + 1)",
             y = "Cell count",
             color = str_replace(string = child_category,
                                 pattern = "_",
                                 replacement = " ")) +
        theme_classic()

    } else {

      df <- data.frame(matrix(nrow = 0, ncol = 3))
      colnames(df) <- c("logCounts", parent_category, child_category)
      parent_child_categories <- vector(mode = "list",
                                        length = meta_data[[parent_category]] %>%
                                          levels() %>%
                                          length())
      names(parent_child_categories) <- meta_data[[parent_category]] %>% levels()

      for (parent_category_level in levels(meta_data[[parent_category]])) {

        df_subset <- meta_data[(meta_data[[parent_category]] %in% parent_category_level), ]
        parent_child_categories[[parent_category_level]] <- df_subset[[child_category]] %>% unique()

      }

      for (parent_category_level in names(parent_child_categories)) {

        for (child_category_level in parent_child_categories[[parent_category_level]]) {

          cell_IDs <- rownames(meta_data[(meta_data[[parent_category]] %in% parent_category_level) &
                                           (meta_data[[child_category]] %in% child_category_level), ])
          gene_abundance <- rowSums(data_obj[, cell_IDs]) %>% log1p()
          df_temp <- data.frame(logCounts = gene_abundance,
                                V1 = parent_category_level,
                                V2 = child_category_level)
          colnames(df_temp)[c(2, 3)] <- c(parent_category, child_category)
          df <- rbind(df, df_temp)

        }

      }

      for (parent_category_level in names(parent_child_categories)) {

        for (child_category_level in parent_child_categories[[parent_category_level]]) {

          df_subset <- df[(df[[parent_category]] %in% parent_category_level) & (df[[child_category]] %in% child_category_level), ]
          density_table <- density(df_subset$logCounts)
          max_y <- ifelse(test = max(density_table$y) > max_y,
                          yes = max(density_table$y),
                          no = max_y)

        }

      }
      rel_dist <- max_y * seq(from = 1.05,
                              by = 0.05,
                              length.out = parent_child_categories %>%
                                unlist() %>%
                                length())
      i <- 1
      for (parent_category_level in names(parent_child_categories)) {

        for (child_category_level in parent_child_categories[[parent_category_level]]) {

          df$y[(df[[parent_category]] %in% parent_category_level) & (df[[child_category]] %in% child_category_level)] <- rel_dist[i]
          i <- i + 1

        }

      }

      p[[child_category]] <- ggplot(data = df,
                                    mapping = aes(x = logCounts,
                                                  color = !!sym(child_category),
                                                  linetype = as.factor(!!sym(parent_category))
                                                  )
                                    ) +
        geom_density(aes(y = after_stat(density)),
                     show.legend = FALSE) +
        stat_density(geom = "line",
                     position = "identity") +
        geom_boxplot(aes(y = y),
                     width = 0.025 * max_y,
                     show.legend = FALSE,
                     outlier.shape = NA) +
        labs(x = "log2(Counts + 1)",
             y = "Cell count",
             color = str_replace(string = child_category,
                                 pattern = "_",
                                 replacement = " "),
             linetype = str_replace(string = parent_category,
                                    pattern = "_",
                                    replacement = " ")
        ) +
        theme_classic()

    }

    parent_category <- child_category

  }

  return(wrap_plots(p))

}

plot_most_expressed <- function(obj) {

  data_obj <- LayerData(obj, assay = "Nanostring", layer = "counts")
  parent_category <- "Slide_name"
  child_category <- "FOV"
  parent_child_categories <- vector(mode = "list",
                                    length = obj@meta.data$Category[[parent_category]] %>%
                                      levels() %>%
                                      length())
  names(parent_child_categories) <- obj@meta.data$Category[[parent_category]] %>%
    levels()
  meta_data <- obj@meta.data$Category

  for (parent_category_level in names(parent_child_categories)) {

    meta_data_subset <- meta_data[meta_data[[parent_category]] %in% parent_category_level, ]
    parent_child_categories[[parent_category_level]] <- meta_data_subset[[child_category]] %>% unique()

  }

  common_top_genes <- c()
  df <- data.frame(matrix(nrow = nrow(data_obj),
                          ncol = 0,
                          dimnames = rownames(data_obj) %>% list()
                          )
                   )
  cell_counts <- list()
  for (parent_category_level in names(parent_child_categories)) {

    for (child_category_level in parent_child_categories[[parent_category_level]]) {

      cell_IDs <- rownames(meta_data)[(meta_data[[parent_category]] %in% parent_category_level) &
                                        (meta_data[[child_category]] %in% child_category_level)]
      gene_abundance <- rowSums(data_obj[, cell_IDs]) %>% sort(decreasing = TRUE)
      df[[paste(parent_category_level, child_category_level, sep = "_")]] <- gene_abundance
      top_genes <- names(gene_abundance)[c(1:20)]
      common_top_genes <- union(top_genes, common_top_genes)
      cell_counts[[paste(parent_category_level, child_category_level, sep = "_")]] <- length(cell_IDs)

    }

  }
  df_subset <- df[common_top_genes, ] %>%
    log1p() %>%
    df_transpose() %>%
    scale() %>%
    df_transpose()
  lim_heatmap_max <- max(df_subset)
  lim_heatmap_min <- min(df_subset)
  col_heatmap <- colorRamp2(c(lim_heatmap_min, 0, lim_heatmap_max),
                            c("blue", "white", "red"))

  anno_df <- obj@meta.data$Category[!(obj@meta.data$Category %>% duplicated()), ]
  anno_df[["FOV"]] <- NULL
  anno_df[["log1p(T.T.C.)"]] <- colSums(df) %>% log1p()
  anno_df[["Cell counts"]] <- unlist(cell_counts)

  col_ha <- HeatmapAnnotation(df = anno_df)

  p <- Heatmap(matrix = as.matrix(df_subset),
               name = "log1p(Counts)",
               row_title = "Top genes",
               row_title_side = "left",
               row_title_gp = gpar(fontsize = 10),
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 5),
               cluster_rows = FALSE,
               column_title = "FOV",
               column_title_side = "bottom",
               column_title_gp = gpar(fontsize = 10),
               column_labels = colnames(df) %>% str_remove(string = ., pattern = ".*_"),
               column_names_side = "bottom",
               column_names_rot = 45,
               column_names_gp = gpar(fontsize = 5),
               cluster_columns = FALSE,
               col = col_heatmap,
               top_annotation = col_ha
               )

  return(p)

}

plot_marker_genes <- function(obj, annotation_column = NA) {

  child_category <- "Slide_name"
  data_obj <- LayerData(obj,
                        assay = "Nanostring",
                        layer = "counts")
  meta_data_image <- obj@meta.data$Image
  meta_data_category <- obj@meta.data$Category
  markers <- c("CD298", "CD45", "PanCK", "CD3", "DAPI")
  protein_markers <- markers %>% paste0("Max.", .)

  data_markers <- data_obj[c("B2M", "PTPRC"), ]
  panck_markers <- rownames(data_obj) %>%
    str_detect(string = ., pattern = "KRT") %>%
    rownames(data_obj)[.]
  panck_data <- data_obj[panck_markers, ] %>% colSums()
  cd3_markers <- c("CD3E", "CD3D", "CD3G")
  cd3_data <- data_obj[cd3_markers, ] %>% colSums()
  data_markers <- rbind(data_markers,
                        panck_data,
                        cd3_data,
                        rep(x = 0, length = ncol(data_obj)))
  rownames(data_markers) <- markers %>% paste0("RNA.", .)

  df <- cbind(meta_data_image[, protein_markers] %>% log1p(),
              t(data_markers) %>% log1p(),
              Slide_name = meta_data_category[, c("Slide_name")])

  p <- list()

  for (marker in markers) {

    df$Anno.CD3 <- obj@meta.data$Annotation[[annotation_column]] %>% as.character()
    df$Anno.CD3[df$Anno.CD3 %in% c("TCell")] <- "TCell"
    df$Anno.CD3[!(df$Anno.CD3 %in% c("TCell"))] <- "Other"
    df$Anno.CD3 <- factor(df$Anno.CD3, levels = c("TCell", "Other"))
    df$Anno.CD45 <- obj@meta.data$Annotation[[annotation_column]] %>% as.character()
    df$Anno.CD45[df$Anno.CD45 %in% c("Myeloid", "TCell", "Bcell", "Mast")] <- "Immune"
    df$Anno.CD45[!(df$Anno.CD45 %in% c("Immune"))] <- "Other"
    df$Anno.CD45 <- factor(df$Anno.CD45, levels = c("Immune", "Other"))
    df$Anno.PanCK <- obj@meta.data$Annotation[[annotation_column]] %>% as.character()
    df$Anno.PanCK[df$Anno.PanCK %in% c("Luminal", "Basal", "Hillock", "Club", "MSMB+")] <- "Epithelial"
    df$Anno.PanCK[!(df$Anno.PanCK %in% c("Epithelial"))] <- "Other"
    df$Anno.PanCK <- factor(df$Anno.PanCK, levels = c("Epithelial", "Other"))
    df$Anno.CD298 <- "Cell" %>% as.factor()
    df$Anno.DAPI <- "Cell" %>% as.factor()

    color_scheme <- c("red", "blue")
    names(color_scheme) <- levels(df[[paste0("Anno.", marker)]])
    alpha_scheme <- c(0.1, 0.1)
    names(alpha_scheme) <- levels(df[[paste0("Anno.", marker)]])

    p1 <-
      ggplot(data = df,
             mapping = aes(x = !!sym(paste0("Max.", marker)),
                           y = !!sym(paste0("RNA.", marker)),
                           color = !!sym(paste0("Anno.", marker)),
                           alpha = !!sym(paste0("Anno.", marker))
                           )
             ) +
      geom_point() +
      scale_color_manual(values = color_scheme) +
      scale_alpha_manual(values = alpha_scheme) +
      theme_classic() +
      theme(legend.position = "bottom",
            panel.background = element_rect(colour = "black", linewidth = 1),
            axis.line = element_line(linewidth = 0.0),
            plot.margin = unit(c(-5, -5, -5, -5), units = "mm"))
    p2 <- ggplot(data = df,
                 mapping = aes(x = !!sym(paste0("RNA.", marker)),
                               fill = !!sym(paste0("Anno.", marker))
                               )
                 ) +
      geom_histogram(aes(y = after_stat(log1p(count))), position = "dodge") +
      scale_fill_manual(values = color_scheme) +
      labs(y = "loge(counts)") +
      coord_flip() +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
      theme_classic() +
      theme(legend.position = "bottom",
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = unit(c(-5, -5, -5, -5), units = "mm"))
    p3 <- ggplot(data = df,
                 mapping = aes(x = !!sym(paste0("Max.", marker)),
                               linetype = !!sym(child_category))) +
      geom_density(aes(y = after_stat(log1p(count)))) +
      labs(y = "log count density") +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "top",
            plot.margin = unit(c(-5, -5, -5, -5), units = "mm"))

    layout <- "
    11112
    33333
    44445"
    p[[paste0("Intensity vs RNA for ", marker)]] <- p3 +
      guide_area() +
      plot_spacer() +
      p1 +
      p2 +
      plot_layout(design = layout,
                  guides = "collect",
                  heights = c(10, 0, 10)) &
      theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))

  }

  return(p)

}

plot_negprobes_hk_genes <- function(obj) {

  data_obj <- LayerData(obj, assay = "Nanostring", layer = "counts")
  parent_category <- "Slide_name"
  child_category <- "FOV"
  parent_child_categories <- vector(mode = "list",
                                    length = obj@meta.data$Category[[parent_category]] %>%
                                      levels() %>%
                                      length())
  names(parent_child_categories) <- obj@meta.data$Category[[parent_category]] %>%
    levels()
  meta_data <- obj@meta.data$Category

  for (parent_category_level in names(parent_child_categories)) {

    meta_data_subset <- meta_data[meta_data[[parent_category]] %in% parent_category_level, ]
    parent_child_categories[[parent_category_level]] <- meta_data_subset[[child_category]] %>% unique()

  }

  neg_probes <- rownames(data_obj) %>%
    str_detect(string = ., pattern = ("NegPrb|Negative")) %>%
    rownames(data_obj)[.]
  df <- data.frame(matrix(nrow = nrow(data_obj),
                          ncol = 0,
                          dimnames = rownames(data_obj) %>% list()
                          )
                   )
  cell_counts <- list()
  for (parent_category_level in names(parent_child_categories)) {

    for (child_category_level in parent_child_categories[[parent_category_level]]) {

      cell_IDs <- rownames(meta_data)[(meta_data[[parent_category]] %in% parent_category_level) &
                                        (meta_data[[child_category]] %in% child_category_level)]
      gene_abundance <- rowSums(data_obj[, cell_IDs]) %>% sort(decreasing = TRUE)
      df[[paste(parent_category_level, child_category_level, sep = "_")]] <- gene_abundance
      cell_counts[[paste(parent_category_level, child_category_level, sep = "_")]] <- length(cell_IDs)

    }

  }
  df_subset <- df[neg_probes, ] %>%
    log1p() %>%
    df_transpose() %>%
    scale() %>%
    df_transpose()
  lim_heatmap_max <- max(df_subset)
  lim_heatmap_min <- min(df_subset)
  col_heatmap <- colorRamp2(c(lim_heatmap_min, 0, lim_heatmap_max),
                            c("blue", "white", "red"))

  anno_df <- obj@meta.data$Category[!(obj@meta.data$Category %>% duplicated()), ]
  anno_df[["FOV"]] <- NULL
  anno_df[["log1p(T.T.C.)"]] <- colSums(df) %>% log1p()
  anno_df[["Cell counts"]] <- unlist(cell_counts)

  col_ha <- HeatmapAnnotation(df = anno_df)

  p <- Heatmap(matrix = as.matrix(df_subset),
               name = "log1p(Counts)",
               row_title = "Top genes",
               row_title_side = "left",
               row_title_gp = gpar(fontsize = 10),
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 5),
               cluster_rows = FALSE,
               column_title = "FOV",
               column_title_side = "bottom",
               column_title_gp = gpar(fontsize = 10),
               column_labels = colnames(df) %>% str_remove(string = ., pattern = ".*_"),
               column_names_side = "bottom",
               column_names_rot = 45,
               column_names_gp = gpar(fontsize = 5),
               cluster_columns = FALSE,
               col = col_heatmap,
               top_annotation = col_ha
  )

  return(p)

}

plot_pca_fovs <- function(obj) {

  data_obj <- LayerData(obj, assay = "Nanostring", layer = "counts")
  norm_data_obj <- t(t(data_obj) / colSums(data_obj))
  parent_category <- "Slide_name"
  child_category <- "FOV"
  parent_child_categories <- vector(mode = "list",
                                    length = obj@meta.data$Category[[parent_category]] %>%
                                      levels() %>%
                                      length())
  names(parent_child_categories) <- obj@meta.data$Category[[parent_category]] %>%
    levels()
  meta_data <- obj@meta.data$Category

  for (parent_category_level in names(parent_child_categories)) {

    meta_data_subset <- meta_data[meta_data[[parent_category]] %in% parent_category_level, ]
    parent_child_categories[[parent_category_level]] <- meta_data_subset[[child_category]] %>% unique()

  }

  df <- data.frame(matrix(nrow = nrow(norm_data_obj),
                          ncol = 0,
                          dimnames = rownames(norm_data_obj) %>% list()
                          )
                   )
  cell_counts <- list()
  for (parent_category_level in names(parent_child_categories)) {

    for (child_category_level in parent_child_categories[[parent_category_level]]) {

      cell_IDs <- rownames(meta_data)[(meta_data[[parent_category]] %in% parent_category_level) &
                                        (meta_data[[child_category]] %in% child_category_level)]
      gene_abundance <- rowSums(norm_data_obj[, cell_IDs]) %>% sort(decreasing = TRUE)
      df[[paste(parent_category_level, child_category_level, sep = "_")]] <- gene_abundance
      cell_counts[[paste(parent_category_level, child_category_level, sep = "_")]] <- length(cell_IDs)

    }

  }
  anno_df <- obj@meta.data$Category[!(obj@meta.data$Category %>% duplicated()), ]
  anno_df[["FOV"]] <- NULL

  pca_df <- PCA(t(df), graph = FALSE, ncp = 10)
  anno_df[, c("PC1", "PC2")] <- pca_df$ind$cos2[, c("Dim.1", "Dim.2")]

  p <- list()

  p[["eig"]] <- fviz_eig(pca_df)
  p[["PCA_plot"]] <- ggplot(data = anno_df,
                            mapping = aes(x = PC1,
                                          y = PC2,
                                          shape = Slide_name,
                                          fill = Patient_ID,
                                          color = Patient_ID)) +
    geom_point(size = 2, stroke = 1.5) +
    scale_shape_manual(values = c(21, 22))

  return(wrap_plots(p, nrow = 2))

}

# plot_expressed_vs_avg_count <- function(obj) {
#
#
#
# }

filter_seurat_obj <- function(
    obj,
    nGenes_min = 0,
    nGenes_max = Inf,
    nTrans_min = 0,
    mito_max = 5,
    ribo.pct = 0,
    hb.pct = 0,
    cell_exp_min,
    regress_cc = TRUE,
    regress_stress = FALSE,
    regress.mito = TRUE,
    regress.ribo = TRUE,
    regress.Hb = TRUE,
    custom.genes = c(),
    gender.check = FALSE) {

  obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^Mt-|^mt-|^MT-")
  counts_data <- LayerData(object = obj, assay  ="RNA", layer = "counts")
  gene_count <- (1*(counts_data > 0)) %>% colSums()
  feature_count <- counts_data %>% colSums()
  selected_c <- which((gene_count > nGenes_min) &
                        (gene_count < nGenes_max) &
                        (feature_count > nTrans_min))
  selected_f <- rownames(obj)[((1*(counts_data > 0)) %>% rowSums()) > cell_exp_min]
  obj <- subset(obj, features = selected_f, cells = selected_c, subset = percent.mt < mito_max)

}

qc_diagnostics <- function(seuratObj,
                           sample_column = NULL,
                           filter_data = TRUE) {

  if (!is.null(sample_column)) {

    sample_list <- vector("list", seuratObj@meta.data[[sample_column]] %>%
                            unique() %>% length())
    names(sample_list) <- seuratObj@meta.data[[sample_column]] %>% unique()

    for (i in names(sample_list)) {
      selected_c <- which(obj@meta.data[[sample_column]] %in% i)
      sample_list[[i]] <- subset(seuratObj, cells = selected_c)
    }

  }

  else {
    sample_list[[1]] <- seuratObj
  }

}

