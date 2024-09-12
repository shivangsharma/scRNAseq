df_transpose <- function(df){
  col_names <- rownames(df)
  row_names <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  df <- as.data.frame(df, row.names = row_names)
  colnames(df) <- col_names
  return(df)
}

not <- function(vec) {

  vec <- as.logical(vec)
  new_vec <- c()
  for (i in c(1:length(vec))) {

    new_vec <- c(new_vec, !vec[i])

  }
  return(new_vec)

}

# convert_to_10X <- function(filepath,
#                          unzip = FALSE,
#                          writeto10Xformat = FALSE,
#                          createseuratobject = FALSE) {
#
#   # dir <- sub("/[^/]+$", "", filepath)
#   dir <- dirname(filepath)
#   orig_threads <- setDTthreads(threads = 40)
#   on.exit(setDTthreads(orig_threads))
#
#   if (unzip == TRUE) {
#     filepath <- gunzip(filepath,
#                        remove = FALSE,
#                        overwrite = TRUE)
#   }
#
#   data_format <- file_ext(filepath)
#
#   if (data_format == "txt") {
#
#     data <- fread(file = filepath)
#
#     if ((data[[1]] %>% typeof) == "character") {
#       gene_names <- data[[1]]
#       cell_ids <- colnames(data)
#       data[[1]] <- NULL
#
#     }
#     else {
#       stop("Gene names not found")
#     }
#
#     data <- sparsify(data)
#     data@Dimnames[[1]] <- gene_names
#     data@Dimnames[[2]] <- cell_ids
#
#   }
#   else if (data_format == "h5") {
#
#     data <- h5read(file = filepath, name = "/matrix")
#     gene_names <- data[["features"]][["name"]]
#     cell_ids <- data[["barcodes"]]
#     data <- sparseMatrix(i = data$indices[] + 1, p = data$indptr[], x = as.numeric(data$data), dims = data$shape, repr = "C")
#     rownames(data) <- gene_names
#     colnames(data) <- cell_ids
#     h5closeAll()
#
#   }
#
#
#   if(writeto10Xformat == TRUE) {
#     write10xCounts(path = file.path(dir, "10X format"),
#                    x = data,
#                    barcodes = data@Dimnames[[2]],
#                    gene.symbol = data@Dimnames[[1]],
#                    version = "3",
#                    overwrite = TRUE,
#                    )
#   }
#
#   return(data)
#
# }
