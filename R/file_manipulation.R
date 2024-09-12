cd.. <- function(dir, base = FALSE) {

  if (base == FALSE) {

    new_dir <- str_remove(string = dir, pattern = paste0("/", basename(dir)))

  } else {

    # Need to find right regex
    new_dir <- str_remove(string = dir, pattern = "")

  }

  new_dir

}
