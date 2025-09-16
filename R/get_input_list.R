#' @title Generating a input list
#' @description Generating a input list from CellChat or Spatalk objects
#' @param sample A vector of sample names
#' @param directory A vector of directories to load the data
#' @param type Type of the input data, either "CellChat" or "Spatalk"
#' @return A list of CellChat or Spatalk objects
#' @export

### Generating a input list
get_input_list <- function(sample, directory, type = "CellChat") {
  if (is.null(directory)) {
    stop("Please provide a directory to load the data...")
  }
  if (length(sample) != length(directory)) {
    stop("Sample names and directory length do not match...")
  }
  if (!type %in% c("CellChat", "Spatalk")) {
    stop("Please provide a valid type: CellChat, Spatalk")
  }
  if (type == "CellChat") {
    message("Loading CellChat data...")
  } else if (type == "Spatalk") {
    message("Loading Spatalk data...")
  }

  input_list <- list()
  for (i in 1:length(sample)) {
    message(paste0("Loading sample: ", sample[i]))
    temp <- readRDS(directory[i])
    input_list[[sample[i]]] <- temp
  }
  return(input_list)
}


