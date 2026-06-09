#' @title loadingLR_cellchat
#' @description Loading CellChat object and extract LR pairs information
#' @param input_list A list of CellChat objects, or NULL if using sample and directory parameters
#' @param sample Optional: character vector of sample names
#' @param directory Optional: character vector of file paths to CellChat RDS files
#' @param value Element to extract, "pval" or "prob"
#' @return A list where each element is a data frame for one sample with columns:
#'   source (sender cell type), target (receiver cell type), ligand, receptor, and value (pval or prob)
#' @importFrom reshape2 melt
#' @importFrom stringr str_split_fixed
#' @export

### Cell Chat Data loading
loadingLR_cellchat <- function(input_list = NULL, sample = NULL, directory = NULL, value = "pval") {

  if (!value %in% c("pval", "prob")) {
    stop("value must be either 'pval' or 'prob'")
  }
  
  # Check if input_list is provided
  if (!is.null(input_list)) {
    # Use provided input_list, ignore sample and directory
    message(paste0("Using provided input_list with ", length(input_list), " samples..."))
    if (length(input_list) < 1) {
      stop("input_list must contain at least 1 sample")
    }
  } else {
    # Load from sample and directory
    if (is.null(sample) || is.null(directory)) {
      stop("Please provide either input_list OR both sample and directory parameters")
    }
    if (length(sample) != length(directory)) {
      stop("Sample names and directory length do not match...")
    }
  }
  
  # Initialize output list
  output_list <- list()
  
  # Process each sample
  if (!is.null(input_list)) {
    # Process from input_list
    message(paste0("Extracting ", value, " from ", length(input_list), " samples..."))
    for (i in seq_along(input_list)) {
      sample_name <- names(input_list)[i]
      if (is.null(sample_name)) sample_name <- paste0("Sample_", i)
      
      message(paste0("Processing sample: ", sample_name))
      cellchat_obj <- input_list[[i]]
      
      # Extract data based on value type
      if (value == "pval") {
        interaction <- reshape2::melt(cellchat_obj@net$pval)
      } else if (value == "prob") {
        interaction <- reshape2::melt(cellchat_obj@net$prob)
      }
      
      # Rename columns
      colnames(interaction) <- c("source", "target", "interaction_name", value)
      
      # Split interaction_name into ligand and receptor
      lr_split <- stringr::str_split_fixed(interaction$interaction_name, "_", 2)
      interaction$ligand <- lr_split[, 1]
      interaction$receptor <- lr_split[, 2]
      
      # Reorder columns
      interaction <- interaction[, c("source", "target", "ligand", "receptor", "interaction_name", value)]
      
      # Store in output list
      output_list[[sample_name]] <- interaction
    }
  } else {
    # Load and process samples one by one
    message("Loading and processing CellChat data...")
    for (i in 1:length(sample)) {
      sample_name <- sample[i]
      message(paste0("Loading and processing sample: ", sample_name))
      
      # Load sample
      cellchat_obj <- readRDS(directory[i])
      
      # Extract data based on value type
      if (value == "pval") {
        interaction <- reshape2::melt(cellchat_obj@net$pval)
      } else if (value == "prob") {
        interaction <- reshape2::melt(cellchat_obj@net$prob)
      }
      
      # Rename columns
      colnames(interaction) <- c("source", "target", "interaction_name", value)
      
      # Split interaction_name into ligand and receptor
      lr_split <- stringr::str_split_fixed(interaction$interaction_name, "_", 2)
      interaction$ligand <- lr_split[, 1]
      interaction$receptor <- lr_split[, 2]
      
      # Reorder columns
      interaction <- interaction[, c("source", "target", "ligand", "receptor", "interaction_name", value)]
      
      # Store in output list
      output_list[[sample_name]] <- interaction
    }
  }
  
  message("Extraction complete!")
  return(output_list)
}


#' @title loadingPWY_cellchat
#' @description Loading CellChat object and extract pathway-level information
#' @param input_list A list of CellChat objects, or NULL if using sample and directory parameters
#' @param sample Optional: character vector of sample names
#' @param directory Optional: character vector of file paths to CellChat RDS files
#' @param value Element to extract, "pval" or "prob"
#' @return A list where each element is a data frame for one sample with columns:
#'   source (sender cell type), target (receiver cell type), pathway, and value (pval or prob)
#' @importFrom reshape2 melt
#' @export

### Cell Chat Pathway Data loading
loadingPWY_cellchat <- function(input_list = NULL, sample = NULL, directory = NULL, value = "pval") {

  if (!value %in% c("pval", "prob")) {
    stop("value must be either 'pval' or 'prob'")
  }
  
  # Check if input_list is provided
  if (!is.null(input_list)) {
    # Use provided input_list, ignore sample and directory
    message(paste0("Using provided input_list with ", length(input_list), " samples..."))
    if (length(input_list) < 1) {
      stop("input_list must contain at least 1 sample")
    }
  } else {
    # Load from sample and directory
    if (is.null(sample) || is.null(directory)) {
      stop("Please provide either input_list OR both sample and directory parameters")
    }
    if (length(sample) != length(directory)) {
      stop("Sample names and directory length do not match...")
    }
  }
  
  # Initialize output list
  output_list <- list()
  
  # Process each sample
  if (!is.null(input_list)) {
    # Process from input_list
    message(paste0("Extracting pathway ", value, " from ", length(input_list), " samples..."))
    for (i in seq_along(input_list)) {
      sample_name <- names(input_list)[i]
      if (is.null(sample_name)) sample_name <- paste0("Sample_", i)
      
      message(paste0("Processing sample: ", sample_name))
      cellchat_obj <- input_list[[i]]
      
      # Check if pathway-level information exists
      if (is.null(cellchat_obj@netP) || length(cellchat_obj@netP) == 0) {
        warning(paste0("Sample ", sample_name, " does not have pathway-level information. Run computeCommunProbPathway() first."))
        next
      }
      
      # Extract data based on value type
      if (value == "pval") {
        if (is.null(cellchat_obj@netP$pval)) {
          warning(paste0("Sample ", sample_name, " does not have pathway pval information."))
          next
        }
        pathway_data <- reshape2::melt(cellchat_obj@netP$pval)
      } else if (value == "prob") {
        if (is.null(cellchat_obj@netP$prob)) {
          warning(paste0("Sample ", sample_name, " does not have pathway prob information."))
          next
        }
        pathway_data <- reshape2::melt(cellchat_obj@netP$prob)
      }
      
      # Rename columns
      colnames(pathway_data) <- c("source", "target", "pathway", value)
      
      # Store in output list
      output_list[[sample_name]] <- pathway_data
    }
  } else {
    # Load and process samples one by one
    message("Loading and processing pathway data...")
    for (i in 1:length(sample)) {
      sample_name <- sample[i]
      message(paste0("Loading and processing sample: ", sample_name))
      
      # Load sample
      cellchat_obj <- readRDS(directory[i])
      
      # Check if pathway-level information exists
      if (is.null(cellchat_obj@netP) || length(cellchat_obj@netP) == 0) {
        warning(paste0("Sample ", sample_name, " does not have pathway-level information. Run computeCommunProbPathway() first."))
        next
      }
      
      # Extract data based on value type
      if (value == "pval") {
        if (is.null(cellchat_obj@netP$pval)) {
          warning(paste0("Sample ", sample_name, " does not have pathway pval information."))
          next
        }
        pathway_data <- reshape2::melt(cellchat_obj@netP$pval)
      } else if (value == "prob") {
        if (is.null(cellchat_obj@netP$prob)) {
          warning(paste0("Sample ", sample_name, " does not have pathway prob information."))
          next
        }
        pathway_data <- reshape2::melt(cellchat_obj@netP$prob)
      }
      
      # Rename columns
      colnames(pathway_data) <- c("source", "target", "pathway", value)
      
      # Store in output list
      output_list[[sample_name]] <- pathway_data
    }
  }
  
  message("Pathway extraction complete!")
  return(output_list)
}






