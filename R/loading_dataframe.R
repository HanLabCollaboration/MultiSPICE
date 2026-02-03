#' @title loadingLR_dataframe
#' @description Load ligand-receptor pair data from user-provided dataframes
#' @param data A dataframe or list of dataframes containing LR interaction data
#' @param sample_col Name of column containing sample IDs (required if data is a single dataframe)
#' @param source_col Name of column containing source/sender cell type (default: "source")
#' @param target_col Name of column containing target/receiver cell type (default: "target")
#' @param ligand_col Name of column containing ligand gene (default: "ligand")
#' @param receptor_col Name of column containing receptor gene (default: "receptor")
#' @param interaction_col Name of column containing interaction name (e.g., "LIGAND_RECEPTOR").
#'   If NULL, will be created from ligand and receptor columns (default: NULL)
#' @param value_col Name of column containing significance values (default: "pval")
#' @param value_type Type of value in value_col: "pval" or "prob" (default: "pval")
#' @return A named list where each element is a dataframe for one sample with columns:
#'   source, target, ligand, receptor, interaction_name, and pval/prob
#' @examples
#' \dontrun{
#' # Example 1: Single dataframe with sample column
#' df <- data.frame(
#'   sample = c("S1", "S1", "S2", "S2"),
#'   source = c("T cell", "B cell", "T cell", "B cell"),
#'   target = c("Tumor", "Tumor", "Tumor", "Tumor"),
#'   ligand = c("IFNG", "TNF", "IFNG", "TNF"),
#'   receptor = c("IFNGR1", "TNFRSF1A", "IFNGR1", "TNFRSF1A"),
#'   pval = c(0.001, 0.01, 0.002, 0.015)
#' )
#' lr_list <- loadingLR_dataframe(df, sample_col = "sample")
#' 
#' # Example 2: List of dataframes (one per sample)
#' df_list <- list(
#'   Sample1 = data.frame(
#'     source = c("T cell", "B cell"),
#'     target = c("Tumor", "Tumor"),
#'     ligand = c("IFNG", "TNF"),
#'     receptor = c("IFNGR1", "TNFRSF1A"),
#'     pval = c(0.001, 0.01)
#'   ),
#'   Sample2 = data.frame(
#'     source = c("T cell", "B cell"),
#'     target = c("Tumor", "Tumor"),
#'     ligand = c("IFNG", "TNF"),
#'     receptor = c("IFNGR1", "TNFRSF1A"),
#'     pval = c(0.002, 0.015)
#'   )
#' )
#' lr_list <- loadingLR_dataframe(df_list)
#' }
#' @export

loadingLR_dataframe <- function(data,
                                sample_col = NULL,
                                source_col = "source",
                                target_col = "target",
                                ligand_col = "ligand",
                                receptor_col = "receptor",
                                interaction_col = NULL,
                                value_col = "pval",
                                value_type = "pval") {
  
  # Validate value_type
  if (!value_type %in% c("pval", "prob")) {
    stop("value_type must be either 'pval' or 'prob'")
  }
  
  # Check if data is a list or dataframe
  if (is.data.frame(data)) {
    # Single dataframe - must have sample column
    if (is.null(sample_col) || !sample_col %in% colnames(data)) {
      stop("For single dataframe input, sample_col must be specified and exist in the dataframe")
    }
    
    # Validate required columns
    required_cols <- c(sample_col, source_col, target_col, ligand_col, receptor_col, value_col)
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    
    # Split by sample
    samples <- unique(data[[sample_col]])
    message(paste0("Processing ", length(samples), " samples from dataframe..."))
    
    output_list <- list()
    for (s in samples) {
      sample_data <- data[data[[sample_col]] == s, ]
      
      # Create interaction name if not provided
      if (is.null(interaction_col) || !interaction_col %in% colnames(sample_data)) {
        sample_data$interaction_name <- paste0(sample_data[[ligand_col]], "_", 
                                               sample_data[[receptor_col]])
      } else {
        sample_data$interaction_name <- sample_data[[interaction_col]]
      }
      
      # Create standardized output
      output_df <- data.frame(
        source = sample_data[[source_col]],
        target = sample_data[[target_col]],
        ligand = sample_data[[ligand_col]],
        receptor = sample_data[[receptor_col]],
        interaction_name = sample_data$interaction_name,
        value = sample_data[[value_col]],
        stringsAsFactors = FALSE
      )
      
      # Rename value column based on value_type
      colnames(output_df)[colnames(output_df) == "value"] <- value_type
      
      output_list[[as.character(s)]] <- output_df
      message(paste0("  Sample ", s, ": ", nrow(output_df), " interactions"))
    }
    
  } else if (is.list(data)) {
    # List of dataframes
    message(paste0("Processing ", length(data), " samples from list..."))
    
    output_list <- list()
    for (i in seq_along(data)) {
      sample_name <- names(data)[i]
      if (is.null(sample_name) || sample_name == "") {
        sample_name <- paste0("Sample_", i)
      }
      
      sample_data <- data[[i]]
      
      if (!is.data.frame(sample_data)) {
        stop(paste("Element", i, "is not a dataframe"))
      }
      
      # Validate required columns
      required_cols <- c(source_col, target_col, ligand_col, receptor_col, value_col)
      missing_cols <- setdiff(required_cols, colnames(sample_data))
      if (length(missing_cols) > 0) {
        stop(paste("Sample", sample_name, "missing required columns:", 
                   paste(missing_cols, collapse = ", ")))
      }
      
      # Create interaction name if not provided
      if (is.null(interaction_col) || !interaction_col %in% colnames(sample_data)) {
        sample_data$interaction_name <- paste0(sample_data[[ligand_col]], "_", 
                                               sample_data[[receptor_col]])
      } else {
        sample_data$interaction_name <- sample_data[[interaction_col]]
      }
      
      # Create standardized output
      output_df <- data.frame(
        source = sample_data[[source_col]],
        target = sample_data[[target_col]],
        ligand = sample_data[[ligand_col]],
        receptor = sample_data[[receptor_col]],
        interaction_name = sample_data$interaction_name,
        value = sample_data[[value_col]],
        stringsAsFactors = FALSE
      )
      
      # Rename value column based on value_type
      colnames(output_df)[colnames(output_df) == "value"] <- value_type
      
      output_list[[sample_name]] <- output_df
      message(paste0("  Sample ", sample_name, ": ", nrow(output_df), " interactions"))
    }
    
  } else {
    stop("data must be either a dataframe or a list of dataframes")
  }
  
  message("Data loading complete!")
  return(output_list)
}


#' @title loadingPWY_dataframe
#' @description Load pathway-level interaction data from user-provided dataframes
#' @param data A dataframe or list of dataframes containing pathway interaction data
#' @param sample_col Name of column containing sample IDs (required if data is a single dataframe)
#' @param source_col Name of column containing source/sender cell type (default: "source")
#' @param target_col Name of column containing target/receiver cell type (default: "target")
#' @param pathway_col Name of column containing pathway name (default: "pathway")
#' @param value_col Name of column containing significance values (default: "pval")
#' @param value_type Type of value in value_col: "pval" or "prob" (default: "pval")
#' @return A named list where each element is a dataframe for one sample with columns:
#'   source, target, pathway, and pval/prob
#' @examples
#' \dontrun{
#' # Example 1: Single dataframe with sample column
#' df <- data.frame(
#'   sample = c("S1", "S1", "S2", "S2"),
#'   source = c("T cell", "B cell", "T cell", "B cell"),
#'   target = c("Tumor", "Tumor", "Tumor", "Tumor"),
#'   pathway = c("IFN-II", "TNF", "IFN-II", "TNF"),
#'   pval = c(0.001, 0.01, 0.002, 0.015)
#' )
#' pwy_list <- loadingPWY_dataframe(df, sample_col = "sample")
#' 
#' # Example 2: List of dataframes (one per sample)
#' df_list <- list(
#'   Sample1 = data.frame(
#'     source = c("T cell", "B cell"),
#'     target = c("Tumor", "Tumor"),
#'     pathway = c("IFN-II", "TNF"),
#'     pval = c(0.001, 0.01)
#'   ),
#'   Sample2 = data.frame(
#'     source = c("T cell", "B cell"),
#'     target = c("Tumor", "Tumor"),
#'     pathway = c("IFN-II", "TNF"),
#'     pval = c(0.002, 0.015)
#'   )
#' )
#' pwy_list <- loadingPWY_dataframe(df_list)
#' }
#' @export

loadingPWY_dataframe <- function(data,
                                 sample_col = NULL,
                                 source_col = "source",
                                 target_col = "target",
                                 pathway_col = "pathway",
                                 value_col = "pval",
                                 value_type = "pval") {
  
  # Validate value_type
  if (!value_type %in% c("pval", "prob")) {
    stop("value_type must be either 'pval' or 'prob'")
  }
  
  # Check if data is a list or dataframe
  if (is.data.frame(data)) {
    # Single dataframe - must have sample column
    if (is.null(sample_col) || !sample_col %in% colnames(data)) {
      stop("For single dataframe input, sample_col must be specified and exist in the dataframe")
    }
    
    # Validate required columns
    required_cols <- c(sample_col, source_col, target_col, pathway_col, value_col)
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    
    # Split by sample
    samples <- unique(data[[sample_col]])
    message(paste0("Processing ", length(samples), " samples from dataframe..."))
    
    output_list <- list()
    for (s in samples) {
      sample_data <- data[data[[sample_col]] == s, ]
      
      # Create standardized output
      output_df <- data.frame(
        source = sample_data[[source_col]],
        target = sample_data[[target_col]],
        pathway = sample_data[[pathway_col]],
        value = sample_data[[value_col]],
        stringsAsFactors = FALSE
      )
      
      # Rename value column based on value_type
      colnames(output_df)[colnames(output_df) == "value"] <- value_type
      
      output_list[[as.character(s)]] <- output_df
      message(paste0("  Sample ", s, ": ", nrow(output_df), " pathway interactions"))
    }
    
  } else if (is.list(data)) {
    # List of dataframes
    message(paste0("Processing ", length(data), " samples from list..."))
    
    output_list <- list()
    for (i in seq_along(data)) {
      sample_name <- names(data)[i]
      if (is.null(sample_name) || sample_name == "") {
        sample_name <- paste0("Sample_", i)
      }
      
      sample_data <- data[[i]]
      
      if (!is.data.frame(sample_data)) {
        stop(paste("Element", i, "is not a dataframe"))
      }
      
      # Validate required columns
      required_cols <- c(source_col, target_col, pathway_col, value_col)
      missing_cols <- setdiff(required_cols, colnames(sample_data))
      if (length(missing_cols) > 0) {
        stop(paste("Sample", sample_name, "missing required columns:", 
                   paste(missing_cols, collapse = ", ")))
      }
      
      # Create standardized output
      output_df <- data.frame(
        source = sample_data[[source_col]],
        target = sample_data[[target_col]],
        pathway = sample_data[[pathway_col]],
        value = sample_data[[value_col]],
        stringsAsFactors = FALSE
      )
      
      # Rename value column based on value_type
      colnames(output_df)[colnames(output_df) == "value"] <- value_type
      
      output_list[[sample_name]] <- output_df
      message(paste0("  Sample ", sample_name, ": ", nrow(output_df), " pathway interactions"))
    }
    
  } else {
    stop("data must be either a dataframe or a list of dataframes")
  }
  
  message("Pathway data loading complete!")
  return(output_list)
}
