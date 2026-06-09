#' Create LR Feature Matrix
#'
#' Creates a feature matrix where rows represent ligand-receptor pairs and columns 
#' represent samples. Each cell contains the count of source-target cell type
#' pairs with significant interactions for that LR pair in that sample.
#' Output is formatted for downstream analysis with limma.
#'
#' @param lr_list A named list of data frames from loadingLR_cellchat, where each element 
#'   is a sample with columns: source, target, interaction_name, pval/prob
#' @param rank_table A data frame containing ranked LR pairs (output from screenLR)
#' @param cutoff Significance cutoff threshold (default: 0.05)
#' @param value Type of values in lr_list, "pval" or "prob" (default: "pval")
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A matrix with LR pairs as rows and samples as columns, containing 
#'         counts of source-target cell type pairs with significant interactions
#'
#' @examples
#' \dontrun{
#' # Load LR data from CellChat objects
#' LR_pairs_list <- loadingLR_cellchat(sample = sample, directory = directory, value = "pval")
#' 
#' # Screen for top LR pairs
#' rank_table <- screenLR(lr_list = LR_pairs_list, type = "rank", threshold = 0.8,
#'                        cutoff = 0.01, value = "pval", all_pairs = FALSE)
#' 
#' # Create LR feature matrix
#' lr_matrix <- LR_matrix(lr_list = LR_pairs_list, rank_table = rank_table, 
#'                        cutoff = 0.05, value = "pval")
#' }
#'
#' @export

LR_matrix <- function(lr_list, 
                      rank_table, 
                      cutoff = 0.05,
                      value = "pval",
                      verbose = TRUE) {
  
  # Validate inputs
  if (!is.list(lr_list) || length(lr_list) == 0) {
    stop("lr_list must be a non-empty list from loadingLR_cellchat")
  }
  
  if (!is.data.frame(rank_table) || !"interaction_name" %in% colnames(rank_table)) {
    stop("rank_table must be a data frame with an 'interaction_name' column")
  }
  
  if (!value %in% c("pval", "prob")) {
    stop("value must be either 'pval' or 'prob'")
  }
  
  if (!is.numeric(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("cutoff must be a numeric value between 0 and 1")
  }
  
  # Extract unique LR pairs and sample names
  candidate_lr <- unique(rank_table$interaction_name)
  sample_names <- names(lr_list)
  
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_along(lr_list))
    names(lr_list) <- sample_names
  }
  
  if (verbose) {
    cat(sprintf("Creating LR feature matrix for %d samples and %d LR pairs...\n", 
                length(sample_names), length(candidate_lr)))
  }
  
  # Initialize matrix
  mat <- matrix(0, nrow = length(sample_names), ncol = length(candidate_lr))
  rownames(mat) <- sample_names
  colnames(mat) <- candidate_lr
  
  # Fill matrix with interaction counts
  for (i in seq_along(sample_names)) {
    s <- sample_names[i]
    sample_data <- lr_list[[s]]
    
    # Filter for significant interactions based on value type
    if (value == "pval") {
      sample_data <- sample_data[sample_data$pval < cutoff, ]
    } else {
      sample_data <- sample_data[sample_data$prob > cutoff, ]
    }
    
    # For each LR pair, count unique source-target combinations
    for (lr in candidate_lr) {
      lr_data <- sample_data[sample_data$interaction_name == lr, ]
      
      if (nrow(lr_data) > 0) {
        # Count unique source-target pairs
        unique_pairs <- unique(paste0(lr_data$source, "--", lr_data$target))
        mat[s, lr] <- length(unique_pairs)
      }
    }
    
    if (verbose && i %% 10 == 0) {
      cat(sprintf("Processed %d/%d samples...\n", i, length(sample_names)))
    }
  }
  
  if (verbose) {
    cat("LR feature matrix creation complete!\n")
    cat(sprintf("Matrix dimensions: %d LR pairs x %d samples\n", 
                ncol(mat), nrow(mat)))
  }
  
  # Transpose so rows are LR pairs and columns are samples (standard for limma)
  return(t(mat))
}


