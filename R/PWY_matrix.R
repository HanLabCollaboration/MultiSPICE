#' Create Pathway Feature Matrix
#'
#' Creates a feature matrix where rows represent pathways and columns represent samples.
#' Each cell contains the count of source-target cell type pairs with significant 
#' pathway interactions in that sample. Output is formatted for downstream analysis with limma.
#'
#' @param pwy_list A named list of data frames from loadingPWY_cellchat, where each element 
#'   is a sample with columns: source, target, pathway, pval/prob
#' @param rank_table A data frame containing ranked pathways (output from screenPWY)
#' @param cutoff Significance cutoff threshold (default: 0.05)
#' @param value Type of values in pwy_list, "pval" or "prob" (default: "prob")
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A matrix with pathways as rows and samples as columns, containing 
#'         counts of source-target cell type pairs with significant pathway interactions
#'
#' @examples
#' \dontrun{
#' # Load pathway data from CellChat objects
#' PWY_list <- loadingPWY_cellchat(sample = sample, directory = directory, value = "prob")
#' 
#' # Screen for top pathways
#' pwy_rank_table <- screenPWY(pwy_list = PWY_list, type = "rank", threshold = 0.8,
#'                             cutoff = 0.05, value = "prob", all_pairs = FALSE)
#' 
#' # Create pathway feature matrix
#' pwy_matrix <- PWY_matrix(pwy_list = PWY_list, rank_table = pwy_rank_table, 
#'                          cutoff = 0.05, value = "prob")
#' }
#'
#' @export

PWY_matrix <- function(pwy_list, 
                       rank_table, 
                       cutoff = 0.05,
                       value = "prob",
                       verbose = TRUE) {
  
  # Validate inputs
  if (!is.list(pwy_list) || length(pwy_list) == 0) {
    stop("pwy_list must be a non-empty list from loadingPWY_cellchat")
  }
  
  if (!is.data.frame(rank_table) || !"pathway" %in% colnames(rank_table)) {
    stop("rank_table must be a data frame with a 'pathway' column")
  }
  
  if (!value %in% c("pval", "prob")) {
    stop("value must be either 'pval' or 'prob'")
  }
  
  if (!is.numeric(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("cutoff must be a numeric value between 0 and 1")
  }
  
  # Extract unique pathways and sample names
  candidate_pwy <- unique(rank_table$pathway)
  sample_names <- names(pwy_list)
  
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_along(pwy_list))
    names(pwy_list) <- sample_names
  }
  
  if (verbose) {
    cat(sprintf("Creating pathway feature matrix for %d samples and %d pathways...\n", 
                length(sample_names), length(candidate_pwy)))
  }
  
  # Initialize matrix
  mat <- matrix(0, nrow = length(sample_names), ncol = length(candidate_pwy))
  rownames(mat) <- sample_names
  colnames(mat) <- candidate_pwy
  
  # Fill matrix with pathway counts
  for (i in seq_along(sample_names)) {
    s <- sample_names[i]
    sample_data <- pwy_list[[s]]
    
    # Filter for significant pathway interactions based on value type
    if (value == "pval") {
      sample_data <- sample_data[sample_data$pval < cutoff, ]
    } else {
      sample_data <- sample_data[sample_data$prob > cutoff, ]
    }
    
    # For each pathway, count unique source-target combinations
    for (pwy in candidate_pwy) {
      pwy_data <- sample_data[sample_data$pathway == pwy, ]
      
      if (nrow(pwy_data) > 0) {
        # Count unique source-target pairs
        unique_pairs <- unique(paste0(pwy_data$source, "--", pwy_data$target))
        mat[s, pwy] <- length(unique_pairs)
      }
    }
    
    if (verbose && i %% 10 == 0) {
      cat(sprintf("Processed %d/%d samples...\n", i, length(sample_names)))
    }
  }
  
  if (verbose) {
    cat("Pathway feature matrix creation complete!\n")
    cat(sprintf("Matrix dimensions: %d pathways x %d samples\n", 
                ncol(mat), nrow(mat)))
  }
  
  # Transpose so rows are pathways and columns are samples (standard for limma)
  return(t(mat))
}
