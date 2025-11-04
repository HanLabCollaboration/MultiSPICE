#' Create LR Feature Matrix
#'
#' Creates a feature matrix where rows represent samples and columns represent 
#' ligand-receptor pairs. Each cell contains the count of significant interactions
#' for that LR pair in that sample.
#'
#' @param input_list A named list of CellChat objects, where names are sample IDs
#' @param rank_table A data frame containing ranked LR pairs (output from screenLR)
#' @param p_threshold P-value threshold for significant interactions (default: 0.05)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A matrix with samples as rows and LR pairs as columns, containing 
#'         interaction counts
#'
#' @examples
#' \dontrun{
#' # After loading CellChat objects
#' input_list <- get_input_list(sample = sample, directory = directory, type = "CellChat")
#' LR_pairs <- loadingLR_cellchat(input_list = input_list, element = "pval")
#' rank_table <- screenLR(pairs_table = LR_pairs, type = "rank", threshold = 0.8,
#'                        p_value = 0.01, value = "pvalue", all_pairs = FALSE)
#' 
#' # Create LR feature matrix
#' lr_matrix <- create_LR_matrix(input_list, rank_table, p_threshold = 0.05)
#' }
#'
#' @export
create_LR_matrix <- function(input_list, 
                              rank_table, 
                              p_threshold = 0.05, 
                              verbose = TRUE) {
  
  # Validate inputs
  if (!is.list(input_list) || length(input_list) == 0) {
    stop("input_list must be a non-empty list of CellChat objects")
  }
  
  if (!is.data.frame(rank_table) || !"LR_pairs" %in% colnames(rank_table)) {
    stop("rank_table must be a data frame with an 'LR_pairs' column")
  }
  
  if (!is.numeric(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("p_threshold must be a numeric value between 0 and 1")
  }
  
  # Extract unique LR pairs and sample names
  candidate_lr <- unique(rank_table$LR_pairs)
  sample_names <- names(input_list)
  
  if (verbose) {
    cat(sprintf("Creating LR feature matrix for %d samples and %d LR pairs...\n", 
                length(sample_names), length(candidate_lr)))
  }
  
  # Initialize matrix
  mat <- matrix(0, nrow = length(sample_names), ncol = length(candidate_lr))
  rownames(mat) <- sample_names
  colnames(mat) <- candidate_lr
  
  # Fill matrix with interaction counts
  for (lr in candidate_lr) {
    counts <- c()
    
    for (s in sample_names) {
      # Get CellChat object for current sample
      obj <- input_list[[s]]
      
      # Extract p-values
      pval <- get_pvalCC(obj)
      
      # Filter for current LR pair
      pval <- pval[pval$lr_pairs %in% lr, ]
      
      # Filter for significant interactions
      pval <- pval[pval$pval < p_threshold, ]
      
      # Count significant interactions
      counts <- c(counts, nrow(pval))
    }
    
    # Assign counts to matrix
    mat[, lr] <- counts
    
    if (verbose) {
      cat(sprintf("Processing %s... done\n", lr))
    }
  }
  
  if (verbose) {
    cat("LR feature matrix creation complete!\n")
    cat(sprintf("Matrix dimensions: %d samples x %d LR pairs\n", 
                nrow(mat), ncol(mat)))
  }
  
  return(mat)
}


#' Get Cell Type Pairs from Rank Table
#'
#' Extracts unique sender-receiver cell type pairs from a rank table
#'
#' @param rank_table A data frame containing ranked LR pairs with sender_celltype 
#'                   and receiver_celltype columns
#'
#' @return A character vector of cell type pairs in format "sender--receiver"
#'
#' @export
get_celltype_pairs <- function(rank_table) {
  
  if (!all(c("sender_celltype", "receiver_celltype") %in% colnames(rank_table))) {
    stop("rank_table must contain 'sender_celltype' and 'receiver_celltype' columns")
  }
  
  celltype_pairs <- unique(paste0(rank_table$sender_celltype, "--", 
                                   rank_table$receiver_celltype))
  
  return(celltype_pairs)
}


