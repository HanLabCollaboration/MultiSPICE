#' @title screenLR
#' @description Screening ligand-receptor pairs based on significance across multiple samples.
#' @param lr_list A list of data frames from loadingLR_cellchat, where each element is a sample
#'   with columns: source, target, ligand, receptor, interaction_name, and pval/prob
#' @param type Screening method, currently only "rank" is supported.
#' @param threshold Proportion of samples in which a ligand-receptor pair must be significant 
#'   to be retained (default is 0.8).
#' @param cutoff Significance cutoff threshold (default is 0.01). For pval: values < cutoff are significant. For prob: values > cutoff are significant.
#' @param value Type of values in the input, "pval" or "prob"
#' @param all_pairs Logical, if TRUE returns all pairs with their ranks, if FALSE returns only 
#'   pairs passing the threshold.
#' @return A data frame of screened ligand-receptor pairs with columns:
#'   source, target, ligand, receptor, interaction_name, rank
#' @importFrom stringr str_split_fixed
#' @examples
#' \dontrun{
#' # Load LR data from CellChat objects
#' lr_data <- loadingLR_cellchat(sample = samples, directory = dirs, value = "pval")
#' # Screen for consistent LR pairs
#' screened <- screenLR(lr_data, threshold = 0.8, cutoff = 0.01)
#' }
#' @export

screenLR <- function(lr_list, type = "rank", threshold = 0.8, cutoff = 0.01,
                     value = "pval", all_pairs = FALSE) {
  
  if (!is.list(lr_list) || length(lr_list) < 1) {
    stop("lr_list must be a non-empty list from loadingLR_cellchat")
  }
  
  if (!value %in% c("pval", "prob")) {
    stop("value must be either 'pval' or 'prob'")
  }
  
  message(paste0("Screening LR pairs across ", length(lr_list), " samples..."))
  
  # Combine all samples into one data frame with sample ID
  combined_data <- list()
  for (i in seq_along(lr_list)) {
    sample_name <- names(lr_list)[i]
    if (is.null(sample_name)) sample_name <- paste0("Sample_", i)
    
    temp <- lr_list[[i]]
    temp$sample <- sample_name
    combined_data[[i]] <- temp
  }
  combined_data <- do.call(rbind, combined_data)
  
  # Create unique identifier for each interaction
  combined_data$unique_id <- paste0(
    combined_data$source, "--",
    combined_data$target, "--",
    combined_data$interaction_name
  )
  
  # Count significant interactions per unique_id
  if (value == "pval") {
    message("Screening based on p-value threshold...")
    combined_data$is_significant <- combined_data$pval < cutoff
  } else {
    message("Screening based on probability threshold...")
    combined_data$is_significant <- combined_data$prob > cutoff
  }
  
  # Calculate rank (number of samples with significant interaction)
  rank_summary <- aggregate(
    is_significant ~ unique_id + source + target + ligand + receptor + interaction_name,
    data = combined_data,
    FUN = sum
  )
  colnames(rank_summary)[colnames(rank_summary) == "is_significant"] <- "rank"
  
  # Calculate minimum samples threshold
  min_sample <- threshold * length(lr_list)
  message(paste0("Total sample size: ", length(lr_list)))
  message(paste0("Minimum samples required: ", min_sample))
  
  # Filter based on threshold
  n_pass <- sum(rank_summary$rank > min_sample)
  message(paste0("Number of unique pairs passing threshold: ", n_pass))
  
  # Sort by rank
  rank_summary <- rank_summary[order(-rank_summary$rank), ]
  rownames(rank_summary) <- NULL
  
  # Remove unique_id column
  rank_summary$unique_id <- NULL
  
  # Return based on all_pairs parameter
  if (all_pairs == FALSE) {
    rank_summary <- rank_summary[rank_summary$rank > min_sample, ]
    return(rank_summary)
  } else {
    return(rank_summary)
  }
}


#' @title screenPWY
#' @description Screening pathway-level communication based on significance across multiple samples.
#' @param pwy_list A list of data frames from loadingPWY_cellchat, where each element is a sample
#'   with columns: source, target, pathway, and pval/prob
#' @param type Screening method, currently only "rank" is supported.
#' @param threshold Proportion of samples in which a pathway must be significant 
#'   to be retained (default is 0.8).
#' @param cutoff Significance cutoff threshold (default is 0.01). For pval: values < cutoff are significant. For prob: values > cutoff are significant.
#' @param value Type of values in the input, "pval" or "prob"
#' @param all_pairs Logical, if TRUE returns all pathways with their ranks, if FALSE returns only 
#'   pathways passing the threshold.
#' @return A data frame of screened pathways with columns:
#'   source, target, pathway, rank
#' @examples
#' \dontrun{
#' # Load pathway data from CellChat objects
#' pwy_data <- loadingPWY_cellchat(sample = samples, directory = dirs, value = "pval")
#' # Screen for consistent pathways
#' screened <- screenPWY(pwy_data, threshold = 0.8, cutoff = 0.01)
#' }
#' @export

screenPWY <- function(pwy_list, type = "rank", threshold = 0.8, cutoff = 0.01,
                      value = "pval", all_pairs = FALSE) {
  
  if (!is.list(pwy_list) || length(pwy_list) < 1) {
    stop("pwy_list must be a non-empty list from loadingPWY_cellchat")
  }
  
  if (!value %in% c("pval", "prob")) {
    stop("value must be either 'pval' or 'prob'")
  }
  
  message(paste0("Screening pathways across ", length(pwy_list), " samples..."))
  
  # Combine all samples into one data frame with sample ID
  combined_data <- list()
  for (i in seq_along(pwy_list)) {
    sample_name <- names(pwy_list)[i]
    if (is.null(sample_name)) sample_name <- paste0("Sample_", i)
    
    temp <- pwy_list[[i]]
    temp$sample <- sample_name
    combined_data[[i]] <- temp
  }
  combined_data <- do.call(rbind, combined_data)
  
  # Create unique identifier for each pathway interaction
  combined_data$unique_id <- paste0(
    combined_data$source, "--",
    combined_data$target, "--",
    combined_data$pathway
  )
  
  # Count significant interactions per unique_id
  if (value == "pval") {
    message("Screening based on p-value threshold...")
    combined_data$is_significant <- combined_data$pval < cutoff
  } else {
    message("Screening based on probability threshold...")
    combined_data$is_significant <- combined_data$prob > cutoff
  }
  
  # Calculate rank (number of samples with significant pathway)
  rank_summary <- aggregate(
    is_significant ~ unique_id + source + target + pathway,
    data = combined_data,
    FUN = sum
  )
  colnames(rank_summary)[colnames(rank_summary) == "is_significant"] <- "rank"
  
  # Calculate minimum samples threshold
  min_sample <- threshold * length(pwy_list)
  message(paste0("Total sample size: ", length(pwy_list)))
  message(paste0("Minimum samples required: ", min_sample))
  
  # Filter based on threshold
  n_pass <- sum(rank_summary$rank > min_sample)
  message(paste0("Number of unique pathways passing threshold: ", n_pass))
  
  # Sort by rank
  rank_summary <- rank_summary[order(-rank_summary$rank), ]
  rownames(rank_summary) <- NULL
  
  # Remove unique_id column
  rank_summary$unique_id <- NULL
  
  # Return based on all_pairs parameter
  if (all_pairs == FALSE) {
    rank_summary <- rank_summary[rank_summary$rank > min_sample, ]
    return(rank_summary)
  } else {
    return(rank_summary)
  }
}


