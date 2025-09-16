#' @title screenLR
#' @description Screening ligand-receptor pairs based on p-value matrix.
#' @param pairs_table A matrix of p-values with rows as unique ligand-receptor pairs and columns as samples.
#' @param type Screening method, currently only "rank" is supported.
#' @param threshold Proportion of samples in which a ligand-receptor pair must be significant to be retained (default is 0.8).
#' @param p_value Significance threshold for p-values (default is 0.01).
#' @param value Type of values in the input matrix, currently only "pvalue" is supported.
#' @param all_pairs Logical, if TRUE returns all pairs with their ranks, if FALSE returns only pairs passing the threshold.
#' @return A data frame of screened ligand-receptor pairs with their ranks.
#' @examples
#' # Example usage:
#' # screened_pairs <- screenLR(pairs_table, type = "rank", threshold = 0.8, p_value = 0.01, value = "pvalue", all_pairs = FALSE)
#' @export


screenLR <- function(pairs_table, type = "rank", threshold = 0.8, p_value = 0.01,
                     value = "pvalue", all_pairs = FALSE) {
  pairs_table <- as.matrix(pairs_table)
  if (value == "pvalue"){
    message("The input is p-value matrix..")
    dat <- apply(pairs_table, 2, pvalue_to_binary, p_value = p_value)
    row_sum <- rowSums(dat)
    min_sample <- threshold*ncol(pairs_table)
    message(paste0("Total sample size is ", ncol(pairs_table)))
    message(paste0("Mininum sample is ", min_sample))
    top_index <- which(row_sum > min_sample)
    message(paste0("There are ", length(top_index), " unique pairs pass the threshold.."))
  }

  rank_table <- data.frame(unique_pair = rownames(dat), rank = row_sum)
  string <- str_split_fixed(rank_table$unique_pair, pattern = "--", n = 3)
  colnames(string) <- c("sender_celltype", "receiver_celltype", "LR_pairs")
  rank_table <- cbind(rank_table, string)
  rank_table <- rank_table[ ,c("sender_celltype", "receiver_celltype", "LR_pairs", "rank")]
  rownames(rank_table) <- NULL

  if (all_pairs == FALSE) {
    rank_table <- rank_table[top_index, ]
    rank_table <- rank_table[order(-rank_table$rank), ]
    return(rank_table)
  }

  if (all_pairs == TRUE) {
    rank_table <- rank_table[order(-rank_table$rank), ]
    return(rank_table)
  }
}


### Short Functions
pvalue_to_binary <- function(vector, p_value) {
  index1 <- which(vector < p_value)
  index2 <- which(vector >= p_value)
  vector[index1] <- 1
  vector[index2] <- 0
  return(vector)
}






