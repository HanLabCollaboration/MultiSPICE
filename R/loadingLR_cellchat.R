#' @title loadingLR_cellchat
#' @description Loading CellChat object and extract LR pairs information
#' @param input_list A list of CellChat objects
#' @param element Element to extract, "pval" or "prob"
#' @return A data frame of LR pairs with samples as columns and LR pairs as rows
#' @importFrom reshape2 melt
#' @export

### Cell Chat Data loading
loadingLR_cellchat <- function(input_list, element = "pval") {

  if (length(input_list) < 2) {
    stop("Need at least 2 samples for muti-sample analysis...")
  }

  message("Gathering LR pair names...")
  all_pairs <- lapply(input_list, get_pairsCC)
  all_pairs <- unique(unlist(all_pairs))

  message("Making unique LR pair names...")
  pairs_table <- lapply(input_list, get_pvalCC)
  pairs_table <- lapply(pairs_table, make_unique_pairsCC)
  pairs_table <- lapply(pairs_table, order_uniqur_pairsCC, all_pairs = all_pairs)

  LR_pairs <- data.frame(pairs = all_pairs)
  message("Extract LR Pair siginificance from muti-sample...")
  for (n in 1:length(input_list)) {
    message(paste0("Extract LR-Pairs from sample ", n))
    temp <- pairs_table[[n]]
    LR_pairs <- cbind(LR_pairs, temp$pval)
  }

  LR_pairs <- data.frame(LR_pairs, row.names = 1)
  colnames(LR_pairs) <- names(input_list)

  if (element == "pval") {
    LR_pairs[is.na(LR_pairs)] <- 1
  }
  return(LR_pairs)
}


### Short Functions (Cell Chat)
get_pairsCC <- function(data){
  interaction <- melt(data@net$prob)
  interaction <- interaction[ ,c("Var1", "Var2", "Var3")]
  pairs <- paste0(interaction$Var1, "--", interaction$Var2, "--", interaction$Var3)
  return(pairs)
}

get_pvalCC <- function(data) {
  temp_p <- melt(data@net$pval)
  colnames(temp_p) <- c("type1", "type2", "lr_pairs", "pval")
  return(temp_p)
}

get_probCC <- function(data) {
  temp_s <- melt(data@net$prob)
  colnames(temp_s) <- c("type1", "type2", "lr_pairs", "prob")
  return(temp_s)
}

# data <- temp_p
make_unique_pairsCC <- function(data) {
  data$unique_pairs <- paste0(data$type1, "--", data$type2, "--", data$lr_pairs)
  return(data)
}

order_uniqur_pairsCC <- function(data, all_pairs) {
  data <- data[match(all_pairs, data$unique_pairs), ]
  data$unique_pairs <- all_pairs
  return(data)
}






