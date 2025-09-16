#' @title Load Spatalk LR results for Multi-sample analysis
#' @description Load Spatalk LR results for Multi-sample analysis
#' @param input_list A list of Spatalk LR results
#' @param element The element to extract from the Spatalk LR results. Default is "pval".
#' @return A data frame of LR pairs and their significance across samples.
#' @importFrom dplyr %>%
#' @export

### Spatalk Data loading
# element = "pval"
loadingLR_spatalk <- function(input_list, element = "pval") {

  if (length(input_list) < 2) {
    stop("Need at least 2 samples for muti-sample analysis...")
  }

  message("Gathering LR pair names...")
  all_pairs <- lapply(input_list, get_pairsSP)
  all_pairs <- unique(unlist(all_pairs))

  message("Making unique LR pair names...")
  pairs_table <- lapply(input_list, get_pvalSP)
  pairs_table <- lapply(pairs_table, make_unique_pairsSP)
  pairs_table <- lapply(pairs_table, order_uniqur_pairsSP, all_pairs = all_pairs)

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

###########################################################################
### Short Functions (Spatalk)
get_pairsSP <- function(data){
  interaction <- data@lrpair
  pairs <- paste0(interaction$celltype_sender, "--", interaction$celltype_receiver,
                  "--", interaction$ligand, "_", interaction$receptor)
  return(pairs)
}


get_pvalSP <- function(data) {
  interaction <- data@lrpair
  temp_p <- interaction[, c("celltype_sender", "celltype_receiver", "ligand", "receptor", "lr_co_ratio_pvalue")]
  colnames(temp_p) <- c("celltype_sender", "celltype_receiver", "ligand", "receptor", "pval")
  return(temp_p)
}

make_unique_pairsSP <- function(data) {
  data$unique_pairs <- paste0(data[["celltype_sender"]], "--", data[["celltype_receiver"]], "--",
                              data[["ligand"]], "_", data[["receptor"]])
  return(data)
}

order_uniqur_pairsSP <- function(data, all_pairs) {
  data <- data[match(all_pairs, data$unique_pairs), ]
  data$unique_pairs <- all_pairs
  return(data)
}

