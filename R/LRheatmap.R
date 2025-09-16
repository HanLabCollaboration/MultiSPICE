#' @title LRheatmap
#' @description Plot heatmap for cell-cell interaction counts
#' @param dat Dataframe containing cell-cell interaction information
#' @param col1 Column name for the first cell type (default: "celltype1")
#' @param col2 Column name for the second cell type (default: "celltype2")
#' @param col Color vector for the heatmap (default: col_vector)
#' @param name Title for the heatmap legend (default: "Counts")
#' @return A heatmap object visualizing cell-cell interaction counts


############ Plot Heatmap function ###################
LRheatmap <- function(dat, col1 = "celltype1", col2 = "celltype2", col = col_vector, name = "Counts") {
  counts <- tibble(var1 = rank_table[[col1]], var2 = rank_table[[col2]])
  counts <- table(counts)
  mat <- matrix(counts, ncol = ncol(counts), dimnames = dimnames(counts))
  f1 = colorRamp2(seq(min(mat), max(mat), length = 3), c("white", "#2a60d9", "darkblue"))
  Heatmap(mat, col = f1, cluster_columns = FALSE, cluster_rows = FALSE, name = name)
}


