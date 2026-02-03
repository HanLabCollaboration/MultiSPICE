#' @title LRheatmap
#' @description Plot heatmap for cell-cell interaction counts. Rows represent source (col1) and columns represent target (col2).
#' @param dat Dataframe containing cell-cell interaction information
#' @param col1 Column name for the first cell type (rows/source, default: "celltype1")
#' @param col2 Column name for the second cell type (columns/target, default: "celltype2")
#' @param name Title for the heatmap legend (default: "Counts")
#' @param count_col Column name containing items to count (e.g., "pathway" or "interaction_name"). If NULL, counts rows.
#' @param cell_colors Named vector of colors for cell types. If NULL, no annotation is added.
#' @param row_title Title for rows (default: uses col1 value)
#' @param column_title Title for columns (default: uses col2 value)
#' @return A heatmap object visualizing cell-cell interaction counts


############ Plot Heatmap function ###################
LRheatmap <- function(dat, col1 = "celltype1", col2 = "celltype2", name = "Counts", 
                      count_col = NULL, cell_colors = NULL, 
                      row_title = NULL, column_title = NULL) {
  
  # Set default titles if not provided
  if (is.null(row_title)) row_title <- col1
  if (is.null(column_title)) column_title <- col2
  
  if (!is.null(count_col)) {
    # Count unique values in count_col for each source-target pair
    count_df <- aggregate(
      dat[[count_col]] ~ dat[[col1]] + dat[[col2]],
      FUN = function(x) length(unique(x))
    )
    colnames(count_df) <- c("var1", "var2", "count")
    
    # Get all unique cell types
    all_var1 <- sort(unique(dat[[col1]]))
    all_var2 <- sort(unique(dat[[col2]]))
    
    # Create matrix with all combinations
    mat <- matrix(0, nrow = length(all_var1), ncol = length(all_var2),
                  dimnames = list(all_var1, all_var2))
    
    # Fill in the counts
    for (i in 1:nrow(count_df)) {
      mat[as.character(count_df$var1[i]), as.character(count_df$var2[i])] <- count_df$count[i]
    }
  } else {
    # Original behavior: count rows
    counts <- data.frame(var1 = dat[[col1]], var2 = dat[[col2]])
    counts <- table(counts)
    mat <- matrix(counts, ncol = ncol(counts), dimnames = dimnames(counts))
  }
  
  f1 = colorRamp2(seq(min(mat), max(mat), length = 3), c("white", "#2a60d9", "darkblue"))
  
  # Create annotations if colors are provided
  if (!is.null(cell_colors)) {
    # Row annotation (source/col1)
    row_ha = rowAnnotation(
      Source = rownames(mat),
      col = list(Source = cell_colors),
      show_legend = FALSE,
      show_annotation_name = FALSE
    )
    
    # Column annotation (target/col2)
    col_ha = HeatmapAnnotation(
      Target = colnames(mat),
      col = list(Target = cell_colors),
      show_legend = FALSE,
      show_annotation_name = FALSE
    )
    
    Heatmap(mat, col = f1, cluster_columns = FALSE, cluster_rows = FALSE, name = name,
            left_annotation = row_ha, top_annotation = col_ha,
            row_names_side = "left", column_names_side = "bottom",
            row_title = row_title, column_title = column_title,
            row_title_side = "left", column_title_side = "bottom")
  } else {
    Heatmap(mat, col = f1, cluster_columns = FALSE, cluster_rows = FALSE, name = name,
            row_title = row_title, column_title = column_title,
            row_title_side = "left", column_title_side = "bottom")
  }
}


