#' Create Heatmap Visualization of Feature Matrix
#'
#' Creates a heatmap visualization of LR pair or pathway feature matrices.
#' Features (LR pairs/pathways) are shown as rows and samples as columns.
#'
#' @param matrix A feature matrix with features as rows and samples as columns
#'   (output from LR_matrix or PWY_matrix)
#' @param features Character vector of feature names to plot. If NULL, plots all features (default: NULL)
#' @param scale_rows Logical, whether to scale rows (Z-score normalization) (default: FALSE)
#' @param column_annotation Named vector or data frame for column annotations (e.g., cluster assignments).
#'   Names should match column names of the matrix. If NULL, no annotation is added (default: NULL)
#' @param annotation_name Name for the annotation legend title when column_annotation is a vector (default: "Cluster")
#' @param annotation_colors Named list of colors for annotations. List names should match 
#'   annotation column names. If NULL, uses default colors (default: NULL)
#' @param name Title for the heatmap legend (default: "Counts")
#' @param cluster_rows Logical, whether to cluster rows (default: TRUE)
#' @param cluster_columns Logical, whether to cluster columns (default: TRUE)
#' @param show_row_names Logical, whether to show row names (default: TRUE for <50 rows, FALSE otherwise)
#' @param show_column_names Logical, whether to show column names (default: TRUE)
#' @param row_names_gp Text properties for row names (default: fontsize = 8)
#' @param column_names_gp Text properties for column names (default: fontsize = 10)
#' @param colors Color scheme for the heatmap. If NULL, uses a blue gradient (default: NULL)
#' @param main Title for the heatmap (default: NULL)
#'
#' @return A ComplexHeatmap object
#'
#' @examples
#' \dontrun{
#' # Create LR feature matrix
#' lr_matrix <- LR_matrix(lr_list = LR_pairs_list, rank_table = rank_table, 
#'                        cutoff = 0.01, value = "pval")
#' 
#' # Visualize all features
#' matrix_heatmap(lr_matrix, name = "LR Counts", main = "LR Pair Feature Matrix")
#' 
#' # Perform clustering
#' pca_results <- dimensionality_clustering(lr_matrix, clustering_method = "louvain")
#' 
#' # Visualize with cluster annotations
#' matrix_heatmap(lr_matrix, column_annotation = pca_results$clusters,
#'                name = "LR Counts", main = "LR Pairs with Cluster Annotations")
#' 
#' # Visualize with row scaling
#' matrix_heatmap(lr_matrix, scale_rows = TRUE, name = "Z-score", 
#'                main = "LR Pair Feature Matrix (Scaled)")
#' 
#' # Visualize selected features only
#' top_features <- c("CXCL12_CXCR4", "VEGFA_VEGFR2", "TNF_TNFRSF1A")
#' matrix_heatmap(lr_matrix, features = top_features, name = "LR Counts", 
#'                main = "Top LR Pairs")
#' }
#'
#' @export

matrix_heatmap <- function(matrix, 
                           features = NULL,
                           scale_rows = FALSE,
                           column_annotation = NULL,
                           annotation_name = "Cluster",
                           annotation_colors = NULL,
                           name = "Counts",
                           cluster_rows = TRUE,
                           cluster_columns = TRUE,
                           show_row_names = NULL,
                           show_column_names = TRUE,
                           row_names_gp = gpar(fontsize = 8),
                           column_names_gp = gpar(fontsize = 10),
                           colors = NULL,
                           main = NULL) {
  
  # Validate input
  if (!is.matrix(matrix) && !is.data.frame(matrix)) {
    stop("Input must be a matrix or data frame")
  }
  
  # Convert to matrix if data frame
  if (is.data.frame(matrix)) {
    matrix <- as.matrix(matrix)
  }
  
  # Filter for selected features if provided
  if (!is.null(features)) {
    # Check if requested features exist in the matrix
    missing_features <- setdiff(features, rownames(matrix))
    if (length(missing_features) > 0) {
      warning(paste("The following features are not found in the matrix and will be skipped:",
                    paste(missing_features, collapse = ", ")))
    }
    
    # Get features that exist
    features_to_plot <- intersect(features, rownames(matrix))
    
    if (length(features_to_plot) == 0) {
      stop("None of the specified features are found in the matrix")
    }
    
    # Subset matrix
    matrix <- matrix[features_to_plot, , drop = FALSE]
  }
  
  # Scale rows if requested (Z-score normalization)
  if (scale_rows) {
    matrix <- t(scale(t(matrix)))
    # Update legend name if using default
    if (name == "Counts") {
      name <- "Z-score"
    }
  }
  
  # Auto-determine whether to show row names based on number of rows
  if (is.null(show_row_names)) {
    show_row_names <- nrow(matrix) < 50
  }
  
  # Set default color scheme if not provided
  if (is.null(colors)) {
    if (scale_rows) {
      # Use red-white-blue for scaled data
      colors <- colorRamp2(
        c(min(matrix, na.rm = TRUE), 0, max(matrix, na.rm = TRUE)), 
        c("blue", "white", "red")
      )
    } else {
      # Use white-blue gradient for count data
      colors <- colorRamp2(
        seq(min(matrix, na.rm = TRUE), max(matrix, na.rm = TRUE), length = 3), 
        c("white", "#2a60d9", "darkblue")
      )
    }
  }
  
  # Create column annotation if provided
  col_ha <- NULL
  if (!is.null(column_annotation)) {
    # Convert to data frame if it's a vector
    if (is.vector(column_annotation)) {
      # Match annotation to matrix columns
      annot_vec <- as.factor(column_annotation[colnames(matrix)])
      annot_df <- as.data.frame(setNames(list(annot_vec), annotation_name))
      rownames(annot_df) <- colnames(matrix)
    } else {
      annot_df <- as.data.frame(column_annotation)
      # Ensure row names match matrix columns
      annot_df <- annot_df[colnames(matrix), , drop = FALSE]
      # If single column data frame, use annotation_name
      if (ncol(annot_df) == 1) {
        colnames(annot_df) <- annotation_name
      }
    }
    
    # Set up annotation colors
    if (is.null(annotation_colors)) {
      # Auto-generate colors for each annotation column
      annot_colors <- list()
      for (col in colnames(annot_df)) {
        n_levels <- length(unique(annot_df[[col]]))
        if (n_levels <= 11) {
          # Use RColorBrewer Set3 palette
          annot_colors[[col]] <- setNames(
            RColorBrewer::brewer.pal(max(3, n_levels), "Set3")[1:n_levels],
            unique(as.character(annot_df[[col]]))
          )
        } else {
          # Use rainbow colors for many levels
          annot_colors[[col]] <- setNames(
            rainbow(n_levels),
            unique(as.character(annot_df[[col]]))
          )
        }
      }
    } else {
      annot_colors <- annotation_colors
    }
    
    # Create HeatmapAnnotation using do.call to properly set names
    annot_args <- list()
    for (col in colnames(annot_df)) {
      annot_args[[col]] <- annot_df[[col]]
    }
    annot_args$col <- annot_colors
    annot_args$show_legend <- TRUE
    annot_args$show_annotation_name <- TRUE
    annot_args$annotation_name_side <- "left"
    annot_args$annotation_name_gp <- gpar(fontsize = 10)
    
    col_ha <- do.call(HeatmapAnnotation, annot_args)
  }
  
  # Create the heatmap
  ht <- Heatmap(
    matrix,
    name = name,
    col = colors,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    column_title = main,
    row_title = "Features",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title_gp = gpar(fontsize = 12),
    top_annotation = col_ha,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  return(ht)
}
