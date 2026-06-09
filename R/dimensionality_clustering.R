#' PCA Dimension Reduction and Clustering
#'
#' Performs PCA dimension reduction on feature matrix and applies graph-based clustering
#' using Louvain, Leiden, or Walktrap algorithms.
#'
#' @param matrix A feature matrix with features as rows and samples as columns
#'   (output from LR_matrix or PWY_matrix)
#' @param scale Logical, whether to scale features before PCA (default: TRUE)
#' @param center Logical, whether to center features before PCA (default: TRUE)
#' @param n_pcs Number of principal components to compute (default: 10)
#' @param n_neighbors Number of nearest neighbors for graph construction (default: 15)
#' @param clustering_method Clustering algorithm: "louvain", "leiden", or "walktrap" (default: "louvain")
#' @param resolution Resolution parameter for Louvain/Leiden clustering (default: 1.0)
#' @param walktrap_steps Number of steps for Walktrap algorithm (default: 4)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{pca}{PCA object from prcomp}
#'   \item{pca_coords}{Matrix of PCA coordinates (samples x PCs)}
#'   \item{clusters}{Vector of cluster assignments}
#'   \item{graph}{igraph object of the k-nearest neighbor graph}
#'   \item{variance_explained}{Proportion of variance explained by each PC}
#'
#' @examples
#' \dontrun{
#' # Create LR feature matrix
#' lr_matrix <- LR_matrix(lr_list = LR_pairs_list, rank_table = rank_table, 
#'                        cutoff = 0.01, value = "pval")
#' 
#' # Perform PCA and Louvain clustering
#' results <- dimensionality_clustering(lr_matrix, clustering_method = "louvain")
#' 
#' # Perform PCA and Leiden clustering
#' results <- dimensionality_clustering(lr_matrix, clustering_method = "leiden", 
#'                                      resolution = 0.8)
#' 
#' # Perform PCA and Walktrap clustering
#' results <- dimensionality_clustering(lr_matrix, clustering_method = "walktrap",
#'                                      walktrap_steps = 5)
#' 
#' # Access results
#' pca_coords <- results$pca_coords
#' clusters <- results$clusters
#' variance <- results$variance_explained
#' }
#'
#' @export

dimensionality_clustering <- function(matrix,
                                     scale = TRUE,
                                     center = TRUE,
                                     n_pcs = 10,
                                     n_neighbors = 15,
                                     clustering_method = "louvain",
                                     resolution = 1.0,
                                     walktrap_steps = 4,
                                     seed = 42,
                                     verbose = TRUE) {
  
  # Load required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it with: install.packages('igraph')")
  }
  
  # Check for leiden if needed
  if (clustering_method == "leiden") {
    if (!requireNamespace("leiden", quietly = TRUE)) {
      stop("Package 'leiden' is required for leiden clustering. Please install it with: install.packages('leiden')")
    }
  }
  
  # Validate inputs
  if (!is.matrix(matrix) && !is.data.frame(matrix)) {
    stop("Input must be a matrix or data frame")
  }
  
  if (!clustering_method %in% c("louvain", "leiden", "walktrap")) {
    stop("clustering_method must be one of: 'louvain', 'leiden', 'walktrap'")
  }
  
  # Convert to matrix if data frame
  if (is.data.frame(matrix)) {
    matrix <- as.matrix(matrix)
  }
  
  # Transpose so samples are rows (required for PCA)
  mat_t <- t(matrix)
  
  if (verbose) {
    cat(sprintf("Performing PCA on %d samples with %d features...\n", 
                nrow(mat_t), ncol(mat_t)))
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Perform PCA
  pca_result <- prcomp(mat_t, scale. = scale, center = center)
  
  # Get PCA coordinates
  pca_coords <- pca_result$x[, 1:min(n_pcs, ncol(pca_result$x)), drop = FALSE]
  
  # Calculate variance explained
  variance_explained <- summary(pca_result)$importance[2, ]
  
  if (verbose) {
    cat(sprintf("PCA complete. Using %d principal components.\n", ncol(pca_coords)))
    cat(sprintf("Cumulative variance explained by PC1-%d: %.2f%%\n", 
                ncol(pca_coords), 
                sum(variance_explained[1:ncol(pca_coords)]) * 100))
  }
  
  # Build k-nearest neighbor graph
  if (verbose) {
    cat(sprintf("Building k-nearest neighbor graph (k=%d)...\n", n_neighbors))
  }
  
  # Compute pairwise distances in PCA space
  dist_matrix <- as.matrix(dist(pca_coords))
  
  # Build KNN graph
  knn_graph <- matrix(0, nrow = nrow(pca_coords), ncol = nrow(pca_coords))
  rownames(knn_graph) <- rownames(pca_coords)
  colnames(knn_graph) <- rownames(pca_coords)
  
  for (i in 1:nrow(dist_matrix)) {
    # Find k nearest neighbors (excluding self)
    neighbors <- order(dist_matrix[i, ])[2:(n_neighbors + 1)]
    knn_graph[i, neighbors] <- 1
  }
  
  # Make graph symmetric (mutual nearest neighbors)
  knn_graph <- pmax(knn_graph, t(knn_graph))
  
  # Create igraph object
  graph <- igraph::graph_from_adjacency_matrix(knn_graph, mode = "undirected", 
                                                weighted = NULL, diag = FALSE)
  
  # Perform clustering
  if (verbose) {
    cat(sprintf("Performing %s clustering...\n", clustering_method))
  }
  
  if (clustering_method == "louvain") {
    cluster_result <- igraph::cluster_louvain(graph, resolution = resolution)
    clusters <- igraph::membership(cluster_result)
    
  } else if (clustering_method == "leiden") {
    # Convert to adjacency matrix for leiden
    adj_matrix <- as.matrix(igraph::as_adjacency_matrix(graph))
    clusters <- leiden::leiden(adj_matrix, resolution_parameter = resolution, seed = seed)
    names(clusters) <- rownames(pca_coords)
    
  } else if (clustering_method == "walktrap") {
    cluster_result <- igraph::cluster_walktrap(graph, steps = walktrap_steps)
    clusters <- igraph::membership(cluster_result)
  }
  
  if (verbose) {
    cat(sprintf("Clustering complete. Identified %d clusters.\n", length(unique(clusters))))
    cat(sprintf("Cluster sizes: %s\n", paste(table(clusters), collapse = ", ")))
  }
  
  # Return results
  return(list(
    pca = pca_result,
    pca_coords = pca_coords,
    clusters = clusters,
    graph = graph,
    variance_explained = variance_explained
  ))
}


#' Plot PCA Results
#'
#' Creates scatter plots of PCA results colored by clusters or other metadata.
#'
#' @param pca_coords Matrix of PCA coordinates (samples x PCs)
#' @param clusters Vector of cluster assignments
#' @param pc_x Principal component for x-axis (default: 1)
#' @param pc_y Principal component for y-axis (default: 2)
#' @param color_by Vector to color points by. If NULL, uses clusters (default: NULL)
#' @param colors Color palette. If NULL, uses default colors (default: NULL)
#' @param point_size Size of points (default: 3)
#' @param main Title for the plot (default: NULL)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Perform clustering
#' results <- dimensionality_clustering(lr_matrix)
#' 
#' # Plot PCA colored by clusters
#' plot_pca(results$pca_coords, results$clusters, pc_x = 1, pc_y = 2)
#' }
#'
#' @export

plot_pca <- function(pca_coords,
                     clusters,
                     pc_x = 1,
                     pc_y = 2,
                     color_by = NULL,
                     colors = NULL,
                     point_size = 3,
                     main = NULL) {
  
  # Validate inputs
  if (!is.matrix(pca_coords) && !is.data.frame(pca_coords)) {
    stop("pca_coords must be a matrix or data frame")
  }
  
  if (max(pc_x, pc_y) > ncol(pca_coords)) {
    stop(sprintf("Requested PCs exceed available dimensions. Maximum PC: %d", ncol(pca_coords)))
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    PC_X = pca_coords[, pc_x],
    PC_Y = pca_coords[, pc_y],
    Sample = rownames(pca_coords)
  )
  
  # Determine what to color by
  if (is.null(color_by)) {
    plot_data$Color <- as.factor(clusters)
    legend_title <- "Cluster"
  } else {
    plot_data$Color <- as.factor(color_by)
    legend_title <- "Group"
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = PC_X, y = PC_Y, color = Color)) +
    geom_point(size = point_size) +
    theme_bw() +
    xlab(sprintf("PC%d", pc_x)) +
    ylab(sprintf("PC%d", pc_y)) +
    labs(color = legend_title) +
    ggtitle(main) +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  return(p)
}
