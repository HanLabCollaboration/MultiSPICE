#' Create Boxplot for Individual Features
#'
#' Creates boxplots comparing feature counts (LR pairs or pathways) across groups.
#' Includes statistical testing and customizable aesthetics.
#'
#' @param matrix A feature matrix with features as rows and samples as columns
#'   (output from LR_matrix or PWY_matrix)
#' @param metadata Data frame containing sample metadata with group information.
#'   Must have a 'sample' column matching matrix column names.
#' @param features Character vector of feature names to plot. If NULL, plots all features (default: NULL)
#' @param group_col Name of the column in metadata containing group labels (default: "Group")
#' @param colors Named vector of colors for groups. If NULL, uses default colors (default: NULL)
#' @param test_method Statistical test method: "wilcox_test", "t_test", or "anova" (default: "wilcox_test")
#' @param p_adjust_method P-value adjustment method: "none", "bonferroni", "BH", etc. (default: "none")
#' @param hide_ns Logical, whether to hide non-significant comparisons (default: FALSE)
#' @param ncol Number of columns for faceting when plotting multiple features (default: 2)
#' @param point_size Size of jittered points (default: 2)
#' @param title_size Font size for plot titles (default: 15)
#' @param axis_title_size Font size for axis titles (default: 14)
#' @param axis_text_size Font size for axis text (default: 12)
#' @param legend_text_size Font size for legend text (default: 14)
#' @param show_x_axis Logical, whether to show x-axis labels (default: FALSE)
#' @param y_label Label for y-axis (default: "Interaction Counts")
#'
#' @return A ggplot object or list of ggplot objects (if multiple features)
#'
#' @examples
#' \dontrun{
#' # Create LR feature matrix
#' lr_matrix <- LR_matrix(lr_list = LR_pairs_list, rank_table = rank_table, 
#'                        cutoff = 0.01, value = "pval")
#' 
#' # Create metadata
#' metadata <- data.frame(
#'   sample = colnames(lr_matrix),
#'   Group = c("Control", "Treatment", "Control", "Treatment")
#' )
#' 
#' # Plot single feature
#' feature_boxplot(lr_matrix, metadata, features = "CXCL12_CXCR4")
#' 
#' # Plot multiple features
#' features <- c("CXCL12_CXCR4", "VEGFA_VEGFR2", "TNF_TNFRSF1A")
#' feature_boxplot(lr_matrix, metadata, features = features, ncol = 3)
#' 
#' # Plot with custom colors
#' colors <- c(Control = "blue", Treatment = "red")
#' feature_boxplot(lr_matrix, metadata, features = "CXCL12_CXCR4", colors = colors)
#' }
#'
#' @export

feature_boxplot <- function(matrix,
                           metadata,
                           features = NULL,
                           group_col = "Group",
                           colors = NULL,
                           test_method = "wilcox_test",
                           p_adjust_method = "none",
                           hide_ns = FALSE,
                           ncol = 2,
                           point_size = 2,
                           title_size = 15,
                           axis_title_size = 14,
                           axis_text_size = 12,
                           legend_text_size = 14,
                           show_x_axis = FALSE,
                           y_label = "Interaction Counts") {
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required. Please install it with: install.packages('ggpubr')")
  }
  
  # Validate inputs
  if (!is.matrix(matrix) && !is.data.frame(matrix)) {
    stop("matrix must be a matrix or data frame")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  
  if (!group_col %in% colnames(metadata)) {
    stop(paste("Column", group_col, "not found in metadata"))
  }
  
  # Convert to matrix if data frame
  if (is.data.frame(matrix)) {
    matrix <- as.matrix(matrix)
  }
  
  # Transpose matrix if needed (ensure samples are columns)
  # Matrix should have features as rows and samples as columns
  
  # Select features to plot
  if (is.null(features)) {
    features <- rownames(matrix)
    message(paste("Plotting all", length(features), "features..."))
  } else {
    # Check if requested features exist
    missing_features <- setdiff(features, rownames(matrix))
    if (length(missing_features) > 0) {
      warning(paste("The following features are not found in the matrix:",
                    paste(missing_features, collapse = ", ")))
    }
    features <- intersect(features, rownames(matrix))
    
    if (length(features) == 0) {
      stop("None of the specified features are found in the matrix")
    }
  }
  
  # Prepare data for plotting
  plot_data_list <- list()
  for (feature in features) {
    # Extract feature values
    feature_values <- matrix[feature, ]
    
    # Match with metadata
    plot_df <- data.frame(
      sample = names(feature_values),
      value = as.numeric(feature_values),
      stringsAsFactors = FALSE
    )
    
    # Merge with metadata
    plot_df <- merge(plot_df, metadata, by = "sample")
    plot_df$feature <- feature
    
    plot_data_list[[feature]] <- plot_df
  }
  
  # Combine all data
  plot_data <- do.call(rbind, plot_data_list)
  
  # Set up colors
  if (is.null(colors)) {
    # Use default ggplot colors or generate colors
    n_groups <- length(unique(plot_data[[group_col]]))
    if (n_groups <= 8) {
      colors <- scales::hue_pal()(n_groups)
      names(colors) <- unique(plot_data[[group_col]])
    }
  }
  
  # Create plot function
  create_plot <- function(feature_name) {
    data_subset <- plot_data[plot_data$feature == feature_name, ]
    
    p <- ggplot(data_subset, aes_string(x = group_col, y = "value", fill = group_col)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = point_size, alpha = 0.6) +
      theme_classic() +
      labs(y = y_label, title = feature_name) +
      theme(
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank()
      )
    
    # Add colors if provided
    if (!is.null(colors)) {
      p <- p + scale_fill_manual(values = colors)
    }
    
    # Hide x-axis if requested
    if (!show_x_axis) {
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    
    # Add statistical test
    p <- p + ggpubr::geom_pwc(
      aes_string(group = group_col),
      tip.length = 0,
      method = test_method,
      label = "p.format",
      p.adjust.method = p_adjust_method,
      hide.ns = hide_ns
    )
    
    return(p)
  }
  
  # Create plots
  if (length(features) == 1) {
    # Single plot
    return(create_plot(features[1]))
  } else {
    # Multiple plots - create faceted plot or list
    if (length(features) <= 6) {
      # Create faceted plot for small number of features
      p <- ggplot(plot_data, aes_string(x = group_col, y = "value", fill = group_col)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, size = point_size, alpha = 0.6) +
        facet_wrap(~ feature, scales = "free_y", ncol = ncol) +
        theme_classic() +
        labs(y = y_label) +
        theme(
          strip.text = element_text(size = title_size - 2, face = "bold"),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size - 1),
          legend.text = element_text(size = legend_text_size),
          legend.title = element_blank()
        )
      
      # Add colors if provided
      if (!is.null(colors)) {
        p <- p + scale_fill_manual(values = colors)
      }
      
      # Hide x-axis if requested
      if (!show_x_axis) {
        p <- p + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
      }
      
      # Add statistical test
      p <- p + ggpubr::geom_pwc(
        aes_string(group = group_col),
        tip.length = 0,
        method = test_method,
        label = "p.format",
        p.adjust.method = p_adjust_method,
        hide.ns = hide_ns
      )
      
      return(p)
    } else {
      # Return list of individual plots for many features
      message(paste("Creating individual plots for", length(features), "features..."))
      plots <- lapply(features, create_plot)
      names(plots) <- features
      return(plots)
    }
  }
}
