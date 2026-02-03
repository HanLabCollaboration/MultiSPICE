#' Calculate Heterogeneity of LR Pairs or Pathways Across Samples
#'
#' Identify LR pairs or pathways that show high variability across samples, indicating
#' sample-specific or heterogeneous communication patterns. This function takes the
#' output from loading functions (e.g., loadingLR_cellchat, loadingPWY_cellchat) and
#' calculates heterogeneity metrics to understand consistency across samples.
#'
#' @param data_list A named list of data frames, where each element represents one sample.
#'   This is the output from loading functions like loadingLR_cellchat() or loadingPWY_cellchat().
#'   Each data frame should have columns:
#'   \itemize{
#'     \item For LR data: interaction_name, pval (and optionally prob, count)
#'     \item For pathway data: pathway_name, pval (and optionally prob, count)
#'   }
#' @param feature_col Character string specifying the column name for features.
#'   Default is "interaction_name" for LR pairs. Use "pathway_name" for pathways.
#' @param value_col Character string specifying which value to use for calculating
#'   heterogeneity. Options: "pval" (default), "prob", or "count".
#' @param pval_threshold Numeric value for p-value threshold. Features with pval >
#'   threshold are considered "absent" in that sample. Default is 0.05.
#' @param method Character string specifying the heterogeneity metric. Options:
#'   \itemize{
#'     \item "frequency" (default): Calculate how many samples each feature appears in
#'     \item "cv": Coefficient of variation of values across samples
#'     \item "range": Range (max - min) of values
#'     \item "all": Calculate all metrics
#'   }
#' @param min_samples Integer specifying minimum number of samples where a feature
#'   must be present to be included. Default is 2.
#' @param top_n Integer specifying how many features to return, sorted by heterogeneity.
#'   If NULL (default), returns all features.
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item feature: Feature name (LR pair or pathway)
#'     \item n_samples: Number of samples where feature is present (significant)
#'     \item frequency: Proportion of samples where feature is present
#'     \item mean_value: Mean value (pval/prob/count) across samples where present
#'     \item Additional columns based on method selected
#'   }
#'
#' @details
#' This function helps identify which features (LR pairs or pathways) are:
#' - **Consistent**: Present in most/all samples (high frequency, low heterogeneity)
#' - **Variable**: Present in some samples but not others (medium frequency)
#' - **Rare**: Present in few samples (low frequency)
#'
#' Understanding heterogeneity is crucial before performing downstream analyses,
#' as it reveals whether you're working with consistent biological patterns or
#' sample-specific noise.
#'
#' @examples
#' \dontrun{
#' # Load LR pair data
#' lr_data <- loadingLR_cellchat(sample = sample_names,
#'                               directory = data_dir,
#'                               value = "pval")
#'
#' # Calculate frequency-based heterogeneity
#' lr_het <- calculate_heterogeneity(
#'   lr_data,
#'   feature_col = "interaction_name",
#'   value_col = "pval",
#'   pval_threshold = 0.05,
#'   method = "frequency"
#' )
#'
#' # View most heterogeneous features
#' head(lr_het)
#'
#' # Calculate all metrics for pathway data
#' pwy_data <- loadingPWY_cellchat(sample = sample_names,
#'                                 directory = data_dir,
#'                                 value = "pval")
#'
#' pwy_het <- calculate_heterogeneity(
#'   pwy_data,
#'   feature_col = "pathway_name",
#'   method = "all"
#' )
#' }
#'
#' @export
calculate_heterogeneity <- function(data_list,
                                    feature_col = "interaction_name",
                                    value_col = "pval",
                                    pval_threshold = 0.05,
                                    method = "frequency",
                                    min_samples = 2,
                                    top_n = NULL) {
  
  # Input validation
  if (!is.list(data_list)) {
    stop("data_list must be a list of data frames")
  }
  
  if (length(data_list) == 0) {
    stop("data_list cannot be empty")
  }
  
  # Check that all elements are data frames
  if (!all(sapply(data_list, is.data.frame))) {
    stop("All elements in data_list must be data frames")
  }
  
  valid_methods <- c("frequency", "cv", "range", "all")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Check that feature_col exists in all data frames
  has_feature_col <- sapply(data_list, function(df) feature_col %in% colnames(df))
  if (!all(has_feature_col)) {
    stop("Column '", feature_col, "' not found in all data frames")
  }
  
  # Get all unique features across all samples
  all_features <- unique(unlist(lapply(data_list, function(df) {
    if ("pval" %in% colnames(df)) {
      as.character(df[[feature_col]][df$pval <= pval_threshold])
    } else {
      as.character(df[[feature_col]])
    }
  })))
  
  if (length(all_features) == 0) {
    warning("No features found meeting the criteria")
    return(data.frame())
  }
  
  n_samples_total <- length(data_list)
  
  # Calculate frequency: how many samples each feature appears in
  feature_frequency <- sapply(all_features, function(feat) {
    sum(sapply(data_list, function(df) {
      # Convert to character to avoid factor level comparison issues
      feat_col_char <- as.character(df[[feature_col]])
      if ("pval" %in% colnames(df)) {
        any(feat_col_char == feat & df$pval <= pval_threshold)
      } else {
        any(feat_col_char == feat)
      }
    }))
  })
  
  # Initialize results data frame
  results <- data.frame(
    feature = all_features,
    n_samples = feature_frequency,
    frequency = feature_frequency / n_samples_total,
    stringsAsFactors = FALSE
  )
  
  # Filter by minimum samples
  keep_idx <- results$n_samples >= min_samples
  if (sum(keep_idx) == 0) {
    warning("No features meet the min_samples criterion")
    return(data.frame())
  }
  results <- results[keep_idx, ]
  
  # Calculate additional metrics based on method
  if (method %in% c("cv", "range", "all") && value_col %in% c("pval", "prob", "count")) {
    
    # Extract values for each feature across samples
    feature_values <- lapply(results$feature, function(feat) {
      vals <- sapply(data_list, function(df) {
        # Convert to character to avoid factor comparison issues
        feat_col_char <- as.character(df[[feature_col]])
        row_idx <- which(feat_col_char == feat)
        if (length(row_idx) > 0 && "pval" %in% colnames(df)) {
          if (df$pval[row_idx[1]] <= pval_threshold) {
            return(df[[value_col]][row_idx[1]])
          }
        } else if (length(row_idx) > 0) {
          return(df[[value_col]][row_idx[1]])
        }
        return(NA)
      })
      return(vals[!is.na(vals)])
    })
    
    # Calculate mean values
    results$mean_value <- sapply(feature_values, function(x) {
      if (length(x) > 0) mean(x) else NA
    })
    
    if (method == "cv" || method == "all") {
      results$sd_value <- sapply(feature_values, function(x) {
        if (length(x) > 1) sd(x) else 0
      })
      results$cv <- ifelse(results$mean_value > 0,
                          results$sd_value / results$mean_value,
                          0)
    }
    
    if (method == "range" || method == "all") {
      results$range <- sapply(feature_values, function(x) {
        if (length(x) > 0) max(x) - min(x) else 0
      })
    }
  }
  
  # Sort results
  if (method == "frequency") {
    # Sort by frequency (descending) then by n_samples
    results <- results[order(-results$frequency, -results$n_samples), ]
  } else if (method == "cv") {
    results <- results[order(-results$cv), ]
  } else if (method == "range") {
    results <- results[order(-results$range), ]
  } else if (method == "all") {
    # Sort by frequency for "all" method
    results <- results[order(-results$frequency), ]
  }
  
  # Select top N if specified
  if (!is.null(top_n) && top_n < nrow(results)) {
    results <- results[1:top_n, ]
  }
  
  rownames(results) <- NULL
  
  return(results)
}


#' Plot Heterogeneity Results
#'
#' Visualize heterogeneity of features across samples. Creates visualizations
#' to show which features are consistent vs. variable across samples.
#'
#' @param heterogeneity_results Data frame output from calculate_heterogeneity()
#' @param plot_type Character string specifying plot type:
#'   \itemize{
#'     \item "frequency" (default): Barplot showing how many samples each feature appears in
#'     \item "scatter": Scatter plot showing mean value vs. frequency
#'   }
#' @param top_n Integer specifying number of features to display. Default is 20.
#' @param color_by Character string specifying column to use for coloring.
#'   Default is "frequency". Set to NULL for uniform color.
#' @param title Character string for plot title. Default is "Feature Heterogeneity Across Samples"
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Load data
#' lr_data <- loadingLR_cellchat(sample = sample_names,
#'                               directory = data_dir,
#'                               value = "pval")
#'
#' # Calculate heterogeneity
#' het_results <- calculate_heterogeneity(lr_data, method = "frequency")
#'
#' # Create frequency barplot
#' p1 <- plot_heterogeneity(het_results, plot_type = "frequency", top_n = 20)
#' print(p1)
#'
#' # Create scatter plot if CV was calculated
#' het_results_cv <- calculate_heterogeneity(lr_data, method = "all")
#' p2 <- plot_heterogeneity(het_results_cv, plot_type = "scatter")
#' print(p2)
#' }
#'
#' @export
#' @import ggplot2
plot_heterogeneity <- function(heterogeneity_results,
                               plot_type = "frequency",
                               top_n = 20,
                               color_by = "frequency",
                               title = "Feature Heterogeneity Across Samples") {
  
  # Input validation
  if (!is.data.frame(heterogeneity_results)) {
    stop("heterogeneity_results must be a data frame")
  }
  
  if (nrow(heterogeneity_results) == 0) {
    stop("heterogeneity_results is empty")
  }
  
  # Select top N features
  plot_data <- head(heterogeneity_results, top_n)
  
  if (plot_type == "frequency") {
    # Bar plot: feature frequency across samples
    plot_data$feature <- factor(plot_data$feature, levels = rev(plot_data$feature))
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = feature, y = n_samples)) +
      ggplot2::geom_bar(stat = "identity",
                       ggplot2::aes(fill = if (!is.null(color_by) && color_by %in% colnames(plot_data))
                                           .data[[color_by]] else NULL)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Feature",
                   y = "Number of Samples",
                   title = title,
                   fill = color_by) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 10),
        legend.position = "right"
      )
    
    if (!is.null(color_by) && color_by %in% colnames(plot_data)) {
      p <- p + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue")
    }
    
  } else if (plot_type == "scatter") {
    # Scatter plot: mean value vs frequency (or CV if available)
    
    if (!"mean_value" %in% colnames(plot_data)) {
      stop("Scatter plot requires 'mean_value' column. Use method='cv' or method='all' in calculate_heterogeneity()")
    }
    
    # Determine y-axis variable
    y_var <- if ("cv" %in% colnames(plot_data)) {
      "cv"
    } else {
      "frequency"
    }
    
    y_label <- if (y_var == "cv") {
      "Coefficient of Variation"
    } else {
      "Frequency (proportion of samples)"
    }
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mean_value, y = .data[[y_var]])) +
      ggplot2::geom_point(size = 3, alpha = 0.7,
                         ggplot2::aes(color = if (!is.null(color_by) && color_by %in% colnames(plot_data))
                                              .data[[color_by]] else NULL)) +
      ggplot2::labs(x = "Mean Value",
                   y = y_label,
                   title = title,
                   color = color_by) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(size = 12),
        legend.position = "right"
      )
    
    if (!is.null(color_by) && color_by %in% colnames(plot_data)) {
      p <- p + ggplot2::scale_color_gradient(low = "blue", high = "red")
    }
    
  } else {
    stop("plot_type must be 'frequency' or 'scatter'")
  }
  
  return(p)
}
