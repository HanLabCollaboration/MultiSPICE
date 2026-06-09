#' Feature labels for heterogeneity (interaction only or source-target-LR composite)
#'
#' @param df A data frame.
#' @param feature_col Column name for the feature identifier.
#' @param feature_key \code{"interaction_name"} uses \code{feature_col} only;
#'   \code{"source_target"} uses \code{paste(source, target, feature, sep = "|")}.
#' @noRd
.hetero_labels_from_df <- function(df, feature_col, feature_key) {
  if (feature_key == "interaction_name") {
    as.character(df[[feature_col]])
  } else {
    paste(
      as.character(df$source),
      as.character(df$target),
      as.character(df[[feature_col]]),
      sep = "|"
    )
  }
}


#' Binary Shannon entropy in log2 for a Bernoulli(p) presence pattern across samples
#'
#' @noRd
.bernoulli_entropy_log2 <- function(p) {
  p <- pmax(pmin(p, 1), 0)
  ifelse(
    p <= 0 | p >= 1,
    0,
    -(p * log2(p) + (1 - p) * log2(1 - p))
  )
}


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
#' @param top_n Integer specifying how many features to return, after sorting by
#'   \code{sort_by}. If NULL (default), returns all features.
#' @param sort_by Character. Column used to order results before optional \code{top_n}
#'   trimming: \code{"frequency"} (default, most shared first), \code{"specificity"}
#'   (rare across cohort first; uses \code{1 - frequency}), \code{"n_samples"},
#'   \code{"cv"}, \code{"range"}, or \code{"presence_entropy"}.
#' @param feature_key For LR tables with \code{source} and \code{target} columns,
#'   \code{"source_target"} defines a feature as sender, receiver, and interaction
#'   (preserves cell-type context). Default \code{"interaction_name"} matches prior
#'   behaviour (ligand-receptor name only).
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item feature: Feature name (LR pair or pathway, or composite key)
#'     \item n_samples: Number of samples where feature is present (significant)
#'     \item frequency: Proportion of samples where feature is present
#'     \item specificity_index: \code{1 - frequency} (high = sample-specific)
#'     \item presence_entropy: Shannon entropy (log2) of the binary presence vector
#'     \item mean_value: Mean value (pval/prob/count) across samples where present
#'     \item Additional columns based on method selected (cv, range, sd_value)
#'   }
#'
#' @details
#' This function helps identify which features (LR pairs or pathways) are:
#' \itemize{
#'   \item **Consistent**: Present in most/all samples (high frequency, low specificity_index)
#'   \item **Variable**: Present in some samples but not others (medium frequency)
#'   \item **Rare**: Present in few samples (low frequency, high specificity_index)
#' }
#'
#' **CV and range:** For \code{method} \code{"cv"}, \code{"range"}, or \code{"all"},
#' values are taken only from samples where the feature is significant (or from all
#' samples if no \code{pval} column). Samples where the feature is absent are omitted
#' from the variance calculation, so CV/range can be inflated for sparse features
#' that are significant in only a few samples. Use \code{sort_by = "specificity"} or
#' \code{lr_presence_matrix()} / \code{sample_lr_dissimilarity()} for a sample-centric view.
#'
#' **Downstream:** \code{LR_matrix()} returns features as rows and samples as columns,
#' matching \code{matrix_heatmap()} and \code{dimensionality_clustering()}.
#'
#' @examples
#' \dontrun{
#' # Load LR pair data
#' lr_data <- loadingLR_cellchat(sample = sample_names,
#'                               directory = data_dir,
#'                               value = "pval")
#'
#' lr_het <- calculate_heterogeneity(
#'   lr_data,
#'   feature_col = "interaction_name",
#'   value_col = "pval",
#'   pval_threshold = 0.05,
#'   method = "frequency",
#'   sort_by = "specificity"
#' )
#'
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
                                    top_n = NULL,
                                    sort_by = "frequency",
                                    feature_key = c("interaction_name", "source_target")) {
  feature_key <- match.arg(feature_key)

  if (!is.list(data_list)) {
    stop("data_list must be a list of data frames")
  }

  if (length(data_list) == 0) {
    stop("data_list cannot be empty")
  }

  if (!all(sapply(data_list, is.data.frame))) {
    stop("All elements in data_list must be data frames")
  }

  valid_methods <- c("frequency", "cv", "range", "all")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  valid_sort <- c(
    "frequency", "specificity", "n_samples", "cv", "range", "presence_entropy"
  )
  if (!sort_by %in% valid_sort) {
    stop("sort_by must be one of: ", paste(valid_sort, collapse = ", "))
  }

  has_feature_col <- sapply(data_list, function(df) feature_col %in% colnames(df))
  if (!all(has_feature_col)) {
    stop("Column '", feature_col, "' not found in all data frames")
  }

  if (feature_key == "source_target") {
    has_st <- sapply(data_list, function(df) {
      all(c("source", "target") %in% colnames(df))
    })
    if (!all(has_st)) {
      stop("feature_key = 'source_target' requires 'source' and 'target' in all data frames")
    }
  }

  all_features <- unique(unlist(lapply(data_list, function(df) {
    if ("pval" %in% colnames(df)) {
      ok <- df$pval <= pval_threshold
    } else {
      ok <- rep(TRUE, nrow(df))
    }
    if (!any(ok)) {
      return(character(0))
    }
    unique(.hetero_labels_from_df(df[ok, , drop = FALSE], feature_col, feature_key))
  })))

  if (length(all_features) == 0) {
    warning("No features found meeting the criteria")
    return(data.frame())
  }

  n_samples_total <- length(data_list)

  feature_frequency <- sapply(all_features, function(feat) {
    sum(sapply(data_list, function(df) {
      if ("pval" %in% colnames(df)) {
        ok <- df$pval <= pval_threshold
      } else {
        ok <- rep(TRUE, nrow(df))
      }
      if (!any(ok)) {
        return(FALSE)
      }
      lab <- .hetero_labels_from_df(df[ok, , drop = FALSE], feature_col, feature_key)
      any(lab == feat)
    }))
  })

  results <- data.frame(
    feature = all_features,
    n_samples = feature_frequency,
    frequency = feature_frequency / n_samples_total,
    stringsAsFactors = FALSE
  )

  results$specificity_index <- 1 - results$frequency
  results$presence_entropy <- .bernoulli_entropy_log2(results$frequency)

  keep_idx <- results$n_samples >= min_samples
  if (sum(keep_idx) == 0) {
    warning("No features meet the min_samples criterion")
    return(data.frame())
  }
  results <- results[keep_idx, , drop = FALSE]

  if (method %in% c("cv", "range", "all") && value_col %in% c("pval", "prob", "count")) {
    if (!all(sapply(data_list, function(df) value_col %in% colnames(df)))) {
      stop("value_col '", value_col, "' not found in all data frames")
    }

    feature_values <- lapply(results$feature, function(feat) {
      vals <- sapply(data_list, function(df) {
        if ("pval" %in% colnames(df)) {
          ok <- df$pval <= pval_threshold
        } else {
          ok <- rep(TRUE, nrow(df))
        }
        if (!any(ok)) {
          return(NA_real_)
        }
        sub <- df[ok, , drop = FALSE]
        lab <- .hetero_labels_from_df(sub, feature_col, feature_key)
        row_idx <- which(lab == feat)
        if (length(row_idx) == 0) {
          return(NA_real_)
        }
        sub[[value_col]][row_idx[1]]
      })
      vals[!is.na(vals)]
    })

    results$mean_value <- sapply(feature_values, function(x) {
      if (length(x) > 0) mean(x) else NA_real_
    })

    if (method == "cv" || method == "all") {
      results$sd_value <- sapply(feature_values, function(x) {
        if (length(x) > 1) stats::sd(x) else 0
      })
      results$cv <- ifelse(
        results$mean_value > 0,
        results$sd_value / results$mean_value,
        0
      )
    }

    if (method == "range" || method == "all") {
      results$range <- sapply(feature_values, function(x) {
        if (length(x) > 0) max(x) - min(x) else 0
      })
    }
  }

  ord <- switch(
    sort_by,
    frequency = order(-results$frequency, -results$n_samples, results$feature),
    specificity = order(-results$specificity_index, results$frequency, results$feature),
    n_samples = order(-results$n_samples, -results$frequency, results$feature),
    cv = {
      if (!"cv" %in% names(results)) {
        stop("sort_by = 'cv' requires method 'cv' or 'all'")
      }
      order(-results$cv, -results$specificity_index, results$feature)
    },
    range = {
      if (!"range" %in% names(results)) {
        stop("sort_by = 'range' requires method 'range' or 'all'")
      }
      order(-results$range, -results$specificity_index, results$feature)
    },
    presence_entropy = order(-results$presence_entropy, -results$specificity_index, results$feature)
  )
  results <- results[ord, , drop = FALSE]
  rownames(results) <- NULL

  if (!is.null(top_n) && top_n < nrow(results)) {
    results <- results[seq_len(top_n), , drop = FALSE]
  }

  results
}


#' Ligand-Receptor presence matrix across samples
#'
#' Builds a binary (0/1) matrix: rows are features, columns are samples. Each cell is 1
#' if the feature is significant in that sample. Compatible with \code{matrix_heatmap()}
#' and \code{dimensionality_clustering()} (features as rows, samples as columns), same
#' orientation as \code{LR_matrix()}.
#'
#' @param data_list Named list of per-sample data frames (e.g. from loadingLR_cellchat).
#' @param feature_col Feature identifier column (default \code{"interaction_name"}).
#' @param pval_threshold Significance threshold on \code{pval} when present.
#' @param feature_key See \code{\link{calculate_heterogeneity}}.
#'
#' @return A numeric matrix with values 0 or 1, rows named by feature, columns by sample.
#'
#' @seealso \code{\link{matrix_heatmap}}, \code{\link{LR_matrix}}, \code{\link{sample_lr_dissimilarity}}
#'
#' @export
lr_presence_matrix <- function(data_list,
                               feature_col = "interaction_name",
                               pval_threshold = 0.05,
                               feature_key = c("interaction_name", "source_target")) {
  feature_key <- match.arg(feature_key)

  if (!is.list(data_list) || !length(data_list)) {
    stop("data_list must be a non-empty list of data frames")
  }
  if (!all(sapply(data_list, is.data.frame))) {
    stop("All elements in data_list must be data frames")
  }

  has_feature_col <- sapply(data_list, function(df) feature_col %in% colnames(df))
  if (!all(has_feature_col)) {
    stop("Column '", feature_col, "' not found in all data frames")
  }

  if (feature_key == "source_target") {
    has_st <- sapply(data_list, function(df) {
      all(c("source", "target") %in% colnames(df))
    })
    if (!all(has_st)) {
      stop("feature_key = 'source_target' requires 'source' and 'target' in all data frames")
    }
  }

  sample_names <- names(data_list)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_along(data_list))
  }

  all_features <- unique(unlist(lapply(data_list, function(df) {
    if ("pval" %in% colnames(df)) {
      ok <- df$pval <= pval_threshold
    } else {
      ok <- rep(TRUE, nrow(df))
    }
    if (!any(ok)) {
      return(character(0))
    }
    unique(.hetero_labels_from_df(df[ok, , drop = FALSE], feature_col, feature_key))
  })))

  if (length(all_features) == 0) {
    warning("No significant features found; returning empty matrix")
    return(matrix(numeric(0), nrow = 0, ncol = length(sample_names),
                  dimnames = list(character(0), sample_names)))
  }

  mat <- matrix(0, nrow = length(all_features), ncol = length(data_list),
                dimnames = list(all_features, sample_names))

  for (j in seq_along(data_list)) {
    df <- data_list[[j]]
    if ("pval" %in% colnames(df)) {
      ok <- df$pval <= pval_threshold
    } else {
      ok <- rep(TRUE, nrow(df))
    }
    if (!any(ok)) {
      next
    }
    sub <- df[ok, , drop = FALSE]
    lab <- .hetero_labels_from_df(sub, feature_col, feature_key)
    u <- unique(lab)
    mat[u, j] <- 1
  }

  mat
}


#' Sample-level LR dissimilarity and overlap summaries
#'
#' Computes pairwise Jaccard similarity on significant feature sets per sample,
#' optional Spearman correlation of within-sample ranks on the feature intersection
#' (requires \code{rank_value_col} in all data frames), and cohort-level exclusivity
#' summaries.
#'
#' @param data_list Named list of per-sample data frames.
#' @param feature_col Feature column (default \code{"interaction_name"}).
#' @param pval_threshold Significance cutoff for membership in each sample's set.
#' @param feature_key See \code{\link{calculate_heterogeneity}}.
#' @param rank_value_col If non-NULL (e.g. \code{"prob"}), Spearman correlation between
#'   samples is computed on features significant in both samples, using this column
#'   (first matching row per feature label). If NULL, Spearman output is omitted.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{jaccard}: symmetric matrix of Jaccard indices between sample sets
#'     \item \code{jaccard_dissimilarity}: \code{1 - jaccard}
#'     \item \code{pairwise_spearman}: symmetric matrix of Spearman correlations, or NULL
#'     \item \code{sample_exclusive_counts}: integer vector, per sample count of features
#'       significant only in that sample
#'     \item \code{summary}: list with \code{mean_offdiag_jaccard},
#'       \code{frac_features_singleton} (fraction of union features present in exactly one sample)
#'   }
#'
#' @seealso \code{\link{lr_presence_matrix}}, \code{\link{plot_lr_sample_comparison}}
#'
#' @export
sample_lr_dissimilarity <- function(data_list,
                                    feature_col = "interaction_name",
                                    pval_threshold = 0.05,
                                    feature_key = c("interaction_name", "source_target"),
                                    rank_value_col = NULL) {
  feature_key <- match.arg(feature_key)

  pres <- lr_presence_matrix(
    data_list,
    feature_col = feature_col,
    pval_threshold = pval_threshold,
    feature_key = feature_key
  )

  sample_names <- colnames(pres)
  n <- length(sample_names)
  if (n == 0) {
    stop("No samples in presence matrix")
  }

  jacc <- matrix(0, n, n, dimnames = list(sample_names, sample_names))
  if (nrow(pres) == 0) {
    return(list(
      jaccard = jacc,
      jaccard_dissimilarity = 1 - jacc,
      pairwise_spearman = NULL,
      sample_exclusive_counts = setNames(integer(n), sample_names),
      summary = list(
        mean_offdiag_jaccard = NA_real_,
        frac_features_singleton = NA_real_
      )
    ))
  }

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      a <- pres[, i] == 1
      b <- pres[, j] == 1
      inter <- sum(a & b)
      uni <- sum(a | b)
      jacc[i, j] <- if (uni > 0) inter / uni else 0
    }
  }

  singleton_mask <- rowSums(pres) == 1
  frac_singleton <- if (nrow(pres) > 0) {
    sum(singleton_mask) / nrow(pres)
  } else {
    NA_real_
  }

  off <- jacc[lower.tri(jacc) | upper.tri(jacc)]
  mean_off <- if (length(off)) mean(off) else NA_real_

  exclusive <- integer(n)
  names(exclusive) <- sample_names
  for (i in seq_len(n)) {
    exclusive[i] <- sum(pres[, i] == 1 & rowSums(pres) == 1)
  }

  spearman_mat <- NULL
  if (!is.null(rank_value_col)) {
    if (!all(sapply(data_list, function(df) rank_value_col %in% colnames(df)))) {
      stop("rank_value_col not found in all data frames")
    }
    spearman_mat <- matrix(NA_real_, n, n, dimnames = list(sample_names, sample_names))
    diag(spearman_mat) <- 1

    .value_for_label <- function(df, lab) {
      if ("pval" %in% colnames(df)) {
        ok <- df$pval <= pval_threshold
      } else {
        ok <- rep(TRUE, nrow(df))
      }
      if (!any(ok)) {
        return(NA_real_)
      }
      sub <- df[ok, , drop = FALSE]
      labs <- .hetero_labels_from_df(sub, feature_col, feature_key)
      hit <- which(labs == lab)
      if (!length(hit)) {
        return(NA_real_)
      }
      sub[[rank_value_col]][hit[1]]
    }

    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        fi <- rownames(pres)[pres[, i] == 1 & pres[, j] == 1]
        if (length(fi) < 2) {
          spearman_mat[i, j] <- spearman_mat[j, i] <- NA_real_
          next
        }
        v1 <- sapply(fi, function(lab) {
          .value_for_label(data_list[[sample_names[i]]], lab)
        })
        v2 <- sapply(fi, function(lab) {
          .value_for_label(data_list[[sample_names[j]]], lab)
        })
        okp <- is.finite(v1) & is.finite(v2)
        if (sum(okp) < 2) {
          spearman_mat[i, j] <- spearman_mat[j, i] <- NA_real_
        } else {
          r <- suppressWarnings(stats::cor(v1[okp], v2[okp], method = "spearman"))
          spearman_mat[i, j] <- spearman_mat[j, i] <- r
        }
      }
    }
  }

  list(
    jaccard = jacc,
    jaccard_dissimilarity = 1 - jacc,
    pairwise_spearman = spearman_mat,
    sample_exclusive_counts = exclusive,
    summary = list(
      mean_offdiag_jaccard = mean_off,
      frac_features_singleton = frac_singleton
    )
  )
}


#' Plot sample-centric LR heterogeneity
#'
#' @param data_list Named list of per-sample LR data frames.
#' @param diss Pre-computed output from \code{sample_lr_dissimilarity()}. If NULL,
#'   it is computed from \code{data_list} and other arguments.
#' @param plot_type One of \code{"jaccard"} (tile heatmap of Jaccard similarity) or
#'   \code{"pca_samples"} (PCA of sample-by-LR presence; features with zero variance
#'   across samples are dropped first).
#' @param feature_col,pval_threshold,feature_key,rank_value_col Passed to
#'   \code{sample_lr_dissimilarity()} when \code{diss} is NULL.
#' @param title Plot title.
#' @param label_samples If TRUE (default) and \code{plot_type = "pca_samples"},
#'   draw sample name labels next to points. Ignored for \code{plot_type = "jaccard"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @import ggplot2
#' @seealso \code{\link{sample_lr_dissimilarity}}, \code{\link{lr_presence_matrix}},
#'   \code{\link{plot_lr_mean_jaccard_boxplot}}
#'
#' @export
plot_lr_sample_comparison <- function(data_list,
                                      diss = NULL,
                                      plot_type = c("jaccard", "pca_samples"),
                                      feature_col = "interaction_name",
                                      pval_threshold = 0.05,
                                      feature_key = c("interaction_name", "source_target"),
                                      rank_value_col = NULL,
                                      title = NULL,
                                      label_samples = TRUE) {
  plot_type <- match.arg(plot_type)
  feature_key <- match.arg(feature_key)

  if (is.null(diss)) {
    diss <- sample_lr_dissimilarity(
      data_list,
      feature_col = feature_col,
      pval_threshold = pval_threshold,
      feature_key = feature_key,
      rank_value_col = rank_value_col
    )
  }

  if (plot_type == "jaccard") {
    J <- diss$jaccard
    sn <- rownames(J)
    df <- expand.grid(sample_i = sn, sample_j = sn, stringsAsFactors = FALSE)
    df$jaccard <- as.vector(J)
    if (is.null(title)) {
      title <- "Pairwise Jaccard similarity (significant LR sets)"
    }
    ggplot2::ggplot(df, ggplot2::aes(x = .data$sample_i, y = .data$sample_j, fill = .data$jaccard)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Jaccard", title = title) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  } else if (plot_type == "pca_samples") {
    pres <- lr_presence_matrix(
      data_list,
      feature_col = feature_col,
      pval_threshold = pval_threshold,
      feature_key = feature_key
    )
    if (!ncol(pres) || !nrow(pres)) {
      stop("Presence matrix is empty; cannot run PCA")
    }
    X <- t(pres)
    col_sd <- apply(X, 2, stats::sd, na.rm = TRUE)
    keep <- is.finite(col_sd) & col_sd > .Machine$double.eps
    if (sum(keep) < 2L) {
      stop("Fewer than two LR features with variation across samples; cannot run PCA")
    }
    X <- X[, keep, drop = FALSE]
    pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
    vx <- summary(pca)$importance[2, ]
    pc_df <- data.frame(
      Sample = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      stringsAsFactors = FALSE
    )
    if (is.null(title)) {
      title <- sprintf(
        "PCA of LR presence (PC1=%.0f%%, PC2=%.0f%% variance)",
        100 * vx[1], 100 * vx[2]
      )
    }
    p <- ggplot2::ggplot(pc_df, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
      ggplot2::geom_point(size = 3, color = "darkred") +
      ggplot2::labs(
        x = sprintf("PC1 (%.0f%%)", 100 * vx[1]),
        y = sprintf("PC2 (%.0f%%)", 100 * vx[2]),
        title = title
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))
    if (isTRUE(label_samples)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$Sample),
        vjust = -0.5,
        size = 3
      )
    }
    p
  }
}


#' Heatmap of LR presence for the most variable features
#'
#' Subsets \code{lr_presence_matrix()} to the top features by row variance (binary
#' variance = p(1-p) with p = row mean) and passes the result to \code{matrix_heatmap()}.
#'
#' @param data_list Named list of per-sample data frames.
#' @param top_n Number of features (rows) to display.
#' @param feature_col,pval_threshold,feature_key Passed to \code{lr_presence_matrix()}.
#' @param ... Additional arguments to \code{matrix_heatmap()}.
#'
#' @return A \code{ComplexHeatmap} object from \code{matrix_heatmap()}.
#'
#' @seealso \code{\link{matrix_heatmap}}, \code{\link{lr_presence_matrix}}
#'
#' @export
plot_lr_presence_heatmap <- function(data_list,
                                     top_n = 50,
                                     feature_col = "interaction_name",
                                     pval_threshold = 0.05,
                                     feature_key = c("interaction_name", "source_target"),
                                     ...) {
  feature_key <- match.arg(feature_key)
  pres <- lr_presence_matrix(
    data_list,
    feature_col = feature_col,
    pval_threshold = pval_threshold,
    feature_key = feature_key
  )
  if (!nrow(pres)) {
    stop("No features to plot")
  }
  rv <- apply(pres, 1, function(r) mean(r) * (1 - mean(r)))
  ord <- order(-rv, names(rv))
  keep <- head(ord, min(top_n, length(ord)))
  matrix_heatmap(pres[keep, , drop = FALSE], ...)
}


#' Boxplot of mean pairwise Jaccard per sample
#'
#' For each sample, computes the mean of off-diagonal \code{jaccard} entries from
#' \code{\link{sample_lr_dissimilarity}} (same summary as the former bar chart).
#' Plots a single cohort-level boxplot of those means, with jittered points and
#' optional sample labels to show each sample's value.
#'
#' @param diss List returned by \code{sample_lr_dissimilarity()} (must contain
#'   a square \code{jaccard} matrix with at least two samples).
#' @param title Plot title.
#' @param label_samples If TRUE (default), draw sample names with \code{geom_text}
#'   at the same jitter positions as the points. If FALSE, jittered points are
#'   still drawn; only the text labels are omitted.
#' @param box_fill Fill color for the box (default \code{"#B8D4F0"}, light blue).
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{sample_lr_dissimilarity}}, \code{\link{plot_lr_sample_comparison}}
#'
#' @export
#' @import ggplot2
plot_lr_mean_jaccard_boxplot <- function(diss,
                                         title = NULL,
                                         label_samples = TRUE,
                                         box_fill = "#B8D4F0") {
  if (!is.list(diss) || !"jaccard" %in% names(diss)) {
    stop("diss must be output from sample_lr_dissimilarity() with element 'jaccard'")
  }
  J <- diss$jaccard
  if (!is.matrix(J) || length(dim(J)) != 2L || nrow(J) != ncol(J)) {
    stop("diss$jaccard must be a square matrix")
  }
  n <- ncol(J)
  if (n < 2L) {
    stop("At least two samples are required for mean off-diagonal Jaccard")
  }
  sn <- colnames(J)
  if (is.null(sn)) {
    sn <- paste0("Sample_", seq_len(n))
  }
  mj <- vapply(seq_len(n), function(i) mean(J[i, -i, drop = TRUE]), NA_real_)
  df <- data.frame(
    sample = sn,
    mean_jaccard = mj,
    stringsAsFactors = FALSE
  )

  if (is.null(title)) {
    title <- "Mean Jaccard by sample"
  }

  pos <- ggplot2::position_jitter(width = 0.12, seed = 1L)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = "", y = .data$mean_jaccard)) +
    ggplot2::geom_boxplot(
      width = 0.35,
      fill = box_fill,
      outlier.shape = NA
    ) +
    ggplot2::geom_point(
      size = 2.5,
      alpha = 0.85,
      color = "steelblue",
      position = pos
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Mean Jaccard vs other samples",
      title = title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  if (isTRUE(label_samples)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$sample),
      position = pos,
      vjust = -1.1,
      size = 2.8
    )
  }

  p
}


#' Distribution of heterogeneity metrics across features
#'
#' Histogram (or boxplot) of a column from \code{\link{calculate_heterogeneity}}
#' output, summarizing how many features are rare vs shared across samples.
#'
#' @param heterogeneity_results Data frame from \code{calculate_heterogeneity()}.
#' @param x Column to plot: \code{frequency} or \code{specificity_index}.
#' @param geom \code{"histogram"} (default) or \code{"boxplot"}.
#' @param bins Number of bins for continuous \code{x} when \code{geom = "histogram"}.
#' @param title Plot title; short defaults use the column name when \code{NULL}.
#' @param point_jitter If TRUE (default) and \code{geom = "boxplot"}, overlay
#'   horizontally jittered points (same idea as \code{\link{plot_lr_mean_jaccard_boxplot}}).
#'   Ignored when \code{geom = "histogram"}.
#' @param jitter_width Horizontal jitter width for \code{geom = "boxplot"} (default 0.05).
#' @param point_size Point size for jittered points when \code{point_jitter} is TRUE (default 1.35).
#' @param point_color Color for jittered points when \code{point_jitter} is TRUE.
#' @param box_fill Fill color for \code{geom = "boxplot"} (default \code{"#B8D4F0"}, light blue).
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{calculate_heterogeneity}}, \code{\link{plot_heterogeneity}}
#'
#' @export
#' @import ggplot2
plot_heterogeneity_distribution <- function(heterogeneity_results,
                                            x = c("frequency", "specificity_index"),
                                            geom = c("histogram", "boxplot"),
                                            bins = 30L,
                                            title = NULL,
                                            point_jitter = TRUE,
                                            jitter_width = 0.05,
                                            point_size = 1.35,
                                            point_color = "steelblue",
                                            box_fill = "#B8D4F0") {
  x <- match.arg(x)
  geom <- match.arg(geom)

  if (!is.data.frame(heterogeneity_results)) {
    stop("heterogeneity_results must be a data frame")
  }
  if (nrow(heterogeneity_results) == 0L) {
    stop("heterogeneity_results is empty")
  }
  if (!x %in% colnames(heterogeneity_results)) {
    stop("Column '", x, "' not found in heterogeneity_results")
  }

  plot_data <- heterogeneity_results
  plot_data[[x]] <- as.numeric(plot_data[[x]])

  if (is.null(title)) {
    title <- switch(
      x,
      frequency = "Frequency",
      specificity_index = "Specificity",
      paste0("Distribution: ", x)
    )
    if (geom == "histogram") {
      title <- paste0(title, " (features)")
    }
  }

  if (geom == "boxplot") {
    pos <- ggplot2::position_jitter(width = jitter_width, height = 0, seed = 1L)
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data[[x]], x = "")) +
      ggplot2::labs(x = NULL, y = x, title = title) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    if (isTRUE(point_jitter)) {
      p <- p + ggplot2::geom_point(
        color = point_color,
        alpha = 0.65,
        size = point_size,
        position = pos
      )
    }
    p <- p + ggplot2::geom_boxplot(
      fill = box_fill,
      outlier.shape = NA,
      width = 0.35
    )
    return(p)
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x]])) +
    ggplot2::geom_histogram(fill = "steelblue", color = "white", bins = bins) +
    ggplot2::labs(x = x, y = "Count", title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
}


#' Plot Heterogeneity Results (scatter)
#'
#' Scatter plot of mean \code{value_col} (e.g. mean p-value) versus coefficient of
#' variation or cohort frequency across samples where the feature is significant.
#' Requires a \code{mean_value} column from \code{\link{calculate_heterogeneity}} with
#' \code{method = "cv"}, \code{method = "range"}, or \code{method = "all"}.
#'
#' @param heterogeneity_results Data frame output from \code{calculate_heterogeneity()}
#'   with \code{mean_value} (and typically \code{cv}) columns.
#' @param top_n Integer specifying number of features to display (first \code{top_n} rows
#'   of \code{heterogeneity_results} in their given order). Re-run
#'   \code{\link{calculate_heterogeneity}} with a different \code{sort_by} to change order.
#' @param color_by Character string specifying column to use for point color.
#'   Default is \code{"frequency"}. Set to \code{NULL} for uniform color.
#' @param title Character string for plot title.
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' lr_data <- loadingLR_cellchat(sample = sample_names,
#'                               directory = data_dir,
#'                               value = "pval")
#' het_results <- calculate_heterogeneity(lr_data, method = "all", sort_by = "specificity")
#' p <- plot_heterogeneity(het_results, top_n = 50, color_by = "specificity_index")
#' }
#'
#' @export
#' @import ggplot2
plot_heterogeneity <- function(heterogeneity_results,
                               top_n = 20,
                               color_by = "frequency",
                               title = "Feature Heterogeneity Across Samples") {
  if (!is.data.frame(heterogeneity_results)) {
    stop("heterogeneity_results must be a data frame")
  }

  if (nrow(heterogeneity_results) == 0) {
    stop("heterogeneity_results is empty")
  }

  plot_data <- head(heterogeneity_results, top_n)

  if (!"mean_value" %in% colnames(plot_data)) {
    stop(
      "Scatter plot requires 'mean_value'. Use method = 'cv', 'range', or 'all' in calculate_heterogeneity()"
    )
  }

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

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$mean_value, y = .data[[y_var]])) +
    ggplot2::geom_point(
      size = 3,
      alpha = 0.7,
      ggplot2::aes(color = if (!is.null(color_by) && color_by %in% colnames(plot_data)) {
        .data[[color_by]]
      } else {
        NULL
      })
    ) +
    ggplot2::labs(
      x = "Mean Value",
      y = y_label,
      title = title,
      color = color_by
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "right"
    )

  if (!is.null(color_by) && color_by %in% colnames(plot_data)) {
    p <- p + ggplot2::scale_color_gradient(low = "blue", high = "red")
  }

  p
}
