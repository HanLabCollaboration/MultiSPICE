#' Create Pathway Gene Sets for Limma Analysis
#'
#' Generate pathway gene sets from CellChat database for use in limma gene set enrichment analysis.
#' Each pathway contains a list of ligand-receptor (LR) interaction pairs that can be used
#' with limma's fry() or camera() functions.
#'
#' @param species Character string specifying the species. Options are "human" (default) or "mouse".
#'
#' @return A named list where each element corresponds to a pathway and contains a character vector
#'   of interaction names (LR pairs) belonging to that pathway.
#'
#' @details
#' This function extracts the CellChat interaction database and organizes ligand-receptor pairs
#' by their pathway annotations. The resulting gene sets can be used directly with limma's
#' gene set testing functions (fry, camera) to test for coordinated changes in pathway activity.
#'
#' @examples
#' \dontrun{
#' # Create pathway gene sets for human
#' pathway_sets <- create_pathway_geneset(species = "human")
#'
#' # Use with limma fry
#' library(limma)
#' pathway_index <- lapply(pathway_sets, function(lr_pairs) {
#'   which(rownames(expression_matrix) %in% lr_pairs)
#' })
#' fry_results <- fry(expression_matrix, index = pathway_index, design = design, contrast = 2)
#' }
#'
#' @export
#' @importFrom CellChat CellChatDB.human CellChatDB.mouse

create_pathway_geneset <- function(species = "human") {
  # Load CellChat database
  if (species == "human") {
    CellChatDB <- CellChat::CellChatDB.human
  } else if (species == "mouse") {
    CellChatDB <- CellChat::CellChatDB.mouse
  } else {
    stop("Species must be 'human' or 'mouse'")
  }

  # Extract interaction database
  interaction_db <- CellChatDB$interaction

  # Create gene sets for each pathway using interaction_name (LR pairs)
  pathway_gene_sets <- lapply(split(interaction_db, interaction_db$pathway_name), function(pathway_data) {
    # Use interaction_name as the "genes" for each pathway
    lr_pairs <- pathway_data$interaction_name

    # Get unique LR pairs
    lr_pairs <- unique(lr_pairs)

    return(lr_pairs)
  })

  # Print summary
  cat("Created pathway gene sets for", length(pathway_gene_sets), "pathways\n")
  cat("Number of LR pairs per pathway:\n")
  print(summary(sapply(pathway_gene_sets, length)))

  return(pathway_gene_sets)
}


#' Plot Pathway Enrichment Results
#'
#' Create a barplot to visualize pathway enrichment results from limma's fry() or camera() analysis.
#' The plot displays pathways ranked by significance with bars colored by direction of enrichment.
#'
#' @param enrichment_results A data frame containing pathway enrichment results from limma's fry() or camera().
#'   Must contain columns: "FDR" (or "adj.P.Val"), "PValue", and "Direction".
#'   Row names should be pathway names.
#' @param fdr_threshold Numeric value specifying the FDR threshold for significance. Default is 0.05.
#' @param top_n Integer specifying the number of top pathways to show if no significant pathways are found.
#'   Default is 20.
#' @param colors Named character vector of colors for "Up", "Down", and "Mixed" directions.
#'   Default is c("Up" = "#377EB8", "Down" = "#E41A1C", "Mixed" = "purple").
#' @param title Character string for the plot title. Default is "Pathway Enrichment Analysis".
#' @param subtitle Character string for the plot subtitle. Default is "Up = enriched in treatment, Down = enriched in control".
#' @param x_label Character string for the x-axis label. Default is "-log10(FDR)".
#' @param y_label Character string for the y-axis label. Default is "Pathway".
#' @param text_size Numeric value for the pathway label text size. Default is 10.
#'
#' @return A ggplot2 object.
#'
#' @details
#' This function visualizes pathway enrichment analysis results by creating a horizontal barplot.
#' Pathways are ranked by significance (FDR), and bar colors indicate the direction of enrichment:
#' - "Up": Pathway is upregulated in the treatment/contrast group
#' - "Down": Pathway is downregulated in the treatment/contrast group
#' - "Mixed": Pathway contains both up and down regulated LR pairs
#'
#' If no pathways meet the FDR threshold, the function will display the top N pathways by p-value.
#'
#' @examples
#' \dontrun{
#' # Run pathway enrichment analysis
#' pathway_sets <- create_pathway_geneset(species = "human")
#' pathway_index <- lapply(pathway_sets, function(lr_pairs) {
#'   which(rownames(expr_matrix) %in% lr_pairs)
#' })
#' fry_results <- fry(expr_matrix, index = pathway_index, design = design, contrast = 2)
#'
#' # Plot results
#' p <- plot_pathway_enrichment(fry_results)
#' print(p)
#'
#' # Customize plot
#' p_custom <- plot_pathway_enrichment(
#'   fry_results,
#'   fdr_threshold = 0.01,
#'   title = "Pathway Enrichment: mCRC vs priCRC",
#'   subtitle = "Up = enriched in mCRC, Down = enriched in priCRC"
#' )
#' }
#'
#' @export
#' @import ggplot2
plot_pathway_enrichment <- function(enrichment_results,
                                    fdr_threshold = 0.05,
                                    top_n = 20,
                                    colors = c("Up" = "#377EB8", "Down" = "#E41A1C", "Mixed" = "purple"),
                                    title = "Pathway Enrichment Analysis",
                                    subtitle = "Up = enriched in treatment, Down = enriched in control",
                                    x_label = "-log10(FDR)",
                                    y_label = "Pathway",
                                    text_size = 10) {

  # Check if required columns exist
  if (!("Direction" %in% colnames(enrichment_results))) {
    stop("enrichment_results must contain a 'Direction' column")
  }

  # Handle FDR column name (could be FDR or adj.P.Val)
  if ("FDR" %in% colnames(enrichment_results)) {
    fdr_col <- "FDR"
  } else if ("adj.P.Val" %in% colnames(enrichment_results)) {
    fdr_col <- "adj.P.Val"
  } else {
    stop("enrichment_results must contain either 'FDR' or 'adj.P.Val' column")
  }

  # Prepare data for plotting
  plot_data <- enrichment_results
  plot_data$pathway <- rownames(plot_data)
  plot_data$log10FDR <- -log10(plot_data[[fdr_col]])

  # Filter significant pathways
  plot_data_sig <- plot_data[plot_data[[fdr_col]] < fdr_threshold, ]

  # If no significant pathways, show top N
  if (nrow(plot_data_sig) == 0) {
    warning(paste("No pathways with", fdr_col, "<", fdr_threshold,
                  ". Showing top", top_n, "pathways instead."))
    plot_data_sig <- plot_data[order(plot_data[[fdr_col]]), ][1:min(top_n, nrow(plot_data)), ]
    subtitle <- paste(subtitle, "\n(No significant pathways at FDR < ", fdr_threshold, ")", sep = "")
  }

  # Order by significance for plotting
  plot_data_sig <- plot_data_sig[order(plot_data_sig$log10FDR), ]
  plot_data_sig$pathway <- factor(plot_data_sig$pathway, levels = plot_data_sig$pathway)

  # Create barplot
  p <- ggplot2::ggplot(plot_data_sig, ggplot2::aes(x = log10FDR, y = pathway, fill = Direction)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_vline(xintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = x_label,
                  y = y_label,
                  title = title,
                  subtitle = subtitle) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = text_size),
      axis.title = ggplot2::element_text(size = text_size + 2),
      plot.title = ggplot2::element_text(size = text_size + 4, face = "bold"),
      legend.position = "right"
    )

  return(p)
}
