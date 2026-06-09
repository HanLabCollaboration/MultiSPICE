#' @title Plot rank plot for LR interactions
#' @description Plot rank plot for LR interactions
#' @param dat Data frame containing the LR interaction information
#' @param colors Color vector for cell types
#' @param type Column name to use for coloring points (e.g., "sender_celltype")
#' @param label_threshold Threshold for labeling points on the plot
#' @param do_raster Logical indicating whether to use rasterized points for large datasets
#' @param main Title for the plot
#' @param label_col Column name containing interaction names for labels (default: "interaction_name")
#' @return A ggplot object representing the rank plot

############ Plot rank plot function ###################
LRrankplot <- function(dat, colors = col_vector, type = "sender_celltype",
                       label_threshold = 8, do_raster = TRUE, main = NULL, 
                       label_col = "interaction_name") {
  plot_data <- dat
  plot_data <- plot_data[order(plot_data[[type]], plot_data[[label_col]]), ]
  plot_data$position <- seq(1:nrow(plot_data))

  if (do_raster == TRUE){
    gra <- ggplot(plot_data, aes(x=position, y=rank, color=.data[[type]])) + geom_point_rast()
  } else {
    gra <- ggplot(plot_data, aes(x=position, y=rank, color=.data[[type]])) + geom_point()
  }

  gra + scale_color_manual(values=colors) + theme_bw() + xlab("Celltype") + ggtitle(main) +
    geom_hline(yintercept=label_threshold, linetype="dashed", color="darkred", linewidth=1) +
    geom_label_repel(data=subset(plot_data, rank > label_threshold), aes(x=position, y=rank, label=.data[[label_col]]), size = 2,
                     nudge_y = 0.1, max.overlaps = 100, show.legend=FALSE) + ylab("Number of siginificant samples") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13),
          legend.text = element_text(size = 13), legend.title = element_text(size = 14),
          panel.grid = element_blank(), title = element_text(face="bold"), axis.text.x = element_blank())
}





