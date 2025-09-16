#' @title LRpiechart
#' @description Plot a pie chart for a given column in the data frame.
#' @param dat A data frame containing the data to be plotted.
#' @param colnum A string specifying the column name to be used for the pie chart.
#' @param colors A vector of colors to be used for the pie chart segments.
#' @param main A string specifying the title of the pie chart.
#' @return A ggplot object representing the pie chart.


LRpiechart <- function(dat, colnum = "celltype1", colors = col_vector, main = NULL){
  plot_data <- data.frame(table(dat[[colnum]]))
  # Get the positions
  plot_data <- plot_data %>%
    mutate(csum = rev(cumsum(rev(Freq))),
           pos = Freq/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Freq/2, pos))

  ggplot(plot_data, aes(x = "" , y = Freq, fill = Var1)) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_label_repel(data = plot_data,
                     aes(y = pos, label = Freq),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = colnum)) +
    theme_void() + scale_fill_manual(values = colors) + ggtitle(main)
}



