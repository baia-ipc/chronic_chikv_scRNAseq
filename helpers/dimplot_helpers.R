
#
# Create a dimplot with two groups highlighted
#
# Mandatory arguments:
#   so: Seurat object
#   grp1: first group to highlight
#   grp2: second group to highlight
#   grpcol: column in metadata to use for grouping
#
# Optional arguments:
#   group_names: names of the groups to highlight
#                [default: use grp1 and grp2]
#   colors_highlight: colors to use for highlighting
#                     [default: darkred and darkblue]
#   size_highlight: size of the highlighted points
#                   [default: 0.1]
#   base_color: color of the non-highlighted points
#               [default: grey]
#   pt_size: size of the non-highlighted points
#            [default: 0.1]
#   split_by: column in metadata to use for splitting
#             [default: no splitting]
#   title: title of the plot
#          [default: "UMAP plot"]
#   reduction_method: method used for dimensionality reduction
#                     [default: "umapharm"]
#   x_label: label for the x-axis
#            [default: "UMAP Harmony 1"]
#   y_label: label for the y-axis
#            [default: "UMAP Harmony 2"]
#   font_size: font size for the plot
#              [default: 10]
#   legend_position: position of the legend
#                    [default: "bottom"]
#
two_groups_dimplot <- function(so, grp1, grp2, grpcol,
                               group_names = NULL,
                               colors_highlight = c("darkred", "darkblue"),
                               size_highlight = 0.1,
                               base_color = "grey",
                               pt_size = 0.1,
                               split_by = NULL,
                               title = "UMAP plot",
                               reduction_method = "umapharm",
                               x_label = "UMAP Harmony 1",
                               y_label = "UMAP Harmony 2",
                               font_size = 10,
                               legend_position = "bottom") {

  highl1 <- WhichCells(so, cells = which(so@meta.data[[grpcol]] == grp1))
  highl2 <- WhichCells(so, cells = which(so@meta.data[[grpcol]] == grp2))

  if (is.null(group_names)) {
    group_names <- c(grp1, grp2)
  }

  cells_highlight <- setNames(list(highl1, highl2), group_names)

  plt <- DimPlot(so, reduction = reduction_method,
                 cells.highlight = cells_highlight,
                 sizes.highlight = size_highlight, pt.size = pt_size,
                 cols.highlight = colors_highlight, split.by = split_by,
                 cols = base_color) +
         ggtitle(title) +
         labs(x = x_label, y = y_label) +
         theme(legend.position = legend_position) +
         theme(text = element_text(size = font_size))

  return(plt)
}

