library(ggplot2)
library(grid)
library(RColorBrewer)

split_print_legend <- function(plt) {
  print(plt + theme(legend.position = "none"))
  tmp <- ggplot_gtable(ggplot_build(plt)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  grid.newpage()
  grid.draw(legend) 
}

no_grays_palette <- c(
  brewer.pal(9, "Set1")[c(1:8)],
  brewer.pal(8, "Set2")[c(1:7)],
  brewer.pal(12, "Set3")[c(1:8)],
  brewer.pal(12, "Set3")[c(10:12)],
  brewer.pal(8, "Dark2")[1:7]
)
