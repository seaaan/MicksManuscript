smallPointSize <- function() 1
pointSize <- function() 2
lineSize <- function() 0.5
thinLineSize <- function() 0.25
axisLineColor <- function() "gray50"
colorBlindPalette <- function() c("#009E73", "#D55E00", "#0072B2", "#CC79A7")

highlight_line_color <- function() "grey50"

hline <- function(x) geom_hline(yintercept = x, color = highlight_line_color())
vline <- function(x) geom_vline(xintercept = x, color = highlight_line_color())

theme_pub <- function(...) {
   theme(
      # backgrounds
      strip.background = element_blank(),
      panel.background = element_blank(),
      legend.key = element_blank(),
      
      # lines
      axis.line = element_line(color = axisLineColor(), size = lineSize()),
      axis.ticks = element_line(color = axisLineColor(), size = lineSize()),
      panel.grid.major = element_line(color = "grey90", size = thinLineSize()),
      panel.grid.minor = element_blank(),
      
      # fonts
      title = element_text(size = 12),
      axis.title = element_text(size = 10),
      strip.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "plain"),
      legend.text = element_text(size = 8, color = "grey50"),
      axis.text = element_text(size = 8),
      
      # remove whitespace at edges of plots
      plot.margin = grid::unit(c(0, 0, 0, 0), "mm"), 
      
      # center plot titles
      plot.title = element_text(hjust = 0.5),
      ...
   ) 
}

y_axis_all_facets <- function() {
   annotate("segment", color = axisLineColor(), 
      size = lineSize(),
      x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
}

x_axis_all_facets <- function() {
   annotate("segment", color = axisLineColor(), 
      size = lineSize(),
      x = -Inf, xend = Inf, y = -Inf, yend = -Inf)
}

add_white_space <- function(top = 0, right = 0, bottom = 0, left = 0) {
   theme(plot.margin = grid::unit(c(top, right, bottom, left), "mm"))
}
