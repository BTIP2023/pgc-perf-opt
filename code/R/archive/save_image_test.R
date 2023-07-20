library(plotly)
# volcano is a numeric matrix that ships with R
fig <- plot_ly(z = ~volcano)
fig <- fig %>% add_surface()

fig

plotly::export(p = fig, #the graph to export
               file = "test.png") #the name and type of file (can be .png, .jpeg, etc.)
htmlwidgets::saveWidget(
  widget = fig, #the plotly object
  file = "figure.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)

library(plotly)

if (!require("processx")) install.packages("processx")
fig <- plot_ly(z = ~volcano) %>% add_surface()
orca(fig, "surface-plot.svg")
orca(fig, "surface-plot.png")

library(plotly)
save_image(fig, "results/test.png")

