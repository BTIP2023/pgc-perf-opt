library(tidyverse)
library(umap)
library(plotly)
library(htmlwidgets)
library(webshot)

# Loading all the data files
data3 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_3.csv')
data5 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_5.csv')
data7 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_7.csv')

umap_fn <- function(data, k) {
  # Sort data so figure labels are sorted
  data <- data[order(data$variant), ]
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col-1)]
  target <- data$variant
  
  umap_reduce <- umap(X, n_neighbors = 15, metric = "euclidean", min_dist = 0.1)
  emb <- umap_reduce$layout
  
  X_o <- emb[, 1]
  Y_o <- emb[, 2]
  
  umap_2d <- plot_ly(data = as.data.frame(emb), x = ~X_o, y = ~Y_o, color = ~target, type = "scatter", mode = "markers") %>%
    layout(xaxis = list(title = "UMAP_1"),
           yaxis = list(title = "UMAP_2"),
           colorway = "Dark2")
  
  # Save as PNG using webshot and htmlwidgets
  pngFile <- paste0("C:/Users/Brylle/Desktop/R Port/dim-reduce/umap-", k, ".png")
  htmlwidgets::saveWidget(umap_2d, "temp.html")
  webshot::webshot("temp.html", pngFile)
  file.remove("temp.html")
  
  # Save as HTML
  htmlFile <- paste0("C:/Users/Brylle/Desktop/R Port/dim-reduce/umap-", k, ".html")
  htmlwidgets::saveWidget(umap_2d, file = htmlFile, selfcontained = TRUE)
}

umap_fn(data3, 3)
umap_fn(data5, 5)
umap_fn(data7, 7)