library(tsne)
library(plotly)
library(htmlwidgets)
library(webshot)

# Loading all the data files
data3 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_3.csv')
data5 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_5.csv')
data7 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_7.csv')

tsne_fn <- function(data, k) {
  data <- data[order(data$variant), ]
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col - 1)]
  target <- data$variant
  set.seed(10)
  tsne_results <- tsne(X, initial_dims = 2, perplexity = 40, max_iter = 300)
  df <- data.frame(X1 = tsne_results[, 1], X2 = tsne_results[, 2], target = target)
  p <- plot_ly(df, x = ~X1, y = ~X2, color = ~target, type = "scatter", mode = "markers") %>%
    layout(xaxis = list(title = "TSNE-2D-1"),
           yaxis = list(title = "TSNE-2D-2"),
           colorway = "Dark2")
  
  # Save as PNG using webshot and htmlwidgets
  pngFile <- paste0("C:/Users/Brylle/Desktop/R Port/dim-reduce/tsne-", k, ".png")
  htmlwidgets::saveWidget(p, "temp.html")
  webshot::webshot("temp.html", pngFile)
  file.remove("temp.html")
  
  # Save as HTML
  htmlFile <- paste0("C:/Users/Brylle/Desktop/R Port/dim-reduce/tsne-", k, ".html")
  htmlwidgets::saveWidget(p, file = htmlFile, selfcontained = TRUE)
}

tsne_fn(data3, 3)
tsne_fn(data5, 5)
# tsne_fn(data7, 7) # takes too long/too much memory