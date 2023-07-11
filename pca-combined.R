library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(plotly)
library(htmlwidgets)
library(webshot)
library(cowplot)

# Loading all the data files
data3 <- read.csv('./data/kmers/kmer_3.csv')
data5 <- read.csv('./data/kmers/kmer_5.csv')
data7 <- read.csv('./data/kmers/kmer_7.csv')

pcaPlot <- function(data, color, label, symbol) {
  slice_col <- which(names(data) == 'strain')
  X <- data[, 2:(slice_col-1)]
  target <- data$variant
  
  # Exclude columns with zero variance
  non_zero_var_cols <- apply(X, 2, var) > 0
  X <- X[, non_zero_var_cols]
  
  if (ncol(X) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  x <- scale(X)
  
  pca <- prcomp(x, center = TRUE, scale = TRUE)
  PC <- as.data.frame(pca$x[, 1:10])
  
  explained_variances <- summary(pca)$importance[2, 1:10] / sum(summary(pca)$importance[2, ])
  
  cat('Variance of each component:', explained_variances, '\n')
  
  PC_values <- seq_along(explained_variances)
  prop_var <- explained_variances
  
  df <- data.frame(PC_values, prop_var, color = I(color), label = label, symbol = symbol)
  
  return(df)
}

df1 <- pcaPlot(data3, '#9400D3', 'k=3', 'circle')
df2 <- pcaPlot(data5, '#FF0000', 'k=5', 'star')
df3 <- pcaPlot(data7, '#0000FF', 'k=7', 'x')

combined_df <- bind_rows(df1, df2, df3)

fig <- plot_ly(data = combined_df, x = ~PC_values, y = ~prop_var, type = "scatter", mode = "lines+markers",
               color = ~color, name = ~label, symbol = ~symbol, symbols = c("circle", "star", "x"),
               marker = list(size = 10))

fig <- fig %>% layout(xaxis = list(title = "Number of Principal Components"),
                      yaxis = list(title = "Proportion of Variance Explained"))

# Save as PNG using webshot and htmlwidgets
pngFile <- "./results/dim-reduce/R/combined-PCA.png"
htmlwidgets::saveWidget(fig, "temp.html")
webshot::webshot("temp.html", pngFile)
file.remove("temp.html")

# Save as HTML
htmlFile <- "./results/dim-reduce/R/combined-PCA.html"
htmlwidgets::saveWidget(fig, file = htmlFile, selfcontained = TRUE)