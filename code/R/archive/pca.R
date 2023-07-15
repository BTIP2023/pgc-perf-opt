library(tidyverse)
library(plotly)
library(factoextra)
library(scales)
library(webshot)

# Loading all the data files
data3 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_3.csv')
data5 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_5.csv')
data7 <- read.csv('C:/Users/Brylle/Desktop/R Port/kmer_7.csv')

pca_plotly <- function(data, k) {
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col-1)]
  target <- data$variant
  
  # Exclude columns with zero variance
  non_zero_var_cols <- apply(X, 2, var) > 0
  X <- X[, non_zero_var_cols]
  
  if (ncol(X) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  x <- scale(X)
  
  pca <- prcomp(x, center = TRUE, scale. = TRUE)
  PC <- as.data.frame(pca$x[, 1:2])
  colnames(PC) <- c('Principal Component 1', 'Principal Component 2')
  FDF <- cbind(PC, variant = target)
  
  cat('Explained variance:', summary(pca)$importance[2, ] / sum(summary(pca)$importance[2, ]), '\n')
  
  fig <- plot_ly(FDF, x = ~`Principal Component 1`, y = ~`Principal Component 2`, color = ~variant) %>%
    add_markers(size = 4, alpha = 0.5) %>%
    layout(
      xaxis = list(title = paste("PC 1 (", percent(summary(pca)$importance[2, 1] / sum(summary(pca)$importance[2, ])), " explained variance)")),
      yaxis = list(title = paste("PC 2 (", percent(summary(pca)$importance[2, 2] / sum(summary(pca)$importance[2, ])), " explained variance)")),
      legend = list(orientation = "h", x = 0.5, y = 1)
    )
  
  # Save as PNG using webshot and htmlwidgets
  pngFile <- paste0("C:/Users/Brylle/Desktop/R Port/dim-reduce/pca-", k, ".png")
  htmlwidgets::saveWidget(fig, "temp.html")
  webshot::webshot("temp.html", pngFile)
  file.remove("temp.html")
  
  # Save as HTML
  htmlFile <- paste0("C:/Users/Brylle/Desktop/R Port/dim-reduce/pca-", k, ".html")
  htmlwidgets::saveWidget(fig, file = htmlFile, selfcontained = TRUE)
}

pca_plotly(data3, 3)
pca_plotly(data5, 5)
pca_plotly(data7, 7)