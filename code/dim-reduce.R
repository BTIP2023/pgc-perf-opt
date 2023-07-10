library(tidyverse)
library(umap)
library(plotly)
library(htmlwidgets)
library(factoextra)
library(scales)
library(Rtsne)
library(webshot)

# Loading all the data files
data3 <- read.csv('./data/kmers/kmer_3.csv')
data5 <- read.csv('./data/kmers/kmer_5.csv')
data7 <- read.csv('./data/kmers/kmer_7.csv')

save_dim_reduction_plot <- function(method, k, p) {
  # Save as PNG using webshot and htmlwidgets
  pngFile <- paste0("./results/dim-reduce/R/", method, "-", k, ".png")
  htmlwidgets::saveWidget(p, "temp.html")
  webshot::webshot("temp.html", pngFile)
  file.remove("temp.html")
  
  # Save as HTML
  htmlFile <- paste0("./results/dim-reduce/R/", method, "-", k, ".html")
  htmlwidgets::saveWidget(p, file = htmlFile, selfcontained = TRUE)
}

pca_fn <- function(data, k) {
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
  
  pca <- prcomp(x, center = TRUE, scale = TRUE)
  PC <- as.data.frame(pca$x[, 1:2])
  colnames(PC) <- c('Principal Component 1', 'Principal Component 2')
  FDF <- cbind(PC, variant = target)
  
  explained_variance_1 <- round(summary(pca)$importance[2, 1] / sum(summary(pca)$importance[2, ]), 4)
  explained_variance_2 <- round(summary(pca)$importance[2, 2] / sum(summary(pca)$importance[2, ]), 4)
  
  cat('Explained variance:', explained_variance_1, explained_variance_2, '\n')
  
  fig <- plot_ly(FDF, x = ~`Principal Component 1`, y = ~`Principal Component 2`, color = ~variant) %>%
    add_markers(size = 4, alpha = 0.5) %>%
    layout(
      xaxis = list(title = paste("PC 1 (", round(explained_variance_1*100, 2), "% explained variance)")),
      yaxis = list(
        title = paste("PC 2 (", round(explained_variance_2*100, 2), "% explained variance)"),
        autorange = "reversed"
      ),
      legend = list(orientation = "h", x = 0.5, y = 1)
    )
  
  save_dim_reduction_plot("pca", k, fig)
}

tsne_fn <- function(data, k) {
  data <- data[order(data$variant), ]
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col - 1)]
  target <- data$variant
  
  set.seed(10)
  tsne_results <- Rtsne(X, dims = 2, perplexity = 40, max_iter = 300)
  df <- data.frame(X1 = tsne_results$Y[, 1], X2 = tsne_results$Y[, 2], target = target)
  p <- plot_ly(df, x = ~X1, y = ~X2, color = ~target, type = "scatter", mode = "markers") %>%
    layout(xaxis = list(title = "TSNE-2D-1"),
           yaxis = list(title = "TSNE-2D-2"),
           colorway = "Dark2")
  
  save_dim_reduction_plot("tsne", k, p)
}

umap_fn <- function(data, k) {
  # Sort data so figure labels are sorted
  data <- data[order(data$variant), ]
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col-1)]
  target <- data$variant
  
  # Exclude columns with zero variance
  non_zero_var_cols <- apply(X, 2, var) > 0
  X <- X[, non_zero_var_cols]
  
  x <- scale(X)
  
  umap_reduce <- umap(x, n_neighbors = 15, metric = "euclidean", min_dist = 0.1, seed = 10)
  emb <- umap_reduce$layout
  
  X_o <- emb[, 1]
  Y_o <- emb[, 2]
  
  umap_2d <- plot_ly(data = as.data.frame(emb), x = ~X_o, y = ~Y_o, color = ~target, type = "scatter", mode = "markers") %>%
    layout(xaxis = list(title = "UMAP_1"),
           yaxis = list(title = "UMAP_2"),
           colorway = "Dark2")
  
  save_dim_reduction_plot("umap", k, umap_2d)
}

data_list <- list(data3, data5, data7)
dimensions <- c(3, 5, 7)

for (i in 1:length(data_list)) {
  pca_fn(data_list[[i]], dimensions[i])
  umap_fn(data_list[[i]], dimensions[i])
  tsne_fn(data_list[[i]], dimensions[i])
}