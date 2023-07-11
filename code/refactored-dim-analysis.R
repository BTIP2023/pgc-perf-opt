library(tidyverse)
library(umap)
library(plotly)
library(htmlwidgets)
library(factoextra)
library(scales)
library(Rtsne)
library(webshot)

# Set paths
data_path <- "./data/kmers/"
results_path <- "./results/dim-reduce/R/"

# Loading all the data files
data3 <- read.csv(paste0(data_path, "kmer_3.csv"))
data5 <- read.csv(paste0(data_path, "kmer_5.csv"))
data7 <- read.csv(paste0(data_path, "kmer_7.csv"))

# Define t-SNE and UMAP parameters
tsne_seed <- 10
tsne_dims <- 2
tsne_perplexity <- 40
tsne_max_iter <- 300
umap_seed <- 10
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1

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

process_data <- function(data) {
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col - 1)]
  target <- data$variant
  
  non_zero_var_cols <- apply(X, 2, var) > 0
  X <- X[, non_zero_var_cols]
  
  if (ncol(X) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  x <- scale(X)
  
  return(list(x = x, target = target))
}

pc_exp <- function(data, color, label, symbol) {
  slice_col <- which(colnames(data) == 'strain')
  X <- data[, 2:(slice_col - 1)]
  target <- data$variant
  
  non_zero_var_cols <- apply(X, 2, var) > 0
  X <- X[, non_zero_var_cols]
  
  if (ncol(X) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  x <- scale(X)
  
  pca <- prcomp(x, center = TRUE, scale = TRUE)
  PC <- as.data.frame(pca$x[, 1:10])
  
  explained_variances <- summary(pca)$importance[2, 1:10] / sum(summary(pca)$importance[2, ])
  
  PC_values <- seq_along(explained_variances)
  prop_var <- explained_variances
  
  df <- data.frame(PC_values, prop_var, color = I(color), label = label, symbol = symbol)
  
  return(df)
}

pc_exp_plot <- function(x_lst) {
  df1 <- pc_exp(data3, '#9400D3', 'k=3', 'circle')
  df2 <- pc_exp(data5, '#FF0000', 'k=5', 'star')
  df3 <- pc_exp(data7, '#0000FF', 'k=7', 'x')
  
  combined_df <- bind_rows(df1, df2, df3)
  
  fig <- plot_ly(data = combined_df, x = ~PC_values, y = ~prop_var, type = "scatter", mode = "lines+markers",
                 color = ~color, name = ~label, symbol = ~symbol, symbols = c("circle", "star", "x"),
                 marker = list(size = 10)) %>% layout(xaxis = list(title = "Number of Principal Components"),
                                                      yaxis = list(title = "Proportion of Variance Explained"))
  
  save_dim_reduction_plot("pca", "combined", fig)
}

pca_fn <- function(x, target) {
  pca_reduce <- prcomp(x, center = TRUE, scale. = TRUE)
  retPC <- as.data.frame(pca_reduce$x[, 1:50])
  
  PC <- as.data.frame(pca_reduce$x[, 1:2])
  colnames(PC) <- c('Principal Component 1', 'Principal Component 2')
  FDF <- cbind(PC, variant = target)
  
  explained_variance_1 <- round(summary(pca_reduce)$importance[2, 1] / sum(summary(pca_reduce)$importance[2, ]), 4)
  explained_variance_2 <- round(summary(pca_reduce)$importance[2, 2] / sum(summary(pca_reduce)$importance[2, ]), 4)
  
  return(list(pca_results = retPC, FDF = FDF, variant = target, explained_variance_1 = explained_variance_1, explained_variance_2 = explained_variance_2))
}

pca_plot <- function(FDF, variant, k, explained_variance_1, explained_variance_2) {
  fig <- plot_ly(FDF, x = ~`Principal Component 1`, y = ~`Principal Component 2`, color = ~variant) %>%
    add_markers(size = 4, alpha = 0.5) %>%
    layout(
      xaxis = list(title = paste("PC 1 (", round(explained_variance_1*100, 2), "% explained variance)")),
      yaxis = list(title = paste("PC 2 (", round(explained_variance_2*100, 2), "% explained variance)")),
      legend = list(orientation = "h", x = 0.5, y = 1)
    )
  
  save_dim_reduction_plot("pca", k, fig)
}

tsne_fn <- function(X, pca_results) {
  set.seed(tsne_seed)
  tsne_reduce <- Rtsne(pca_results, dims = tsne_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter, pca = FALSE)
  
  return(tsne_reduce)
}

tsne_plot <- function(tsne_reduce, target, k) {
  df <- data.frame(X1 = tsne_reduce$Y[, 1], X2 = tsne_reduce$Y[, 2], target = target)
  fig <- plot_ly(df, x = ~X1, y = ~X2, color = ~target, type = "scatter", mode = "markers") %>%
    layout(xaxis = list(title = "TSNE-2D-1"),
           yaxis = list(title = "TSNE-2D-2"),
           colorway = "Dark2")
  
  save_dim_reduction_plot("tsne", k, fig)
} 

umap_fn <- function(x) {
  umap_reduce <- umap(x, n_neighbors = umap_n_neighbors, metric = umap_metric, min_dist = umap_min_dist, seed = umap_seed)
  return(umap_reduce)
}

umap_plot <- function(umap_reduce, target, k) {
  emb <- umap_reduce$layout
  
  X_o <- emb[, 1]
  Y_o <- emb[, 2]
  
  fig <- plot_ly(data = as.data.frame(emb), x = ~X_o, y = ~Y_o, color = ~target, type = "scatter", mode = "markers") %>%
    layout(xaxis = list(title = "UMAP_1"),
           yaxis = list(title = "UMAP_2"),
           colorway = "Dark2")
  
  save_dim_reduction_plot("umap", k, fig)
}

main_fn <- function(data3, data5, data7) {
  data_list <- list(data3, data5, data7)
  dimensions <- c(3, 5, 7)
  
  for (i in 1:length(data_list)) {
    data <- data_list[[i]]
    k <- dimensions[i]
    results <- process_data(data)
    
    pc_exp_plot()
    
    pca_reduce <- pca_fn(results$x, results$target)
    pca_plot(pca_reduce$FDF, pca_reduce$target, k, results$explained_variance_1, results$explained_variance_2)
    
    tsne_reduce <- tsne_fn(data, pca_reduce$pca_results)
    tsne_plot(tsne_reduce, results$target, k)
    
    umap_reduce <- umap_fn(results$x)
    umap_plot(umap_reduce, results$target, k)
  }
}

main_fn(data3, data5, data7)