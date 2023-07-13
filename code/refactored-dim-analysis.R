# Install and load the required packages
if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(
  tidyverse,
  umap,
  plotly,
  htmlwidgets,
  factoextra,
  scales,
  Rtsne,
  RColorBrewer
)

# -----Functions-----

# Function for saving plots as PNG and HTML
save_plot <- function(method, k, p) {
  # Save as PNG
  filename <- paste0(method, "-", k, ".png")
  path <- "./results/dim-reduce/R/"
  ggsave(filename, p, path, device = "png", dpi = 300)
  
  # Convert ggplot object to ggplotly
  p <- ggplotly(p)
  
  # Save as HTML
  htmlFile <- paste0("./results/dim-reduce/R/", method, "-", k, ".html")
  htmlwidgets::saveWidget(p, file = htmlFile, selfcontained = TRUE)
}

# Function for pre-processing and scaling of data 
process_data <- function(data) {
  # Drop metadata, NAs, etc.
  # Include k-mers
  
  # Determine the columns to use
  slice_col <- which(colnames(data) == "strain")
  X <- data[, 2:(slice_col - 1)]
  target <- data$variant
  
  # Check for columns that have zero variance for PCA
  non_zero_var_cols <- apply(X, 2, var) > 0
  new_X <- X[, non_zero_var_cols]
  
  if (ncol(new_X) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  # Scale data
  x <- scale(new_X) # can be optional
  
  return(list(X = X, x = x, target = target))
}

# Function that performs PCA
pca_fn <- function(x) {
  pca_df <- prcomp(x, center = TRUE, scale. = TRUE)
  return(pca_df)
}

# Function that plots PCA results
pca_plot <- function(pca_df, target, k) {
  # Plot two principal components
  PC <- as.data.frame(pca_df$x[, 1:2])
  colnames(PC) <- c("Principal Component 1", "Principal Component 2")
  FDF <- cbind(PC, variant = target)
  
  # Create ggplot object
  p <- ggplot(FDF, aes(x = `Principal Component 1`, y = `Principal Component 2`, color = variant)) +
    geom_point(size = 4, alpha = 0.5) +
    xlab(paste("PC 1 (", round(summary(pca_df)$importance[2,1]  * 100, 2), "% explained variance)")) +
    ylab(paste("PC 2 (", round(summary(pca_df)$importance[2,2] * 100, 2), "% explained variance)")) +
    labs(color = "Label") +
    scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("pca", k, p)
}

# Function that performs t-SNE
tsne_fn <- function(pca_results) {
  set.seed(tsne_seed)
  tsne_df <- Rtsne(pca_results, dims = tsne_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter, check_duplicates = FALSE, pca = FALSE)
  return(tsne_df)
}

# Function that plots t-SNE results
tsne_plot <- function(tsne_df, target, k) {
  df <- data.frame(X1 = tsne_df$Y[, 1], X2 = tsne_df$Y[, 2], target = target)
  
  # Create ggplot object
  p <- ggplot(df, aes(x = X1, y = X2, color = target)) +
    geom_point() +
    xlab("TSNE-2D-1") +
    ylab("TSNE-2D-2") +
    labs(color = "Label") +
    scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("tsne", k, p)
}

# Function that performs UMAP
umap_fn <- function(x) {
  umap_df <- umap(x, n_neighbors = umap_n_neighbors, metric = umap_metric, min_dist = umap_min_dist, seed = umap_seed)
  return(umap_df)
}

# Function that plots UMAP results
umap_plot <- function(umap_df, target, k) {
  emb <- umap_df$layout
  
  X_o <- emb[, 1]
  Y_o <- emb[, 2]
  
  # Create ggplot object
  p <- ggplot(data = as.data.frame(emb), aes(x = X_o, y = Y_o, color = target)) +
    geom_point() +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    labs(color = "Label") +
    scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("umap", k, p)
}

# -----END of Functions-----

# -----START-----

# Set k for k-mer analysis
k <- 7

# Set paths
data_path <- "./data/archive/"
results_path <- "./results/dim-reduce/R/"

# Loading the data file
data <- read.csv(paste0(data_path, "kmer_", k, ".csv"))

# Define t-SNE and UMAP parameters
tsne_seed <- 10
tsne_dims <- 2
tsne_perplexity <- 40
tsne_max_iter <- 300
umap_seed <- 10
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1

# Process the data
data <- process_data(data)

# Perform PCA
pca_df <- pca_fn(data$x)

# Plot PCA results
pca_plot(pca_df, data$target, k)

# Perform t-SNE using PCA results
tsne_df <- tsne_fn(pca_df$x)

# Plot t-SNE results
tsne_plot(tsne_df, data$target, k)

# Perform UMAP
umap_df <- umap_fn(data$x)

# Plot UMAP results
umap_plot(umap_df, data$target, k)

# -----END-----