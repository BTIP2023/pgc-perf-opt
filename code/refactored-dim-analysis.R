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
  tsne,
  RColorBrewer,
  ggfortify
)

# -----Functions-----

# Function for saving plots as PNG and HTML
save_plot <- function(method, k, p) {
  # Save as PNG
  filename <- paste0(method, "-", k, ".png")
  ggsave(filename, p, results_path, device = "png", width=5, height=5, dpi = 300, bg = "white")
  
  # Convert ggplot object to ggplotly
  p <- ggplotly(p)
  
  # Save as HTML
  htmlFile <- paste0(results_path, method, "-", k, ".html")
  htmlwidgets::saveWidget(p, file = htmlFile, selfcontained = TRUE)
}

# Function for pre-processing and scaling of data 
process_data <- function(data) {
  # Drop metadata, NAs, etc.
  # Include k-mers
  
  # Determine the columns to use
  slice_col <- which(colnames(data) == "strain")
  x <- data[, 2:(slice_col - 1)]
  target <- data$variant
  
  # Check for columns that have zero variance for PCA
  non_zero_var_cols <- apply(x, 2, var) > 0
  x <- x[, non_zero_var_cols]

  if (ncol(x) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }

  # # Scale data
  x <- scale(x) # can be optional
  
  return(list(x = x, target = target))
}

# Function that performs PCA
pca_fn <- function(x) {
  pca_df <- prcomp(x, center = TRUE)
  return(pca_df)
}

# Function that plots PCA results
pca_plot <- function(pca_df, data, k) {
  p <- autoplot(pca_df, data = data, colour = 'variant') + scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("pca", k, p)
}

# Function that performs t-SNE (Rtsne library)
Rtsne_fn <- function(pca_results) {
  set.seed(tsne_seed)
  tsne_df <- Rtsne(pca_results, dims = tsne_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter, check_duplicate = FALSE, pca = FALSE)
  return(tsne_df)
}

# Function that performs t-SNE (tsne library)
tsne_fn <- function(pca_results) {
  set.seed(tsne_seed)
  tsne_df <- tsne(pca_results, k = tsne_dims, initial_dims = tsne_initial_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter)
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
k <- 3

# Set paths
data_path <- "./data/archive/"

results_path <- "./results/dim-reduce/R/"

# Check if the directory already exists
if (!dir.exists(results_path)) {
  # Create the directory if it doesn't exist
  dir.create(results_path, recursive = TRUE)
}

# Loading the data file
raw_data <- read.csv(paste0(data_path, "kmer_", k, ".csv"))

# Define t-SNE and UMAP parameters
tsne_seed <- 10
tsne_dims <- 2
tsne_perplexity <- 40
tsne_max_iter <- 300
tsne_initial_dims <- 50
umap_seed <- 10
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1

# Process the data
data <- process_data(raw_data)

# Perform PCA
pca_df <- pca_fn(data$x)

# Plot PCA results
pca_plot(pca_df, raw_data, k)

# Create screeplot
p <- fviz_eig(pca_df)
save_plot("screeplot", k, p)

# Create graph of individuals
p <- fviz_pca_ind(pca_df, 
                  col.ind = "cos2", # Color by the quality of representation
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
                  )
save_plot("indivgraph", k, p)

# Create graph of variables
p <- fviz_pca_var(pca_df,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
save_plot("vargraph", k, p)

# Create biplot of individuals and variables
p <- fviz_pca_biplot(pca_df,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
save_plot("biplot", k, p)

# Perform t-SNE(Rtsne library) using PCA results
tsne_df <- Rtsne_fn(pca_df$x)

# # Perform t-SNE(tsne Library) using PCA results
# TO DO: Check/Fix params needeed (Read docu)
# tsne_df <- tsne_fn(pca_df$x)

# Plot t-SNE results
tsne_plot(tsne_df, data$target, k)

# Perform UMAP
umap_df <- umap_fn(data$x)

# Plot UMAP results
umap_plot(umap_df, data$target, k)

# -----END-----