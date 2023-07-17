# Load Source
source("code/R/kmer-analysis.R")

# Install and load the required packages
if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(
  tidyverse,
  lubridate,
  umap,
  plotly,
  htmlwidgets,
  factoextra,
  scales,
  Rtsne,
  tsne,
  RColorBrewer,
  ggfortify,
  ggbiplot
)

# Hide warnings
options(warn = -1) # Address open issue in plot_ly: warning 'bar' objects don't have these attributes: 'mode'

# -----Functions-----

# Function to extract the timestamp
get_time <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  as.numeric(gsub(".csv", "", parts[3]))
}

# Function to search and read the CSV file
read_kmer_csv <- function(data_path, k) {
  # Get the list of files matching the pattern
  file_pattern <- paste0("kmer_", k, "_", ".*\\.csv$")
  file_list <- list.files(path = data_path, pattern = file_pattern, full.names = FALSE)

  # Check if any files are found
  if (length(file_list) == 0) {
    message("No files found for k = ", k)
    return(NULL)
  }
  
  # Sort the strings based on the timestamp in descending order
  sorted_strings <- file_list[order(sapply(file_list, get_time), decreasing = TRUE)]
  
  df <- read.csv(paste0(data_path, sorted_strings[1]))

  return(df)
}

# Function for saving 2D plots as PNG and HTML
save_plot <- function(method, k, p, is3D=FALSE) {
  # File name for saving
  filename <- paste0(method, "-", k, ".png")
  
  # Note: 3D plots are plot_ly objects, 2D plots are ggplot objects.
  if (is3D) {
    # Save plot_ly obj. as PNG
    save_image(p, paste0(results_path, filename))
  } else {
    # Save as PNG
    ggsave(filename, p, results_path, device = "png", width=5, height=5, dpi = 300, bg = "white")
    # Convert ggplot object to ggplotly
    p <- ggplotly(p)
  }
  
  # Save as HTML
  htmlFile <- paste0(results_path, method, "-", k, ".html")
  htmlwidgets::saveWidget(p, file = htmlFile, selfcontained = TRUE)
}

# Function for pre-processing and scaling of data 
pre_process <- function(data, col_name) {
  # Extract year from date column (This is needed for labeling of points in 3D plots)
  df$year <- format(as.Date(df$date), "%Y")
  
  # Determine the columns to use (drop metadata, retain k-mers)
  slice_col <- which(colnames(data) == "strain")
  x <- data[, 2:(slice_col - 1)]
  target <- data[[col_name]]
  
  # Check for columns that have zero variance for PCA
  non_zero_var_cols <- apply(x, 2, var) > 0
  x <- x[, non_zero_var_cols]
  
  if (ncol(x) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  # Scale data
  x <- scale(x)
  
  return(list(x = x, target = target))
}

# Function that performs PCA
pca_fn <- function(x) {
  pca_df <- prcomp(x, center = TRUE)
  return(pca_df)
}

# Function for 2D PCA plot
pca_plot <- function(pca_df, data, k) {
  p <- autoplot(pca_df, data = data, colour = col_name) + scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("pca", k, p)
}

# Function for 3D PCA plot
pca_3d <- function(pca_df, df, col_name) {
  pc <- as.data.frame(pca_df$x[, 1:3])
  
  # Use plot_ly for 3D visualization
  p <- plot_ly(pc, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = df[[col_name]], 
               text = paste("Variant: ", df$variant, "<br>",
                            "Sex: ", df$sex, "<br>",
                            "Division Exposure: ", df$division_exposure, "<br>",
                            "Year: ", format(as.Date(df$date), "%Y"), "<br>",
                            "Strain: ", df$strain, "<br>",
                            "Pangolin Lineage: ", df$pangolin_lineage))
  
  # Save plot as PNG and HTML
  save_plot("3d-pca", k, p, is3D=TRUE)
}

screeplot <- function(pca_df) {
  p <- fviz_eig(pca_df,
                xlab = "Number of Principal Components")
  
  # Save plot as PNG and HTML
  save_plot("screeplot", k, p)
}

indiv <- function(pca_df) {
  p <- fviz_pca_ind(pca_df,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,     # Avoid text overlapping
               label = list(ind = list(label = "Quality of Representation")),
               xlab = "PC1",
               ylab = "PC2"
  )
  
  # Save plot as PNG and HTML
  save_plot("indivgraph", k, p)
}

vars <- function(pca_df) {
  p <- fviz_pca_ind(pca_df,
                    col.ind = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE,        # Avoid text overlapping
                    label = list(ind = list(label = "Contribution to PC")),
                    xlab = "PC1",
                    ylab = "PC2"
  )
  
  # Save plot as PNG and HTML
  save_plot("vargraph", k, p)
}

biplot <- function(pca_df) {
  # # [Old method] Create biplot of individuals and variables
  # p <- fviz_pca_biplot(pca_df,
  #                 col.var = "#2E9FDF", # Variables color
  #                 col.ind = "#696969"  # Individuals color
  # )
  
  # Create biplot of individuals and variables (using ggbiplot)
  p <- ggbiplot(pca_df, obs.scale = 1, var.scale = 1,
                groups = target, ellipse = TRUE, circle = TRUE) +
                scale_color_discrete(name = '') +
                theme(legend.direction = 'horizontal', legend.position = 'top')
  
  # Save plot as PNG and HTML
  save_plot("biplot", k, p)
}

# Function that performs t-SNE (Rtsne library)
Rtsne_fn <- function(pca_results, tsne_dims) {
  set.seed(tsne_seed)
  tsne_df <- Rtsne(pca_results, dims = tsne_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter, check_duplicate = FALSE, pca = FALSE)
  return(tsne_df)
}

# Function that includes visualization for each t-SNE iteration
ecb <- function(x) {
  epoc_df <- data.frame(x,target = df[[col_name]])
  
  plt <- ggplot(epoc_df,aes(x = X1, y = X2,label = target,color = target)) + geom_text()
  
  print(plt)
}

# Function that performs t-SNE (tsne library)
tsne_fn <- function(pca_results, tsne_dims) {
  set.seed(tsne_seed)
  if (tsne_dims == 2) {
    tsne_df <- tsne(pca_results, k = tsne_dims, initial_dims = tsne_initial_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter, epoch_callback = ecb)
  } else {
    tsne_df <- tsne(pca_results, k = tsne_dims, initial_dims = tsne_initial_dims, perplexity = tsne_perplexity, max_iter = tsne_max_iter)
  }
  
  return(tsne_df)
}

# Function for 2D t-SNE plot
tsne_plot <- function(tsne_df, target, k, is_tsne) {
  if (is_tsne) {
    df <- data.frame(X1 = tsne_df[, 1], X2 = tsne_df[, 2], target = target)
  } else {
    df <- data.frame(X1 = tsne_df$Y[, 1], X2 = tsne_df$Y[, 2], target = target)
  }
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

# Function for 3D t-SNE plot
tsne_3d <- function(tsne_df, df, col_name) {
  final <- cbind(data.frame(tsne_df), df[[col_name]])
  p <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers", color = ~df[[col_name]], 
               text = paste("Variant: ", df$variant, "<br>",
                            "Sex: ", df$sex, "<br>",
                            "Division Exposure: ", df$division_exposure, "<br>",
                            "Year: ", format(as.Date(df$date), "%Y"), "<br>",
                            "Strain: ", df$strain, "<br>",
                            "Pangolin Lineage: ", df$pangolin_lineage))
  
  # Save plot as PNG and HTML
  save_plot("3d-tsne", k, p, is3D=TRUE)
}

# Function that performs UMAP
umap_fn <- function(x, umap_dims) {
  umap_df <- umap(x, n_components = umap_dims, n_neighbors = umap_n_neighbors, metric = umap_metric, min_dist = umap_min_dist, random_state = umap_seed)
  return(umap_df)
}

# Function for 2D UMAP plot
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

# Function for 3D UMAP plot
umap_3d <- function(umap_df, df, col_name) {
  final <- cbind(data.frame(umap_df[["layout"]]), df[[col_name]])
  p <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers", color = ~df[[col_name]], 
               text = paste("Variant: ", df$variant, "<br>",
                            "Sex: ", df$sex, "<br>",
                            "Division Exposure: ", df$division_exposure, "<br>",
                            "Year: ", format(as.Date(df$date), "%Y"), "<br>",
                            "Strain: ", df$strain, "<br>",
                            "Pangolin Lineage: ", df$pangolin_lineage))
  
  p <- p %>% add_markers() 
  p <- p %>% layout(scene = list(xaxis = list(title = '0'), 
                                 yaxis = list(title = '1'), 
                                 zaxis = list(title = '2'))) 
  
  # Save plot as PNG and HTML
  save_plot("3d-umap", k, p, is3D=TRUE)
}

# -----END of Functions-----



# -----START-----

# Set k for k-mer analysis (Valid values: 3, 5, 7)
k <- 7

# Set paths
data_path <- "./data/kmers/"

results_path <- "./results/dim-reduce/R/"

# Check if the directory already exists
if (!dir.exists(results_path)) {
  # Create the directory if it doesn't exist
  dir.create(results_path, recursive = TRUE)
}

# Search and read the CSV file
df <- read_kmer_csv(data_path, k)

# Define t-SNE and UMAP parameters
tsne_seed <- 0
tsne_perplexity <- 40
tsne_max_iter <- 1000
tsne_initial_dims <- 50
umap_seed <- 0
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1

# Define Target Column Name
col_name <- "variant"

# Pre-process the data
data <- pre_process(df, col_name)
x <- data$x
target <- data$target

# Perform PCA
pca_df <- pca_fn(x)

# Generate 2D PCA plot
pca_plot(pca_df, df, k)

# Generate 3D PCA plot (does not run PCA again)
pca_3d(pca_df, df, col_name)

# Generate screeplot
screeplot(pca_df)

# Generate graph of individuals
indiv(pca_df)

# Generate graph of variables
vars(pca_df)

# Generate biplot
biplot(pca_df)

# Perform t-SNE via 'tsne' library using PCA results (in 2 dimensions)
# # Note: Uncomment the next two line to use tsne; otherwise, comment them
# is_tsne <- TRUE
# tsne_df <- tsne_fn(pca_df$x, 2)

# Perform t-SNE via 'Rtsne' library using PCA results (in 2 dimensions)
# # Note: Uncomment the next two line to use Rtsne; otherwise, comment them
is_tsne <- FALSE
tsne_df <- Rtsne_fn(pca_df$x, 2)

# Generate 2D t-SNE plot
tsne_plot(tsne_df, target, k, is_tsne)

# Generate 3D t-SNE plot (runs t-SNE again in 3 dimensions)
# # Note: Uncomment the two succeeding lines to use tsne; otherwise, comment them
# tsne_df <- tsne_fn(pca_df$x, 3)
# tsne_3d(tsne_df, df, col_name)

# # Note: Uncomment the two succeeding lines to use Rtsne; otherwise, comment them
tsne_df <- Rtsne_fn(pca_df$x, 3)
tsne_3d(tsne_df$Y, df, col_name)

# Perform UMAP (in 2 dimensions)
umap_df <- umap_fn(x, 2)

# Generate 2D UMAP plot
umap_plot(umap_df, target, k)

# Generate 3D UMAP plot (runs t-SNE again in 3 dimensions)
umap_df <- umap_fn(x, 3)
umap_3d(umap_df, df, col_name)

# -----END-----