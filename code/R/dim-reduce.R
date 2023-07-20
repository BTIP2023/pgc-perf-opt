# TO DO: Fix readability and documentation of this file

# -----Functions-----

# Function to extract the timestamp from the kmer files
extract_time <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  as.numeric(gsub(".csv", "", parts[3]))
}

# Function to search and read the CSV file
read_kmer_csv <- function(data_path, k) {
  print("Reading CSV file...")
  # Get the list of files matching the pattern
  file_pattern <- paste0("kmer_", k, "_", ".*\\.csv$")
  file_list <- list.files(
    path = data_path, pattern = file_pattern,
    full.names = FALSE
  )
  
  # Check if any files are found
  if (length(file_list) == 0) {
    message("No files found for k = ", k)
    return(NULL)
  }
  
  # Sort the strings based on the timestamp in descending order
  sorted_strings <- file_list[order(sapply(file_list, extract_time),
                                    decreasing = TRUE
  )]
  
  df <- read.csv(paste(data_path, sorted_strings[1], sep = "/"))
  
  return(df)
}

# Function for saving 2D plots as PNG and HTML
save_plot <- function(method, k, p, is_3d = FALSE) {
  print("Saving plot...")
  # File name for saving
  filename <- paste0(method, "-", k, ".png")
  
  # Note: 3D plots are plot_ly objects, 2D plots are ggplot objects.
  if (is_3d) {
    print("ISSUE")
    # Save plot_ly obj. as PNG
    # save_image(p, paste(results_path, filename, sep = "/"))
  } else {
    # Save as PNG
    ggsave(filename, p, results_path,
           device = "png", width = 5, height = 5,
           dpi = 300, bg = "white"
    )
    # Convert ggplot object to ggplotly
    p <- ggplotly(p) 
  }
  
  # Save as HTML
  html_file <- paste0(results_path, "/", method, "-", k, ".html")
  htmlwidgets::saveWidget(p, file = html_file, selfcontained = TRUE)
}

# Function for pre-processing and scaling of data
pre_process <- function(df, col_name) {
  print("Pre-processing and scaling data...")
  # Extract year from date column
  
  # Determine the columns to use (drop metadata, retain k-mers)
  slice_col <- which(colnames(df) == "strain")
  x <- df[, 2:(slice_col - 1)]
  target <- df[[col_name]]
  
  # Check for columns that have zero variance for PCA
  non_zero_var_cols <- apply(x, 2, var) > 0
  x <- x[, non_zero_var_cols]
  
  if (ncol(x) < 2) {
    stop("Insufficient columns with non-zero variance for PCA.")
  }
  
  # Scale data
  x <- scale(x)
  
  # (This is needed for labeling of points in 3D plots)
  df$year <- format(as.Date(df$date), "%Y")

  return(list(x = x, target = target))
}

# Function to execute before main dim-reduce codes
pre_reduce <- function(results_path, data_path, k, col_name) {
  # Check if the directory already exists
  if (!dir.exists(results_path)) {
    # Create the directory if it doesn't exist
    dir.create(results_path, recursive = TRUE)
  }
  # Search and read the CSV file
  df <- read_kmer_csv(data_path, k)
  # Pre-process the data
  data <- pre_process(df, col_name)
  return(list(df = df, data = data))
}

# Function that performs PCA
pca_fn <- function(x) {
  print("Performing PCA...")
  pca_df <- prcomp(x, center = TRUE)
  return(pca_df)
}

# Function for 2D PCA plot
pca_plot <- function(pca_df, df, k, target) {
  variant <- df[[target]]
  print("Generating 2D PCA plot...")
  p <- autoplot(pca_df, data = df) +
    geom_point(aes(color = variant, text = paste(
      "Identifier: ", df$gisaid_epi_isl, "\n",
      "Variant: ", df$variant, "\n",
      "Sex: ", df$sex, "\n",
      "Division Exposure: ", df$division_exposure, "\n",
      "Year: ", format(as.Date(df$date), "%Y"), "\n",
      "Strain: ", df$strain, "\n",
      "Pangolin Lineage: ", df$pangolin_lineage
    ))) +
    scale_color_brewer(palette = "Set1")
  # Save plot as PNG and HTML
  save_plot("pca", k, p)
}

# Function for 3D PCA plot
pca_3d <- function(pca_df, df, target, k) {
  print("Generating 3D PCA plot...")
  pc <- as.data.frame(pca_df$x[, 1:3])
  
  # Use plot_ly for 3D visualization
  p <- plot_ly(pc,
               x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
               mode = "markers", color = df[[target]],
               text = paste(
                 "Identifier: ", df$gisaid_epi_isl, "\n",
                 "Variant: ", df$variant, "\n",
                 "Sex: ", df$sex, "\n",
                 "Division Exposure: ", df$division_exposure, "\n",
                 "Year: ", format(as.Date(df$date), "%Y"), "\n",
                 "Strain: ", df$strain, "\n",
                 "Pangolin Lineage: ", df$pangolin_lineage
               )
  )
  
  # Save plot as PNG and HTML
  save_plot("3d-pca", k, p, is_3d = TRUE)
}

# Function that generates scree plot from PCA results
screeplot <- function(pca_df, k) {
  print("Generating scree plot...")
  p <- fviz_eig(pca_df,
                xlab = "Number of Principal Components"
  )
  
  # Save plot as PNG and HTML
  save_plot("screeplot", k, p)
}

# Function that generates factor loadings of first
# n_components Principal Components
factor_loadings <- function(pca_df, x, k, n_components) {
  print(paste("Generating factor loadings plot of first", n_components,
              "PCs...",
              sep = " "
  ))
  # Extract factor loadings
  loadings <- pca_df$rotation
  
  # Plot bar plots for factor loadings of n principal components
  for (i in 1:n_components) {
    # Create a data frame for the factor loadings
    loadings_df <- data.frame(
      variable = colnames(x),
      loading = loadings[, i]
    )
    
    # Create a bar plot using ggplot2
    p <- ggplot(loadings_df, aes(x = variable, y = loading)) +
      geom_bar(stat = "identity", fill = "blue") +
      labs(
        title = paste("Principal Component", i),
        x = "Variables", y = "Factor Loadings"
      )
    
    # Save plot as PNG and HTML
    save_plot(paste("loadings", i, sep = "-"), k, p)
  }
}

# Function that generates graph of individuals from PCA results
indiv <- function(pca_df, k) {
  print("Generating graph of individuals...")
  p <- fviz_pca_ind(pca_df,
                    col.ind = "cos2", # Color by the quality of representation
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE, # Avoid text overlapping
                    xlab = "PC1",
                    ylab = "PC2",
  )
  
  # Save plot as PNG and HTML
  save_plot("indivgraph", k, p)
}

# Function that generates graph of variables from PCA results
vars <- function(pca_df, k) {
  print("Generating graph of variables...")
  p <- fviz_pca_var(pca_df,
                    col.var = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE, # Avoid text overlapping
                    xlab = "PC1",
                    ylab = "PC2"
  )
  
  # Save plot as PNG and HTML
  save_plot("vargraph", k, p)
}

# Function that generates biplot from PCA results
biplot <- function(pca_df, k) {
  print("Generating biplot...")
  # # [Old method] Create biplot of individuals and variables
  p <- fviz_pca_biplot(pca_df,
                       col.var = "#2E9FDF", # Variables color
                       col.ind = "#696969", # Individuals color
                       addEllipses = TRUE,
                       xlab = "PC1",
                       ylab = "PC2"
  )
  
  # # Create biplot of individuals and variables (using ggbiplot)
  # p <- ggbiplot(pca_df,
  #   obs.scale = 1, var.scale = 1,
  #   groups = target, ellipse = TRUE, circle = TRUE
  # ) +
  #   scale_color_discrete(name = "") +
  #   theme(legend.direction = "horizontal", legend.position = "top")
  
  # Save plot as PNG and HTML
  save_plot("biplot", k, p)
}

# Function that performs t-SNE (Rtsne library)
rtsne_fn <- function(pca_results, tsne_dims, tsne_perplexity, tsne_max_iter, tsne_seed) {
  print("Performing t-SNE...")
  set.seed(tsne_seed)
  tsne_df <- Rtsne(pca_results,
                   dims = tsne_dims, perplexity = tsne_perplexity,
                   max_iter = tsne_max_iter, check_duplicate = FALSE,
                   pca = FALSE
  )
  return(tsne_df)
}

# Function that includes visualization for each t-SNE iteration
ecb <- function(x) {
  epoc_df <- data.frame(x, target = df[[col_name]])
  
  plt <- ggplot(epoc_df, aes(
    x = X1, y = X2,
    label = df[[target]], color = df[[target]]
  )) +
    geom_text()
  
  print(plt)
}

# Function that performs t-SNE (tsne library)
tsne_fn <- function(pca_results, tsne_dims, tsne_initial_dims, tsne_perplexity, tsne_max_iter, tsne_seed) {
  print("Performing t-SNE...")
  set.seed(tsne_seed)
  if (tsne_dims == 2) {
    tsne_df <- tsne(pca_results,
                    k = tsne_dims,
                    initial_dims = tsne_initial_dims,
                    perplexity = tsne_perplexity,
                    max_iter = tsne_max_iter,
                    # epoch_callback = ecb
    )
  } else {
    tsne_df <- tsne(pca_results,
                    k = tsne_dims,
                    initial_dims = tsne_initial_dims,
                    perplexity = tsne_perplexity,
                    max_iter = tsne_max_iter
    )
  }
  
  return(tsne_df)
}

# Function for 2D t-SNE plot
tsne_plot <- function(tsne_df, df, target, k, is_tsne) {
  print("Generating 2D t-SNE plot...")
  if (is_tsne) {
    tsne_df <- data.frame(X1 = tsne_df[, 1], X2 = tsne_df[, 2], target = df[[target]])
  } else {
    tsne_df <- data.frame(X1 = tsne_df$Y[, 1], X2 = tsne_df$Y[, 2], target = df[[target]])
  }
  # Create ggplot object
  p <- ggplot(tsne_df, aes(x = X1, y = X2, color = target)) +
    geom_point(aes(text = paste(
      "Identifier: ", df$gisaid_epi_isl, "\n",
      "Variant: ", df$variant, "\n",
      "Sex: ", df$sex, "\n",
      "Division Exposure: ", df$division_exposure, "\n",
      "Year: ", format(as.Date(df$date), "%Y"), "\n",
      "Strain: ", df$strain, "\n",
      "Pangolin Lineage: ", df$pangolin_lineage
    ))) +
    xlab("TSNE-2D-1") +
    ylab("TSNE-2D-2") +
    labs(color = "Label") +
    scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("tsne", k, p)
}

# Function for 3D t-SNE plot
tsne_3d <- function(tsne_df, df, target, k) {
  print("Generating 3D t-SNE plot...")
  final <- cbind(data.frame(tsne_df), df[[target]])
  p <- plot_ly(final,
               x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers",
               color = ~target,
               text = paste(
                 "Identifier: ", df$gisaid_epi_isl, "<br>",
                 "Variant: ", df$variant, "<br>",
                 "Sex: ", df$sex, "<br>",
                 "Division Exposure: ", df$division_exposure, "<br>",
                 "Year: ", format(as.Date(df$date), "%Y"), "<br>",
                 "Strain: ", df$strain, "<br>",
                 "Pangolin Lineage: ", df$pangolin_lineage
               )
  )
  
  # Save plot as PNG and HTML
  save_plot("3d-tsne", k, p, is_3d = TRUE)
}

# Function that performs UMAP
umap_fn <- function(x, umap_dims, umap_n_neighbors, umap_metric, umap_min_dist, umap_seed) {
  print("Performing UMAP...")
  umap_df <- umap(x,
                  n_components = umap_dims, n_neighbors = umap_n_neighbors,
                  metric = umap_metric, min_dist = umap_min_dist,
                  random_state = umap_seed
  )
  return(umap_df)
}

# Function for 2D UMAP plot
umap_plot <- function(umap_df, df, target, k) {
  print("Generating 2D UMAP plot...")
  emb <- umap_df$layout
  
  x_o <- emb[, 1]
  y_o <- emb[, 2]
  
  # Create ggplot object
  p <- ggplot(
    data = as.data.frame(emb),
    aes(x = x_o, y = y_o, color = df[[target]])
  ) +
    geom_point(aes(text = paste(
      "Identifier: ", df$gisaid_epi_isl, "\n",
      "Variant: ", df$variant, "\n",
      "Sex: ", df$sex, "\n",
      "Division Exposure: ", df$division_exposure, "\n",
      "Year: ", format(as.Date(df$date), "%Y"), "\n",
      "Strain: ", df$strain, "\n",
      "Pangolin Lineage: ", df$pangolin_lineage
    ))) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    labs(color = "Label") +
    scale_color_brewer(palette = "Set1")
  
  # Save plot as PNG and HTML
  save_plot("umap", k, p)
}

# Function for 3D UMAP plot
umap_3d <- function(umap_df, df, target, k) {
  print("Generating 3D UMAP plot...")
  final <- cbind(data.frame(umap_df[["layout"]]), target)
  p <- plot_ly(final,
               x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers",
               color = ~target,
               text = paste(
                 "Identifier: ", df$gisaid_epi_isl, "<br>",
                 "Variant: ", df$variant, "<br>",
                 "Sex: ", df$sex, "<br>",
                 "Division Exposure: ", df$division_exposure, "<br>",
                 "Year: ", format(as.Date(df$date), "%Y"), "<br>",
                 "Strain: ", df$strain, "<br>",
                 "Pangolin Lineage: ", df$pangolin_lineage
               )
  )
  
  p <- p %>% add_markers()
  p <- p %>% layout(scene = list(
    xaxis = list(title = "0"),
    yaxis = list(title = "1"),
    zaxis = list(title = "2")
  ))
  
  # Save plot as PNG and HTML
  save_plot("3d-umap", k, p, is_3d = TRUE)
}

# Main Function
dim_reduce <- function(k, data_path, results_path, tsne_seed, tsne_perplexity,
                       tsne_max_iter, tsne_initial_dims, umap_seed,
                       umap_n_neighbors, umap_metric, umap_min_dist, col_name) {
  # -----START-----

  pre_reduce_res <- pre_reduce(results_path, data_path, k, col_name)
  
  df <- pre_reduce_res$df                # df is the original dataset
  x <- pre_reduce_res$data$x             # x is the scaled data
  target <- col_name                     # target is the column used for 
                                         # clustering

  # Perform PCA
  pca_df <- pca_fn(x)

  # Generate 2D PCA plot
  pca_plot(pca_df, df, k, target)

  # Generate 3D PCA plot (does not run PCA again)
  pca_3d(pca_df, df, target, k)

  # Generate screeplot
  screeplot(pca_df, k)

  # Generate factor loadings plot of first 3 principal components
  factor_loadings(pca_df, x, k, 3)

  # Generate graph of individuals
  indiv(pca_df, k)

  # Generate graph of variables
  vars(pca_df, k)

  # Generate biplot
  biplot(pca_df, k)

  # Perform t-SNE via 'tsne' library using PCA results (in 2 dimensions)
  # # Note: Uncomment the next two line to use tsne; otherwise, comment them
  is_tsne <- TRUE
  tsne_df <- tsne_fn(pca_df$x, 2, tsne_initial_dims, tsne_perplexity, 
                     tsne_max_iter, tsne_seed)

  # Perform t-SNE via 'Rtsne' library using PCA results (in 2 dimensions)
  # # Note: Uncomment the next two line to use Rtsne; otherwise, comment them
  # is_tsne <- FALSE
  # tsne_df <- rtsne_fn(pca_df$x, 2)

  # Generate 2D t-SNE plot
  tsne_plot(tsne_df, df, target, k, is_tsne)

  # Generate 3D t-SNE plot (runs t-SNE again in 3 dimensions)
  # # Note: Uncomment the two succeeding lines to use tsne;
  # # otherwise, comment them
  tsne_df <- tsne_fn(pca_df$x, 3, tsne_initial_dims, tsne_perplexity, 
                     tsne_max_iter, tsne_seed)
  tsne_3d(tsne_df, df, target, k)

  # # Note: Uncomment the two succeeding lines to use Rtsne;
  # # otherwise, comment them
  # tsne_df <- rtsne_fn(pca_df$x, 3, tsne_perplexity, tsne_max_iter, tsne_seed)
  # tsne_3d(tsne_df$Y, df, target, k)

  # Perform UMAP (in 2 dimensions)
  umap_df <- umap_fn(x, 2, umap_n_neighbors, umap_metric, 
                     umap_min_dist, umap_seed)

  # Generate 2D UMAP plot
  umap_plot(umap_df, df, target, k)

  # Generate 3D UMAP plot (runs t-SNE again in 3 dimensions)
  umap_df <- umap_fn(x, 3, umap_n_neighbors, umap_metric, 
                     umap_min_dist, umap_seed)
  umap_3d(umap_df, df, target, k)

  # -----END-----
}
