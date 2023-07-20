if (!require("pacman"))
  install.packages("pacman")
library(pacman)

if (pacman::p_detectOS() == "Linux" && !pacman::p_exists(xml2, local = TRUE)) {
  install.packages("xml2", dependencies = TRUE, INSTALL_opts = c("--no-lock"))
  pacman::p_load(xml2)
}

# Load required packages
pacman::p_load(plyr, dplyr, GGally, ggplot2, ggthemes, ggvis,
               httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse,
               ape, kmer, readr, validate, gsubfn, seqinr,
               umap, htmlwidgets, factoextra, scales,
               Rtsne, tsne, RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark)
install_github("vqv/ggbiplot", upgrade = FALSE, quiet = TRUE)
pacman::p_load(ggbiplot)

# Load Source
source("code/R/dim-reduce.R")

# Function for setting backend before performing operation(function)
set_backend <- function(operation, args_list, backend) {
  flexiblas_load_backend(backend)
  do.call(operation, args_list)
}

# Benchmarking function
benchmark_backends <- function(operation, args_list, backends, times, unit) {
  results <- list()
  for (backend in backends) {
    # Run the benchmark for the current backend
    bm_result <- microbenchmark(set_backend(operation, args_list, backend),
                                times = times, unit = unit)
    # Store the results
    results[[backend]] <- bm_result$time
  }
  results
}

# Benchmark plotting function
plot_results <- function(method_benchmark, method) {
  # Combine the results into a data frame
  results_df <- data.frame(Backend = rep(selected_backends, each = 1), 
                           Time = unlist(method_benchmark))
  # Create the bar plot
  p <- ggplot(results_df, aes(x = Backend, y = Time, fill = Backend)) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Benchmark: ", method),
         x = "Backend",
         y = "Execution Time (ms)") +
    theme_minimal()
  
  # Create directory for storing plots
  results_path <- "results/benchmark/R/FlexiBLAS"
  
  # Check if the directory already exists
  if (!dir.exists(results_path)) {
    # Create the directory if it doesn't exist
    dir.create(results_path, recursive = TRUE)
  }
  
  # Save the bar plot
  ggsave(filename=paste0(method, "-", k, "-benchmark.png"), 
         plot = p, path = results_path, 
         device = "png", width = 5, height = 5,
         dpi = 300, bg = "white")
         
  # Convert ggplot object to ggplotly
  p <- ggplotly(p) 
  
  # Save as HTML
  html_file <- paste0(results_path, "/", method, "-", k, ".html")
  htmlwidgets::saveWidget(p, file = html_file, selfcontained = TRUE)
}

# Backends to compare
selected_backends <- c("NETLIB", "ATLAS", "OPENBLASOPENMP",
                       "MKLOPENMP", "BLISOPENMP", "__FALLBACK__")

# dim-reduce.R::dim_reduce() parameters
seed <- 1234
k_vals <- c(3, 5, 7)
data_path_kmers <- "data/kmers"
results_path_dimreduce <- "results/dim-reduce/R"
tsne_perplexity <- 40
tsne_max_iter <- 1000
tsne_initial_dims <- 50
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1
target_col <- "variant"

# Benchmarking parameters
bm_times <- 1 
bm_unit <- "seconds"

for (k in k_vals) {
  pre_reduce_res <- pre_reduce(results_path_dimreduce, data_path_kmers, k, target_col)
  
  df <- pre_reduce_res$df                # df is the original dataset
  x <- pre_reduce_res$data$x             # x is the scaled data
  target <- target_col                   # target is the column used for 
                                         # clustering
  
  # Run an iteration of PCA for t-SNE use
  pca_df <- pca_fn(x)
  
  # Benchmark backends on pca_fn
  pca_benchmark <- benchmark_backends(pca_fn, list(x),
                                      selected_backends,
                                      bm_times,
                                      bm_unit)
  
  # Plot PCA benchmark results
  plot_results(pca_benchmark, "pca")
  
  # Benchmark backends on tsne_fn (2-dimensions)
  tsne_benchmark <- benchmark_backends(tsne_fn,
                                       list(pca_df$x, 2, tsne_initial_dims,
                                            tsne_perplexity, tsne_max_iter,
                                            tsne_seed = seed),
                                       selected_backends,
                                       bm_times,
                                       bm_unit)

  # Plot t-SNE benchmark results (2D)
  plot_results(tsne_benchmark, "tsne-2d")
  
  # Benchmark backends on tsne_fn (3-dimensions)
  tsne_benchmark <- benchmark_backends(tsne_fn,
                                       list(pca_df$x, 3, tsne_initial_dims,
                                            tsne_perplexity, tsne_max_iter,
                                            tsne_seed = seed),
                                       selected_backends,
                                       bm_times,
                                       bm_unit)
  
  # Plot t-SNE benchmark results (3D)
  plot_results(tsne_benchmark, "tsne-3d")
  
  # Benchmark backends on umap_fn (2-dimensions)
  umap_benchmark <- benchmark_backends(umap_fn,
                                       list(x, 2, umap_n_neighbors,
                                            umap_metric, umap_min_dist,
                                            umap_seed = seed),
                                       selected_backends,
                                       bm_times,
                                       bm_unit)
  
  # Plot UMAP benchmark results (2D)
  plot_results(umap_benchmark, "umap-2d")
  
  # Benchmark backends on umap_fn (3-dimensions)
  umap_benchmark <- benchmark_backends(umap_fn,
                                       list(x, 3, umap_n_neighbors,
                                            umap_metric, umap_min_dist,
                                            umap_seed = seed),
                                       selected_backends,
                                       bm_times,
                                       bm_unit)
  
  # Plot UMAP benchmark results (3D)
  plot_results(umap_benchmark, "umap-3d")
  
  # Benchmark backends on dim_reduce
  dr_benchmark <- benchmark_backends(dim_reduce, 
                                     list(3, data_path_kmers, 
                                          results_path_dimreduce,
                                          tsne_seed = seed, tsne_perplexity,
                                          tsne_max_iter, tsne_initial_dims,
                                          umap_seed = seed, umap_n_neighbors,
                                          umap_metric, umap_min_dist, 
                                          col_name = target_col), 
                                     selected_backends, bm_times, bm_unit)
  
  # Plot dim_reduce benchmark results
  plot_results(dr_benchmark, "dim-red")
  
  print(paste0("Benchmarking ", k, "-mer Dimensionality Reduction Done :>"))
}