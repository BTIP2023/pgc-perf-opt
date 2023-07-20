# File: pipeline.R
# Contains the complete bioinformatics pipeline
# for convenient runs and benchmarks.
# This will source functions from code/R/

# INSTALL AND LOAD PACKAGES ################################
options(repos = "https://cloud.r-project.org/")

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")
library(pacman)

### ATTN: IN LINUX SYSTEMS, CONSULT README FOR ADDITIONAL PREREQUISITES
### BEFORE RUNNING ANY SCRIPT. ISSUE: tidyverse installation.
### This is a non-issue if code is ran within cpu.Dockerfile.
### This cannot be scripted because this requires sudo priveleges.

# Install xml2 in advance to prep for tidyverse installation in Linux.
# Note that in Windows RStudio, this is installed by default.
# If you're getting xml2 errors on Windows, you broke something lol.
if (pacman::p_detectOS() == "Linux" && !pacman::p_exists(xml2, local = TRUE)) {
  install.packages("xml2", dependencies = TRUE, INSTALL_opts = c("--no-lock"))
  pacman::p_load(xml2)
}

# Use pacman to load add-on packages as desired.
# TODO: Remove redundancies in dependencies. E.g., dplyr and ggplot2
# are already dependencies of tidyverse.
pacman::p_load(dplyr, GGally, ggplot2, ggthemes, ggvis,
               httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse,
               ape, kmer, readr, validate, gsubfn, seqinr,
               umap, htmlwidgets, factoextra, scales, ggbiplot,
               Rtsne, tsne, RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark, flexiblas)
install_github("vqv/ggbiplot", upgrade = FALSE, quiet = TRUE)

# LOAD SOURCES #############################################
source("code/R/helper.R")
source("code/R/preprocess.R")
source("code/R/kmer-analysis.R")
source("code/R/dim-reduce.R")
source("code/R/clustering-variant.R")
source("code/R/clustering-region.R")

# SET PARAMETERS ###########################################
# pipeline.R general parameters
seed <- 1234
stamp <- get_time()
kmer_list <- c(3, 5, 7)

# microbenchmark parameters
bm_times <- 1L   # how many times should routine be evaluated
bm_units <- "seconds"

# preprocess.R::preprocess() parameters
data_path_gisaid <- "data/GISAID"
extract_path <- "data/GISAID/datasets"
strat_size <- 100
country_exposure <- "Philippines"
write_fastacsv <- FALSE

# kmer-analysis.R::get_kmers() parameters

# dim-reduce.R::dim_reduce() parameters
data_path_kmers <- "data/kmers"
results_path_dimreduce <- "results/dim-reduce/R"
tsne_perplexity <- 40
tsne_max_iter <- 1000
tsne_initial_dims <- 50
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1
target_col <- "variant"

# AGNES Clustering Parameters :: dendogram_create_x()
results_path_agnes <- "results/dendrogram"

# RUN PIPELINE #############################################
# Step 1: preprocess()
list[fasta_all, metadata_all] <- preprocess(data_path_gisaid, extract_path, seed,
                                            strat_size, country_exposure,
                                            write_fastacsv, stamp)
# Step 2: get_kmers()
for (k in kmer_list) {
  get_kmers(fasta_all, metadata_all, k, stamp)
}

# Step 3: dim_reduce()
for (k in kmer_list) {
  dim_reduce(k, data_path_kmers, results_path_dimreduce,
             tsne_seed = seed, tsne_perplexity,
             tsne_max_iter, tsne_initial_dims,
             umap_seed = seed, umap_n_neighbors,
             umap_metric, umap_min_dist, col_name = target_col)
}

# Function for setting backend before performing operation(function)
set_backend <- function(operation, args_list, backend) {
  flexiblas_load_backend(backend)
  do.call(operation, args_list)
}

# Benchmarking function
benchmark_backends <- function(operation, args_list, backends, 
                               bm_times, bm_units) {
  results <- list()
  for (backend in backends) {
    # Run the benchmark for the current backend
    bm_result <- microbenchmark(set_backend(operation, args_list, backend),
                                times = bm_times, unit = bm_unit)
    # Store the results
    results[[backend]] <- bm_result$time
  }
  results
}

# Benchmark plotting function
plot_results <- function(method_benchmark, method) {
  # Combine the results into a data frame
  results_df <- data.frame(Backend = rep(selected_backends, each = 100), 
                           Time = unlist(method_benchmark))
  # Create the bar plot
  fig <- ggplot(results_df, aes(x = Backend, y = Time, fill = Backend)) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Benchmark: ", method),
         x = "Backend",
         y = "Execution Time (ms)") +
    theme_minimal()
  
  # Create directory for storing plots
  results_path <- "results/benchmark/R"
  
  # Check if the directory already exists
  if (!dir.exists(results_path)) {
    # Create the directory if it doesn't exist
    dir.create(results_path, recursive = TRUE)
  }
  
  # Save the bar plot
  ggsave(filename=paste0(method, "-benchmark.png"), 
         plot = fig, path = results_path, 
         device = "png", width = 5, height = 5,
         dpi = 300, bg = "white"
  )
}

# microbenchmark parameters
bm_times <- 100
bm_unit <- "seconds"

# Backends to compare
selected_backends <- c("NETLIB", "ATLAS", "OPENBLASOPENMP",
                       "MKLOPENMP", "BLISOPENMP", "__FALLBACK__")

data <- pre_reduce(results_path_dimreduce, data_path_kmers, k, df, col_name)

x <- data$x
target <- data$target

# Run an iteration of PCA for t-SNE use
pca_df <- pca_fn(x)

# Benchmark backends on pca_fn
pca_benchmark <- benchmark_backends(pca_fn, list(x), 
                                    selected_backends, bm_times, bm_unit)

# Plot PCA benchmark results
plot_results(pca_benchmark, "t-SNE")

# Benchmark backends on tsne_fn (2-dimensions)
tsne_benchmark <- benchmark_backends(tsne_fn, list(pca_df, 2), 
                                     selected_backends, bm_times, bm_unit)

# Plot t-SNE benchmark results (2D)
plot_results(tsne_benchmark, "t-SNE (2D)")

# Benchmark backends on tsne_fn (3-dimensions)
tsne_benchmark <- benchmark_backends(tsne_fn, list(pca_df, 3), 
                                     selected_backends, bm_times, bm_unit)

# Plot t-SNE benchmark results (3D)
plot_results(tsne_benchmark, "t-SNE (3D)")

# Benchmark backends on umap_fn (2-dimensions)
umap_benchmark <- benchmark_backends(umap_fn, list(x, 2), 
                                     selected_backends, bm_times, bm_unit)

# Plot UMAP benchmark results (2D)
plot_results(umap_benchmark, "t-SNE (2D)")

# Benchmark backends on umap_fn (3-dimensions)
umap_benchmark <- benchmark_backends(umap_fn, list(x, 3), 
                                     selected_backends, bm_times, bm_unit)

# Plot UMAP benchmark results (3D)
plot_results(umap_benchmark, "t-SNE (3D)")

print("All operations completed successfully!")