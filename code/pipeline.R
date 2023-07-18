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
### This is a non-issue if code is ran within cpu-gpu.Dockerfile.
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
# is already a dependency of tidyverse.
pacman::p_load(dplyr, GGally, ggplot2, ggthemes, ggvis,
               httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse,
               ape, kmer, readr, validate, gsubfn, seqinr,
               umap, htmlwidgets, factoextra, scales, ggbiplot,
               Rtsne, tsne, RColorBrewer, ggfortify, devtools)
install_github("vqv/ggbiplot", upgrade = FALSE)

# validate used for %vin% operator 
# gsubfn used to destructure more than one return value
# devtools supports install_github for installing ggbiplot
# Note: Divisive k-means clustering available via kmer::cluster

# LOAD SOURCES #############################################
source("code/R/helper.R")
source("code/R/preprocess.R")
source("code/R/kmer-analysis.R")
source("code/R/dim-reduce.R")

# SET PARAMETERS ###########################################
# pipeline.R general parameters
seed <- 1234
stamp <- get_time()
kmer_list <- c(3, 5, 7)

## Critical routines:
# [preprocess, get_kmers per k (kcount), get_kmers for all k,
#  ...]
benchmark_mode <- TRUE
benchmark_times <- 1L   # how many times should routine be evaluated

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

# RUN PIPELINE #############################################

if (!benchmark_mode) {
  # Step 1: preprocess()
  list[fasta_all, metadata_all] <- preprocess(data_path, extract_path, seed,
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
} else {
  # Initialize benchmark results collector, write to log later
  benchmark_results = list()

  # Benchmark:
  # To benchmark preprocess starting from extraction, delete:
  # data/GISAID/datasets/
  # To benchmark while generating data, set write_fastacsv = TRUE.
  microbenchmark(
    preprocess = list[fasta_all, metadata_all] <-
      preprocess(data_path, extract_path, seed,
                 strat_size, country_exposure,
                 write_fastacsv, stamp),
    get_kmers_loop = for (k in kmer_list) {
      get_kmers(fasta_all, metadata_all, k, stamp)
    },
    list = list(test = 1+1),
    times = benchmark_times,
    unit = "seconds",
    control = list(order = "inorder", warmup = 2L)
  )
}

print("All operations completed successfully!")
# Step 1: preprocess()
list[fasta_all, metadata_all] <- preprocess(data_path, extract_path, seed,
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

print("All operations completed successfully!")