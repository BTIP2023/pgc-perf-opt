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
pacman::p_load(plyr, dplyr, GGally, ggplot2, ggthemes, ggvis,
               httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse,
               ape, kmer, readr, validate, gsubfn, seqinr,
               umap, htmlwidgets, factoextra, scales,
               Rtsne, tsne, RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark,
               highcharter, glue)
if (!require(ggbiplot))
  install_github("vqv/ggbiplot", upgrade = FALSE, quiet = TRUE)
pacman::p_load(ggbiplot)

# validate used for %vin% operator 
# gsubfn used to destructure more than one return value
# devtools supports install_github for installing ggbiplot
# Note: Divisive k-means clustering available via kmer::cluster

# LOAD SOURCES #############################################
source("code/R/helper.R")
source("code/R/preprocess.R")
source("code/R/kmer-analysis.R")
source("code/R/dim-reduce.R")
source("code/R/clustering-variant.R")
source("code/R/clustering-region.R")

# SET PARAMETERS ###########################################
# pipeline.R general parameters
# stamp <- [get_time():str|NULL]
# if stamp = "", then generated files won't be timestamped
seed <- 1234
stamp <- get_time()
write_fastacsv <- TRUE
kmer_list <- c(3, 5, 7)

# preprocess.R::get_sample() parameters
# strat_size: no. of samples per stratum. Current nrow(data) = 24671.
# Also consider using sample_frac for proportionate allocation.
gisaid_data_path <- "data/GISAID"
gisaid_extract_path <- "data/GISAID/datasets"
strat_size <- 25000
country_exposure <- "Philippines"

# preprocess.R::auxiliary parameters
interm_write_path <- "data/interm"
compile_write_path <- "data/overview"
treemaps_write_path <- "data/overview/treemaps"

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

# dim-reduce.R::dim_reduce() filtering parameters - OPTIONAL
factor1 <- "variant"
values1 <- c("Omicron", "Omicron Sub")
factor2 <- "year"
values2 <- c("2023")

# clusterting-x.R::dendogram_create_x() parameters
results_path_agnes <- "results/dendrogram"

# RUN BENCHMARK #############################################
# Benchmark parameters
bm_times <- 2L   # how many times should routine be evaluated
bm_units <- "seconds"
bm_log_path <- "benchmarks"
OS <- pacman::p_detectOS()
# valid values: ["ALL"|"SOME" (Linux only)|"NONE"]
mitigations <- "ALL"

# Benchmark Notes:
# get_sample: to start from extraction, delete data/GISAID/datasets/
# generate_interm always happens for this benchmark
# TODO: Add descriptions for my functions to the plot results
results <- microbenchmark(
  get_sample = list[fasta_all, metadata_all] <-
    get_sample(gisaid_data_path,
               gisaid_extract_path,
               seed, strat_size,
               country_exposure),
  sanitize_sample = metadata_all <- sanitize_sample(metadata_all),
  generate_interm = generate_interm(fasta_all, metadata_all,
                                    interm_write_path, stamp),
  compile_overview = metadata_all <- compile_overview(metadata_all,
                                                      compile_write_path,
                                                      stamp),
  make_treemaps = make_treemaps(metadata_all, treemaps_write_path, stamp),
  get_kmers_all = for (k in kmer_list) {
    get_kmers(fasta_all, metadata_all, k, stamp)
  },
  get_kmers_3 = get_kmers(fasta_all, metadata_all, 3, stamp),
  get_kmers_5 = get_kmers(fasta_all, metadata_all, 5, stamp),
  get_kmers_7 = get_kmers(fasta_all, metadata_all, 7, stamp),
  # dim_reduce = for (k in kmer_list) {
  #   dim_reduce(k, data_path_kmers, results_path_dimreduce,
  #              tsne_seed = seed, tsne_perplexity,
  #              tsne_max_iter, tsne_initial_dims,
  #              umap_seed = seed, umap_n_neighbors,
  #              umap_metric, umap_min_dist, col_name = target_col)
  # },
  # dim_reduce_3 = dim_reduce(3, data_path_kmers, results_path_dimreduce,
  #                           tsne_seed = seed, tsne_perplexity,
  #                           tsne_max_iter, tsne_initial_dims,
  #                           umap_seed = seed, umap_n_neighbors,
  #                           umap_metric, umap_min_dist, col_name = target_col),
  # dim_reduce_5 = dim_reduce(5, data_path_kmers, results_path_dimreduce,
  #                           tsne_seed = seed, tsne_perplexity,
  #                           tsne_max_iter, tsne_initial_dims,
  #                           umap_seed = seed, umap_n_neighbors,
  #                           umap_metric, umap_min_dist, col_name = target_col),
  # dim_reduce_7 = dim_reduce(7, data_path_kmers, results_path_dimreduce,
  #                           tsne_seed = seed, tsne_perplexity,
  #                           tsne_max_iter, tsne_initial_dims,
  #                           umap_seed = seed, umap_n_neighbors,
  #                           umap_metric, umap_min_dist, col_name = target_col),
  times = bm_times,
  unit = bm_units,
  control = list(order = "inorder", warmup = 2L)
)

print("All operations completed successfully!")

# Write hardware specs and parameters used to log.txt
message("Writing logs... ", appendLF = FALSE)
write_to_log(output_dir = "results", filename = "log.txt",
             log_string = sprintf("timestamp = %s\nseed = %d, strat_size = %d, k-value = %d",
                                  stamp, seed, strat_size, k))
message("DONE.")

# CLEAN UP #################################################

# Clear environment
# rm(list = ls()) 

# Clear packages (unloading them before another adds another compat check)
p_unload(all)  # Remove all add-ons

# Clear plots but only if there IS a plot
while (!is.null(dev.list())) dev.off()

# Clear console
# cat("\014")  # ctrl+L

# Clear mind :)