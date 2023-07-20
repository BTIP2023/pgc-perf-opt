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
               highcharter)
install_github("vqv/ggbiplot", upgrade = FALSE, quiet = TRUE)
pacman::p_load(ggbiplot)

# Set highcharter options
options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)))

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
seed <- 1234
stamp <- get_time()
kmer_list <- c(3, 5, 7)

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

#Step 4: AGNES Clustering by Variant
for (k in kmer_list) {
  dendrogram_create_variant(k, data_path_kmers, results_path_agnes)
}

#Step 5: AGNES Clustering by Region
for (k in kmer_list){
  dendrogram_create_region(k, data_path_kmers, results_path_agnes)
}

print("All operations completed successfully!")

# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages (unloading them before another adds another compat check)
p_unload(all)  # Remove all add-ons

# Clear plots but only if there IS a plot
while (!is.null(dev.list())) dev.off()

# Clear console
# cat("\014")  # ctrl+L

# Clear mind :)