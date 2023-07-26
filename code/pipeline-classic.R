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
               microbenchmark, data.table,
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
strat_size <- 100
country_exposure <- "Philippines"

# preprocess.R::auxiliary parameters
interm_write_path <- "data/interm"
compile_write_path <- "data/overview"
treemaps_write_path <- "data/overview/treemaps"
heatmaps_write_path <- "data/overview/heatmaps"

# dim-reduce.R::dim_reduce() parameters
kmers_data_path <- "data/kmers"
dimreduce_write_path <- "results/dim-reduce/R"
tsne_perplexity <- 40
tsne_max_iter <- 1000
tsne_initial_dims <- 50
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1
color <- "variant"
shape <- "sex"

# dim-reduce.R::dim_reduce() filtering parameters - OPTIONAL
#factor1 <- "variant"
#values1 <- c("Omicron", "Omicron Sub")
#factor2 <- "year"
#values2 <- c("2023")

# clustering-x.R::dendogram_create_x() parameters
results_path_agnes <- "results/dendrogram"

# RUN PIPELINE #############################################

# Step 1: get_sample()
list[fasta_all, metadata_all] <- get_sample(gisaid_data_path,
                                            gisaid_extract_path,
                                            seed, strat_size,
                                            country_exposure)

# Step 1.5A: sanitize_sample()
metadata_all <- sanitize_sample(metadata_all)

# Step 1.5B: generate_interm()
# Note that at strat_size > nrow(Omicron), you'll be writing around 700MB
# of fasta_all_stamp.csv, so be cautious of generate_interm's space usage.
if (write_fastacsv)
  generate_interm(fasta_all, metadata_all, interm_write_path, stamp)

# Step 1.5C: compile_overview()
# compile_overview drops the submitting_lab and authors column
# after compilation, hence the reassignment to metadata_all.
metadata_all <- compile_overview(metadata_all, compile_write_path, stamp)

# Step 1.5D: make_treemaps()
# NOTE: The treemap() function in helper.R
# can generate any treemap you can think of, yeah!
make_treemaps(metadata_all, treemaps_write_path, stamp)

# Step 2: get_kmers()
for (k in kmer_list) {
  get_kmers(fasta_all, metadata_all, k, stamp)
}

# Step 2.5: generate_heatmap()
for (k in kmer_list){
  generate_heatmap(kmers_data_path, heatmaps_write_path, k)
}


# Step 3: dim_reduce()
for (k in kmer_list) {
  dim_reduce(k, kmers_data_path, dimreduce_write_path,
             tsne_seed = seed, tsne_perplexity,
             tsne_max_iter, tsne_initial_dims,
             umap_seed = seed, umap_n_neighbors,
             umap_metric, umap_min_dist, color = color, shape = shape,
             filter1_factor = factor1, filter1_values = values1, # OPTIONAL
             filter2_factor = factor2, filter2_values = values2) # OPTIONAL
}

#Step 4: AGNES Clustering by Variant
for (k in kmer_list) {
  dendrogram_create_variant(k, kmers_data_path, results_path_agnes)
}

#Step 5: AGNES Clustering by Region
for (k in kmer_list){
  dendrogram_create_region(k, kmers_data_path, results_path_agnes)
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