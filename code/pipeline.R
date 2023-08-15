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
pacman::p_load(plyr, GGally, ggthemes, ggvis, plotly, psych,
               htmlwidgets, rio, markdown, shiny, tidyverse,
               ape, seqinr, kmer, validate, gsubfn,
               Rtsne, tsne, umap, factoextra, scales,
               RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark, data.table, highcharter)
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
source("code/R/clustering.R")

# SET PARAMETERS ###########################################
# pipeline.R general parameters
# stamp <- [get_time()|str(strat_size)|NULL]
# if stamp = NULL, then generated files won't be timestamped
seed <- 1234
stamp <- get_time()
write_fastacsv <- TRUE
kmer_list <- c(3, 5, 7)
# strat_size: no. of samples per stratum. Current nrow(data) = 24671.
# Also consider using sample_frac for proportionate allocation.
# Note that valid strat_size will only be those with corresponding
# files in `data/interm` and `data/kmers`
strat_size <- 100
# include_plots: whether to include plots in execution
include_plots <- TRUE

# preprocess.R::get_sample() parameters
gisaid_data_path <- "data/GISAID"
gisaid_extract_path <- "data/GISAID/datasets"
country_exposure <- "Philippines"

# preprocess.R::auxiliary parameters
interm_write_path <- "data/interm"
compile_write_path <- "data/overview"

# dim-reduce.R::dim_reduce() parameters
kmers_data_path <- "data/kmers"
dr_write_path <- "results/dim-reduce/R"
tsne_perplexity <- 40
tsne_max_iter <- 1000
tsne_initial_dims <- 50
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1
dr_color <- "variant"
dr_shape <- "year"

# dim-reduce.R::dim_reduce() filters
# Origs: (variant, c("Omicron", "Omicron Sub"), year, c(2023))
# Other options for factors: date, year, region_exposure, country_exposure,
# division_exposure, division_code, ph_region, age, age_group,
# sex, pangolin_lineage, variant
# Values for values can be scalar or vectors.
dr_factor1 <- NULL
dr_values1 <- NULL
dr_factor2 <- NULL
dr_values2 <- NULL

# clustering-x.R::dendogram_create_x() parameters
agnes_write_path <- "results/dendrogram"

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
make_treemaps(metadata_all, compile_write_path, stamp)

# Step 2: get_kmers()
for (k in kmer_list) {
  get_kmers(fasta_all, metadata_all, k, stamp)
}

# GET KMERS FROM PRE-WRITTEN FILES (depends on strat_size)
# kmers is list of kmer dataframes
kmers <- list()
for (i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  k_path <- sprintf("%s/kmer_%d_%d.csv", kmers_data_path, k, strat_size)
  message(sprintf("Reading %s for later... ", k_path), appendLF = FALSE)
  kmers[[i]] <- utils::read.csv(k_path)
  message("DONE!")
}

# Step 2.5: generate_kmer_heatmap and generate_kmer_wordcloud
# for (i in 1:length(kmer_list)) {
#   k <- kmer_list[i]
#   generate_kmer_heatmap(kmers[[i]], kmers_data_path, k)
#   generate_kmer_wordcloud(kmers[[i]], kmers_data_path, k)
# }

for (i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  generate_kmer_wordcloud(kmers[[i]], kmers_data_path, k)
}

# Step 3: dim_reduce()
for (i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  dim_reduce(k, kmers[[i]], dr_write_path,
             tsne_seed = seed, tsne_perplexity,
             tsne_max_iter, tsne_initial_dims,
             umap_seed = seed, umap_n_neighbors,
             umap_metric, umap_min_dist, color = dr_color, shape = dr_shape,
             filter1_factor = dr_factor1, filter1_values = dr_values1,
             filter2_factor = dr_factor2, filter2_values = dr_values2,
             include_plots = include_plots)
}

#Step 4: AGNES Clustering by Variant
for (i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  dendrogram_create_variant(k, kmers[[i]], agnes_write_path, include_plots)
}

#Step 5: AGNES Clustering by Region
for (i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  dendrogram_create_region(k, kmers[[i]], agnes_write_path, include_plots)
}

message("All operations completed successfully!")

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