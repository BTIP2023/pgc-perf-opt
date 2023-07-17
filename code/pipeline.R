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
               ape, kmer, readr, validate, gsubfn, seqinr)

# validate used for %vin% operator 
# gsubfn used to destructure more than one return value
# Note: Divisive k-means clustering available via kmer::cluster

# LOAD SOURCES #############################################
source("code/R/preprocess.R")
source("code/R/helper.R")
source("code/R/kmer-analysis.R")

# SET PARAMETERS ###########################################
# pipeline.R general parameters
seed <- 1234
stamp <- get_time()

# TODO: Params for until which step of the pipeline should be run.
# TODO: Which pipeline step is to be solely performed (for dev).

# preprocess.R::preprocess() parameters
data_path <- "data/GISAID"
extract_path <- "data/GISAID/datasets"
strat_size <- 100
country_exposure <- "Philippines"
write_fastacsv <- TRUE

# kmer-analysis.R::get_kmers() parameters
kmer_list <- c(3, 5, 7)

# RUN PIPELINE #############################################

# Step 1: preprocess()
list[fasta_all, metadata_all] <- preprocess(data_path, extract_path, seed,
                                            strat_size, country_exposure,
                                            write_fastacsv, stamp)

# Step 2: get_kmers()
# Suggestion: Move get_kmers() for k loop to pipeline.R so as to
# integrate get_kmers() with the for k loops of dim-reduce and clustering.
# get_kmers will therefore be returning the kmer data frames, minimizing
# the need for repeat reading of csvs.
get_kmers(fasta_all, metadata_all, kmer_list, stamp)

# Step 3: ...

print("Operations completed successfully!")