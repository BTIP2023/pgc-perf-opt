# File: kmer-analysis.R

# INSTALL AND LOAD PACKAGES ################################

library(datasets)  # Load base packages manually

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman", repos = "https://cran.case.edu")
library(pacman)

# Use pacman to load add-on packages as desired
# First call are standard packages for the project
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, psych,
               rio, rmarkdown, shiny,
               stringr, tidyr, tidyverse)
# Second call are file-specific packages
pacman::p_load(ape, kmer, readr, lubridate, stringr, validate, gsubfn)

# Note: gsubfn is used to destructure more than one return value

# Load preprocess() function
source('code/preprocess.R') 

# LOAD DATA ################################################
# Assumption: tar filename format is country-variant-etc.
# Extract GISAID tars to data/GISAID/datasets/country-variant
# If data/GISAID/datasets already exist, do not do this routine
# Note: Each tar = {tsv, fasta}

# Note: All paths are relative to project root. 

# Get fastaAll and metaDataAll from sourced preprocess() function,
# preprocess() defaults to seed = 10 or stratSize = 100 if not provided.
# This function performs:
# 1. Data parsing and augmentation
# 2. Stratified random sampling
# 3. Sanitation
list[fastaAll, metaDataAll] <- preprocess('data/GISAID', 'data/GISAID/datasets',
                                          seed = 10, stratSize = 100,
                                          country_exposure = 'Philippines')

# At this point, fastaAll and metaDataAll are SR sampled, sanitized, and 1:1

# Consider line below if we still want an intermediate fasta file
# write.FASTA(fastaAll, 'data/fastaAll.fasta')

# DEF LOCAL FUNCTIONS ###########################################

# Returns a data frame of kmers
kmer <- function(fasta, metaData, k){
  kmers <- kcount(fasta, k = k)
  kmer_df <- data.frame(kmers)
  
  print(dim(kmer_df))
  print(dim(metaData))
  
  # Append meta
  kmer_df <- cbind(kmer_df, metaData)
  
  return(kmer_df)
}

# WORK WITH DATA ###########################################
kmer_list = list(3,5,7)

for (k in kmer_list) {
  print(sprintf("Performing %d-mer analysis...", k))
  kmers <- kcount(fastaAll, k = k)
  kmer_df <- data.frame(kmers)
  
  # Append meta
  kmer_df <- cbind(kmer_df, metaDataAll)
  
  # Write to a csv file in data/kmers
  # Rewrites file if it already exists
  outputDir <- paste('data/kmers', sprintf("kmer_%d.csv", k), sep='/')
  print(paste("Writing kmer data to", outputDir))
  write.csv(kmer_df, outputDir)
}

# CLEAN UP #################################################

# Clear environment
rm(list = ls())

# Clear packages
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base

# Clear plots but only if there IS a plot
# while (!is.null(dev.list())) dev.off()

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)