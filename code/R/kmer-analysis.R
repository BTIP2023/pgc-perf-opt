# File: kmer-analysis.R

# INSTALL AND LOAD PACKAGES ################################

library(datasets)  # Sample datasets for testing

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman", repos = "https://cran.case.edu")
library(pacman)

### ATTN: IN LINUX SYSTEMS, CONSULT README FOR ADDITIONAL PREREQUISITES
### BEFORE RUNNING ANY SCRIPT. ISSUE: tidyverse installation.
### This cannot be scripted because this requires sudo priveleges.

# Install xml2 in advance to prep for tidyverse installation in Linux.
# Note that in Windows RStudio, this is installed by default.
# If you're getting xml2 errors on Windows, you broke something lol.
# Ref: https://medium.com/@jamie84mclaughlin/installing-r-and-the-tidyverse-on-ubuntu-20-04-60170020649b
if (pacman::p_detectOS() == 'Linux' & !pacman::p_exists(xml2, local=T)) {
  install.packages('xml2', dependencies = T, INSTALL_opts = c("--no-lock"))
  pacman::p_load(xml2)
}

# Use pacman to load add-on packages as desired
# First p_load call are standard packages for the project
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse)

# Second p_load call are file-specific packages
pacman::p_load(ape, kmer, readr, lubridate, stringr, validate, gsubfn, seqinr)

# Note: gsubfn is used to destructure more than one return value

# LOAD SOURCES #############################################
source('code/R/preprocess.R')
source('code/R/helper.R')
stamp <- timeString()     # get timestamp for file suffix

# LOAD DATA ################################################
# Note: All paths are relative to project root. 

### SET THESE PARAMETERS ###
seed <- 10
stratSize <- 100
country_exposure <- 'Philippines'
write_fastacsv <- TRUE

# Get fastaAll and metaDataAll from sourced preprocess() function
list[fastaAll, metaDataAll] <- preprocess('data/GISAID', 'data/GISAID/datasets',
                                          seed = seed, stratSize = stratSize,
                                          country_exposure = country_exposure,
                                          write_fastacsv = write_fastacsv,
                                          stamp = stamp)

# At this point, fastaAll and metaDataAll are SR sampled, sanitized, and 1:1

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
  outputDir <- paste('data/kmers',
                     sprintf("kmer_%d_%s.csv", k, stamp),
                     sep='/')
  print(paste("Writing kmer data to", outputDir))
  write.csv(kmer_df, outputDir)
}

# Write parameters used to text file
# TODO: Add runtime log capability to paramsLog for benchmarking
paramsLog(outputDir = 'data/kmers', filename = 'log.txt',
          paramString = sprintf("timestamp = %s\nseed = %d, stratSize = %d",
                                stamp, seed, stratSize))

# CLEAN UP #################################################

# Clear environment
rm(list = ls())

# Clear packages
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base

# Clear plots but only if there IS a plot
# while (!is.null(dev.list())) dev.off()

# Clear console
# cat("\014")  # ctrl+L

# Clear memory
gc()

print("Operation completed successfully!")