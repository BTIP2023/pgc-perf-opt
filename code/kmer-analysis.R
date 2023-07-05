# Reference: https://github.com/ai-covariants/analysis-mutations
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
pacman::p_load(ape, kmer, readr, lubridate, stringr)

# LOAD DATA ################################################
# Assumption: tar filename format is country-variant-etc.
# Extract GISAID tars to data/GISAID/datasets/country-variant
# If data/GISAID/datasets already exist, do not do this routine
# Note: Each tar = {tsv, fasta}

# GISAID data path from root
datapath <- 'data/GISAID'
# GISAID data extraction path after getting untarred
extractpath <- 'data/GISAID/datasets'

if (dir.exists(extractpath)) {
  print("GISAID data already extracted from tar archives.")
} else {
  tars <- list.files(datapath, pattern = '.+\\.tar')
  for (filename in tars) {
    subdir <- str_match(filename, pattern='[^-]+-[^-]+')
    print(paste("Extracting to:", subdir))
    untar(paste(datapath, filename, sep='/'),
          exdir = paste(extractpath, subdir, sep='/'),
          verbose = FALSE)
  }
}

# DEFINE FUNCTIONS ###########################################

# Function computes the kmer of given length, returns a data frame
kmer_df <- function(fastaPath,tsvPath,k){
  fastaPath <- paste('data/GISAID/datasets', fastaPath, sep='/')
  tsvPath <- paste('data/GISAID/datasets', tsvPath, sep='/')
  fasta <- read.FASTA(fastaPath)
  tsv <- read_tsv(tsvPath)
  print("Calling problems")
  problems(fasta)
  problems(tsv)
  # kmer3 = kcount(kmer, k = k)
  # target = rep((variant),length(kmer))
  # kmer_df1 = data.frame(kmer3, target)
  # return(kmer_df1)
}

# WORK WITH DATA ###########################################

# Different lengths of kmers to be used 
kmer_list <- list(3,5,7)

for(k in kmer_list) {
  # fasta contains sequence, while tsv contains metadata
  fastas <- list.files('data/GISAID/datasets', recursive = TRUE, pattern = '.+\\.fasta')
  tsvs <- list.files('data/GISAID/datasets', recursive = TRUE, pattern = '.+\\.tsv')
  nfiles <- length(fastas)
  for (i in 1:nfiles) {
    fastapath <- paste(extractpath, fastas[i], sep='/')
    tsvpath <- paste(extractpath, tsvs[i], sep='/')
    print(paste("Reading", fastapath))
    print(paste("Reading", tsvpath))
    # kmer_df(fastas[i],tsvs[i],k)
  }

  # data = bind_rows(alpha, beta, delta, gamma, omicron)
  # outputName = sprintf("covid_kmer_%d.csv",k)
  # #store the combined data in a csv file 
  # write.csv(data, outputName)
}

# CLEAN UP #################################################

# Clear environment
# rm(list = ls()) 

# Clear packages
# p_unload(all)  # Remove all add-ons
# detach("package:datasets", unload = TRUE)  # For base

# Clear plots but only if there IS a plot
# while (!is.null(dev.list())) dev.off()

# Clear console
# cat("\014")  # ctrl+L

# Clear mind :)