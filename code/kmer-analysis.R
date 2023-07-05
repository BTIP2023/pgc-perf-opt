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
pacman::p_load(ape, kmer, readr, lubridate, stringr, validate)

# LOAD DATA ################################################
# Assumption: tar filename format is country-variant-etc.
# Extract GISAID tars to data/GISAID/datasets/country-variant
# If data/GISAID/datasets already exist, do not do this routine
# Note: Each tar = {tsv, fasta}

# GISAID data path from root
dataPath <- 'data/GISAID'
# GISAID data extraction path after getting untarred
extractPath <- 'data/GISAID/datasets'

if (dir.exists(extractPath)) {
  print("GISAID data already extracted from tar archives.")
} else {
  tars <- list.files(dataPath, pattern = '.+\\.tar')
  for (fileName in tars) {
    subdir <- str_match(fileName, pattern='[^-]+-[^-]+')
    print(paste("Extracting to:", subdir))
    untar(paste(dataPath, fileName, sep='/'),
          exdir = paste(extractPath, subdir, sep='/'),
          verbose = FALSE)
  }
}

# DEFINE FUNCTIONS ###########################################
omicron_sub = c('ba275', 'xbb', 'xbb_1.5', 'xbb_1.16', 'xbb1.91')

# Function computes the kmer of given length, returns a data frame
kmer_df <- function(fastaPath, tsvPath, variant, k){
  fasta <- read.FASTA(fastaPath)
  meta <- as.data.frame(read_tsv(tsvPath, col_select = c(1,2,5,10,11,12,14,16,17,19)))
  # Not removing raw date as I believe it is useful for sorting or can be parsed on an as needed basis
  meta <- meta %>%
    dplyr::mutate(year = lubridate::year(date),
                  month = lubridate::month(date),
                  day = lubridate::day(date))
  kmers <- kcount(fasta, k = k)
  target <- rep((variant),length(kmers))
  kmer_df <- data.frame(kmers, target)
  
  # Append meta
  kmer_df <- cbind(kmer_df, meta)
  
  return(kmer_df)
}

# WORK WITH DATA ###########################################

# Different lengths of kmers to be used 
kmer_list <- list(3)

# Instantiate alpha accumulators

for(k in kmer_list) {
  # fasta contains sequence, while tsv contains meta
  fastas <- list.files('data/GISAID/datasets', recursive = TRUE, pattern = '.+\\.fasta')
  tsvs <- list.files('data/GISAID/datasets', recursive = TRUE, pattern = '.+\\.tsv')
  nfiles <- length(fastas)
  # instantiate accumulator
  data <- data.frame()
  for (i in 1:nfiles) {
    fastaPath <- paste(extractPath, fastas[i], sep='/')
    tsvPath <- paste(extractPath, tsvs[i], sep='/')
    variant <- str_match(fastaPath, pattern = '(?<=-).*(?=\\/)')
    if (variant %vin% omicron_sub) {
      variant <- 'omicron_sub'
    }
    print(paste("Reading", fastaPath))
    print(paste("Reading", tsvPath))
    temp <- kmer_df(fastaPath, tsvPath, variant, k)
    data <- bind_rows(data, temp)
  }
  outputName = sprintf("kmer_%d.csv", k)
  # #store the combined data in a csv file 
  write.csv(paste('data/kmers'), outputName)
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