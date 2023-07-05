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
  # get list of tar files in data, tar = {tsv, csv}

  # fasta contains sequence, while tsv contains metadata
  # fastas <- list.files('data/GISAID/datasets', pattern = '.+\\.fasta')
  # tsvs <- list.files('data/GISAID/datasets', pattern = '.+\\.tsv')
  # nfiles <- length(fastas)  # equivalent to length(tsvs)
  # # Problems 4, 12
  # i <- 10
  # print(paste("Reading", fastas[i]))
  # print(paste("Reading", tsvs[i]))
  # fastaPath <- paste('data/GISAID/datasets', fastas[i], sep='/')
  # tsvPath <- paste('data/GISAID/datasets', tsvs[i], sep='/')
  # fasta <- read.FASTA(fastaPath)
  # metadata <- as.data.frame(read_tsv(tsvPath, col_select = c(1,2,5,10,11,12,19)))
  # metadata <- metadata %>%
  #   dplyr::mutate(year = lubridate::year(date), 
  #                 month = lubridate::month(date), 
  #                 day = lubridate::day(date))
  # 
  # print("Calling problems")
  # problems(fasta)
  # problems(metadata)
  # for (i in 1:nfiles) {
    # print(paste("Reading", fastas[i]))
    # print(paste("Reading", tsvs[i]))
    # kmer_df(fastas[i],tsvs[i],k)
  # }
  # alpha = kmer_df('alpha.fasta','alpha',k)
  # beta = kmer_df('beta.fasta','beta',k)
  # gamma = kmer_df('gamma.fasta','gamma',k)
  # omicron = kmer_df('omicron.fasta','omicron',k)
  # delta = kmer_df('delta.fasta','delta',k)
  # 
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