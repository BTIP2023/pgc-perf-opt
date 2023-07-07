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

# Merge all extracted fasta and tsv files
# fasta contains sequence, while tsv contains meta
omicron_sub = c('ba275', 'xbb', 'xbb_1.5', 'xbb_1.16', 'xbb1.91')

fastas <- list.files(extractPath, recursive = TRUE, pattern = '.+\\.fasta')
tsvs <- list.files(extractPath, recursive = TRUE, pattern = '.+\\.tsv')
nfiles <- length(fastas)

fastaAll <- data.frame()
metaDataAll <- data.frame()

# TODO: Alternative approach: append fasta list as column to metadata
for (i in 1:nfiles) {
  fastaPath <- paste(extractPath, fastas[i], sep='/')
  tsvPath <- paste(extractPath, tsvs[i], sep='/')
  variant <- str_match(fastaPath, pattern = '(?<=-).*(?=\\/)')
  if (variant %vin% omicron_sub) {
    variant <- 'Omicron Sub'
  }
  variant <- str_to_title(variant)
  print(paste("Reading", fastaPath))
  print(paste("Reading", tsvPath))
  
  # Merge fasta file with accumulator
  fasta <- read.FASTA(fastaPath)
  fastaAll <- c(fastaAll, fasta)
  
  # Merge metaData file with accumulator
  # Defer sanitation after random sampling so fasta and metaData maintains 1:1
  metaData <- as.data.frame(read_tsv(tsvPath,
                                     col_select = c(1,5,10,11,12,14,16,17,19)))
  # Not removing raw date as I believe it is useful for sorting or can be parsed on an as needed basis
  metaData <- metaData %>%
    dplyr::mutate(year = lubridate::year(date),
                  month = lubridate::month(date),
                  day = lubridate::day(date),
                  variant = as.character(variant))
  
  metaData$length <- as.integer(metaData$length)
  metaData$age <- as.integer(metaData$age)
  metaData$year <- as.integer(metaData$year)
  metaData$month <- as.integer(metaData$month)
  metaData$day <- as.integer(metaData$day)
  
  metaDataAll <- bind_rows(metaDataAll, metaData)
}

rm(fasta)
rm(metaData)

# May put merging of fastaAll and metaDataAll here

# At this point, fastaAll and metaDataAll contains the needed data
# Now do stratified random sampling of <sampleSize> samples
seed = 10         # seed for random number generator
sampleSize = 100  # sample size per stratum
set.seed(seed)

# Note: append row names to column for fasta sampling
# Drop before exporting!
metaGrouped <- metaDataAll %>%
  dplyr::group_by(variant) %>%
  tibble::rownames_to_column()

# Do not preserve grouping structure (below only) to avoid NULL groups
droppedVariants <- filter(metaGrouped, n() < sampleSize)
metaGrouped <- filter(metaGrouped, n() >= sampleSize) %>%
  sample_n(sampleSize)
metaDataAll <- bind_rows(metaGrouped, droppedVariants)

rm(droppedVariants)
rm(metaGrouped)

set.seed(NULL)  # reset seed (rest of code is true random)

idxs <- as.integer(metaDataAll$rowname)
fastaAll <- fastaAll[idxs]

metaDataAll = subset(metaDataAll, select = -c(rowname) )  # drop rowname column

# Drop rows with NA values and type mismatches
# For dropping, get the idxs of the dropped rows and also drop them in fastaAll
drop_idxs1 <- which(is.na(metaDataAll), arr.ind=TRUE)[,1]
drop_idxs2 <- c(which(is.numeric(metaDataAll$sex)),
                which(!(metaDataAll$sex %vin% list("Male", "Female"))))
drop_idxs3 <- which(lengths(fastaAll) == 0)
drop_idxs <- unique(c(drop_idxs1, drop_idxs2, drop_idxs3))

fastaAll <- fastaAll[is.na(pmatch(1:nrow(fastaAll), drop_idxs))]
metaDataAll <- metaDataAll[is.na(pmatch(1:nrow(metaDataAll), drop_idxs)),]

# At this point, fastaAll and metaDataAll are sanitized and 1:1
# write.FASTA(fastaAll, 'data/fastaAll.fasta')

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
  
  # Write to a csv file in data/kmers-updated
  # Rewrites file if it already exists
  outputDir <- paste('data/kmers-updated', sprintf("kmer_%d.csv", k), sep='/')
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