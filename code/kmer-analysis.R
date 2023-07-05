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
pacman::p_load(ape, kmer)

# LOAD DATA ################################################
#Function computes the kmer of given length 
kmer_df = function(filePath,variant,k){
  kmer = read.FASTA(filePath)
  kmer3 = kcount(kmer, k = k)
  target = rep((variant),length(kmer))
  kmer_df1 = data.frame(kmer3, target)
  return(kmer_df1)
}

# WORK WITH DATA ###########################################




#Different lengths of kmers to be used 
kmer_list = list(3,5,7)

for(k in kmer_list)
{
  alpha = kmer_df('alpha.fasta','alpha',k)
  beta = kmer_df('beta.fasta','beta',k)
  gamma = kmer_df('gamma.fasta','gamma',k)
  omicron = kmer_df('omicron.fasta','omicron',k)
  delta = kmer_df('delta.fasta','delta',k)
  
  data = bind_rows(alpha, beta, delta, gamma, omicron)
  outputName = sprintf("covid_kmer_%d.csv",k)
  #store the combined data in a csv file 
  write.csv(data, outputName)
}

# CLEAN UP #################################################

# Reset wd to project root
setwd(getwd())

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base

# Clear plots but only if there IS a plot
while (!is.null(dev.list())) dev.off()

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)