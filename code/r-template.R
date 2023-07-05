# File: r-template.R

# INSTALL AND LOAD PACKAGES ################################

library(datasets)  # Load base packages manually

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman", repos = "https://cran.case.edu")
library(pacman)

# Use pacman to load add-on packages as desired
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, psych,
               rio, rmarkdown, shiny, 
               stringr, tidyr, tidyverse)

# LOAD DATA ################################################



# WORK WITH DATA ###########################################



# CLEAN UP #################################################

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

