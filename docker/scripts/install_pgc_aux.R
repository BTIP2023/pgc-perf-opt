# Install other R tweaks

# ggbiplot: https://github.com/vqv/ggbiplot#installation
# Note that this uses libcurl, not curl (avoid conflict by not installing curl)
library(devtools)
devtools::install_github("vqv/ggbiplot")

# knitr 1.42 for presentations
library(remotes)
remotes::install_version("knitr", version = "1.42",
                          repos = "https://cran.r-project.org")

# save_image
# In Linux, so no need for
# `reticulate::conda_install("r-reticulate", "python-kaleido==0.1.*"")
# Can't install to /home/rstudio/.local/share/r-miniconda even on closing all connections.
# Seems like a bug in rocker/tidyverse source code wherein it won't install to ~/.local/share
# if it already has contents
install.packages("reticulate")
library(reticulate)
Sys.setenv(RETICULATE_MINICONDA_PATH = "/home/rstudio/r-miniconda")
reticulate::install_miniconda()
reticulate::conda_install("r-reticulate", "python-kaleido")
reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
