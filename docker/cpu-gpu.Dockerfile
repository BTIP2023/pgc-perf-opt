# escape=`

# Dockerfile for GPU performance side objective.

# Base Image: rocker/ml:4.3.1
# Adds CUDA support to rocker/tidyverse:4.3.1
# Image stack:
# - rocker/ml:4.3.1
# - rocker/cuda:4.3.1
# - rocker/tidyverse:4.3.1 (implicit)
# - rocker/rstudio:4.3.1 (implicit)
# - rocker/r-ver:4.3.1
FROM rocker/ml:4.3.1

### Python
# This also already includes Python 3.10.6 with base packages via python3.
# Install the rest of the dependencies manually via pip.
# In particular, install plotly.

### R and RStudio
# Bundled with tidyverse.
# RStudio Server is available via exposed local port, configure creds.

# Usage:
# https://rocker-project.org/images/versioned/rstudio.html
# https://rocker-project.org/images/versioned/cuda.html

# Ref: https://rocker-project.org/images/
