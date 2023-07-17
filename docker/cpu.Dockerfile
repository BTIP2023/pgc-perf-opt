# escape=`

# Main objective: Perf eval/opt wrt Spectre/Meltdown CPU patches.

# Base Image: rocker/tidyverse:4.3.1
# Slim image for benchmarking and developing R workflows.
# Image stack:
# - rocker/tidyverse:4.3.1
# - rocker/rstudio:4.3.1
# - rocker/r-ver:4.3.1
FROM rocker/tidyverse:4.3.1

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
