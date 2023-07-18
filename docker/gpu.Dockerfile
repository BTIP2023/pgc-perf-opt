# escape=`

# Side objective: Perf eval/opt wrt NVIDIA CUDA implementation.

# Base Image: rocker/ml:4.3.1
# Adds CUDA support to rocker/tidyverse:4.3.1
# Image stack:
# - rocker/ml:4.3.1
# - rocker/cuda:4.3.1
# - rocker/tidyverse:4.3.1 (implicit)
# - rocker/rstudio:4.3.1 (implicit)
# - rocker/r-ver:4.3.1
FROM rocker/ml:4.3.1

# Copy repository (might have to adjust for diff. obj)
COPY . /home/src/pgc-perf-opt

# Bind a volume for development (might have to adjust for diff. obj)
# Notes: Container has a /home/rstudio directory.

# Building (from project root):
# docker build -t pgc-perf-opt/gpu . -f ./docker/gpu.Dockerfile

### Python
# Comes with Python 3.10.6 with base packages via python3.
# Install the rest of the dependencies manually via pip.
# In particular, install plotly.

### R and RStudio
# Comes with R 4.3.1 and RStudio Server, bundled with tidyverse.
# RStudio Server is available via exposed local port, configure creds.

# Building: docker build -f 

### Python
# This also already includes Python 3.10.6 with base packages via python3.
# Install the rest of the dependencies manually via pip.
# In particular, install plotly.

### R and RStudio
# Bundled with tidyverse.
# RStudio Server is available via exposed local port, configure creds.

# Usage:
# https://rocker-project.org/images/versioned/cuda.html
# Same as that of rocker/tidyverse, with a few exceptions:

### R console
# CPU-only
# docker run --rm -ti rocker/cuda
# Machines with nvidia-docker and GPU support
# docker run --gpus all --rm -ti rocker/cuda

### RStudio Server instance
# CPU-only
# docker run -p 8787:8787 rocker/ml
# Machines with nvidia-docker and GPU support
# docker run --gpus all -p 8787:8787 rocker/ml

### Switching Python versions through R
# reticulate::install_miniconda()
# reticulate::conda_install(packages = "python=3.7")

# Ref: https://rocker-project.org/images/
