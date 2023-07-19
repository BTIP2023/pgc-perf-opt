# escape=`

# Side objective: Perf eval/opt wrt NVIDIA CUDA implementation.
# Also serves as: Base GPU compute Image for PGC Performance Optimization Group.

# Base Image: rocker/ml:4.3.1
# Adds CUDA support to rocker/tidyverse:4.3.1
# Image stack:
# - rocker/ml:4.3.1
# - rocker/cuda:4.3.1
# - rocker/tidyverse:4.3.1 (implicit)
# - rocker/rstudio:4.3.1 (implicit)
# - rocker/r-ver:4.3.1
FROM rocker/ml:4.3.1

LABEL organization="Philippine Genome Center - Core Facility for Bioinformatics" `
      description="For bioinformatics performance evaluation and `
        optimization with respect to Spectre/Meltdown CPU patches. `
        Also serves as Base GPU Image for PGC Performance Optimization Group"

# Copy local repository snapshot
# Notes: Container has a /home/rstudio directory.
#   - Comment out presentations/ in .dockerignore if you wish
#   - to work on presentations in the container.
COPY . /home/rstudio/pgc-perf-opt

# Change working directory to project root
WORKDIR /home/rstudio/pgc-perf-opt

# Bind a volume for development
VOLUME ["/home/rstudio/pgc-perf-opt"]

# Install project base R, Python, and system-level dependencies
RUN ./docker/scripts/install_pgc_base.sh

### Python
# Comes with Python 3.10.6 with base packages via python3.
# Install the rest of the dependencies manually via pip.
# In particular, install plotly.

### R and RStudio
# Comes with R 4.3.1 and RStudio Server, bundled with tidyverse.
# RStudio Server is available via exposed local port, configure creds.

# Building (from project root):
# docker build -t pgc-perf-opt/gpu . -f ./docker/cpu.Dockerfile

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
