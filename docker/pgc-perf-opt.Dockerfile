# Base Image for PGC Performance Optimization Group.
# Natively supports: Performance eval/opt wrt Spectre/Meltdown CPU patches. 

# Base Image: rocker/tidyverse:4.3.1
# Slim image for benchmarking and developing R workflows.
# By slim, we mean a minimized Ubuntu 22.04 (no man pages, etc.).
# Image stack:
# - rocker/tidyverse:4.3.1
# - rocker/rstudio:4.3.1
# - rocker/r-ver:4.3.1
FROM rocker/tidyverse:4.3.1

LABEL organization="Philippine Genome Center - Core Facility for Bioinformatics" \
      description="Base image for PGC Performance Optimization Group. \
        Natively supports: Performance eval/opt wrt Spectre/Meltdown CPU patches." \
      maintainer="Yenzy Urson S. Hebron <yshebron@up.edu.ph>"

# Copy local repository snapshot (see .dockerignore)
# Notes: Container has a /home/rstudio directory.
#   - Comment out presentations/ in .dockerignore if you wish
#   - to work on presentations in the container.
COPY . /home/rstudio/pgc-perf-opt

RUN --mount=type=cache,target=/var/cache/apt apt-get update
RUN apt-get install -y libpython2-dev
RUN apt-get install -y libpython3-dev

# Change working directory to project root
WORKDIR /home/rstudio/pgc-perf-opt

# Install project base R, Python, and system-level dependencies
RUN ./docker/scripts/install_pgc_base.sh

# Auxiliary Installations of R packages (and Python by necessity)
RUN Rscript ./docker/scripts/install_pgc_aux.R

# Set RStudio Server working directory to workspace
RUN echo "setwd('/home/rstudio/pgc-perf-opt')" > /home/rstudio/.Rprofile

# Make RStudio Server use reticulated environment on startup
# /home/rstudio/.local/share/r-miniconda/envs/r-reticulate
RUN echo "library(reticulate)" >> /home/rstudio/.Rprofile
RUN echo "Sys.setenv(RETICULATE_MINICONDA_PATH = '/home/rstudio/r-miniconda')" >> /home/rstudio/.Rprofile
RUN echo "Sys.setenv(RETICULATE_MINICONDA_ENABLED = 'true')" >> /home/rstudio/.Rprofile
RUN echo "reticulate::py_config()" >> /home/rstudio/.Rprofile
RUN echo "reticulate::use_miniconda('reticulate')" >> /home/rstudio/.Rprofile

# Note: Currently using bind mounts instead of volumes for dev convenience
# and because of skill issue.

### Python
# Comes with Python 3.10.6 with base packages via python3.
# However, we only use Python 3.9.16, installed via reticulate to ~/r-miniconda.

### R and RStudio
# Comes with R 4.3.1 and RStudio Server, bundled with tidyverse.
# RStudio Server is available via exposed local port, configure creds.

# Building (from project root):
# docker build -t pgc-perf-opt/cpu . -f ./docker/cpu.Dockerfile

# Usage:
# https://rocker-project.org/images/versioned/rstudio.html
# Password: docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 pgc-perf-opt/cpu
# Root: docker run --rm -ti -e ROOT=true -p 8787:8787 pgc-perf-opt/cpu
# Disable Authentication: docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 pgc-perf-opt/cpu
# User ID and Group ID: docker run --rm -ti -e USERID=1001 -e GROUPID=1001 -p 8787:8787 pgc-perf-opt/cpu
# To not expose: -p 127.0.0.1:8787:8787

# Ref: https://rocker-project.org/images/
