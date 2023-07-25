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

ARG NCPUS=-1

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

# R packages installation, removed redundancies from install_tidyverse.sh
RUN install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    pacman \
    plyr \
    GGally \
    ggthemes \
    ggvis \
    plotly \
    psych \
    rio \
    markdown \
    shiny \
    devtools \
    microbenchmark \
    reticulate \
    highcharter

# For kmer-analysis.R and sources
RUN install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    ape \
    kmer \
    validate \
    gsubfn \
    seqinr

# For dim-reduce.R and sources
RUN install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    umap \
    htmlwidgets \
    factoextra \
    scales \
    Rtsne \
    tsne \
    RColorBrewer \
    ggfortify

# For clustering and sources
RUN install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    ggdendro \
    dendextend \
    cluster \
    colorspace

# Upgrade for curl compat with R remotes, then install deferred packs
RUN apt-get upgrade -y
RUN apt-get install curl -y
# Install for factoextra
RUN apt-get install cmake -y

# Flexiblas installation
RUN ./docker/scripts/install_flexiblas.sh
RUN curl https://csc.mpi-magdeburg.mpg.de/mpcsc/software/flexiblas/flexiblas-3.3.0.tar.gz | tar -xz
WORKDIR /home/rstudio/pgc-perf-opt/flexiblas-3.3.0
RUN fakeroot dpkg-buildpackage -us -uc
RUN dpkg -i ../libflexiblas-*.deb
RUN debian/rules clean
WORKDIR /home/rstudio/pgc-perf-opt

# Clean up of install temps
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/downloaded_packages

## Strip binary installed lybraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
RUN strip /usr/local/lib/R/site-library/*/libs/*.so

# Auxiliary Installations of R packages (and Python by necessity)
# Done in root R, so no need for clean-up
RUN Rscript ./docker/scripts/install_pgc_aux.R

# Set RStudio Server working directory to workspace
RUN echo "setwd('/home/rstudio/pgc-perf-opt')" > /home/rstudio/.Rprofile

# Make RStudio Server use reticulated environment on startup
# /home/rstudio/.local/share/r-miniconda/envs/r-reticulate
RUN echo "library(reticulate)" >> /home/rstudio/.Rprofile
RUN echo "Sys.setenv(RETICULATE_MINICONDA_PATH = '/home/rstudio/r-miniconda')" >> /home/rstudio/.Rprofile
RUN echo "Sys.setenv(RETICULATE_MINICONDA_ENABLED = 'true')" >> /home/rstudio/.Rprofile
RUN echo "reticulate::py_config()" >> /home/rstudio/.Rprofile
RUN echo "reticulate::use_miniconda('r-reticulate')" >> /home/rstudio/.Rprofile

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
