#!/bin/bash

# Bash script to install project's base dependencies

set -e

## build ARGs
NCPUS=${NCPUS:--1}

# a function to install apt packages only if they are not installed
function apt_install() {
    if ! dpkg -s "$@" >/dev/null 2>&1; then
        if [ "$(find /var/lib/apt/lists/* | wc -l)" = "0" ]; then
            apt-get update
        fi
        apt-get install -y --no-install-recommends "$@"
    fi
}

# Install Linux utility libraries
# Defer curl installation for compatibility with install_github()
apt_install \
    tree \
    jq \
    htop \
    texinfo \
    nano \
    less \
    python3-pip

# R packages installation, removed redundancies from install_tidyverse.sh
install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
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
    reticulate

# Auxiliary R packages (more complicated installs)
Rscript ./docker/scripts/install_pgc_aux.R
python3 -m pip install plotly

# Install deferred curl
apt-get update && apt-get upgrade && apt-get install curl -y
apt install cmake -y # install for factoextra

# For kmer-analysis.R and sources
install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    ape \
    kmer \
    validate \
    gsubfn \
    seqinr

# For dim-reduce.R and sources
install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    umap \
    htmlwidgets \
    factoextra \
    scales \
    Rtsne \
    tsne \
    RColorBrewer \
    ggfortify

# For clustering and sources
install2.r --error --skipmissing --skipinstalled -n "$NCPUS" \
    ggdendro \
    dendextend \
    cluster \
    colorspace

# Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## Strip binary installed lybraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
strip /usr/local/lib/R/site-library/*/libs/*.so

echo -e "\nInstall PGC Perf Opt base packages, done!"