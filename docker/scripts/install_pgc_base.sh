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
# libglpk-dev is required for highcharter installation in Ubuntu
apt_install \
    tree \
    jq \
    htop \
    texinfo \
    nano \
    less \
    libglpk-dev

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
    reticulate \
    highcharter

echo -e "\nInstall PGC Perf Opt base apt packages, done!"
