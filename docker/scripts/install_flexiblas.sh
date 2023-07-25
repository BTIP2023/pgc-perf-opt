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

# For Flexiblas
apt_install \
    debhelper \
    cmake \
    gfortran \
    gcc \
    libatlas-base-dev \
    libopenblas-dev \
    build-essential \
    fakeroot

# More lib dependencies
apt_install libopenblas-openmp-dev \
    libopenblas-serial-dev \
    libopenblas64-openmp-dev \
    libopenblas64-pthread-dev \
    libopenblas64-serial-dev \
    libblis-openmp-dev \
    libblis-pthread-dev \
    libblis-serial-dev \
    libblis64-openmp-dev \
    libblis64-pthread-dev \
    libblis64-serial-dev \
    libblas64-dev \
    liblapack64-dev \
    libmkl-dev

# update-alternatives --config libblas.so
# update-alternatives --config libblas.so.3

# update-alternatives --config liblapack.so
# update-alternatives --config liblapack.so.3

echo -e "\nInstall Flexiblas packages, done!"
