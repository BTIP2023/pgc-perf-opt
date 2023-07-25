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

echo -e "\nInstall Flexiblas packages, done!"
