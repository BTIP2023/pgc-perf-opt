# escape=`

# Main objective: Perf eval/opt wrt Spectre/Meltdown CPU patches.

# Base Image: rocker/tidyverse:4.3.1
# Slim image for benchmarking and developing R workflows.
# Image stack:
# - rocker/tidyverse:4.3.1
# - rocker/rstudio:4.3.1
# - rocker/r-ver:4.3.1
FROM rocker/tidyverse:4.3.1

# Copy repository (might have to adjust for diff. obj)
COPY . /home/src/pgc-perf-opt

# Bind a volume for development (might have to adjust for diff. obj)
# Notes: Container has a /home/rstudio directory.


### Python
# Comes with Python 3.10.6 with base packages via python3.
# Install the rest of the dependencies manually via pip.
# In particular, install plotly.

### R and RStudio
# Comes with R 4.3.1 and RStudio Server, bundled with tidyverse.
# RStudio Server is available via exposed local port, configure creds.

# Building (from project root):
# docker build -t pgc-perf-opt/cpu . -f ./docker/cpu.Dockerfile

# Usage:
# https://rocker-project.org/images/versioned/rstudio.html
# Password: docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 rocker/rstudio
# Root: docker run --rm -ti -e ROOT=true -p 8787:8787 rocker/rstudio
# Disable Authentication: docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 rocker/rstudio
# User ID and Group ID: docker run --rm -ti -e USERID=1001 -e GROUPID=1001 -p 8787:8787 rocker/rstudio
# To not expose: -p 127.0.0.1:8787:8787

# Ref: https://rocker-project.org/images/
