# Use rocker/tidyverse:4.3.1 as base image.
# Layers:
# - rocker/rstudio:4.3.1
# - rocker/r-ver:4.3.1

FROM rocker/tidyverse:4.3.1

# Add CUDA support to rocker/tidyverse:4.3.1
# Layers:
# - rocker/ml:4.3.1
# - rocker/cuda:4.3.1
# - rocker/r-ver:4.3.1
FROM rocker/ml:4.3.1

RUN apt update && apt install -y 

COPY . /app
WORKDIR /app/code/R
CMD RScript kmer-analysis.R

# Usage:
# https://rocker-project.org/images/versioned/rstudio.html
# https://rocker-project.org/images/versioned/cuda.html

# Ref: https://rocker-project.org/images/
