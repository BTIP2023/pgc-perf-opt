# Use r-base as base image.
# This is the official Dockerhub image for R.
# We use this instead of Rocker because Rocker
# uses RStudio Server. I'm not sure, still testing this out.

FROM r-base
COPY . /app
WORKDIR /app/code/R
CMD RScript/kmer-analysis.R