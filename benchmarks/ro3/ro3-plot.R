# File: ro3-plot.R
# Plot results summary from ro3.csv or ro3-fs.csv,
# whichever is more complete.

# IDEA: Run pipeline to see which needs to be profiled and which
# needs to be microbenchmarked, then set those in a list and make a tuple
# that indicates whether the result has been profiled or microbenchmarked.
# Combine results with the same benchmark approach in the same chart.
# Will still make individual charts.

# INSTALL AND LOAD PACKAGES ################################
options(repos = "https://cloud.r-project.org/")

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")
library(pacman)

### ATTN: IN LINUX SYSTEMS, CONSULT README FOR ADDITIONAL PREREQUISITES
### BEFORE RUNNING ANY SCRIPT. ISSUE: tidyverse installation.
### This is a non-issue if code is ran within cpu.Dockerfile.
### This cannot be scripted because this requires sudo priveleges.

# Install xml2 in advance to prep for tidyverse installation in Linux.
# Note that in Windows RStudio, this is installed by default.
# If you're getting xml2 errors on Windows, you broke something lol.
if (pacman::p_detectOS() == "Linux" && !pacman::p_exists(xml2, local = TRUE)) {
  install.packages("xml2", dependencies = TRUE, INSTALL_opts = c("--no-lock"))
  pacman::p_load(xml2)
}

# Use pacman to load add-on packages as desired.
# TODO: Remove redundancies in dependencies. E.g., dplyr and ggplot2
# are already dependencies of tidyverse.
pacman::p_load(plyr, dplyr, GGally, ggplot2, ggthemes, ggvis,
               httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse,
               ape, kmer, readr, validate, gsubfn, seqinr,
               umap, htmlwidgets, factoextra, scales,
               Rtsne, tsne, RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark, data.table,
               highcharter, glue)
if (!require(ggbiplot))
  install_github("vqv/ggbiplot", upgrade = FALSE, quiet = TRUE)
pacman::p_load(ggbiplot)

# validate used for %vin% operator 
# gsubfn used to destructure more than one return value
# devtools supports install_github for installing ggbiplot
# Note: Divisive k-means clustering available via kmer::cluster

# LOAD SOURCES #############################################
source("code/R/helper.R")
source("code/R/preprocess.R")
source("code/R/kmer-analysis.R")
source("code/R/dim-reduce.R")
source("code/R/clustering.R")


# ggplot(data = df, mapping = aes(x = operation, y = time)) +
#   geom_boxplot() +
#   grey_theme

# LOAD THEN WORK WITH DATA
data1 <- readr::read_csv("benchmarks/ro3/ro3-1000-amd-all.csv")
data2 <- readr::read_csv("benchmarks/ro3/ro3-1000-amd-none.csv")
data3 <- readr::read_csv("benchmarks/ro3/ro3-1000-intel-all.csv")
data4 <- readr::read_csv("benchmarks/ro3/ro3-1000-intel-none.csv")

# data1 is mitigations=all, data2 is mitigations=off
plotter <- function(data1, data2, processor) {
  opnames <- colnames(data1)
  # df is accumulator
  df <- tibble()
  for (i in 1:length(opnames)) {
    rows <- data1[i] %>%
      dplyr::rename("time" = opnames[[i]]) %>%
      mutate(operation = paste0(opnames[[i]], "_all"), mitigations = "all", .before = time)
    df <- df %>%
      bind_rows(rows) %>%
      drop_na()
  }
  for (i in 1:length(opnames)) {
    rows <- data2[i] %>%
      dplyr::rename("time" = opnames[[i]]) %>%
      mutate(operation = paste0(opnames[[i]], "_none"), mitigations = "none", .before = time)
    df <- df %>%
      bind_rows(rows) %>%
      drop_na()
  }
  
  p <- plot_ly(df, y = ~time, color = ~operation, type = "box")
  p <- p %>% layout(title = list(text = sprintf("%s Mitigations All vs None", processor)),
                    x = 0.5,
                    xref = "paper",
                    yaxis = list(type = "log", showgrid=T,
                                 ticks="outside",
                                 tickvals=c(0.001,0.01,0.1,1,10,100,1000), autorange=TRUE))
}

html1 <- plotter(data1, data2, "AMD")
html2 <- plotter(data3, data4, "Intel")

saveWidget2(html1, "benchmarks/ro3/amd.html")
saveWidget2(html2, "benchmarks/ro3/intel.html")

