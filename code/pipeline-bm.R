# File: pipeline-bm.R
# Benchmark for Research Objective 3 (ro3) vulnerability mitigations.
# Results are stored in benchmarks/ro3.

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

# SET PARAMETERS ###########################################
# pipeline.R general parameters
# stamp <- [get_time():str|NULL]
# if stamp = "", then generated files won't be timestamped
seed <- 1234
stamp <- get_time()
write_fastacsv <- TRUE
kmer_list <- c(3, 5, 7)
# strat_size: no. of samples per stratum. Current nrow(data) = 24671.
# Also consider using sample_frac for proportionate allocation.
# Note that valid strat_size will only be those with corresponding
# files in `data/interm` and `data/kmers`
strat_size <- 100

# preprocess.R::get_sample() parameters
gisaid_data_path <- "data/GISAID"
gisaid_extract_path <- "data/GISAID/datasets"
country_exposure <- "Philippines"

# preprocess.R::auxiliary parameters
interm_write_path <- "data/interm"
compile_write_path <- "data/overview"
treemaps_write_path <- "data/overview/treemaps"
heatmaps_write_path <- "data/overview/heatmaps"

# dim-reduce.R::dim_reduce() parameters
kmers_data_path <- "data/kmers"
dimreduce_write_path <- "results/dim-reduce/R"
tsne_perplexity <- 40
tsne_max_iter <- 1000
tsne_initial_dims <- 50
umap_n_neighbors <- 15
umap_metric <- "euclidean"
umap_min_dist <- 0.1
color <- "variant"
shape <- "sex"

# dim-reduce.R::dim_reduce() filtering parameters - OPTIONAL
factor1 <- "variant"
values1 <- c("Omicron", "Omicron Sub")
factor2 <- "year"
values2 <- c("2023")

# clustering-x.R::dendogram_create_x() parameters
agnes_write_path <- "results/dendrogram"

# Benchmark parameters
bm_times <- 1L   # how many times should routine be evaluated
bm_write_path <- "benchmarks/ro3"
OS <- pacman::p_detectOS()
if (OS == "Windows") {
  processor <- stringr::str_extract(system("systeminfo",
                                           intern = TRUE)[17], "AMD|Linux")
} else {
  processor <- stringr::str_extract(system("lscpu | grep 'Model name'",
                                           intern = TRUE), "AMD|Linux")
}
# valid values: ["ALL"|"SOME" (Linux only)|"NONE"]
# Mitigations are automatically detected in Linux by write_to_log.
# Don't forget to set this to current runtime configuration of the benchmark.
mitigations <- "ALL"

# HELPER FUNCTIONS ##########################################
# Benchmark passed operation.
# Returns list of benchmark results.
# [ expr unit min lq mean median uq max neval ]
# Note that microbenchmark always returns nanoseconds
bm_cpu <- function(op, args, use_profiling, unit,
                   times = 3L, warmup = 100000L) {
  # system.time: better for longer running code chunks
  # microbenchmark: better for fast-running code chunks
  if (use_profiling) {
    # Warmup before actual benchmarking operation
    for (j in 1:warmup) {
      A <- matrix(c(15,20,25,15,20,25,15,20,25), ncol=3, nrow=3)
      B <- matrix(c(35,26,18,30,25,17,37,28,20), ncol=3, nrow=3)
      warmer <- A %*% B
    }
    # all_times stores all the elapsed times in each system.time
    all_times <- c()
    for (j in 1:times) {
      profile <- system.time(do.call(op, args))
      all_times <- c(all_times, profile["elapsed"])
    }
    # Compute min, lq, mean, median, uq, and max
    summ <- validate::summary(all_times)
    # Perform necessary conversion of data from seconds to desired unit
    if (unit == "nanoseconds") {
      summ[-7] <- lapply(summ[-7], function(x){round(x*1e+09,3)})
    } else if (unit == "microseconds") {
      summ[-7] <- lapply(summ[-7], function(x){round(x*1e+06,3)})
    } else if (unit == "milliseconds") {
      summ[-7] <- lapply(summ[-7], function(x){round(x*1e+03,3)})
    } else if (unit == "minutes") {
      summ[-7] <- lapply(summ[-7], function(x){round(x/60,3)})
    }
    # Switch mean and median to get proper ordering
    tmp <- summ[3]
    summ[3] <- summ[4]
    summ[4] <- tmp
    names(summ) <- c("min", "lq", "mean", "median", "uq", "max")
    res <- c(summ, times)
  } else {
    # Start benchmark
    summ <- summary(microbenchmark(do.call(op, args),
                                   times = times, unit = unit,
                                   control = list(order = "inorder",
                                                  warmup = warmup)))[-1]
    # Perform necessary conversion of data from ns to desired unit
    if (unit == "seconds") {
      summ[-7] <- lapply(summ[-7], function(x){round(x*1e-09,3)})
    } else if (unit == "milliseconds") {
      summ[-7] <- lapply(summ[-7], function(x){round(x*1e-06,3)})
    } else if (unit == "microseconds") {
      summ[-7] <- lapply(summ, function(x){round(x*1e-03,3)})
    } else if (unit == "minutes") {
      summ[-7] <- lapply(summ[-7], function(x){round(x*1e-09/60,3)})
    }
    res <- summ
  }
}

# TODO!
plot_bm <- function(results) {
  
}

# Looper for get_kmers
get_kmers_all <- function(kmer_list, fasta_all, metadata_all, stamp) {
  for (k in kmer_list) {
    get_kmers(fasta_all, metadata_all, k, stamp)
  }
}

# Helpers for dim-reduce algorithms
# Looper for pca_fn
pca_fn_all <- function(draux) {
  for (i in 1:length(draux)) {
    pca_fn(draux[[i]][[2]])
  }
}

# Looper for tsne_fn
tsne_fn_all <- function(draux, D, tsne_initial_dims, tsne_perplexity,
                        tsne_max_iter, tsne_seed) {
  for (i in 1:length(draux)) {
    tsne_fn(draux[[i]][[3]], D, tsne_initial_dims,
            tsne_perplexity, tsne_max_iter,
            tsne_seed = seed)
  }
}

# Looper for umap_fn
umap_fn_all <- function(draux, D, umap_n_neighbors,
                        umap_metric, umap_min_dist, umap_seed) {
  for (i in 1:length(draux)) {
    umap_fn(draux[[i]][[3]], D, umap_n_neighbors,
            umap_metric, umap_min_dist,
            umap_seed = seed)
  }
}

# Loopers for clustering functions
dendro_var_all <- function(kmer_list, kmers, agnes_write_path) {
  for (i in 1:length(kmer_list)) {
    k <- kmer_list[i]
    dendrogram_create_variant(k, kmers[[i]], agnes_write_path)
  }
}

dendro_reg_all <- function(kmer_list, kmers, agnes_write_path) {
  for (i in 1:length(kmer_list)) {
    k <- kmer_list[i]
    dendrogram_create_region(k, kmers[[i]], agnes_write_path)
  }
}

# PARSE DATA FILES ##########################################
f_path <- sprintf("%s/fasta_all_%s.fasta", interm_write_path, strat_size)
m_path <- sprintf("%s/metadata_all_%s.csv", interm_write_path, strat_size)
message(sprintf("Reading %s... ", f_path), appendLF = FALSE)
fasta_all <- ape::read.FASTA(f_path)
message("DONE!")
message(sprintf("Reading %s... ", m_path), appendLF = FALSE)
metadata_all  <- readr::read_csv(m_path, show_col_types = FALSE)
message("DONE!")
# kmers is list of kmer dataframes
kmers <- list()
for (i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  k_path <- sprintf("%s/kmer_%d_%d.csv", kmers_data_path, k, strat_size)
  message(sprintf("Reading %s for later... ", k_path), appendLF = FALSE)
  kmers[[i]] <- utils::read.csv(k_path)
  message("DONE!")
}

# RUN BENCHMARK #############################################
NROWS <- nrow(metadata_all)
message(sprintf("Running pipeline-bm.R benchmark on %s with mitigations: %s", OS, mitigations))
message(sprintf("Number of selected samples are: %d", NROWS))

# dim-reduce-aux: Prepare for dim-reduce algorithms
# List format: {(df_k, x_k, pca_df_k), ...}
draux <- list()
for (i in 1:length(kmer_list)) {
  pre_reduce_res <- pre_reduce(dimreduce_write_path,
                               kmers[[i]], kmer_list[i],
                               factor1, values1, factor2, values2)
  df <- pre_reduce_res$df                # df is the original dataset
  x <- pre_reduce_res$x                  # x is the scaled data
  
  # Run iterations of pca for each k for tsne to use
  pca_df <- pca_fn(x)
  
  # Store for later
  draux[[i]] <- list(df, x, pca_df$x)
}

# Benchmark Notes:
# get_sample: to start from extraction, delete data/GISAID/datasets/
# generate_interm always happens for this benchmark
# TODO: Add descriptions for my functions to the plot results

# Initialize list of operations to benchmark and their arguments
# Format: {operation:function, args:list}
ops <- list(
            # preprocess.R
            list(get_sample,
                 list(gisaid_data_path, gisaid_extract_path,
                      seed, strat_size, country_exposure)),
            list(sanitize_sample,
                 list(metadata_all)),
            list(generate_interm,
                 list(fasta_all, metadata_all, interm_write_path, stamp)),
            list(compile_overview,
                 list(metadata_all, compile_write_path, stamp)),
            list(make_treemaps,
                 list(metadata_all, treemaps_write_path, stamp)),
            # kmer-analysis.R
            list(get_kmers,
                 list(fasta_all, metadata_all, 3, stamp)),
            list(get_kmers,
                 list(fasta_all, metadata_all, 5, stamp)),
            list(get_kmers,
                 list(fasta_all, metadata_all, 7, stamp)),
            list(get_kmers_all,
                 list(kmer_list, fasta_all, metadata_all, stamp)),
            # dim-reduce.R
            ## PCA
            list(pca_fn,   # k = 3
                 list(draux[[1]][[2]])),
            list(pca_fn,   # k = 5
                 list(draux[[2]][[2]])),
            list(pca_fn,   # k = 7
                 list(draux[[3]][[2]])),
            list(pca_fn_all,
                 list(draux)),
            ## TSNE 2D
            list(tsne_fn,                 # k = 3
                 list(draux[[1]][[3]], 2, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            list(tsne_fn,                  # k = 5
                 list(draux[[2]][[3]], 2, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            list(tsne_fn,                   # k = 7
                 list(draux[[3]][[3]], 2, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            list(tsne_fn_all,
                 list(draux, 2, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            ## TSNE 3D
            list(tsne_fn,                   # k = 3
                 list(draux[[1]][[3]], 3, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            list(tsne_fn,                   # k = 5
                 list(draux[[2]][[3]], 3, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            list(tsne_fn,                   # k = 7
                 list(draux[[3]][[3]], 3, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            list(tsne_fn_all,
                 list(draux, 3, tsne_initial_dims,
                      tsne_perplexity, tsne_max_iter,
                      tsne_seed = seed)),
            ## UMAP 2D
            list(umap_fn,                 # k = 3
                 list(draux[[1]][[3]], 2, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            list(umap_fn,                 # k = 5
                 list(draux[[2]][[3]], 2, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            list(umap_fn,                 # k = 7
                 list(draux[[3]][[3]], 2, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            list(umap_fn_all,
                 list(draux, 2, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            ## UMAP 3D
            list(umap_fn,                 # k = 3
                 list(draux[[1]][[3]], 3, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            list(umap_fn,                 # k = 5
                 list(draux[[2]][[3]], 3, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            list(umap_fn,                 # k = 7
                 list(draux[[3]][[3]], 3, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            list(umap_fn_all,
                 list(draux, 3, umap_n_neighbors,
                      umap_metric, umap_min_dist,
                      umap_seed = seed)),
            # Clustering AGNES
            list(dendrogram_create_variant,
                 list(3, kmers[[1]], agnes_write_path)),
            list(dendrogram_create_variant,
                 list(5, kmers[[2]], agnes_write_path)),
            list(dendrogram_create_variant,
                 list(7, kmers[[3]], agnes_write_path)),
            list(dendro_var_all,
                 list(kmer_list, kmers, agnes_write_path)),
            list(dendrogram_create_region,
                 list(3, kmers[[1]], agnes_write_path)),
            list(dendrogram_create_region,
                 list(5, kmers[[2]], agnes_write_path)),
            list(dendrogram_create_region,
                 list(7, kmers[[3]], agnes_write_path)),
            list(dendro_reg_all,
                 list(kmer_list, kmers, agnes_write_path))
            )

# Also initialize names of the functions (can't get it programmatically)
names <- list(
              "get_sample",
              "sanitize_sample",
              "generate_interm",
              "compile_overview",
              "make_treemaps",
              "get_kmers_3",
              "get_kmers_5",
              "get_kmers_7",
              "get_kmers_all",
              "pca_3",
              "pca_5",
              "pca_7",
              "pca_all",
              "tsne_2d_3",
              "tsne_2d_5",
              "tsne_2d_7",
              "tsne_2d_all",
              "tsne_3d_3",
              "tsne_3d_5",
              "tsne_3d_7",
              "tsne_3d_all",
              "umap_2d_3",
              "umap_2d_5",
              "umap_2d_7",
              "umap_2d_all",
              "umap_3d_3",
              "umap_3d_5",
              "umap_3d_7",
              "umap_3d_all",
              "agnes_var_3",
              "agnes_var_5",
              "agnes_var_7",
              "agnes_var_all",
              "agnes_reg_3",
              "agnes_reg_5",
              "agnes_reg_7",
              "agnes_reg_all"
              )

# Addon: profiling boolean list and units char list for finer control
control <- list(
                # preprocess.R
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # kmer-analysis.R
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # PCA
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # TSNE 2D
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # TSNE 3D
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # UMAP 2D
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # UMAP 3D
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # Clustering AGNES variant
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                # Clustering AGNES region
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds"),
                list(TRUE, "seconds")
                )

# Initialize results dataframe
cols <-  c("op", "profiler", "unit",
           "min", "lq", "mean", "median", "uq", "max", "neval",
           "strat_size", "processor", "mitigations", "timestamp")
results <- data.frame(matrix(nrow = 0, ncol = length(cols)))
colnames(results) <- cols
results[, 1:3] <- sapply(results[, 1:3], as.character)
results[, 3:9] <- sapply(results[, 3:9], as.numeric)
results[, 10:13] <- sapply(results[, 10:13], as.character)
results[14] <- sapply(results[14], as.numeric)

# BENCHMARKER
# Get results and append to dataframe (actual benchmarking part)
res <- list()
for (i in 1:length(ops)) {
  # Get arguments
  op <- ops[[i]][[1]]
  opname <- names[[i]]
  args <- ops[[i]][[2]]
  use_profiling <- control[[i]][[1]]
  unit <- control[[i]][[2]]
  
  profiler <- character(0)
  if (use_profiling)
    profiler <- "system.time"
  else
    profiler <- "microbenchmark"
  
  # BENCHMARKER
  res <- bm_cpu(op, args, use_profiling, unit, times = bm_times)
  
  # Append results to results
  results[nrow(results)+1, 1:3] <- c(opname, profiler, unit)
  results[nrow(results), 4:9] <- res[1:6]
  results[nrow(results), 10:11] <- c(res[[7]], strat_size)
  results[nrow(results), 12:14] <- c(processor, mitigations, timestamp)
}

# Write (append) results to accumulator file in bm_write_path
filepath <- paste(bm_write_path, "ro3.csv", sep = "/")
if (!file.exists(filepath)) {
  readr::write_csv(results, filepath)
} else {
  readr::write_csv(results, filepath, append = TRUE)
}

# Write hardware specs and parameters used to log.txt

param_string <- paste(c("---------PARAMETERS---------",
  paste0("MITIGATIONS STATUS:\t", mitigations),
  paste0("seed:\t\t\t", seed),
  paste0("kmer_list:\t\t", paste(kmer_list, collapse = ", ")),
  paste0("strat_size:\t\t", strat_size),
  paste0("number_samples:\t\t", NROWS),
  paste0("country_exposure:\t", country_exposure),
  paste0("tsne_perplexity:\t", tsne_perplexity),
  paste0("tsne_max_iter:\t\t", tsne_max_iter),
  paste0("umap_n_neighbors:\t", umap_n_neighbors),
  paste0("umap_metric:\t\t", umap_metric),
  paste0("umap_min_dist:\t\t", umap_min_dist),
  paste0("color:\t\t\t", color),
  paste0("shape:\t\t\t", shape, "\n")),
  collapse = "\n")

message("Writing logs... ", appendLF = FALSE)
write_to_log(bm_write_path, "bm_log.txt", param_string, stamp)
message("Writing logs... DONE!")

message("All operations completed successfully!")

# CLEAN UP #################################################

# Clear environment
# rm(list = ls()) 

# Clear packages (unloading them before another adds another compat check)
p_unload(all)  # Remove all add-ons

# Clear plots but only if there IS a plot
while (!is.null(dev.list())) dev.off()

# Clear console
# cat("\014")  # ctrl+L

# Clear mind :)