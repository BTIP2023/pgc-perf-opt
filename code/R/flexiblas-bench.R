# File: flexiblas-bench.R
# Contains the code for benchmarking BLAS/LAPACK endpoints against the pipeline
# and some of its components.

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
shape <- "year"
include_plots <- TRUE

# dim-reduce.R::dim_reduce() filtering parameters - OPTIONAL
# factor1 <- "variant"
# values1 <- c("Omicron", "Omicron Sub")
# factor2 <- "year"
# values2 <- c("2023")

# clustering-x.R::dendogram_create_x() parameters
agnes_write_path <- "results/dendrogram"

# FUNCTIONS ###############################################
# Benchmarking function
benchmark_backends <- function(operation, args_list, backends, times, unit, 
                               use_profiling=FALSE) {
  message("Performing Benchmark...")
  results <- list()
  for (backend in backends) {
    # Set backend before performing benchmark
    flexiblas_load_backend(backend)
    
    # Run the benchmark for the current backend
    if (use_profiling) {
      # Simulate 2 rounds of warm-up process of microbenchmark
      for (i in 1:2) {
        do.call(operation, args_list)
      }
      # Start benchmark
      all_times <- c()
      for (i in 1:times) {
        bm_result <- system.time(do.call(operation, args_list))
        all_times <- c(all_times, bm_result["elapsed"])
      }
      all_times <- all_times*1e3            # Convert time unit: s to ms
      # elapsed_time <- mean(all_times*1e3) # Convert time unit: s to ms
    }
    else {
      # Start benchmark
      bm_result <- microbenchmark(do.call(operation, args_list),
                                  times = times, unit = unit)
      # Get mean time
      # elapsed_time <- summary(bm_result)$mean
      all_times <- c(bm_result$time)*(1/1e6) # Convert time unit: ns to ms
    }
    # Store the results (both individual times in all_times and mean time)
    results[[backend]] <- all_times
    # results[[backend]]$overall <- elapsed_time
  }
  
  return(results)
}

# Benchmark plotting function
plot_results <- function(method_res, method) {
  message("Plotting Benchmark Results...")
  # Combine the results into a data frame
  # res_df <- data.frame(Backend = as.character(selected_backends), 
  #                          Time = unlist(method_res))
  
  res_df <- method_res$data
  
  # # Create the bar plot
  # p <- ggplot(res_df, aes(y = Backend, x = Time, fill = Backend)) +
  #   geom_bar(stat = "identity", width = 0.3) +
  #   geom_text(aes(label = round(Time, 3), hjust = 1.25)) +
  #   labs(title = paste0("Benchmark: ", toupper(method)),
  #        y = "Backend",
  #        x = "Execution Time (milliseconds)",
  #        subtitle = "Note: Execution Time in ms; lower is better") +
  #   theme_minimal() +
  #   theme(plot.title = element_text(hjust = 0.5))
  
  # Create the box plot
  p <- ggplot(res_df, aes(y = backend, x = time, fill = backend)) +
    geom_boxplot() +
    # geom_text(aes(label = round(time, 3), hjust = 1.25)) +
    labs(title = paste0("Benchmark: ", toupper(method)),
         y = "Backend",
         x = "Execution Time (milliseconds)",
         subtitle = "Note: Execution Time in ms; lower is better") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ggplot(xx1, aes(x = variable, y = value)) + 
  #   geom_boxplot()
  
  
  # Create directory for storing plots
  results_path <- "results/benchmark/R/FlexiBLAS"
  
  # Check if the directory already exists
  if (!dir.exists(results_path)) {
    # Create the directory if it doesn't exist
    dir.create(results_path, recursive = TRUE)
  }
  
  # Save the bar plot
  ggsave(filename=paste0(method, "-", k, "-benchmark.png"), 
         plot = p, path = results_path, 
         device = "png", width = 12, height = 6,
         dpi = 300, bg = "white")
  
  # Convert ggplot object to ggplotly
  p <- ggplotly(p, width = 1500, height = 700)#, tooltip = c("x"))
  p <- p %>% layout(title = list(text = paste0("Benchmark: ", toupper(method)),
                                 x = 0.5,
                                 xref = "paper"))
  p <- p %>% layout(annotations = list(
                      text = "Note: Execution Time in ms; lower is better",
                      x = 0, y = 5.6,
                      showarrow = FALSE,
                      font = list(size = 12)
                    ))
  
  # Save as HTML
  html_file <- paste0(results_path, "/", method, "-", k, "-benchmark.html")
  saveWidget2(p, file = html_file)
}

# Functions for Statistical Tests
# Function that determines whether to use One-way ANOVA or Kruskal-Wallis
should_anova <- function(data, alpha_value) {
  # Check Assumption #1: Normality using Shapiro-Wilk Test
  message("Performing Shapiro-Wilk Test...")
  shapiro_res <- shapiro.test(data$time)
  p_val <- shapiro_res$p.value
  message(paste0("Shapiro-Wilk Test is DONE w/ p-value ", p_val))
  
  if(p_val > alpha_value) {
    # is_normal <- TRUE
    message("The data is normally distributed. We proceed with Bartlett's Test to test homogeneity of variance.")
    # Check Assumption #2: Equal Variance
    message("Performing Bartlett's Test...")
    p_val <- bartlett_res$p.value
    message(paste0("Bartlett's Test is DONE w/ p-value ", p_val))

    if(p_val > alpha_value) {
      message("There is homogeneity of variance in the data. We can use ANOVA.")
      return(TRUE)
    }
    else {
      message("There is no homogeneity of variance in the data. We must use Kruskal-Wallis Test instead of ANOVA.")
      return(FALSE)
    }
  }
  else {
    # is_normal <- FALSE
    message("The data is not normally distributed. We must use Kruskal-Wallis Test instead of ANOVA.")
    return(FALSE)
  }
  
  # Note that Assumption #3: Independence is assumed.
}

# Function that determines whether there is a statistical difference in the data
check_stat_diff <- function(data, use_anova, alpha_val){
  if(use_anova) {
    # Perform ANOVA
    message("Performing ANOVA...")
    aov_res <- aov(time ~ backend, data = data)
    p_val <- summary(aov_res)[[1]]$Pr[1]
    message(paste0("ANOVA is DONE w/ p-value ", p_val))
    
    if(p_val <= alpha_value) {
      message("At least one of the values is statistically different from the others.")
      # Perform TukeyHSD
      message("Performing Tukey's Test...")
      tukey_res <- TukeyHSD(aov_res, conf.level=1-alpha_value)
      message("Tukey's Test is DONE. The summary of results is as follows:")
      print(tukey_res)
      adj_p_vals <- tukey_res$backend[,3]
      # Filter the rows where the adjusted p-value is less than alpha
      significant_pairs <- list()
      for (i in 1:length(adj_p_vals)) {
        val <- as.numeric(adj_p_vals[i])
        if(as.numeric(adj_p_vals[i]) <= alpha_value) {
          significant_pairs <- append(significant_pairs, names(adj_p_vals[i]))
        }
      }
      significant_pairs <- unlist(significant_pairs)
      return(list(is_diff = TRUE, significant_pairs = significant_pairs))
    }
    else {
      message("The values are not statistically different.")
      return(list(is_diff = FALSE, diff_vals = NULL))
    }
  }
  else {
    message("Performing Kruskal-Wallis Test...")
    kruskal_res <- kruskal.test(time ~ backend, data = data)
    p_val <- kruskal_res$p.value 
    message(paste0("Kruskal-Wallis Test is DONE w/ p-value ", p_val))
    
    if(p_val <= alpha_value) {
      message("At least one of the values is statistically different from the others.")
      message("Performing Dunn's Test...")
      # Perform Dunn's Test
      dunn_res <- dunnTest(time ~ backend, data = data, method="bonferroni")
      message("Dunn's Test is DONE. The summary of results is as follows:")
      print(dunn_res)
      # Filter the rows where the adjusted p-value is less than alpha
      significant_pairs <- dunn_res$res$Comparison[which(dunn_res$res$P.adj <= alpha_value)]
      return(list(is_diff = TRUE, significant_pairs = significant_pairs))
    }
    else {
      message("The values are not statistically different.")
      return(list(is_diff = FALSE, diff_vals = NULL))
    }
  }
}

# Main function that assesses the execution times in the benchmark
assess_bm <- function(method_bm, selected_backends, bm_times) {
  data <- data.frame(backend = rep(selected_backends, times = bm_times), 
                     time = unlist(method_bm))
  use_anova <- should_anova(data, alpha_value)
  is_diff <- check_stat_diff(data, use_anova, alpha_value)
  return(list(data = data, is_diff = is_diff))
}

# SET PARAMETERS ###########################################
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
shape <- "year"
include_plots <- FALSE

# dim-reduce.R::dim_reduce() filtering parameters - OPTIONAL
factor1 <- "variant"
values1 <- c("Omicron", "Omicron Sub")
factor2 <- "year"
values2 <- c("2023")

# Benchmarking parameters
bm_times <- 5 
bm_unit <- "milliseconds"
alpha_value <- 0.05

# Backends to compare
selected_backends <- c("NETLIB", "ATLAS", "OPENBLASSERIAL",
                       "MKLSERIAL", "BLISSERIAL")

for(i in 1:length(kmer_list)) {
  k <- kmer_list[i]
  pre_reduce_res <- pre_reduce(results_path_dimreduce,
                               kmers[[i]], k, factor1, 
                               values1, factor2, values2)
  df <- pre_reduce_res$df                # df is the original dataset
  x <- pre_reduce_res$x                  # x is the scaled data
  # Run an iteration of PCA for t-SNE use
  pca_df <- pca_fn(x)
  
  # Benchmark backends on pca_fn
  pca_bm <- benchmark_backends(pca_fn, list(x),
                               selected_backends,
                               bm_times,
                               bm_unit,
                               use_profiling = FALSE)
  
  # Perform Statistical Analysis on PCA results
  pca_res <- assess_bm(pca_bm, selected_backends, bm_times)
  
  # Plot PCA benchmark results
  plot_results(pca_res, "pca")

  # # Benchmark backends on tsne_fn (3-dimensions)
  # tsne_bm <- benchmark_backends(tsne_fn,
  #                               list(pca_df$x, 3, tsne_initial_dims,
  #                                    tsne_perplexity, tsne_max_iter,
  #                                    tsne_seed = seed),
  #                               selected_backends,
  #                               bm_times,
  #                               bm_unit)
  # 
  # # Perform Statistical Analysis on t-SNE results
  # tsne_res <- assess_bm(tsne_bm, selected_backends, bm_times)
  # 
  # # Plot t-SNE benchmark results
  # plot_results(tsne_bm, "tsne")
  # 
  # # Benchmark backends on umap_fn (3-dimensions)
  # umap_bm <- benchmark_backends(umap_fn,
  #                               list(x, 3, umap_n_neighbors,
  #                                    umap_metric, umap_min_dist,
  #                                    umap_seed = seed),
  #                               selected_backends,
  #                               bm_times,
  #                               bm_unit)
  # 
  # # Perform Statistical Analysis on UMAP results
  # umap_res <- assess_bm(umap_bm, selected_backends, bm_times)
  # 
  # # Plot UMAP benchmark results
  # plot_results(umap_bm, "umap")
  # 
  # # Benchmark backends on dim_reduce
  # dimred_bm <- benchmark_backends(dim_reduce,
  #                                 list(k, data_path_kmers,
  #                                      kmers,
  #                                      tsne_seed = seed, tsne_perplexity,
  #                                      tsne_max_iter, tsne_initial_dims,
  #                                      umap_seed = seed, umap_n_neighbors,
  #                                      umap_metric, umap_min_dist,
  #                                      color = color, shape = shape,
  #                                      filter1_factor = factor1,
  #                                      filter1_values = values1,
  #                                      filter2_factor = factor2,
  #                                      filter2_values = values2,
  #                                      include_plots),
  #                                 selected_backends,
  #                                 bm_times,
  #                                 bm_unit,
  #                                 use_profiling = TRUE)
  # 
  # # Perform Statistical Analysis on Dimensionality Reduction results
  # dimred_res <- assess_bm(dimred_bm, selected_backends, bm_times)
  # 
  # # Plot dim_reduce benchmark results
  # plot_results(dimred_bm, "dimred")
  # 
  # # Benchmark backends on entire pipeline
  # pipeline_file <- "code/pipeline-classic.R"
  # pipeline_bm <- benchmark_backends(source,
  #                                   list(pipeline_file),
  #                                   selected_backends,
  #                                   bm_times,
  #                                   bm_unit,
  #                                   use_profiling = TRUE)
  # 
  # # Perform Statistical Analysis on Pipeline results
  # pipeline_res <- assess_bm(pipeline_bm, selected_backends, bm_times)
  # 
  # # Plot pipeline benchmark results
  # plot_results(pipeline_bm, "pipeline")

  print(paste0("Benchmarking ", k, "-mer Done :>"))
}