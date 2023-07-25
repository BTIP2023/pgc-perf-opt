# TO DO: Add loggers, fix bugs, remove unecessary comments

# File: flexiblas-bench.R
# Contains the code for benchmarking BLAS/LAPACK endpoints against the pipeline
# and some of its components.

# INSTALL AND LOAD PACKAGES ################################
# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")
library(pacman)

if (pacman::p_detectOS() == "Linux" && !pacman::p_exists(xml2, local = TRUE)) {
  install.packages("xml2", dependencies = TRUE, INSTALL_opts = c("--no-lock"))
  pacman::p_load(xml2)
}

# Load required packages
pacman::p_load(plyr, dplyr, GGally, ggplot2, ggthemes, ggvis,
               httr, lubridate, plotly, psych,
               rio, markdown, rmarkdown, shiny,
               stringr, tidyr, tidyverse,
               ape, kmer, readr, validate, gsubfn, seqinr,
               umap, htmlwidgets, factoextra, scales,
               Rtsne, tsne, RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark, flexiblas, car, FSA)
# RATE LIMIT HIT <Uncomment Later>
# install_github("vqv/ggbiplot", upgrade = FALSE, quiet = TRUE)
# pacman::p_load(ggbiplot)


# LOAD SOURCE #############################################
source("code/R/dim-reduce.R")
source("code/R/helper.R")

# FUNCTIONS ################################
# Benchmarking function
benchmark_backends <- function(operation, args_list, backends, times, unit, 
                               use_profiling=FALSE) {
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
      # elapsed_time <- mean(all_times*1e3) # Convert unit of time to ms
    }
    else {
      # Start benchmark
      bm_result <- microbenchmark(do.call(operation, args_list),
                                  times = times, unit = unit)
      # Get mean time
      # elapsed_time <- summary(bm_result)$mean
      all_times <- c(bm_result$time)
    }
    # Store the results
    results[[backend]] <- all_times
  }
  
  return(results)
}

# Benchmark plotting function
plot_results <- function(method_benchmark, method) {
  # Combine the results into a data frame
  results_df <- data.frame(Backend = as.character(selected_backends), 
                           Time = unlist(method_benchmark))
  
  # Create the bar plot
  p <- ggplot(results_df, aes(y = Backend, x = Time, fill = Backend)) +
    geom_bar(stat = "identity", width = 0.3) +
    geom_text(aes(label = round(Time, 3), hjust = 1.25)) +
    labs(title = paste0("Benchmark: ", toupper(method)),
         y = "Backend",
         x = "Execution Time (milliseconds)",
         subtitle = "Note: Execution Time in ms; lower is better") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
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
  p <- ggplotly(p, width = 1500, height = 700, tooltip = c("x"))
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
  saveWidget2(p, path = html_file)
}

# Functions for Statistical Tests
should_anova <- function(data, alpha_value) {
  # Check Assumption #1: Normality using Shapiro-Wilk Test
  shapiro_res <- shapiro.test(data$time)
  if(shapiro_res$p.value > alpha_value) {
    is_normal <- TRUE
  }
  else {
    is_normal <- FALSE
  }
  
  # Check Assumption #2: Equal Variance 
  # If normal, use Bartlett's Test. Otherwise, use Levene's Test
  if(is_normal) {
    bartlett_res <- bartlett.test(time ~ backend, data=data)
    if(bartlett_res$p.value > alpha_value) {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }
  else {
    # Compare: [answer: levene.test is defunct]
    # wrong args!!
    # levene_res <- levene.test(time ~ backend, data = data,
    #                           method = "correction.factor")
    
    # library(car)
    levene_res <- leveneTest(time ~ backend, data = data)
    p_val <- levene_res[[as.name("Pr(>F)")]][1]
    # print(levene_res)
    # print(p_val)
    if(p_val > alpha_value) {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
    
  # Note that Assumption #3: Independence is assumed.
  }
}

check_stat_diff <- function(){
  # write this tomorrow
  warning("You ran check_stat_diff but haven't written it yet.")
}

# shapiro_test <- function(method_bm, selected_backends) {
#   shapiro_test_results <- list()
#   for (backend in selected_backends) {
#     shapiro_test_results[[backend]] <- shapiro.test(method_bm[[backend]])
#   }
#   
#   # Check if each backend's data is normally distributed
#   for (backend in backends) {
#     p_value <- shapiro_test_results[[backend]]$p.value
#     if (p_value <= 0.05) {
#       return(FALSE)
#     }
#   }
#   return(TRUE)
# }

# levene_test <- function(method_bm, selected_backends, bm_times) {
#   data <- data.frame(backend = rep(selected_backends, times = bm_times), 
#                      time = unlist(method_bm))
#   leveneTest(time ~ backend, data = data)
#   
# }



# SET PARAMETERS ###########################################
# dim-reduce.R::dim_reduce() parameters
seed <- 1234
k_vals <- c(3)
data_path_kmers <- "data/kmers"
results_path_dimreduce <- "results/dim-reduce/R"
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

# Benchmarking parameters
bm_times <- 3 
bm_unit <- "milliseconds"
alpha_value <- 0.05

# Backends to compare
selected_backends <- c("NETLIB", "ATLAS", "OPENBLASSERIAL",
                       "MKLSERIAL", "BLISSERIAL")

for (k in k_vals) {
  pre_reduce_res <- pre_reduce(results_path_dimreduce,
                               data_path_kmers, k, factor1, 
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
  
  # TO DO: Create a function: create_df -> generalized method_bm
  data <- data.frame(backend = rep(selected_backends, times = bm_times), 
                     time = unlist(pca_bm))
  
  use_anova <- should_anova(data, alpha_value)
  
  # is_diff <- check_stat_diff(use_anova)
  
  # Decide if succeeding lines should be in a function 
  if(!use_anova) {
    # [DONE, not yet checked]
    # TO DO: Perform ANOVA -> Perform Tukey's Test (TukeyHSD)
    # Compute the analysis of variance
    aov_res <- aov(time ~ backend, data = data)
    print("Done ANOVA")
    # print(summary(aov_res))
    # p_val <- try[[as.name("Pr(>F)")]] # TO FIX: bug here in extraction of p.value
    p_val <- summary(aov_res)[[1]]$Pr[1]
    
    # print(p_val)
    if(p_val <= alpha_value) {
      # Perform TukeyHSD
      tukey_res <- TukeyHSD(aov_res, conf.level=1-alpha_value)
      print("Done Tukey")
    }
  }
  else {
    # [DONE, not yet checked]
    # TO DO: Perform Kruskal-Wallis Test -> Perform Dunn's Test
    kruskal_res <- kruskal.test(time ~ backend, data = data)
    # Summary of the analysis
    print(summary(kruskal_res)) # remove in final (?)
    if(kruskal_res$p.value <= alpha_value) {
      # Perform Dunn's Test
      # x.x NOT FOUND x.x
      # dunn_res <- dunn.test(time ~ backend, data = data, method="bonferroni")
      # print(dunn_res$p.value)
      # Compare: (print the two dunn's test results)
      dunn_res <- dunnTest(time ~ backend, data = data, method="bonferroni")
      print(dunn_res)
      print(dunn_res$p.value)
    }
    print("Done Kruskal")
  }
  
  # # Plot PCA benchmark results
  # plot_results(pca_bm, "pca")

  # # Benchmark backends on tsne_fn (2-dimensions)
  # tsne_bm <- benchmark_backends(tsne_fn,
  #                               list(pca_df$x, 2, tsne_initial_dims,
  #                                    tsne_perplexity, tsne_max_iter,
  #                                    tsne_seed = seed),
  #                               selected_backends,
  #                               bm_times,
  #                               bm_unit)
  # 
  # # Plot t-SNE benchmark results (2D)
  # plot_results(tsne_bm, "tsne-2d")
  # 
  # # Benchmark backends on tsne_fn (3-dimensions)
  # tsne_bm <- benchmark_backends(tsne_fn,
  #                               list(pca_df$x, 3, tsne_initial_dims,
  #                                    tsne_perplexity, tsne_max_iter,
  #                                    tsne_seed = seed),
  #                               selected_backends,
  #                               bm_times,
  #                               bm_unit)
  # 
  # # Plot t-SNE benchmark results (3D)
  # plot_results(tsne_bm, "tsne-3d")
  # 
  # # Benchmark backends on umap_fn (2-dimensions)
  # umap_bm <- benchmark_backends(umap_fn,
  #                               list(x, 2, umap_n_neighbors,
  #                                    umap_metric, umap_min_dist,
  #                                    umap_seed = seed),
  #                               selected_backends,
  #                               bm_times,
  #                               bm_unit)
  # 
  # # Plot UMAP benchmark results (2D)
  # plot_results(umap_bm, "umap-2d")
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
  # # Plot UMAP benchmark results (3D)
  # plot_results(umap_bm, "umap-3d")
  # 
  # # Benchmark backends on dim_reduce
  # dimred_bm <- benchmark_backends(dim_reduce,
  #                                 list(k, data_path_kmers,
  #                                      results_path_dimreduce,
  #                                      tsne_seed = seed, tsne_perplexity,
  #                                      tsne_max_iter, tsne_initial_dims,
  #                                      umap_seed = seed, umap_n_neighbors,
  #                                      umap_metric, umap_min_dist,
  #                                      col_name = target_col),
  #                                 selected_backends,
  #                                 bm_times,
  #                                 bm_unit,
  #                                 use_profiling = TRUE)
  # 
  # # Plot dim_reduce benchmark results
  # plot_results(dimred_bm, "dim-red")
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
  # # Plot pipeline benchmark results
  # plot_results(pipeline_bm, "pipeline")
  # 
  # print(paste0("Benchmarking ", k, "-mer Done :>"))
}