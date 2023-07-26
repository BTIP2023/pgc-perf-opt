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
  print("Performing Shapiro-Wilk Test...")
  shapiro_res <- shapiro.test(data$time)
  p_val <- shapiro_res$p.value
  print(paste0("Shapiro-Wilk Test is DONE w/ p-value ", p_val))
  
  if(p_val > alpha_value) {
    # is_normal <- TRUE
    print("The data is normally distributed. We proceed with Bartlett's Test to test homogeneity of variance.")
    # Check Assumption #2: Equal Variance
    print("Performing Bartlett's Test...")
    p_val <- bartlett_res$p.value
    print(paste0("Bartlett's Test is DONE w/ p-value ", p_val))

    if(p_val > alpha_value) {
      print("There is homogeneity of variance in the data. We can use ANOVA.")
      return(TRUE)
    }
    else {
      print("There is no homogeneity of variance in the data. We must use Kruskal-Wallis Test instead of ANOVA.")
      return(FALSE)
    }
  }
  else {
    # is_normal <- FALSE
    print("The data is not normally distributed. We must use Kruskal-Wallis Test instead of ANOVA.")
    return(FALSE)
  }
  
  # # # [TO DISCARD] THIS PART IS WRONG :< No need to perform Levene's Test...
  # # Check Assumption #2: Equal Variance 
  # # If normal, use Bartlett's Test. Otherwise, use Levene's Test
  # if(is_normal) {
  #   print("Performing Bartlett's Test...")
  #   bartlett_res <- bartlett.test(time ~ backend, data=data)
  #   p_val <- bartlett_res$p.value 
  #   print(paste0("Bartlett's Test is DONE w/ p-value ", p_val))
  #   
  #   if(p_val > alpha_value) {
  #     print("There is homogeneity of variance in the data. We can use ANOVA.")
  #     return(TRUE)
  #   }
  #   else {
  #     print("There is no homogeneity of variance in the data. We cannot use ANOVA.")
  #     return(FALSE)
  #   }
  # }
  # else {
  #   print("Performing Levene's Test...")
  #   levene_res <- leveneTest(time ~ backend, data = data)
  #   p_val <- levene_res[[as.name("Pr(>F)")]][1]
  #   print(paste0("Levene's Test is DONE w/ p-value ", p_val))
  #   
  #   if(p_val > alpha_value) {
  #     print("There is homogeneity of variance in the data. We can use ANOVA.")
  #     return(TRUE)
  #   }
  #   else {
  #     print("There is no homogeneity of variance in the data. We cannot use ANOVA.")
  #     return(FALSE)
  #   }
  # }
  
  # Note that Assumption #3: Independence is assumed.
}

check_stat_diff <- function(use_anova, alpha_val){
  if(use_anova) {
    # [DONE, not yet checked]
    # TO DO: Perform ANOVA -> Perform Tukey's Test (TukeyHSD)
    # Compute the analysis of variance
    print("Performing ANOVA...")
    aov_res <- aov(time ~ backend, data = data)
    p_val <- summary(aov_res)[[1]]$Pr[1]
    print(paste0("ANOVA is DONE w/ p-value ", p_val))
    
    if(p_val <= alpha_value) {
      print("At least one of the values is statistically different from the others.")
      # Perform TukeyHSD
      print("Performing Tukey's Test...")
      tukey_res <- TukeyHSD(aov_res, conf.level=1-alpha_value)
      print("Tukey's Test is DONE. The summary of results is as follows:")
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
      print("The values are not statistically different.")
      return(list(is_diff = FALSE, diff_vals = NULL))
    }
  }
  else {
    # [DONE, not yet checked]
    # TO DO: Perform Kruskal-Wallis Test -> Perform Dunn's Test
    print("Performing Kruskal-Wallis Test...")
    kruskal_res <- kruskal.test(time ~ backend, data = data)
    p_val <- kruskal_res$p.value 
    print(paste0("Kruskal-Wallis Test is DONE w/ p-value ", p_val))
    
    if(p_val <= alpha_value) {
      print("At least one of the values is statistically different from the others.")
      print("Performing Dunn's Test...")
      # Perform Dunn's Test
      dunn_res <- dunnTest(time ~ backend, data = data, method="bonferroni")
      print("Dunn's Test is DONE. The summary of results is as follows:")
      print(dunn_res)
      # Filter the rows where the adjusted p-value is less than alpha
      significant_pairs <- dunn_res$res$Comparison[which(dunn_res$res$P.adj <= alpha_value)]
      return(list(is_diff = TRUE, significant_pairs = significant_pairs))
    }
    else {
      print("The values are not statistically different.")
      return(list(is_diff = FALSE, diff_vals = NULL))
    }
  }
}

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
bm_times <- 10 
bm_unit <- "milliseconds"
alpha_value <- 0.05

# Backends to compare
selected_backends <- c("NETLIB", "ATLAS", "OPENBLASSERIAL",
                       "MKLSERIAL", "BLISSERIAL")

for(k in k_vals) {
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
  
  is_diff <- check_stat_diff(use_anova, alpha_value)
  print(is_diff$significant_pairs[1])
  
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