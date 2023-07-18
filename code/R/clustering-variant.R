# File: clustering-variant.R

# INSTALL AND LOAD PACKAGES ################################
library(datasets)  # Load base packages manually

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman", repos = "https://cran.case.edu")
library(pacman)

# Use pacman to load add-on packages as desired
# First call are standard packages for the project
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, psych,
               rio, rmarkdown, shiny, 
               stringr, tidyr, tidyverse)
# Second call are file-specific packages
pacman::p_load(ggdendro, RColorBrewer, readr, cluster,
               dendextend, colorspace)

# WORK WITH DATA ###########################################

dendrogram_create_variant <- function(k, data_path, results_path)
{
  
  get_time <- function(string) {
    parts <- strsplit(string, "_")[[1]]
    as.numeric(gsub(".csv", "", parts[3]))
  }
  
  file_pattern <- paste0("kmer_", k, "_", ".*\\.csv$")
  file_list <- list.files(
    path = data_path, pattern = file_pattern,
    full.names = FALSE
  )
  
  # Check if any files are found
  if (length(file_list) == 0) {
    message("No files found for k = ", k)
    return(NULL)
  }
  
  # Sort the strings based on the timestamp in descending order
  sorted_strings <- file_list[order(sapply(file_list, get_time),
                                    decreasing = TRUE
  )]
  
  data <- read_csv(paste(data_path, sorted_strings[1], sep = "/"))
  
  # read kmer file
  df = subset(data, select = -c(...1))
  dat <- df %>%
    mutate(sample_name = paste('var', seq(1:nrow(data)), sep = '_')) # 
  metadata <- dat %>%
    select(sample_name, variant)
  numeric_data <- dat %>% select(-c(variant))
  
  # normalize data to values from 0 to 1
  slice_col <- which(colnames(data) == "strain")
  numeric_data_norm <- numeric_data %>%
    select(sample_name, everything()) %>%
    pivot_longer(cols = 2:(slice_col - 1), values_to = 'value', names_to = 'type') %>%
    group_by(type) %>%
    mutate(value_norm = (value-min(value))/(max(value)-min(value))) %>% # normalize data to values 0–1
    select(sample_name, value_norm) %>%
    pivot_wider(names_from = 'type', values_from = 'value_norm') %>%
    column_to_rownames('sample_name')
  
  # create dendrogram from AGNES
  model <- agnes(numeric_data_norm, metric = "euclidean",method="single")
  dendrogram <- as.dendrogram(model)
  
  # extract dendrogram segment data
  dendrogram_data <- dendro_data(dendrogram)
  dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data

  # get terminal dendrogram segments
  dendrogram_ends <- dendrogram_segments %>%
    filter(yend == 0) %>% # filter for terminal dendrogram ends
    left_join(dendrogram_data$labels, by = 'x') %>% # .$labels contains the row names from dist_matrix (i.e., sample_name)
    dplyr::rename(sample_name = label) %>%
    left_join(metadata, by = 'sample_name')
  dendrogram_end<-subset(dendrogram_ends,sample_name!="<NA>")

  # generate color variant color palette
  variant_color <- brewer.pal(n = 6, name = 'Paired')

  p <- ggplot() +
    geom_segment(data = dendrogram_segments,
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_segment(data = dendrogram_end,
                 aes(x=x, y=y.x, xend=xend, yend=yend, color = variant, text = paste('sample name: ', sample_name,
                                                                                      '<br>',
                                                                                      'Variant: ', variant))) + # test aes is for plotly
    scale_color_manual(values = variant_color) + scale_y_reverse() + coord_flip() +
    theme_bw() + theme(legend.position = 'right') + ylab('Distance') + xlab('Sequence')
 
   ggp <- ggplotly(p)
  ggp
  
  #Saves dendrogram as an RData file and PNG in the results folder
  save(ggp, file = file.path(results_path, paste0('clustering-variant', "-", k, ".RData")))
  ggsave(p, file = file.path(results_path, paste0('clustering-variant', "-", k, ".PNG")))

}
