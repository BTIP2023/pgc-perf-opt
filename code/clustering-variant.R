# Reference: https://github.com/ai-covariants/analysis-mutations
# File: kmer-analysis.R

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
dendrogram_create = function(filePath)
{
  data <- read_csv(filePath)
  df = subset(data, select = -c(...1))
  dat <- df %>%
    mutate(sample_name = paste('var', seq(1:nrow(data)), sep = '_')) # 
  metadata <- dat %>%
    select(sample_name, variant)
  numeric_data <- dat %>% select(-c(variant))
  
  # normalize data to values from 0 to 1 
  numeric_data_norm <- numeric_data %>%
    select(sample_name, everything()) %>%
    pivot_longer(cols = 2:"TTTTTTT", values_to = 'value', names_to = 'type') %>%
    group_by(type) %>%
    mutate(value_norm = (value-min(value))/(max(value)-min(value))) %>% # normalize data to values 0–1
    select(sample_name, value_norm) %>%
    pivot_wider(names_from = 'type', values_from = 'value_norm') %>%
    column_to_rownames('sample_name')
  
  # create dendrogram from AGNES
  model <- agnes(numeric_data_norm, metric = "euclidean",method="single")
  dendrogram <- as.dendrogram(model)
  plot(dendrogram)
  
  # extract dendrogram segment data
  dendrogram_data <- dendro_data(dendrogram)
  dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data
  print(dendrogram_segments)

  # get terminal dendrogram segments
  dendrogram_ends <- dendrogram_segments %>%
    filter(yend == 0) %>% # filter for terminal dendrogram ends
    left_join(dendrogram_data$labels, by = 'x') %>% # .$labels contains the row names from dist_matrix (i.e., sample_name)
    rename(sample_name = label) %>%
    left_join(metadata, by = 'sample_name')
  dendrogram_end<-subset(dendrogram_ends,sample_name!="<NA>")
  #print(dendrogram_end)

  # Generate custom color palette for dendrogram ends based on metadata variable
  # unique_vars <- levels(factor(dendrogram_ends$target)) %>%
  #   as.data.frame() %>% rownames_to_column("row_id")
  # color_count <- length(unique(unique_vars$.)) # count number of unique variables
  # get_palette <- colorRampPalette(brewer.pal(n = 8, name = "Set1")) # RColorBrewer
  # palette <- get_palette(color_count) %>% # produce RColorBrewer palette based on number of unique variables in metadata
  #   as.data.frame() %>%
  #   rename("color" = ".") %>%
  #   rownames_to_column(var = "row_id")
  # color_list <- left_join(unique_vars, palette, by = "row_id") %>%
  #   select(-row_id)
  # species_color <- as.character(color_list$color)
  # names(species_color) <- color_list$.
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
  
  # Saves dendogram as an RData file in the presentation 1 assets folder
  save(ggp, file = file.path('presentations/research-updates-1/presentation_files/', "clustering-variant.RData"))
  
  #Saves dendogram as an RData file and PNG in the results folder
  save(ggp, file = file.path('results/dendogram/', 'clustering-variant.RData'))
  ggsave(p, file = file.path('results/dendogram/', 'clustering-variant.png'))

}

dendrogram_create('data/kmers/kmer_7.csv')

# CLEAN UP #################################################

# Clear environment
rm(list = ls())

# Clear packages
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base

# Clear plots but only if there IS a plot
# while (!is.null(dev.list())) dev.off()

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
