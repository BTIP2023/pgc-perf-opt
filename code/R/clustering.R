dendrogram_create_variant <- function(k, kmers, results_path, include_plots = FALSE) { 
  # process kmers dataframe, note that this was read outside the function 
  df <- kmers
  dat <- df %>%
    mutate(sample_name = paste('var', seq(1:nrow(df)), sep = '_')) # 
  metadata <- dat %>%
    select(sample_name, variant)
  numeric_data <- dat %>% select(-c(variant))
  
  # normalize data to values from 0 to 1
  print("Pre-processing data...")
  slice_col <- which(colnames(df) == "strain")
  numeric_data_norm <- numeric_data %>%
    select(sample_name, everything()) %>%
    pivot_longer(cols = 2:(slice_col - 1), values_to = 'value', names_to = 'type') %>%
    group_by(type) %>%
    mutate(value_norm = (value-min(value))/(max(value)-min(value))) %>% # normalize data to values 0–1
    select(sample_name, value_norm) %>%
    pivot_wider(names_from = 'type', values_from = 'value_norm') %>%
    column_to_rownames('sample_name')
  
  # create dendrogram from AGNES
  print("Creating dendrogram...")
  model <- agnes(numeric_data_norm, metric = "euclidean",method="single")
  dendrogram <- as.dendrogram(model)
  
  # extract dendrogram segment data
  print("Plotting dendrogram...")
  dendrogram_data <- dendro_data(dendrogram)
  dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data
  
  # get terminal dendrogram segments
  dendrogram_ends <- dendrogram_segments %>%
    filter(yend == 0) %>% # filter for terminal dendrogram ends
    left_join(dendrogram_data$labels, by = 'x') %>% # .$labels contains the row names from dist_matrix (i.e., sample_name)
    dplyr::rename(sample_name = label) %>%
    left_join(metadata, by = 'sample_name')
  dendrogram_end<-subset(dendrogram_ends,sample_name!="<NA>")
  
  if (!include_plots) {
    return()
  }
  
  # generate color variant color palette
  variant_color <- brewer.pal(n = 6, name = 'Paired')
  
  p <- ggplot() +
    geom_segment(data = dendrogram_segments,
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_segment(data = dendrogram_end,
                 aes(x=x, y=y.x, xend=xend, yend=yend, color = variant, label = paste('sample name: ', sample_name,
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

dendrogram_create_region <- function(k, kmers, results_path, include_plots = FALSE) {
  # process kmers dataframe, note that this was read outside the function 
  df <- kmers
  dat <- df %>%
    mutate(sample_name = paste('var', seq(1:nrow(df)), sep = '_')) # 
  metadata <- dat %>%
    select(sample_name, division_exposure)
  numeric_data <- dat %>% select(-c(division_exposure))
  
  # normalize data to values from 0 to 1
  slice_col <- which(colnames(df) == "strain")
  numeric_data_norm <- numeric_data %>%
    select(sample_name, everything()) %>%
    pivot_longer(cols = 2:(slice_col - 1), values_to = 'value', names_to = 'type') %>%
    group_by(type) %>%
    mutate(value_norm = (value-min(value))/(max(value)-min(value))) %>% # normalize data to values 0–1
    select(sample_name, value_norm) %>%
    pivot_wider(names_from = 'type', values_from = 'value_norm') %>%
    column_to_rownames('sample_name')
  
  # create dendrogram from AGNES
  model <- agnes(numeric_data_norm, metric = "euclidean")
  dendrogram <- as.dendrogram(model)
  plot(dendrogram)
  
  # extract dendrogram segment data
  dendrogram_data <- dendro_data(dendrogram)
  dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data
  head(dendrogram_segments)
  
  # get terminal dendrogram segments
  dendrogram_ends <- dendrogram_segments %>%
    filter(yend == 0) %>% # filter for terminal dendrogram ends
    left_join(dendrogram_data$labels, by = 'x') %>% # .$labels contains the row names from dist_matrix (i.e., sample_name)
    dplyr::rename(sample_name = label) %>%
    left_join(metadata, by = 'sample_name') 
  
  if (!include_plots) {
    return()
  }
  
  # Generate custom color palette for dendrogram ends based on metadata variable
  unique_vars <- levels(factor(dendrogram_ends$division_exposure)) %>% 
    as.data.frame() %>% rownames_to_column('row_id') 
  # count number of unique variables
  color_count <- length(unique(unique_vars$.))
  # get RColorBrewer palette
  get_palette <- colorRampPalette(brewer.pal(n = 12, name = 'Paired'))
  # produce RColorBrewer palette based on number of unique variables in metadata:
  palette <- get_palette(color_count) %>% 
    as.data.frame() %>%
    dplyr::rename('color' = '.') %>%
    rownames_to_column(var = 'row_id')
  color_list <- left_join(unique_vars, palette, by = 'row_id') %>%
    select(-row_id)
  var_color <- as.character(color_list$color)
  names(var_color) <- color_list$.
  
  p <- ggplot() +
    geom_segment(data = dendrogram_segments, 
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_segment(data = dendrogram_ends,
                 aes(x=x, y=y.x, xend=xend, yend=yend, color = division_exposure, label = paste('sample name: ', sample_name,
                                                                                     '<br>',
                                                                                     'Region ', division_exposure))) + # test aes is for plotly
    scale_color_manual(values = var_color) + scale_y_reverse() + coord_flip() +
    theme_bw() + theme(legend.position = 'right') + ylab('Distance') + xlab('Sequence') +
    guides(color = guide_legend(ncol = 1))
  
  ggp <- ggplotly(p)
  ggp
  
  #Saves dendogram as an RData file and PNG in the results folder
  save(ggp, file = file.path(results_path, paste0('clustering-region', "-", k, ".RData")))
  ggsave(p, file = file.path(results_path, paste0('clustering-region', "-", k, ".PNG")), 
         limitsize=FALSE, height=7, width=20)
}

