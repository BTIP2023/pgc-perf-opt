# Load required R packages
library(tidyverse)
library(highcharter) 

# Set highcharter options
options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)))

# Load a demo data
data("mpg", package = "ggplot2")

# Summary table
summary.table <- mpg %>% 
  group_by(manufacturer) %>% 
  summarise(
    nb_cars = n(), 
    nb_model = length(unique(model))
  ) %>% 
  arrange(-nb_cars, -nb_model)
summary.table

hc <- summary.table %>%
  hchart(
    "treemap", 
    hcaes(x = manufacturer, value = nb_cars, color = nb_model)
  )