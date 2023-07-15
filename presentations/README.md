# Presentations Folder
Presentations for research updates and the final presentation

## Creating a new presentation
* create a copy of 'presentation-template' directory and rename to 'research-updates-N'
* modify 'presentation.qmd' file in Rstudio
* add presentation image assets to 'presentation_assets/img' directory

## Dependencies
* 'knitr' package in R must be in version 1.42

To install knitr 1.42 using only the R console, do the following:
```
install.packages('remotes', repos = "https://cran.r-project.org")
require(remotes)
remotes::install_version('knitr', version = '1.42', repos = "https://cran.r-project.org")
detach(remotes, unload=TRUE)
```

