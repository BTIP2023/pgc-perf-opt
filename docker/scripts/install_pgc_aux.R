# Install other R tweaks

# ggbiplot: https://github.com/vqv/ggbiplot#installation
library(devtools)
install_github("vqv/ggbiplot")

# knitr 1.42 for presentations
library(remotes)
remotes::install_version("knitr", version = "1.42",
                          repos = "https://cran.r-project.org")

# save_image
# In Linux, so no need for
# `reticulate::conda_install("r-reticulate", "python-kaleido==0.1.*"")
library(reticulate)
reticulate::install_miniconda()
reticulate::conda_install(packages = "python=3.10.6")
reticulate::conda_install("r-reticulate", "python-kaleido")
reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
reticulate::use_miniconda("r-reticulate")