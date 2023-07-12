# Code Guide

## Prerequisites

### R/Rstudio
The project's R/RStudio codebase imports several packages,
but this is already managed in-line using `pacman`, an R package manager.
Thus, the R code already handles all the prerequisite package installations
by itself.

### Python/Conda/Jupyter
In contrast, a development environment must be manually set up for the project's
Python/Jupyter codebase. Refer to the list below for prerequisite Python packages:
- pandas=2.0.1
- numpy=1.24.3
- scikit-learn=1.2.2
- matplotlib=3.7.1
- umap-learn=0.5.3
- plotly=5.15.0
- python-kaleido=0.1.0
- seaborn=0.12.2
- pip (optional)
  
Needed for jupyter notebook intergration in custom environment
- nb_conda_kernels
- ipykernel
- ipywidgets

The steps below would guide you through the process.
1. Install **Conda** which is bundled with [**Anaconda**](https://www.anaconda.com/download). 
Conda is the recommended virtual environment and package
manager for the Python code.
2. Once Conda has been installed, open the Anaconda Navigator.
3. There, go to `Environments` (can be clicked in the left tab).
4. Click `Create` to create a new Conda environment for this project.
Only the Python environment is required, for which you may use Python version 
3.11.4. Optionally, you may also create an R environment with R version 3.6.1.
5. After creating the environment, search the above-mentioned prerequisite Python packages 
and install them in the newly created environment. `conda install -n <env-name> <pythonpackage=version>`
6. Activate the environment by entering `conda activate <env-name>` in a terminal.
7. [Optional] Change to correct drive letter if working directory is on a different drive.
8. Run jupyter notebook environment `jupyter notebook`
9. You are now ready to begin Python development using Juypyter Notebooks in this project!

For further assistance in setting-up your Conda environment,
you may visit the following:
- https://conda.io/projects/conda/en/latest/user-guide/getting-started.html
- https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

*Remarks:*
Some IDEs do not default to the Conda environment, hence ensure that
your IDE is set to use the new environment.
Note that any new package will have to be added to the environment.

This approach to package management in Python will be maintained
while a container for the project is still being set up.

## Data Parsing, Preprocessing, and k-mer Analysis
Generation of the data to be fed into the dimensionality reduction and results analysis
workflow is handled by `kmer-analysis.R` which sources `preprocess.R` and `helper.R`.

In kmer-analysis.R, the following parameters may be set:
- `seed` for reproducibility of random values.
- `stratSize` for number of samples per variant.
- `write_fastacsv` for turning fasta and csv writeback on or off.

This code outputs intermediate files (`fastaAll.fasta` and metaDataAll.csv`) to `data/interm`, and
finally outputs the generated k-mers augmented with metadata to `data/kmers`. These files have
been .gitignored because their sizes exceed the Github limit, and because of the Github cap on LFS bandwidth.

Note: Warnings that appear (parsing errors due to conflicting column types and values,
and NAs due to type coercion) can be ignored as these are handled in the code.

### Note on Generated Files
The generated files are timestamped and persistent (not overwritten). Because
they can be huge, especially for large stratSize, keep an eye on their space usage
and prune old and unused generated files in `data/interm` and `data/kmers`.

## TODOs
- [ ] Update this README once main pipeline has been established wherein kmer-analysis.R is also just another source.
- [ ] Add auto-pruner for old generated files (might do timestamp arithmetic)
- [ ] Add auto-setter to latest generated files in `dim-reduce.R` and similar files.