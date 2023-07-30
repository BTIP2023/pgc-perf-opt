# pgc-perf-opt
Performance Evaluation and Optimization of an Unsupervised Machine Learning Pipeline for Discriminating Major SARS-CoV-2 Variants in the Philippines

A highly scaleable unsupervised machine learning (UML) pipeline for discriminating major SARS-CoV-2 variants in the Philippines. Also includes performance benchmarks of said workflow according to the applied linear algebra library backend and processor vulnerability mitigations.

**Research Objectives:**
1. Create an improved UML variant discriminator pipeline that can work on any raw GISAID input, esp. for Philippines SARS-CoV-2 data, with enhanced usability, efficiency, and scalability.
2. Understand the impact of an optimized linear algebra library backend on the UML bioinformatics workflow.
3. Understand the impact of processor vulnerability mitigations on the UML bioinformatics workflow.

## Data
This section describes the characteristics of the data such as sources, sizes, and formats, motivating the approach to data wrangling and downstream analysis. Different data configurations may also be used in different benchmarks.

### On GISAID Data
For data obtained from GISAID, only the **accession numbers** are referenced in this remote repository instead of the actual dataset. GISAID (2012) gives the following reason:

> GISAID does not promote the release of data to databases where access to data is anonymous and the rights of the submitter are relinquished.  GISAID already provides the public with open access to data in a transparent way.

Please see `data/README.md` for further instructions.

The accession numbers can be found in `data/overview/accession.all`.

## Code
This section contains the code for the variant discriminator workflow. Said workflow found in `code/pipeline-classic.R` has the following structure:
![code/pipeline-classic.R](presentations/pipeline-flowchart.png)

Raw GISIAD data is placed in `data/GISAID`. The source `code` then does the following:
- Data extraction, wrangling, sanitation, overview compilation, and augmentation of `data/GISAID`.
- k-mer counting which produces augmented k-mer, metadata, and heatmap files in `data/kmers` and `data/overview`.
- UML dimensionality reduction and clustering techniques, the results of which are stored in `results` and presented in `presentations`.

Pilot code may also be written for more novel workflows. For instance, see [this](https://www.frontiersin.org/articles/10.3389/fbioe.2015.00035/full) wheat mutation analysis article.

## Results
This section will summarize the performance evaluations of each benchamrk. Possibly, results from novel bioinformatics workflows or benchmark approaches may also be discussed here.

## Benchmarks
Relevant hardware, software, and data configurations must be explicitly noted for each benchmark.

## Presentations
This section will contain directories for the research proposal presentations, research updates and the final research presentation.

## Custom Plotting
Factor1 = variant
Factor2 = div
filterf1 <- (“omicron” “omicron_sub”)
Filterf2 <- (“all regions from luzon”)

Filter_pca <- (pca_df, data, k, filterf1, filterf2),

 pca_plot <- function(pca_df, data, k, color=factor1, filterf1, shape=factor2, filterf2)

Filtertsne
tsneplot

Filterumap
Umap plot



---
## References
Chandra, R., Bansal, C., Kang, M., Blau, T., Agarwal, V., Singh, P., Wilson, L. O. W., & Vasan, S. (2023). Unsupervised machine learning framework for discriminating major variants of concern during COVID-19. *PLOS ONE, 18(5),* e0285719. https://doi.org/10.1371/journal.pone.0285719. Reference repository at [ai-covariants/analysis-mutations](https://github.com/ai-covariants/analysis-mutations).

GISAID. (2012). *FAQ.* https://gisaid.org/help/faq/

Kessler, N., Bonte, A., Albaum, S. P., Mäder, P., Messmer, M., Goesmann, A., Niehaus, K., Langenkämper, G., & Nattkemper, T. W. (2015). Learning to Classify Organic and Conventional Wheat – A Machine Learning Driven Approach Using the MeltDB 2.0 Metabolomics Analysis Platform. *Frontiers in Bioengineering and Biotechnology, 3.* https://www.frontiersin.org/articles/10.3389/fbioe.2015.00035
