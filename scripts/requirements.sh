#!/bin/bash

echo "Checking and install required packages......"

# Required packages
Rscript -e 'install.packages(c("BiocManager", \
                               "devtools", \
                               "doParallel", \
                               "doRNG", \
                               "doSNOW", \
                               "dplyr", \
                               "foreach", \
                               "igraph", \
                               "Kendall", \
                               "modEvA", \
                               "patchwork", \
                               "progress", \
                               "snow", \
                               "ggplot2", \
                               "RcppArmadillo"), \
                             repos = "https://cloud.r-project.org")'

Rscript -e 'install.packages("digest", \
                             repos = c("https://eddelbuettel.r-universe.dev", \
                                       "https://cloud.r-project.org"))'

Rscript -e 'BiocManager::install(c("ComplexHeatmap", "edgeR"))'

Rscript -e 'devtools::install_github("statOmics/tradeSeq")'

