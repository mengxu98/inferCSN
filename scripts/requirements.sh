#!/bin/bash

# Function to check if R is installed
check_r_installed() {
    if command -v R >/dev/null 2>&1; then
        echo "R is already installed......"
        exit 0
    else
        echo "R is not installed......"
        return 1
    fi
}

# Function to install R
install_r() {
    echo "Installing R......"
    # Add any additional commands required to install R here
    # For example, on Ubuntu-based systems, you can use the following command:
    sudo apt-get update
    sudo apt-get install -y r-base
}

# Check if R is installed
check_r_installed || install_r

# Checking and install required packages......
echo "Checking and install required packages......"

# Required packages
Rscript -e 'install.packages(c("BiocManager", \
                               "devtools", \
                               "doSNOW", \
                               "dplyr", \
                               "foreach", \
                               "igraph", \
                               "magritte", \
                               "patchwork", \
                               "progress", \
                               "purr", \
                               "Rcpp", \
                               "snow", \
                               "ggplot2", \
                               "RcppArmadillo", \
                               "cowplot", \
                               "gtools", \
                               "circlize", \
                               "Kendall", \
                               "precrec", \
                               "ggnetwork"), \
                             repos = "https://cloud.r-project.org")'

Rscript -e 'install.packages("digest", \
                             repos = c("https://eddelbuettel.r-universe.dev", \
                                       "https://cloud.r-project.org"))'

Rscript -e 'BiocManager::install(c("ComplexHeatmap"))'


