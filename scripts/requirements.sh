#!/bin/bash

# Function to check if R is installed
check_r_installed() {
    if command -v R >/dev/null 2>&1; then
        echo "R is already installed......"
        install_packages
    else
        echo "R is not installed......"
        install_r
        install_packages
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

# Function to check and install required packages......
install_packages() {
	echo "Checking and install required packages......"

	# Required packages
	Rscript -e 'install.packages(c("BiocManager", "devtools", "doSNOW", "dplyr", "foreach"), repos = "https://cloud.r-project.org")'
	Rscript -e 'install.packages(c("igraph", "magritte", "patchwork", "progress", "purrr"), repos = "https://cloud.r-project.org")'
	Rscript -e 'install.packages(c("Rcpp", "snow", "ggplot2", "RcppArmadillo", "cowplot"), repos = "https://cloud.r-project.org")'
	Rscript -e 'install.packages(c("gtools", "circlize", "precrec", "ggnetwork"), repos = "https://cloud.r-project.org")'

	# Note: the package 'digest' may cause errors during the installation process on the macOS system
	Rscript -e 'install.packages("digest", repos = c("https://eddelbuettel.r-universe.dev", "https://cloud.r-project.org"))'

	Rscript -e 'BiocManager::install("ComplexHeatmap")'
}

# Check if R is installed
check_r_installed
