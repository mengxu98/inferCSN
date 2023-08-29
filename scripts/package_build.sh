#!/bin/bash

# If add new C++ functions
Rscript -e 'Rcpp::compileAttributes()'

# Build new document
Rscript -e 'devtools::document()'

# Check R package
Rscript -e 'devtools::check()'
