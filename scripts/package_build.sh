#!/bin/bash

# If add new C/C++ function(s)
Rscript -e 'Rcpp::compileAttributes()'

# Build new document
Rscript -e 'devtools::document()'

# Fast check test and example(s)
Rscript -e 'devtools::test()'
Rscript -e 'devtools::run_examples()'

# Check R package
Rscript -e 'devtools::check()'
