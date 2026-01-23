# Changelog

## inferCSN 1.2.3

- **func**:
  - Optimize `check_parameters()` for check input parameters.
- **docs**:
  - Updated documentation across multiple functions.

## inferCSN 1.2.2

- **func**:
  - Removed `network_sift()` and `weight_sift()` functions.
  - Removed S4 method for `data.frame`.
- **docs**:
  - Updated documentation across multiple functions.

## inferCSN 1.2.1

- **func**:
  - Fixed warning on Windows platform for better compatibility.
- **docs**:
  - Updated documentation across multiple functions.

## inferCSN 1.2.0

CRAN release: 2025-10-16

- **func**:
  - Removed `srm` related functions and files.
  - Removed `RcppArmadillo` dependency.
  - Removed `Makevars` files for simplified build process.
  - Reorganized package structure by moving and removing files.
- **refactor**:
  - Updated package template structure.
  - Cleaned up package dependencies and build configuration.
- **docs**:
  - Updated package logo.
  - Improved documentation formatting.

## inferCSN 1.1.9

- **func**:
  - Enhanced logo printing functionality.
- **docs**:
  - Updated documentation across multiple functions.

## inferCSN 1.1.8

- **func**:
  - Added `inferCSN_logo()` function for displaying package logo.
  - Added new parameters for `log_message()` function.
  - Modified `.cores_detect()` function for better CPU core detection.
  - Removed `PKG_LIBS` from `Makevars` and `Makevars.win` files.
  - Moved utility functions to `thisutils` package.
- **refactor**:
  - Updated `pkgdown` webpage theme.
  - Fixed warnings and errors in GCC compilation.
  - Optimized document structure.
- **docs**:
  - Added examples in README.
  - Updated package description.
  - Updated `renv` configuration files.

## inferCSN 1.1.7

CRAN release: 2025-03-30

- **func**:
  - Added `sparse_cov_cor()` function for sparse covariance and
    correlation computation.
  - Removed deprecated `sparseCovCor` function.
  - Updated documentation for `sparse_cor()` function.

## inferCSN 1.1.6

CRAN release: 2025-03-27

- **func**:
  - Fixed error in `table_to_matrix()` function.
  - Optimized performance for matrix operations.
- **bugs**:
  - Fixed error in pkgdown build process.
- **docs**:
  - Added example for
    [`calculate_accuracy()`](https://mengxu98.github.io/inferCSN/reference/calculate_accuracy.md)
    function.
  - Modified examples for
    [`plot_network_heatmap()`](https://mengxu98.github.io/inferCSN/reference/plot_network_heatmap.md)
    function.
  - Updated badge display.
  - Updated ProjectId.

## inferCSN 1.1.3

- **func**:
  - Added new network evaluation functions for assessing network
    quality.
  - Updated import records for better dependency management.
- **docs**:
  - Edited documentation for multiple functions.
  - Modified documentation for
    [`plot_network_heatmap()`](https://mengxu98.github.io/inferCSN/reference/plot_network_heatmap.md)
    function.
  - Updated `test-coverage.yaml` file.
  - Renamed files for better organization.
- **refactor**:
  - Added check workflow for automated testing.
  - Modified import statements.

## inferCSN 1.1.2

- **func**:
  - Renamed `fit_sparse_regression()` to
    [`fit_srm()`](https://mengxu98.github.io/inferCSN/reference/fit_srm.md)
    for shorter function name.
  - Modified parameters for multiple functions.
- **docs**:
  - Updated package logo.
  - Updated badge display.

## inferCSN 1.1.1

- **func**:
  - Renamed function `fit_sparse_regression()` to
    [`fit_srm()`](https://mengxu98.github.io/inferCSN/reference/fit_srm.md).

## inferCSN 1.1.0

- **func**:
  - Modified `algorithm` parameter to non-exported for internal use.

## inferCSN 1.0.9

- **func**:
  - Added new utility functions including `%s%`, `matrix_to_table()`,
    `pearson_correlation()`, `simulate_sparse_matrix()`,
    [`plot_coefficients()`](https://mengxu98.github.io/inferCSN/reference/plot_coefficients.md),
    `split_indices()`, `subsampling_fun()`, and `weight_sift()`.
  - Added new C++ functions including
    [`network_format()`](https://mengxu98.github.io/inferCSN/reference/network_format.md),
    `prepare_calculate_metrics()` for improved performance.
  - Added `CITATION` file for proper citation.
  - Renamed function from `percent_samples` to `subsampling`.
  - Imported native R pipe operator `|>` for better compatibility with R
    \>= 4.1.0.
- **refactor**:
  - Renamed function `%s%` to `%ss%`.
  - Renamed functions to
    [`calculate_metrics()`](https://mengxu98.github.io/inferCSN/reference/calculate_metrics.md).
  - Modified structure of `sparse_regression()` output.
  - Modified normalization method implementation.
  - Updated `log_message()` function.
  - Improved documentation across multiple functions.
- **bugs**:
  - Fixed multiple coding style issues.
  - Fixed error in various functions.
  - Fixed warning when package installation.
  - Fixed download badge display.
- **docs**:
  - Optimized function documentation.
  - Updated `.gitignore`.
  - Updated packages record information.
- **deps**:
  - Deleted suggests packages that are no longer needed.
  - Deleted useless packages from dependencies.

## inferCSN 1.0.8

CRAN release: 2024-08-24

- **func**:
  - Added `log_message()` function for consistent logging across
    package.
  - Added export of `log_message()` for user access.
  - Modified documentation for multiple functions.
- **deps**:
  - Added import of `cli` package for better message formatting.
  - Modified Imports section.
  - Removed exports for internal functions.
- **docs**:
  - Updated README.
  - Improved documentation consistency.
- **refactor**:
  - Coding style improvements.

## inferCSN 1.0.7

- **func**:
  - Removed export of `prepare_network_data()` as internal function.
  - Modified documentation for multiple functions.
- **bugs**:
  - Fixed error for `normalization()` when setting `method` to
    `softmax`.

## inferCSN 1.0.6

- **func**:
  - Added new visualization functions:
    [`plot_embedding()`](https://mengxu98.github.io/inferCSN/reference/plot_embedding.md)
    and `plot_weight_distribution()`.
  - Added new parameters for `plot_weight_distribution()`.
  - Added `asMatrixParallel()` for parallel matrix operations.
  - Added `RcppParallel` support for improved performance.
  - Added `PKG_LIBS` information in Makevars.
  - Added new function `check_sparsity()` for sparsity checking.
  - Added new global variables.
  - Renamed `performance_evaluation` functions.
  - Renamed files and functions for better organization.
  - Support for sparse matrix operations.
  - Updated R version requirement.
  - Reduced computation complexity for better performance.
  - Removed export of `check_parameters()`.
- **bugs**:
  - Fixed error for `.weight_sift()`.
  - Fixed error and improved coding style in
    [`plot_scatter()`](https://mengxu98.github.io/inferCSN/reference/plot_scatter.md).
  - Fixed document issues.
- **docs**:
  - Modified documentation for multiple functions.
  - Modified examples, parameters, and documentation.
  - Updated `renv` files.
  - Updated `Makevars` files.
  - Updated document for some parameters.
- **refactor**:
  - Coding style improvements.
  - Renamed functions and files.

## inferCSN 1.0.5

CRAN release: 2024-06-26

- **func**:
  - Added `entropy` to filter edges in network.
  - Added new function `as_matrix()` for matrix conversion.
  - Added new function
    [`plot_scatter()`](https://mengxu98.github.io/inferCSN/reference/plot_scatter.md)
    for scatter plot visualization.
  - Added new parameter for
    [`plot_scatter()`](https://mengxu98.github.io/inferCSN/reference/plot_scatter.md).
  - Added new function `network_sift()` for network filtering (renamed
    from previous version).
  - Added new functions: `r_square()`,
    [`plot_contrast_networks()`](https://mengxu98.github.io/inferCSN/reference/plot_contrast_networks.md),
    [`plot_dynamic_networks()`](https://mengxu98.github.io/inferCSN/reference/plot_dynamic_networks.md).
  - Added new function `map_parallel()` and `parallelize_fun()` for
    parallel processing.
  - Added `.softmax()` function.
  - Added new parameter `heatmap_size_lock` for `network.heatmap()`.
  - Modified S4 method for `inferCSN`.
  - Removed `dynamic.networks`, `contrast.networks` functions
    (renamed/replaced).
  - Renamed `net.format` to
    [`network_format()`](https://mengxu98.github.io/inferCSN/reference/network_format.md).
  - Renamed `dynamic.networks` to
    [`plot_static_networks()`](https://mengxu98.github.io/inferCSN/reference/plot_static_networks.md).
  - Fixed errors and added new parameters for `network_sift()`.
  - Deleted unnecessary functions.
  - Modified import functions and packages.
  - Using `as_matrix()` instead of
    [`as.matrix()`](https://rdrr.io/r/base/matrix.html) for better
    performance.
- **deps**:
  - Removed `progress` package and added `pbapply` package.
  - Added new packages for enhanced functionality.
  - Updated imports.
- **bugs**:
  - Fixed error in `check.parameters()`.
  - Fixed multiple errors.
  - Fixed warning.
  - Fixed note info.
- **docs**:
  - Updated `badges`.
  - Modified examples.
  - Updated `_pkgdown.yml`.
  - Modified examples and figure style.
  - Updated `renv` file.
  - Modified parameter description of functions.
  - Updated `.gitignore`.
  - Updated README.
- **refactor**:
  - Coding style improvements across multiple files.
  - Improved function organization.

## inferCSN 1.0.3

CRAN release: 2024-04-17

- **func**:
  - Added new function `weight_filter()` for filtering weights.
  - Added new parameters for `network.heatmap()`.
  - Set `crossweight` example to `\dontrun{}` for CRAN compatibility.
- **docs**:
  - Updated NAMESPACE and import records.

## inferCSN 1.0.2

- **func**:
  - Added new function `normalization()` for data normalization.
  - Added new function `crossweight()` for cross-weight analysis.
  - Added new function `crossweight_params()`.
  - Added new function `extract_net()` (later deleted).
  - Added new example data `example_meta_data`.
  - Updated `crossweight` and related functions.
  - Updated example data.
  - Added new parameters for `calculate.gene.rank()`.
  - Renamed function to `filter_sort.matrix` to handle CRAN warning.
  - Renamed multiple functions for consistency.
  - Deleted function `extract_net()`.
- **bugs**:
  - Fixed error when `sample`.
  - Fixed multiple errors.
  - Fixed pkgdown error.
- **docs**:
  - Updated README.
  - Updated `crossweight.Rd`.
  - Updated `_pkgdown.yml`.
  - Updated NAMESPACE.
  - Updated `renv.lock`.
  - Modified examples.
  - Formatted progress information.
- **refactor**:
  - Fixed coding style issues.

## inferCSN 1.0.1

CRAN release: 2024-03-15

- **func**:
  - Added new function `filter.sort.matrix()` for efficient matrix
    filtering and sorting.
  - Optimized memory consumption for handling large datasets.
  - Modified function used to handle data.
  - Modified parameter `k_folds`.
  - Modified parameter `regulators_num`.
  - Renamed multiple functions for consistency: `calculate.gene.rank()`,
    `sub.model.fit()`, `model.fit()`, `single.network()`.
  - Added examples for multiple functions including `network.heatmap()`,
    `dynamic.networks()`.
  - Added new parameter for `net.format()`.
  - Modified some functions for better performance.
  - Fixed error code.
- **bugs**:
  - Fixed multiple errors.
- **docs**:
  - Updated README with improved examples.
  - Updated `_pkgdown.yml`.
  - Updated `.gitignore`.
  - Updated `cran_release.R`.
  - Updated `.Rbuildignore`.
  - Modified parameter information.
  - Added test examples for `table_to_matrix()`.
  - Updated logo.svg.
- **refactor**:
  - Removed test files.
  - Updated NAMESPACE.
  - Updated renv file.
  - Added examples throughout.

## inferCSN 1.0.0

CRAN release: 2024-02-09

- **feat**:
  - First major release of inferCSN package.
  - Core functionality for inferring cell-specific gene regulatory
    networks.
  - Sparse regression-based network inference algorithm.
  - Support for single-cell RNA-seq data.
  - Network visualization and analysis functions.
  - Performance evaluation metrics.
- **docs**:
  - Complete documentation for all exported functions.
  - Comprehensive README with usage examples.
  - Package website with pkgdown.

## inferCSN 0.99.9

CRAN release: 2024-01-10

- **func**:
  - Updated `table_to_matrix.cpp` for improved performance.
  - Updated version to 0.99.9 in preparation for CRAN submission.
- **docs**:
  - Updated `cran_release.R`.
  - Fixed example information.
  - Added new function and renamed files.

## inferCSN 0.99.8

CRAN release: 2023-12-04

- **func**:
  - Updated
    [`inferCSN()`](https://mengxu98.github.io/inferCSN/reference/inferCSN.md)
    function with improvements.
  - Updated C++ functions for better performance.
  - Updated and formatted code across multiple files.
  - Updated example data.
- **bugs**:
  - Fixed pkgdown error.
- **docs**:
  - Fixed error and modified documentation.
  - Updated `renv.lock`.
  - Updated `requirements.sh`.

## inferCSN 0.99.7

CRAN release: 2023-10-30

- **func**:
  - Fixed error when compiling on CentOS system.
  - Updated `requirements.sh`.
- **refactor**:
  - Version update for CRAN submission preparation.

## inferCSN 0.1.1

- **func**:
  - Initial package development.
  - Basic sparse regression implementation.
  - C++ integration for performance.
- **docs**:
  - Basic documentation structure.
  - Initial README.

## inferCSN 0.1.0

- **feat**:
  - Initial version.
