# Calculate Spearman Stability Across Runs

Calculates the median and dispersion of pairwise Spearman correlations
across multiple predicted networks. This metric measures ranking
stability across repeated runs or perturbations.

## Usage

``` r
calculate_stability_spearman(network_tables, ground_truth)
```

## Arguments

- network_tables:

  A list of predicted network tables

- ground_truth:

  A data frame of ground truth network

## Value

A list containing Spearman stability metrics
