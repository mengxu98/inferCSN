# Calculate Jaccard Stability Across Runs

Calculates the median and dispersion of pairwise top-k Jaccard overlap
across multiple predicted networks. This metric measures edge-set
stability across repeated runs or perturbations.

## Usage

``` r
calculate_stability_jaccard(network_tables, ground_truth)
```

## Arguments

- network_tables:

  A list of predicted network tables

- ground_truth:

  A data frame of ground truth network

## Value

A list containing Jaccard stability metrics
