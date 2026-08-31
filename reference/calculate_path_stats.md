# Calculate Path Statistics

Characterises false-positive predicted edges by their shortest-path
distances in the ground truth network. This helps distinguish arbitrary
false positives from edges that are close to the true network topology,
such as shortcut edges across 2-5 step paths.

## Usage

``` r
calculate_path_stats(network_table, ground_truth, top_k = NULL)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame of ground truth network

- top_k:

  Number of top-ranked edges to retain. Defaults to the number of
  ground-truth edges.

## Value

A list containing path statistics
