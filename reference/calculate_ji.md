# Calculate Jaccard Index

Calculates the Jaccard Index (JI) — the size of the intersection divided
by the size of the union of predicted and true edge sets.

## Usage

``` r
calculate_ji(network_table, ground_truth, tf_edges = FALSE)
```

## Arguments

- network_table:

  A data frame of predicted network structure.

- ground_truth:

  A data frame of ground truth network.

- tf_edges:

  Whether to restrict the candidate edge universe to TF-to-gene edges.
  Default is \`FALSE\`.

## Value

A list with a \`metrics\` data frame containing JI.
