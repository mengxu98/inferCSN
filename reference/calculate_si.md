# Calculate Set Intersection

Calculates the Set Intersection (SI) — the number of predicted edges
that exactly match true edges in the ground truth network.

## Usage

``` r
calculate_si(network_table, ground_truth, tf_edges = FALSE)
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

A list with a \`metrics\` data frame containing SI.
