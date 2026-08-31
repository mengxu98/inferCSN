# Calculate Accuracy

Calculates overall classification accuracy for predicted edges against a
ground truth network after optimal thresholding.

## Usage

``` r
calculate_accuracy(network_table, ground_truth, tf_edges = FALSE)
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

A list with a \`metrics\` data frame containing Accuracy.
