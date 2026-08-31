# Calculate Motif Ratios

Calculates mutual interaction, feedforward-loop and feedback-loop ratios
against the ground truth network. These ratios compare local network
motif counts in the predicted top-ranked network against the
corresponding motif counts in the reference network.

## Usage

``` r
calculate_motif_ratios(
  network_table,
  ground_truth,
  top_k = NULL,
  strategy = "beeline"
)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame of ground truth network

- top_k:

  Number of top-ranked edges to retain. Defaults to the number of
  ground-truth edges.

- strategy:

  Edge-selection strategy. Use \`"beeline"\` to match the BEELINE motif
  evaluator, or \`"default"\` to use the package's general top-ranked
  edge selector. Default is \`"beeline"\`.

## Value

A list containing motif ratios
