# Calculate Signed Early Precision Ratio

Calculates activation and inhibitory early precision ratios using a
signed ground truth. These metrics evaluate positive and negative
regulatory edges separately, using the same random-baseline
normalization as EPR.

## Usage

``` r
calculate_signed_epr(network_table, ground_truth)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame with columns regulator, target and type (\`+\` or \`-\`)

## Value

A list containing activation and inhibitory EPR metrics
