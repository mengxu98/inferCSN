# Format a network table

Format a network table

## Usage

``` r
network_format(
  network_table,
  regulators = NULL,
  targets = NULL,
  abs_weight = TRUE
)
```

## Arguments

- network_table:

  Network edge table.

- regulators, targets:

  Nodes to include.

- abs_weight:

  Whether to use absolute weights and add interaction signs.

## Value

A formatted network edge table.
