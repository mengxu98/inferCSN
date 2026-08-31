# Rank TFs and genes in network

Rank TFs and genes in network

## Usage

``` r
calculate_gene_rank(
  network_table,
  regulators = NULL,
  targets = NULL,
  directed = FALSE
)
```

## Arguments

- network_table:

  Network edge table.

- regulators, targets:

  Nodes to include.

- directed:

  Whether the network is directed.

## Value

A table of gene rank.

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2026-08-31 02:39:55] Inferring network for <matrix/array>...
#> ◌ [2026-08-31 02:39:55] Checking parameters...
#> ✔ [2026-08-31 02:39:55] Inferring network done
#> ℹ [2026-08-31 02:39:55] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6
head(calculate_gene_rank(network_table))
#>   gene rank_value regulator
#> 1   g6  0.2148641      TRUE
#> 2   g2  0.1918057      TRUE
#> 3   g4  0.1737057      TRUE
#> 4   g5  0.1715660      TRUE
#> 5   g1  0.1344885      TRUE
#> 6   g3  0.1135699      TRUE
head(calculate_gene_rank(network_table, regulators = "g1"))
#>   gene rank_value regulator
#> 1   g1  0.4864865      TRUE
#> 2   g5  0.3084459     FALSE
#> 3   g6  0.2050676     FALSE
head(calculate_gene_rank(network_table, targets = "g1"))
#>   gene rank_value regulator
#> 1   g1  0.4864865     FALSE
#> 2   g6  0.3118919      TRUE
#> 3   g5  0.2016216      TRUE
```
