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

  The weight data table of network.

- regulators:

  Regulators list.

- targets:

  Targets list.

- directed:

  Whether the network is directed.

## Value

A table of gene rank.

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2026-01-23 02:15:53] Inferring network for <matrix/array>...
#> ◌ [2026-01-23 02:15:53] Checking parameters...
#> ℹ [2026-01-23 02:15:53] Using L0 sparse regression model
#> ℹ [2026-01-23 02:15:53] Using 1 core
#> ⠙ [2026-01-23 02:15:53] Running for g1 [1/18] ■■■                              …
#> ✔ [2026-01-23 02:15:53] Completed 18 tasks in 173ms
#> 
#> ℹ [2026-01-23 02:15:53] Building results
#> ✔ [2026-01-23 02:15:53] Inferring network done
#> ℹ [2026-01-23 02:15:53] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
head(calculate_gene_rank(network_table))
#>   gene rank_value regulator
#> 1  g18 0.05883983      TRUE
#> 2  g17 0.05873281      TRUE
#> 3   g9 0.05869090      TRUE
#> 4   g5 0.05863236      TRUE
#> 5   g7 0.05811175      TRUE
#> 6   g6 0.05756879      TRUE
head(calculate_gene_rank(network_table, regulators = "g1"))
#>   gene rank_value regulator
#> 1   g1 0.46396396      TRUE
#> 2   g2 0.18784178     FALSE
#> 3  g18 0.13390621     FALSE
#> 4  g17 0.02111819     FALSE
#> 5   g5 0.02038973     FALSE
#> 6   g9 0.01973099     FALSE
head(calculate_gene_rank(network_table, targets = "g1"))
#>   gene rank_value regulator
#> 1   g1 0.46396396     FALSE
#> 2  g18 0.22091630      TRUE
#> 3   g2 0.09045695      TRUE
#> 4  g17 0.03444956      TRUE
#> 5   g9 0.01902045      TRUE
#> 6  g16 0.01740632      TRUE
```
