# Filter and sort matrix

Filter and sort matrix

## Usage

``` r
filter_sort_matrix(network_matrix, regulators = NULL, targets = NULL)
```

## Arguments

- network_matrix:

  The matrix of network weight.

- regulators:

  Regulators list.

- targets:

  Targets list.

## Value

Filtered and sorted matrix

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2025-12-17 14:37:08] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:08] Checking parameters...
#> ℹ [2025-12-17 14:37:08] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:08] Using 1 core
#> ℹ [2025-12-17 14:37:08] Building results
#> ✔ [2025-12-17 14:37:08] Building network done
colnames(network_table) <- c("row", "col", "value")
network_matrix <- thisutils::table_to_matrix(network_table)
filter_sort_matrix(network_matrix)[1:6, 1:6]
#>             g1          g2            g3          g4         g5           g6
#> g1  0.00000000  0.66844490  0.0180501564 -0.03658957 0.04489503 -0.030280993
#> g2  0.35630356  0.00000000  0.5491722960  0.02103272 0.01366921 -0.007459899
#> g3  0.01235281  0.70713027  0.0000000000  0.68526346 0.12187906 -0.001181702
#> g4 -0.02896993  0.03201306  0.8103229841  0.00000000 0.65615667  0.124827620
#> g5  0.03866820  0.02176632  0.1539425396  0.70453814 0.00000000  0.655663735
#> g6 -0.02979658 -0.01242983 -0.0008958639  0.14958079 0.72762920  0.000000000

filter_sort_matrix(
  network_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g4")
)
#>            g3          g4
#> g1 0.01805016 -0.03658957
#> g2 0.54917230  0.02103272
```
