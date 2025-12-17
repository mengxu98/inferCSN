# Weight sift

Remove edges with smaller weights in the reverse direction.

## Usage

``` r
weight_sift(table)
```

## Arguments

- table:

  A data frame with three columns: "regulator", "target", and "weight".

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2025-12-17 14:37:36] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:36] Checking parameters...
#> ℹ [2025-12-17 14:37:36] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:36] Using 1 core
#> ⠙ [2025-12-17 14:37:36] Running [1/18] Processing: g1  ETA:  0s
#> ✔ [2025-12-17 14:37:36] Completed 18 tasks in 174ms
#> 
#> ℹ [2025-12-17 14:37:36] Building results
#> ✔ [2025-12-17 14:37:36] Building network done
weight_sift(network_table) |> head()
#>   regulator target     weight
#> 1       g18     g1 -0.9223177
#> 2       g17    g18  0.8770468
#> 3        g4     g3  0.8103230
#> 4       g16    g15  0.7659245
#> 5       g17    g16  0.7558764
#> 6       g12    g11  0.7444053
```
