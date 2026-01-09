# Format network table

Format network table

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

  The weight data table of network.

- regulators:

  Regulators list.

- targets:

  Targets list.

- abs_weight:

  Logical value, default is *`TRUE`*, whether to perform absolute value
  on weights, and when set `abs_weight` to *`TRUE`*, the output of
  weight table will create a new column named `Interaction`.

## Value

Formated network table

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2026-01-09 07:09:33] Inferring network for <dense matrix>...
#> ◌ [2026-01-09 07:09:33] Checking parameters...
#> ℹ [2026-01-09 07:09:33] Using L0 sparse regression model
#> ℹ [2026-01-09 07:09:33] Using 1 core
#> ℹ [2026-01-09 07:09:33] Building results
#> ✔ [2026-01-09 07:09:34] Inferring network done
#> ℹ [2026-01-09 07:09:34] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18

network_format(
  network_table,
  regulators = "g1"
)
#>    regulator target      weight Interaction
#> 1         g1     g2 0.668444900  Activation
#> 2         g1    g18 0.467602224  Repression
#> 3         g1    g17 0.047607627  Activation
#> 4         g1     g5 0.044895033  Activation
#> 5         g1     g9 0.042442058  Repression
#> 6         g1     g4 0.036589574  Repression
#> 7         g1     g6 0.030280993  Repression
#> 8         g1     g7 0.024014462  Activation
#> 9         g1    g16 0.018303456  Repression
#> 10        g1     g3 0.018050156  Activation
#> 11        g1    g15 0.015156339  Activation
#> 12        g1    g14 0.013970876  Repression
#> 13        g1    g12 0.011942835  Activation
#> 14        g1    g10 0.009393328  Activation
#> 15        g1     g8 0.009042402  Activation
#> 16        g1    g11 0.009009966  Activation
#> 17        g1    g13 0.001787438  Repression

network_format(
  network_table,
  regulators = "g1",
  abs_weight = FALSE
)
#>    regulator target       weight
#> 1         g1     g2  0.668444900
#> 2         g1    g18 -0.467602224
#> 3         g1    g17  0.047607627
#> 4         g1     g5  0.044895033
#> 5         g1     g9 -0.042442058
#> 6         g1     g4 -0.036589574
#> 7         g1     g6 -0.030280993
#> 8         g1     g7  0.024014462
#> 9         g1    g16 -0.018303456
#> 10        g1     g3  0.018050156
#> 11        g1    g15  0.015156339
#> 12        g1    g14 -0.013970876
#> 13        g1    g12  0.011942835
#> 14        g1    g10  0.009393328
#> 15        g1     g8  0.009042402
#> 16        g1    g11  0.009009966
#> 17        g1    g13 -0.001787438

network_format(
  network_table,
  targets = "g3"
)
#>    regulator target       weight Interaction
#> 1         g4     g3 0.8103229841  Activation
#> 2         g2     g3 0.5491722960  Activation
#> 3         g5     g3 0.1539425396  Activation
#> 4        g12     g3 0.0560346972  Repression
#> 5        g18     g3 0.0558371348  Activation
#> 6        g17     g3 0.0523816646  Repression
#> 7        g14     g3 0.0504639514  Activation
#> 8         g7     g3 0.0425107240  Activation
#> 9        g15     g3 0.0361683897  Repression
#> 10       g10     g3 0.0313224187  Activation
#> 11       g16     g3 0.0270465616  Activation
#> 12       g13     g3 0.0258582091  Activation
#> 13       g11     g3 0.0229433379  Repression
#> 14        g1     g3 0.0180501564  Activation
#> 15        g8     g3 0.0099732575  Activation
#> 16        g9     g3 0.0095382542  Repression
#> 17        g6     g3 0.0008958639  Repression

network_format(
  network_table,
  regulators = c("g1", "g3"),
  targets = c("g3", "g5")
)
#>   regulator target     weight Interaction
#> 1        g3     g5 0.12187906  Activation
#> 2        g1     g5 0.04489503  Activation
#> 3        g1     g3 0.01805016  Activation
```
