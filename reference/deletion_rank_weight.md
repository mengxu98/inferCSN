# Rank signed edges by deletion evidence

Converts signed coefficients and their deletion delta-BIC into the
signed ordinal weights returned by \[inferCSN()\]. Edges are ranked by
decreasing deletion evidence; evidence within \`1e-12 \* (1 +
abs(maximum))\` of a group's maximum shares one rank, and the sign of
the coefficient is kept.

## Usage

``` r
deletion_rank_weight(coefficient, deletion_delta_bic)
```

## Arguments

- coefficient:

  Signed coefficients.

- deletion_delta_bic:

  Non-negative deletion delta-BIC, one per coefficient.

## Value

Signed weights in (-1, 1); zero, non-finite or negative-evidence entries
receive 0.

## Examples

``` r
deletion_rank_weight(c(0.5, -1.2, 0.3, 0), c(4, 9, 4, 2))
#> [1]  0.25 -0.75  0.25  0.00
```
