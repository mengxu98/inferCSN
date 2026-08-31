# Subsample an expression matrix

Subsample an expression matrix

## Usage

``` r
subsampling(
  matrix,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  seed = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- matrix:

  Input matrix.

- subsampling_method:

  Subsampling strategy.

- subsampling_ratio:

  Fraction retained or aggregated.

- seed:

  Random seed.

- verbose:

  Whether to report progress.

- ...:

  Additional arguments.

## Value

The subsampled matrix.
