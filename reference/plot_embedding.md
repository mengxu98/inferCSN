# Plot an embedding

Plot an embedding

## Usage

``` r
plot_embedding(
  matrix,
  labels = NULL,
  method = "pca",
  colors = RColorBrewer::brewer.pal(length(unique(labels)), "Set1"),
  seed = 1,
  point_size = 1,
  cores = 1
)
```

## Arguments

- matrix:

  Input matrix.

- labels:

  Point labels.

- method:

  Embedding method.

- colors:

  Point colors.

- seed:

  Random seed.

- point_size:

  Point size.

- cores:

  Number of threads.

## Value

An embedding plot.
