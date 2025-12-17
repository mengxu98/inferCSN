# **inferCSN**

## **Introduction**

[inferCSN](https://mengxu98.github.io/inferCSN/) is an R package for
**infer**ring **C**ell-**S**pecific gene regulatory **N**etwork from
single-cell RNA data.

![inferCSN+ workflow diagram](reference/figures/inferCSN.svg)

## **Installation**

You can install the released version from
[CRAN](https://github.com/cran) use:

``` r
install.packages("inferCSN")
# or
if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("inferCSN")
```

You can install the development version from
[GitHub](https://github.com/mengxu98/inferCSN) use
[pak](https://github.com/r-lib/pak):

``` r
if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/inferCSN")
```

## **Usage**

### **Examples**

``` r
library(inferCSN)
data(example_matrix)

network <- inferCSN(
  example_matrix
)
```

More functions and usages about
[inferCSN](https://mengxu98.github.io/inferCSN/)? Please reference
[here](https://mengxu98.github.io/inferCSN/reference/index.html).

## **Citation**

If you use [inferCSN](https://github.com/mengxu98/inferCSN) in your
work, please cite it reference
[here](https://mengxu98.github.io/inferCSN/authors.html#citation).
