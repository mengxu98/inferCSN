# **inferCSN** <img src="man/figures/logo.svg" align="right" width="120"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/inferCSN)](https://github.com/cran/inferCSN) [![R-CMD-check](https://github.com/mengxu98/inferCSN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mengxu98/inferCSN/actions/workflows/R-CMD-check.yaml) [![test-coverage](https://github.com/mengxu98/inferCSN/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/mengxu98/inferCSN/actions/workflows/test-coverage.yaml) [![pkgdown](https://github.com/mengxu98/inferCSN/actions/workflows/pkgdown.yaml/badge.svg)](https://mengxu98.github.io/inferCSN/reference/index.html) [![RStudio CRAN mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/inferCSN)](https://CRAN.R-project.org/package=inferCSN)

<!-- badges: end -->

## **Introduction**

[inferCSN](https://mengxu98.github.io/inferCSN/) is an R package for **infer**ring **C**ell-**S**pecific gene regulatory **N**etwork from single-cell RNA data.

<img src="man/figures/inferCSN.svg" alt="inferCSN+ workflow diagram" width="75%"/>

## **Installation**

You can install the released version from [CRAN](https://github.com/cran) use:

``` r
install.packages("inferCSN")
# or
if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("inferCSN")
```

You can install the development version from [GitHub](https://github.com/mengxu98/inferCSN) use [pak](https://github.com/r-lib/pak):

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

More functions and usages about [inferCSN](https://mengxu98.github.io/inferCSN/)? Please reference [here](https://mengxu98.github.io/inferCSN/reference/index.html).

## **Citation**

If you use [inferCSN](https://github.com/mengxu98/inferCSN) in your work, please cite it reference [here](https://mengxu98.github.io/inferCSN/authors.html#citation).
