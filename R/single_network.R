#' @title Construct network for single target gene
#'
#' @inheritParams inferCSN
#' @param matrix An expression matrix.
#' @param target The target gene.
#'
#' @return A data frame of the single target gene network.
#'  The data frame has three columns:
#'  \item{regulator}{The regulator genes.}
#'  \item{target}{The target gene.}
#'  \item{weight}{The weight of the regulator gene on the target gene.}
#' @export
#' @examples
#' data("example_matrix")
#' head(
#'   single_network(
#'     example_matrix,
#'     regulators = colnames(example_matrix),
#'     target = "g1"
#'   )
#' )
#' head(
#'   single_network(
#'     example_matrix,
#'     regulators = colnames(example_matrix),
#'     target = "g1",
#'     cross_validation = TRUE
#'   )
#' )
#'
#' single_network(
#'   example_matrix,
#'   regulators = c("g1", "g2", "g3"),
#'   target = "g1"
#' )
#' single_network(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   target = "g1"
#' )
single_network <- function(
    matrix,
    regulators,
    target,
    cross_validation = FALSE,
    seed = 1,
    penalty = "L0",
    r_squared_threshold = 0,
    n_folds = 5,
    verbose = TRUE,
    ...) {
  regulators <- setdiff(regulators, target)
  if (length(regulators) < 2) {
    thisutils::log_message(
      "less than 2 regulators found while modeling: ", target,
      message_type = "warning",
      verbose = verbose
    )
    return()
  }
  x <- matrix[, regulators]
  y <- matrix[, target]

  result <- fit_srm(
    x, y,
    cross_validation = cross_validation,
    seed = seed,
    penalty = penalty,
    n_folds = n_folds,
    verbose = verbose,
    ...
  )

  r_squared <- result$metrics$r_squared
  if (r_squared > r_squared_threshold) {
    coefficients <- result$coefficients$coefficient |>
      thisutils::normalization(
        method = "unit_vector",
        ...
      )

    if (length(coefficients) != ncol(x)) {
      coefficients <- rep(0, ncol(x))
    }
  } else {
    coefficients <- rep(0, ncol(x))
  }

  return(
    data.frame(
      regulator = regulators,
      target = target,
      weight = coefficients
    )
  )
}
