#' @title Construct network for single target gene
#'
#' @inheritParams inferCSN
#' @param matrix An expression matrix.
#' @param regulators Candidate regulator genes.
#' @param target The target gene.
#'
#' @param pseudotime Optional pseudotime vector or branch matrix passed to
#' [inferCSN()].
#' @param max_support_size Optional support-size cap passed to [inferCSN()].
#' @param lag_fraction Fractional state lag passed to [inferCSN()].
#' @param lag_steps Optional integer state lag passed to [inferCSN()].
#' @param cores Number of inference workers.
#'
#' @return A data frame containing only selected edges for the requested target.
#' The data frame has three columns: regulator, target, and weight.
#'
#' @export
#' @examples
#' data(example_matrix)
#' head(
#'   single_network(
#'     example_matrix,
#'     regulators = colnames(example_matrix),
#'     target = "g1"
#'   )
#' )
#' single_network(
#'   example_matrix,
#'   regulators = c("g1", "g2", "g3"),
#'   target = "g1"
#' )
single_network <- function(
  matrix,
  regulators,
  target,
  pseudotime = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE
) {
  regulators <- setdiff(regulators, target)
  if (length(regulators) < 1) {
    thisutils::log_message(
      "No candidate regulators found when modeling: {.val {target}}",
      message_type = "warning",
      verbose = verbose
    )
    return(data.frame(
      regulator = character(),
      target = character(),
      weight = numeric()
    ))
  }
  inferCSN(
    object = matrix,
    pseudotime = pseudotime,
    regulators = regulators,
    targets = target,
    max_support_size = max_support_size,
    lag_fraction = lag_fraction,
    lag_steps = lag_steps,
    cores = cores,
    verbose = verbose
  )
}
