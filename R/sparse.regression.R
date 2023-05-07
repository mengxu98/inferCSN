#' @title sparse.regression
#'
#' @param X The data matrix
#' @param y The response vector
#' @inheritParams inferCSN
#'
#' @importFrom stats coef
#'
#' @return The coefficients
#' @export
#'
sparse.regression <- function(X, y,
                              crossValidation = FALSE,
                              penalty = "L0",
                              algorithm = "CD",
                              maxSuppSize = NULL,
                              nFolds = 10,
                              verbose = FALSE) {
  if (crossValidation) {
    tryCatch({
      fit <- inferCSN.cvfit(
        X, y,
        penalty = penalty,
        algorithm = algorithm,
        maxSuppSize = maxSuppSize,
        nFolds = nFolds
      )
      gamma <- fit$fit$gamma[which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))]
      lambdaList <- print(fit) %>% dplyr::filter(gamma == gamma, )
      if (maxSuppSize %in% lambdaList$maxSuppSize) {
        lambda <- lambdaList$maxSuppSize[which(lambdaList$maxSuppSize == maxSuppSize)]
      } else {
        lambda <- min(lambdaList$lambda)
      }
    },
    error = function(e) {
      if (verbose) message("Cross validation error, used fit instead......")
      fit <- inferCSN.fit(
        X, y,
        penalty = penalty,
        algorithm = algorithm,
        maxSuppSize = maxSuppSize
      )
      fitInf <- print(fit)
      lambda <- fitInf$lambda[fitInf$suppSize %>% which.max()]
      gamma <- fitInf$gamma[fitInf$suppSize %>% which.max()]
    }
    )
  } else {
    fit <- inferCSN.fit(
      X, y,
      penalty = penalty,
      algorithm = algorithm,
      maxSuppSize = maxSuppSize
    )
    fitInf <- print(fit)
    lambda <- fitInf$lambda[fitInf$suppSize %>% which.max()]
    gamma <- fitInf$gamma[fitInf$suppSize %>% which.max()]
  }
  return(coef(fit, lambda = lambda, gamma = gamma) %>% as.vector() %>% .[-1])
}

#' sub.inferCSN
#'
#' @param regulatorsMatrix regulatorsMatrix
#' @param targetsMatrix targetsMatrix
#' @param target target
#' @inheritParams inferCSN
#'
#' @return The weight data table of sub-network.
#' @export
#'
sub.inferCSN <- function(regulatorsMatrix = NULL,
                         targetsMatrix = NULL,
                         target = NULL,
                         crossValidation = FALSE,
                         penalty = "L0",
                         algorithm = "CD",
                         maxSuppSize = NULL,
                         nFolds = 10,
                         verbose = FALSE) {
  X <- as.matrix(regulatorsMatrix[, setdiff(colnames(regulatorsMatrix), target)])
  y <- targetsMatrix[, target]

  if (is.null(maxSuppSize)) maxSuppSize <- ncol(X)

  coefficients <- sparse.regression(
    X, y,
    crossValidation = crossValidation,
    penalty = penalty,
    algorithm = algorithm,
    maxSuppSize = maxSuppSize,
    nFolds = nFolds,
    verbose = verbose
  )
  coefficients <- coefficients / sum(abs(coefficients))
  if (length(coefficients) != ncol(X)) coefficients <- 0
  return(data.frame(regulator = colnames(X), target = target, weight = coefficients))
}
