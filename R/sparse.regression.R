globalVariables(c("."))

#' @title sparse.regression
#'
#' @param X The data matrix
#' @param y The response vector
#' @inheritParams inferCSN
#'
#' @importFrom stats coef
#'
#' @return A vector of weight
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
      # fit_inf <- print(fit)
      # lambda_list <- fit_inf[which(fit_inf$gamma == gamma),]
      lambda_list <- print(fit) %>% dplyr::filter(gamma == gamma, )
      if (is.null(maxSuppSize)) {
        lambda <- min(lambda_list$lambda)
      } else {
        if (maxSuppSize %in% lambda_list$maxSuppSize) {
          lambda <- lambda_list$maxSuppSize[which(lambda_list$maxSuppSize == maxSuppSize)]
        } else {
          lambda <- min(lambda_list$lambda)
        }
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
      fit_inf <- print(fit)
      lambda <- fit_inf$lambda[fit_inf$suppSize %>% which.max()]
      gamma <- fit_inf$gamma[fit_inf$suppSize %>% which.max()]
    }
    )
  } else {
    fit <- inferCSN.fit(
      X, y,
      penalty = penalty,
      algorithm = algorithm,
      maxSuppSize = maxSuppSize
    )
    fit_inf <- print(fit)
    lambda <- fit_inf$lambda[fit_inf$suppSize %>% which.max()]
    gamma <- fit_inf$gamma[fit_inf$suppSize %>% which.max()]
  }
  coefficients <- coef(fit, lambda = lambda, gamma = gamma) %>% as.vector() %>% .[-1]
  return(coefficients)
}

#' sub.inferCSN
#'
#' @param regulatorsMatrix regulatorsMatrix
#' @param targetsMatrix targetsMatrix
#' @param target target
#' @inheritParams inferCSN
#'
#' @return A data table
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

  coefficients <- sparse.regression(
    X, y,
    crossValidation = crossValidation,
    penalty = penalty,
    algorithm = algorithm,
    maxSuppSize = maxSuppSize,
    nFolds = nFolds,
    verbose = verbose
  ) %>% abs()
  coefficients <- coefficients / sum(coefficients)
  if (length(coefficients) != ncol(X)) coefficients <- 0
  weightd <- data.frame(regulator = colnames(X), target = target, weight = coefficients)
  return(weightd)
}
