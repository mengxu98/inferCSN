#' @title sparse.regression
#'
#' @param X The data matrix
#' @param y The response vector
#' @inheritParams inferCSN
#'
#' @importFrom stats "coef"
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
      fit_inf <- print(fit)
      optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
      gamma <- fit$fit$gamma[optimalGammaIndex]
      lambda_list <- fit_inf[which(fit_inf$gamma == gamma),]
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
      fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
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
    fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
    lambda <- fit_inf$lambda[1]
    gamma <- fit_inf$gamma[1]
  }
  temp <- coef(fit, lambda = lambda, gamma = gamma)
}
