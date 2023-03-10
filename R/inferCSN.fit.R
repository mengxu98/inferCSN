#' inferCSN.fit
#'
#' @param X The rows are samples and the columns are genes of the matrix
#' @param y Target vector
#' @param penalty penalty = penalty
#' @param crossValidation Cross validation
#' @param nFolds  N folds cross validation
#' @param maxSuppSize maxSuppSize = maxSuppSize
#' @param verbose Print detailed information
#' @param nGamma nGamma = nGamma
#'
#' @return A vector
#' @export
#'
inferCSN.fit <- function(X, y,
                         crossValidation = crossValidation,
                         penalty = penalty,
                         maxSuppSize = maxSuppSize,
                         nFolds = nFolds,
                         nGamma = nGamma,
                         verbose = verbose) {
  if (crossValidation) {
    tryCatch(
      {
        if (verbose) message("---------- Using cross validation ----------")
        fit <- L0Learn::L0Learn.cvfit(X, y,
                                      penalty = penalty,
                                      maxSuppSize = maxSuppSize,
                                      nFolds = nFolds,
                                      nGamma = nGamma
        )
        fit_inf <- print(fit)
        optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
        gamma <- fit$fit$gamma[optimalGammaIndex]
        lambda_list <- fit_inf[which(fit_inf$gamma == gamma), ]
        if (is.null(maxSuppSize)) {
          lambda <- min(lambda_list$lambda)
        } else {
          if (maxSuppSize %in% lambda_list$maxSuppSize) {
            lambda <- lambda_list$maxSuppSize[which(lambda_list$maxSuppSize == maxSuppSize)]
          } else {
            lambda <- min(lambda_list$lambda)
          }
        }
        temp <- coef(fit, lambda = lambda, gamma = gamma)
      },
      error = function(e) {
        if (verbose) message("---------- Cross validation error, used fit instead ----------")
        fit <- L0Learn::L0Learn.fit(X, y,
                                    penalty = penalty,
                                    maxSuppSize = maxSuppSize,
                                    nGamma = nGamma
        )
        fit_inf <- print(fit)
        fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
        lambda <- fit_inf$lambda[1]
        gamma <- fit_inf$gamma[1]
        temp <- coef(fit, lambda = lambda, gamma = gamma)
      }
    )
  } else {
    fit <- L0Learn::L0Learn.fit(X, y,
                                penalty = penalty,
                                maxSuppSize = maxSuppSize,
                                nGamma = nGamma
    )
    fit_inf <- print(fit)
    fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
    lambda <- fit_inf$lambda[1]
    gamma <- fit_inf$gamma[1]
  }
  temp <- coef(fit, lambda = lambda, gamma = gamma)
}
