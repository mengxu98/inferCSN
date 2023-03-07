#' inferCSN.fit
#'
#' @param X The rows are samples and the columns are genes of the matrix
#' @param Y The vector
#' @param penalty penalty = penalty
#' @param crossValidation Cross validation
#' @param nFolds nFolds = 10
#' @param seed seed = 1
#' @param maxSuppSize maxSuppSize = maxSuppSize
#' @param nGamma nGamma = 5
#' @param gammaMin gammaMin = 0.0001
#' @param gammaMax gammaMax = 10
#'
#' @return A vector
#' @export
#'
inferCSN.fit <- function(X, Y,
                         crossValidation = FALSE,
                         penalty = penalty,
                         nFolds = 10,
                         seed = 1,
                         maxSuppSize = maxSuppSize,
                         nGamma = 5,
                         gammaMin = 0.0001,
                         gammaMax = 10) {
  if (crossValidation) {
    tryCatch(
      {
        message("---------- Using cross validation ----------")
        fit <- L0Learn::L0Learn.cvfit(X, Y,
                                      penalty = penalty,
                                      maxSuppSize = maxSuppSize,
                                      nFolds = 10,
                                      seed = 1,
                                      nGamma = 5,
                                      gammaMin = 0.0001,
                                      gammaMax = 10
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
        message("---------- Cross validation error, used fit instead ----------")
        fit <- L0Learn::L0Learn.fit(X, Y,
                                    penalty = penalty,
                                    maxSuppSize = maxSuppSize,
                                    nGamma = 5,
                                    gammaMin = 0.0001,
                                    gammaMax = 10
        )
        fit_inf <- print(fit)
        fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
        lambda <- fit_inf$lambda[1]
        gamma <- fit_inf$gamma[1]
        temp <- coef(fit, lambda = lambda, gamma = gamma)
      }
    )
  } else {
    fit <- L0Learn::L0Learn.fit(X, Y,
                                penalty = penalty,
                                maxSuppSize = maxSuppSize,
                                nGamma = 5,
                                gammaMin = 0.0001,
                                gammaMax = 10
    )
    fit_inf <- print(fit)
    fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
    lambda <- fit_inf$lambda[1]
    gamma <- fit_inf$gamma[1]
  }
  temp <- coef(fit, lambda = lambda, gamma = gamma)
}
