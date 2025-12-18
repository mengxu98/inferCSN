#' @title Sparse regression model
#'
#' @md
#' @inheritParams inferCSN
#' @inheritParams single_network
#' @param x The matrix of regulators.
#' @param y The vector of target.
#' @param regulators_num The number of regulators for target.
#'
#' @return A list of the sparse regression model.
#' The list has the three components: model, metrics, and coefficients.
#'
#' @export
#' @examples
#' data(example_matrix)
#' fit_srm(
#'   x = example_matrix[, -1],
#'   y = example_matrix[, 1]
#' )
fit_srm <- function(
    x, y,
    cross_validation = FALSE,
    seed = 1,
    penalty = "L0",
    regulators_num = ncol(x),
    n_folds = 5,
    verbose = TRUE,
    ...) {
  if (cross_validation) {
    fit <- try(
      L0Learn::L0Learn.cvfit(
        x, y,
        penalty = penalty,
        maxSuppSize = regulators_num,
        nFolds = n_folds,
        seed = seed,
        ...
      )
    )

    if (any(class(fit) == "try-error")) {
      thisutils::log_message(
        "Cross-validation error, setting {.arg cross_validation} to {.pkg FALSE} and re-train",
        message_type = "warning",
        verbose = verbose
      )
      fit <- try(
        L0Learn::L0Learn.fit(
          x, y,
          penalty = penalty,
          maxSuppSize = regulators_num,
          ...
        )
      )
      if (any(class(fit) == "try-error")) {
        return(
          list(
            model = fit,
            metrics = list(r_squared = 0),
            coefficients = list(
              variable = colnames(x),
              coefficient = rep(0, ncol(x))
            )
          )
        )
      }
      fit_inf <- print(fit)
      lambda <- fit_inf$lambda[which.max(fit_inf$suppSize)]
      gamma <- fit_inf$gamma[which.max(fit_inf$suppSize)]
    } else {
      best_gamma_idx <- which.min(sapply(fit$cvMeans, min))
      gamma <- fit$fit$gamma[best_gamma_idx]
      cv_means <- fit$cvMeans[[best_gamma_idx]]
      best_lambda_idx <- which.min(cv_means)
      lambda <- fit$fit$lambda[[best_gamma_idx]][best_lambda_idx]
    }
  } else {
    fit <- try(
      L0Learn::L0Learn.fit(
        x, y,
        penalty = penalty,
        maxSuppSize = regulators_num,
        ...
      )
    )
    if (any(class(fit) == "try-error")) {
      return(
        list(
          model = fit,
          metrics = list(
            r_squared = 0
          ),
          coefficients = list(
            variable = colnames(x),
            coefficient = rep(0, ncol(x))
          )
        )
      )
    }
    fit_inf <- print(fit)
    lambda <- fit_inf$lambda[which.max(fit_inf$suppSize)]
    gamma <- fit_inf$gamma[which.max(fit_inf$suppSize)]
  }

  pred_y <- as.numeric(
    predict(
      fit,
      newx = x,
      lambda = lambda,
      gamma = gamma
    )
  )

  return(
    list(
      model = fit,
      metrics = list(
        r_squared = thisutils::r_square(y, pred_y)
      ),
      coefficients = list(
        variable = colnames(x),
        coefficient = as.vector(
          coef(
            fit,
            lambda = lambda,
            gamma = gamma
          )
        )[-1]
      )
    )
  )
}
