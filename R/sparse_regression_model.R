#' @title Sparse regression model
#'
#' @md
#' @inheritParams inferCSN
#' @inheritParams single_network
#' @param x The matrix of regulators.
#' @param y The vector of target.
#' @param regulators_num The number of non-zore coefficients, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#'
#' @return A list of the sparse regression model.
#'  The list has the following components:
#'  \item{model}{The sparse regression model.}
#'  \item{metrics}{A list of metrics.}
#'  \item{coefficients}{A list of coefficients.}
#'
#' @export
#' @examples
#' data("example_matrix")
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
      sparse_regression(
        x, y,
        penalty = penalty,
        regulators_num = regulators_num,
        cross_validation = cross_validation,
        n_folds = n_folds,
        seed = seed,
        verbose = verbose,
        ...
      )
    )

    if (any(class(fit) == "try-error")) {
      log_message(
        "cross validation error, setting `cross_validation` to `FALSE` and re-train model.",
        message_type = "warning",
        verbose = verbose
      )
      fit <- try(
        sparse_regression(
          x, y,
          penalty = penalty,
          regulators_num = regulators_num,
          cross_validation = FALSE,
          verbose = verbose,
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
      gamma <- fit$fit$gamma[which(
        unlist(lapply(
          fit$cvMeans, min
        )) == min(unlist(lapply(fit$cvMeans, min)))
      )]
      lambda_list <- dplyr::filter(
        print(fit),
        gamma == gamma
      )
      if (regulators_num %in% lambda_list$suppSize) {
        lambda <- lambda_list$lambda[which(
          lambda_list$suppSize == regulators_num
        )]
      } else {
        lambda <- min(lambda_list$lambda)
      }
    }
  } else {
    fit <- try(
      sparse_regression(
        x, y,
        penalty = penalty,
        regulators_num = regulators_num,
        cross_validation = FALSE,
        verbose = verbose,
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
        r_squared = r_square(y, pred_y)
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

#' @title Fit a sparse regression model
#'
#' @description
#'  Computes the regularization path for the specified loss function and penalty function.
#'
#' @md
#' @inheritParams fit_srm
#' @param algorithm The type of algorithm used to minimize the objective function, default is *`CD`*.
#' Currently *`CD`* and *`CDPSI`* are supported.
#' The *`CDPSI`* algorithm may yield better results, but it also increases running time.
#' @param loss The loss function.
#' @param nLambda The number of Lambda values to select.
#' @param nGamma The number of Gamma values to select.
#' @param gammaMax The maximum value of Gamma when using the `L0L2` penalty.
#' For the `L0L1` penalty this is automatically selected.
#' @param gammaMin The minimum value of Gamma when using the `L0L2` penalty.
#' For the `L0L1` penalty, the minimum value of gamma in the grid is set to gammaMin * gammaMax.
#' Note that this should be a strictly positive quantity.
#' @param partialSort If `TRUE`, partial sorting will be used for sorting the coordinates to do greedy cycling.
#'  Otherwise, full sorting is used.
#' @param maxIters The maximum number of iterations (full cycles) for `CD` per grid point.
#' @param rtol The relative tolerance which decides when to terminate optimization,
#' based on the relative change in the objective between iterations.
#' @param atol The absolute tolerance which decides when to terminate optimization,
#' based on the absolute L2 norm of the residuals.
#' @param activeSet If `TRUE`, performs active set updates.
#' @param activeSetNum The number of consecutive times a support should appear before declaring support stabilization.
#' @param maxSwaps The maximum number of swaps used by `CDPSI` for each grid point.
#' @param scaleDownFactor This parameter decides how close the selected Lambda values are.
#' @param screenSize The number of coordinates to cycle over when performing initial correlation screening.
#' @param autoLambda Ignored parameter. Kept for backwards compatibility.
#' @param lambdaGrid A grid of Lambda values to use in computing the regularization path.
#' @param excludeFirstK This parameter takes non-negative integers.
#' @param intercept If `FALSE`, no intercept term is included in the model.
#' @param lows Lower bounds for coefficients.
#' @param highs Upper bounds for coefficients.
#'
#' @references
#'  Hazimeh, Hussein et al.
#'  “L0Learn: A Scalable Package for Sparse Learning using L0 Regularization.”
#'  J. Mach. Learn. Res. 24 (2022): 205:1-205:8.
#'
#'  Hazimeh, Hussein and Rahul Mazumder.
#'  “Fast Best Subset Selection: Coordinate Descent and Local Combinatorial Optimization Algorithms.”
#'  Oper. Res. 68 (2018): 1517-1537.
#'
#'  https://github.com/hazimehh/L0Learn/blob/master/R/fit.R
#'
#' @return An S3 object describing the regularization path
#' @export
#' @examples
#' data("example_matrix")
#' fit <- sparse_regression(
#'   example_matrix[, -1],
#'   example_matrix[, 1]
#' )
#' head(coef(fit))
sparse_regression <- function(
    x, y,
    penalty = "L0",
    algorithm = c("CD", "CDPSI"),
    regulators_num = ncol(x),
    cross_validation = FALSE,
    n_folds = 5,
    seed = 1,
    loss = "SquaredError",
    nLambda = 100,
    nGamma = 5,
    gammaMax = 10,
    gammaMin = 0.0001,
    partialSort = TRUE,
    maxIters = 200,
    rtol = 1e-6,
    atol = 1e-9,
    activeSet = TRUE,
    activeSetNum = 3,
    maxSwaps = 100,
    scaleDownFactor = 0.8,
    screenSize = 1000,
    autoLambda = NULL,
    lambdaGrid = list(),
    excludeFirstK = 0,
    intercept = TRUE,
    lows = -Inf,
    highs = Inf,
    verbose = TRUE,
    ...) {
  algorithm <- match.arg(algorithm)

  if ((rtol < 0) || (rtol >= 1)) {
    stop("The specified rtol parameter must exist in [0, 1).")
  }
  if (atol < 0) {
    stop("The specified atol parameter must exist in [0, INF).")
  }
  if (!(loss %in% c("SquaredError", "Logistic", "SquaredHinge"))) {
    stop("The specified loss function is not supported.")
  }

  # Check binary classification for logistic and squared hinge loss
  if (loss == "Logistic" || loss == "SquaredHinge") {
    if (dim(table(y)) != 2) {
      stop("Only binary classification is supported. Make sure y has only 2 unique values.")
    }
    y <- factor(y, labels = c(-1, 1))
    y <- as.numeric(levels(y))[y]

    if (penalty == "L0") {
      if ((length(lambdaGrid) != 0) && (length(lambdaGrid) != 1)) {
        stop(
          "L0 Penalty requires 'lambdaGrid' to be a list of length 1.
      Where lambdaGrid[[1]] is a list or vector of decreasing positive values."
        )
      }
      penalty <- "L0L2"
      nGamma <- 1
      gammaMax <- 1e-7
      gammaMin <- 1e-7
    }
  }

  # Handle Lambda Grids
  if (length(lambdaGrid) != 0) {
    if (!is.null(autoLambda) && !autoLambda) {
      log_message(
        "'autoLambda' is ignored and inferred if 'lambdaGrid' is supplied.",
        message_type = "warning",
        verbose = verbose
      )
    }
    autoLambda <- FALSE
  } else {
    autoLambda <- TRUE
    lambdaGrid <- list(0)
  }

  # Check lambda grid for L0 penalty
  if (penalty == "L0" && !autoLambda) {
    bad_lambda_grid <- FALSE
    if (length(lambdaGrid) != 1) bad_lambda_grid <- TRUE
    current <- Inf
    for (nxt in lambdaGrid[[1]]) {
      if (nxt > current) {
        bad_lambda_grid <- TRUE
        break
      }
      if (nxt < 0) {
        bad_lambda_grid <- TRUE
        break
      }
      current <- nxt
    }
    if (bad_lambda_grid) {
      stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
       Where 'lambdaGrid[[1]]' is a list or vector of decreasing positive values.")
    }
  }

  # Check lambda grid for L0L2 penalties
  if (penalty != "L0" && !autoLambda) {
    bad_lambda_grid <- FALSE
    if (length(lambdaGrid) != nGamma) {
      log_message(
        "'nGamma' is ignored and replaced with length(lambdaGrid).",
        message_type = "warning",
        verbose = verbose
      )
      nGamma <- length(lambdaGrid)
    }
    for (i in seq_along(lambdaGrid)) {
      current <- Inf
      for (nxt in lambdaGrid[[i]]) {
        if (nxt > current) {
          bad_lambda_grid <- TRUE
          break
        }
        if (nxt < 0) {
          bad_lambda_grid <- TRUE
          break
        }
        current <- nxt
      }
      if (bad_lambda_grid) break
    }
    if (bad_lambda_grid) {
      stop("L0L2 Penalty requires 'lambdaGrid' to be a list of length 'nGamma'.
       Where 'lambdaGrid[[i]]' is a list or vector of decreasing positive values.")
    }
  }

  p <- dim(x)[[2]]
  withBounds <- FALSE

  # Check if bounds are specified
  if ((!identical(lows, -Inf)) || (!identical(highs, Inf))) {
    withBounds <- TRUE

    # Check bounds for CDPSI algorithm
    if (algorithm == "CDPSI") {
      if (any(lows != -Inf) || any(highs != Inf)) {
        stop("Bounds are not YET supported for CDPSI algorithm.")
      }
    }

    # Adjust lows and highs to vectors if scalar is provided
    if (.is_scalar(lows)) {
      lows <- lows * rep(1, p)
    } else if (!all(sapply(lows, .is_scalar)) || length(lows) != p) {
      stop("Lows must be a vector of real values of length p.")
    }

    if (.is_scalar(highs)) {
      highs <- highs * rep(1, p)
    } else if (!all(sapply(highs, .is_scalar)) || length(highs) != p) {
      stop("Highs must be a vector of real values of length p.")
    }

    # Check bounds conditions
    if (any(lows >= highs) || any(lows > 0) || any(highs < 0)) {
      stop("Bounds must conform to the following conditions: Lows <= 0, Highs >= 0, Lows < Highs.")
    }
  }

  m <- list()
  if (!cross_validation) {
    if (methods::is(x, "sparseMatrix")) {
      m <- srm_model_sparse(
        x, y, loss, penalty, algorithm, regulators_num,
        nLambda, nGamma, gammaMax, gammaMin, partialSort,
        maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
        scaleDownFactor, screenSize, !autoLambda, lambdaGrid,
        excludeFirstK, intercept, withBounds, lows, highs
      )
    } else {
      m <- srm_model_dense(
        x, y, loss, penalty, algorithm, regulators_num,
        nLambda, nGamma, gammaMax, gammaMin, partialSort,
        maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
        scaleDownFactor, screenSize, !autoLambda, lambdaGrid,
        excludeFirstK, intercept, withBounds, lows, highs
      )
    }
  } else {
    set.seed(seed)
    if (methods::is(x, "sparseMatrix")) {
      m <- srm_model_cv_sparse(
        x, y, loss, penalty, algorithm, regulators_num,
        nLambda, nGamma, gammaMax, gammaMin, partialSort,
        maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
        scaleDownFactor, screenSize, !autoLambda, lambdaGrid, n_folds,
        seed, excludeFirstK, intercept, withBounds, lows, highs
      )
    } else {
      m <- srm_model_cv_dense(
        x, y, loss, penalty, algorithm, regulators_num,
        nLambda, nGamma, gammaMax, gammaMin, partialSort,
        maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
        scaleDownFactor, screenSize, !autoLambda, lambdaGrid, n_folds,
        seed, excludeFirstK, intercept, withBounds, lows, highs
      )
    }
  }

  settings <- list()
  settings[[1]] <- intercept
  names(settings) <- c("intercept")

  for (i in seq_along(m$SuppSize)) {
    last <- length(m$SuppSize[[i]])
    if (m$SuppSize[[i]][last] > regulators_num) {
      if (last == 1) {
        log_message(
          "only 1 element in path with support size > regulators_num.
            Try increasing 'regulators_num' to resolve the issue.",
          message_type = "warning",
          verbose = verbose
        )
      } else {
        m$SuppSize[[i]] <- m$SuppSize[[i]][-last]
        m$Converged[[i]] <- m$Converged[[i]][-last]
        m$lambda[[i]] <- m$lambda[[i]][-last]
        m$a0[[i]] <- m$a0[[i]][-last]
        m$beta[[i]] <- methods::as(m$beta[[i]][, -last], "sparseMatrix")
        if (!cross_validation) {
          m$CVMeans[[i]] <- m$CVMeans[[i]][-last]
          m$CVSDs[[i]] <- m$CVSDs[[i]][-last]
        }
      }
    }
  }

  fit <- list(
    beta = m$beta,
    lambda = lapply(m$lambda, signif, digits = 6),
    a0 = m$a0,
    converged = m$Converged,
    suppSize = m$SuppSize,
    gamma = m$gamma,
    penalty = penalty,
    loss = loss,
    settings = settings
  )

  if (is.null(colnames(x))) {
    varnames <- seq_len(dim(x)[2])
  } else {
    varnames <- colnames(x)
  }
  fit$varnames <- varnames
  class(fit) <- "srm"
  fit$n <- dim(x)[1]
  fit$p <- dim(x)[2]

  if (!cross_validation) {
    g <- fit
  } else {
    g <- list(fit = fit, cvMeans = m$CVMeans, cvSDs = m$CVSDs)
    class(g) <- "srm_cv"
  }

  return(g)
}
