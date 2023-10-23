#' @title Sparse regression model
#'
#' @param X The data matrix
#' @param y The response vector
#'
#' @inheritParams inferCSN
#'
#' @return The coefficients
#' @export
#'
sparse.regression <- function(X, y,
                              crossValidation = FALSE,
                              seed = 1,
                              penalty = "L0",
                              algorithm = "CD",
                              maxSuppSize = NULL,
                              nFolds = 10,
                              kFolds = NULL,
                              rThreshold = 0,
                              verbose = FALSE) {
  if (!is.null(kFolds)) {
    samples <- sample(kFolds / 10 * nrow(X))
    testX <- X[-samples, ]
    X <- X[samples, ]
    testy <- y[-samples]
    y <- y[samples]
  }

  if (crossValidation) {
    fit <- try(inferCSN.fit(X, y,
                            crossValidation = crossValidation,
                            seed = seed,
                            penalty = penalty,
                            algorithm = algorithm,
                            maxSuppSize = maxSuppSize,
                            nFolds = nFolds))

    if (class(fit)[1] == "try-error") {
      if (verbose) message("Cross validation error, used fit instead.")
      fit <- inferCSN.fit(X, y,
                          penalty = penalty,
                          algorithm = algorithm,
                          maxSuppSize = maxSuppSize)

      fitInf <- print(fit)
      lambda <- fitInf$lambda[which.max(fitInf$suppSize)]
      gamma <- fitInf$gamma[which.max(fitInf$suppSize)]
    } else {
      gamma <- fit$fit$gamma[which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))]
      lambdaList <- dplyr::filter(print(fit), gamma == gamma)
      if (maxSuppSize %in% lambdaList$maxSuppSize) {
        lambda <- lambdaList$maxSuppSize[which(lambdaList$maxSuppSize == maxSuppSize)]
      } else {
        lambda <- min(lambdaList$lambda)
      }
    }
  } else {
    fit <- inferCSN.fit(X, y,
                        penalty = penalty,
                        algorithm = algorithm,
                        maxSuppSize = maxSuppSize)

    fitInf <- print(fit)
    lambda <- fitInf$lambda[which.max(fitInf$suppSize)]
    gamma <- fitInf$gamma[which.max(fitInf$suppSize)]
  }

  if (rThreshold == 0) {
    return(as.vector(coef(fit, lambda = lambda, gamma = gamma))[-1])
  } else {
    r <- 1
    if (!is.null(kFolds)) {
      predy <- as.numeric(predict(fit,
                                  newx = testX,
                                  lambda = lambda,
                                  gamma = gamma))
      if (length(testy) == length(predy)) {
        if (stats::var(testy) != 0 && stats::var(predy) != 0) {
          r <- stats::cor(testy, predy)
        }
      }
    }

    if (r >= rThreshold) {
      return(as.vector(coef(fit, lambda = lambda, gamma = gamma))[-1])
    } else {
      return(0.0001)
    }
  }
}

#' @title Sparse regression model for single gene
#'
#' @param regulatorsMatrix Regulators matrix
#' @param targetsMatrix Targets matrix
#' @param target Target genes
#'
#' @inheritParams inferCSN
#'
#' @return The weight data table of sub-network
#' @export
#'
sub.inferCSN <- function(regulatorsMatrix,
                         targetsMatrix,
                         target = NULL,
                         crossValidation = FALSE,
                         seed = 1,
                         penalty = "L0",
                         algorithm = "CD",
                         maxSuppSize = NULL,
                         nFolds = 10,
                         kFolds = NULL,
                         rThreshold = 0,
                         verbose = FALSE) {
  X <- regulatorsMatrix[, setdiff(colnames(regulatorsMatrix), target)]
  if (is(X, "sparseMatrix")) X <- as.matrix(X)
  y <- targetsMatrix[, target]

  if (is.null(maxSuppSize)) maxSuppSize <- ncol(X)

  coefficients <- sparse.regression(X, y,
                                    crossValidation = crossValidation,
                                    seed = seed,
                                    penalty = penalty,
                                    algorithm = algorithm,
                                    maxSuppSize = maxSuppSize,
                                    nFolds = nFolds,
                                    kFolds = kFolds,
                                    rThreshold = rThreshold,
                                    verbose = verbose)

  coefficients <- coefficients / sum(abs(coefficients))
  if (length(coefficients) != ncol(X)) coefficients <- 0.0001
  return(data.frame(regulator = colnames(X),
                    target = target,
                    weight = coefficients))
}

#' @title Fit a sparse regression model
#'
#' @description Computes the regularization path for the specified loss function and penalty function
#'
#' @param x The data matrix
#' @param y The response vector
#' @param loss The loss function
#' @param nLambda The number of Lambda values to select
#' @param nGamma The number of Gamma values to select
#' @param gammaMax The maximum value of Gamma when using the L0L2 penalty
#' @param gammaMin The minimum value of Gamma when using the L0L2 penalty
#' @param partialSort If TRUE, partial sorting will be used for sorting the coordinates to do greedy cycling. Otherwise, full sorting is used
#' @param maxIters The maximum number of iterations (full cycles) for CD per grid point
#' @param rtol The relative tolerance which decides when to terminate optimization (based on the relative change in the objective between iterations)
#' @param atol The absolute tolerance which decides when to terminate optimization (based on the absolute L2 norm of the residuals)
#' @param activeSet If TRUE, performs active set updates
#' @param activeSetNum The number of consecutive times a support should appear before declaring support stabilization
#' @param maxSwaps The maximum number of swaps used by CDPSI for each grid point
#' @param scaleDownFactor This parameter decides how close the selected Lambda values are
#' @param screenSize The number of coordinates to cycle over when performing initial correlation screening
#' @param autoLambda Ignored parameter. Kept for backwards compatibility
#' @param lambdaGrid A grid of Lambda values to use in computing the regularization path
#' @param excludeFirstK This parameter takes non-negative integers
#' @param intercept If FALSE, no intercept term is included in the model
#' @param lows Lower bounds for coefficients
#' @param highs Upper bounds for coefficients
#'
#' @inheritParams inferCSN
#'
#' @return An S3 object describing the regularization path
#' @export
#'
inferCSN.fit <- function(x, y,
                         penalty = "L0",
                         algorithm = "CD",
                         maxSuppSize = 100,
                         crossValidation = FALSE,
                         nFolds = 10,
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
                         highs = Inf) {
  # Check parameter values
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
  if (loss == "Logistic" | loss == "SquaredHinge") {
    if (dim(table(y)) != 2) {
      stop("Only binary classification is supported. Make sure y has only 2 unique values.")
    }
    y <- factor(y, labels = c(-1, 1)) # Returns a vector of strings
    y <- as.numeric(levels(y))[y]

    # Adjust parameters for L0 penalty
    if (penalty == "L0") {
      if ((length(lambdaGrid) != 0) && (length(lambdaGrid) != 1)) {
        stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
                Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
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
      warning("'autoLambda' is ignored and inferred if 'lambdaGrid' is supplied.")
    }
    autoLambda <- FALSE
  } else {
    autoLambda <- TRUE
    lambdaGrid <- list(0)
  }

  # Check lambda grid for L0 penalty
  if (penalty == "L0" && !autoLambda) {
    bad_lambdaGrid <- FALSE
    if (length(lambdaGrid) != 1) bad_lambdaGrid <- TRUE
    current <- Inf
    for (nxt in lambdaGrid[[1]]) {
      if (nxt > current) {
        bad_lambdaGrid <- TRUE
        break
      }
      if (nxt < 0) {
        bad_lambdaGrid <- TRUE
        break
      }
      current <- nxt
    }
    if (bad_lambdaGrid) {
      stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
           Where 'lambdaGrid[[1]]' is a list or vector of decreasing positive values.")
    }
  }

  # Check lambda grid for L0L2 penalties
  if (penalty != "L0" && !autoLambda) {
    bad_lambdaGrid <- FALSE
    if (length(lambdaGrid) != nGamma) {
      warning("'nGamma' is ignored and replaced with length(lambdaGrid).....")
      nGamma <- length(lambdaGrid)
    }
    for (i in 1:length(lambdaGrid)) {
      current <- Inf
      for (nxt in lambdaGrid[[i]]) {
        if (nxt > current) {
          bad_lambdaGrid <- TRUE
          break
        }
        if (nxt < 0) {
          bad_lambdaGrid <- TRUE
          break
        }
        current <- nxt
      }
      if (bad_lambdaGrid) break
    }
    if (bad_lambdaGrid) {
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
    if (is.scalar(lows)) {
      lows <- lows * rep(1, p)
    } else if (!all(sapply(lows, is.scalar)) || length(lows) != p) {
      stop("Lows must be a vector of real values of length p.")
    }

    if (is.scalar(highs)) {
      highs <- highs * rep(1, p)
    } else if (!all(sapply(highs, is.scalar)) || length(highs) != p) {
      stop("Highs must be a vector of real values of length p.")
    }

    # Check bounds conditions
    if (any(lows >= highs) || any(lows > 0) || any(highs < 0)) {
      stop("Bounds must conform to the following conditions: Lows <= 0, Highs >= 0, Lows < Highs.")
    }
  }

  # Call appropriate C++ function based on matrix type
  M <- list()
  if (!crossValidation) {
    if (is(x, "sparseMatrix")) {
      M <- .Call("_inferCSN_inferCSNFit_sparse",
                 PACKAGE = "inferCSN", x, y, loss, penalty,
                 algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
                 partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
                 scaleDownFactor, screenSize, !autoLambda, lambdaGrid,
                 excludeFirstK, intercept, withBounds, lows, highs)
    } else {
      M <- .Call("_inferCSN_inferCSNFit_dense",
                 PACKAGE = "inferCSN", x, y, loss, penalty,
                 algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
                 partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
                 scaleDownFactor, screenSize, !autoLambda, lambdaGrid,
                 excludeFirstK, intercept, withBounds, lows, highs)
    }
  } else {
    set.seed(seed)
    if (is(x, "sparseMatrix")) {
      M <- .Call("_inferCSN_inferCSNCV_sparse",
                 PACKAGE = "inferCSN", x, y, loss, penalty,
                 algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
                 partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
                 scaleDownFactor, screenSize, !autoLambda, lambdaGrid, nFolds,
                 seed, excludeFirstK, intercept, withBounds, lows, highs)
    } else {
      M <- .Call("_inferCSN_inferCSNCV_dense",
                 PACKAGE = "inferCSN", x, y, loss, penalty,
                 algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
                 partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
                 scaleDownFactor, screenSize, !autoLambda, lambdaGrid, nFolds,
                 seed, excludeFirstK, intercept, withBounds, lows, highs)
    }
  }

  settings <- list()
  settings[[1]] <- intercept
  names(settings) <- c("intercept")

  # Remove potential support sizes exceeding maxSuppSize
  for (i in 1:length(M$SuppSize)) {
    last <- length(M$SuppSize[[i]])
    if (M$SuppSize[[i]][last] > maxSuppSize) {
      if (last == 1) {
        warning("Warning! Only 1 element in path with support size > maxSuppSize.
                Try increasing maxSuppSize to resolve the issue.")
      } else {
        M$SuppSize[[i]] <- M$SuppSize[[i]][-last]
        M$Converged[[i]] <- M$Converged[[i]][-last]
        M$lambda[[i]] <- M$lambda[[i]][-last]
        M$a0[[i]] <- M$a0[[i]][-last]
        M$beta[[i]] <- as(M$beta[[i]][, -last], "sparseMatrix")
        if (!crossValidation) {
          M$CVMeans[[i]] <- M$CVMeans[[i]][-last]
          M$CVSDs[[i]] <- M$CVSDs[[i]][-last]
        }
      }
    }
  }

  fit <- list(beta = M$beta,
              lambda = lapply(M$lambda, signif, digits = 6),
              a0 = M$a0,
              converged = M$Converged,
              suppSize = M$SuppSize,
              gamma = M$gamma,
              penalty = penalty,
              loss = loss,
              settings = settings)

  if (is.null(colnames(x))) {
    varnames <- 1:dim(x)[2]
  } else {
    varnames <- colnames(x)
  }
  fit$varnames <- varnames
  class(fit) <- "inferCSN"
  fit$n <- dim(x)[1]
  fit$p <- dim(x)[2]

  if (!crossValidation) {
    G <- fit
  } else {
    G <- list(fit = fit, cvMeans = M$CVMeans, cvSDs = M$CVSDs)
    class(G) <- "inferCSNCV"
  }
  return(G)
}
