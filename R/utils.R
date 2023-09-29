utils::globalVariables(c("x",
                         "y",
                         "xend",
                         "yend",
                         "weight",
                         "Interaction",
                         "name",
                         "regulator",
                         "target",
                         "degree",
                         "edges",
                         "curvetype"))

#' @title Check input parameters
#'
#' @inheritParams inferCSN
#'
#' @return NULL
#' @export
#'
check.parameters <- function(matrix,
                             penalty,
                             algorithm,
                             crossValidation,
                             seed,
                             nFolds,
                             kFolds,
                             rThreshold,
                             regulators,
                             targets,
                             maxSuppSize,
                             verbose,
                             cores) {
  if (verbose) message.success("Checking input parameters.")

  matrixErrorMessage <- paste0("Parameter matrix must be a two-dimensional matrix,
                               where each column corresponds to a gene and each row corresponds to a sample/cell.")
  if (!is.matrix(matrix) && !is.array(matrix)) {
    cli::cli_abort(c(matrixErrorMessage,
                     x = "You've supplied a {.cls {class(matrix)}}."))
  }

  if (length(dim(matrix)) != 2) {
    cli::cli_abort(c(x = matrixErrorMessage))
  }

  if (is.null(colnames(matrix))) {
    cli::cli_abort(c(x = "Parameter matrix must contain the names of the genes as colnames."))
  }

  # Check the penalty term of the regression model
  if (!any(c("L0", "L0L2") == penalty)) {
    cli::cli_abort(c(x = "inferCSN does not support '{penalty}' penalty regression.\n",
                     i = "Please set penalty item as 'L0' or 'L0L2'."))
  }

  # Check the algorithm of the regression model
  if (!any(c("CD", "CDPSI") == algorithm)) {
    cli::cli_abort(c(x = "inferCSN does not support '{algorithm}' algorithmn.\n",
                     i = "Please set algorithm item as 'CD' or 'CDPSI'."))
  }

  if (!is.numeric(seed)) {
    seed <- 1
    message.warning("Supplied seed is not a valid integer, initialize 'seed' to 1.")
  }

  if (!is.null(kFolds)) {
    if (!(kFolds > 0 && kFolds < 10)) {
      cli::cli_abort(c(x = "Please set 'kFolds' value between: (0, 10)."))
    }
  }

  if (!is.null(targets)) {
    if (!is.vector(targets)) {
      cli::cli_abort(c(x = "Parameter 'targets' must a vector (of indices or gene names)."))
    }

    if (is.numeric(targets)) {
      if(max(targets) > nrow(matrix)) cli::cli_abort(c(x = "At least one index in 'targets' exceeds the number of genes."))
      if(min(targets) <= 0) cli::cli_abort(c(x = "The indexes in 'targets' should be >=1."))
    }

    if(any(table(targets) > 1)) cli::cli_abort(c(x = "Please provide each target only once."))

    if(is.character(targets)){
      targetsInMatrix <- intersect(targets, colnames(matrix))
      if(length(targetsInMatrix) == 0) {
        cli::cli_abort(c(x = "The genes must contain at least one target."))
      }

      if(length(targetsInMatrix) < length(targets)) {
        message.warning("Only {length(targetsInMatrix)} out of {length(targets)} target genes are in the expression matrix.")
      }
    }
  }

  if (!is.null(regulators)) {
    if(is.list(regulators)) {
      if(!all(names(regulators) %in% targets)) {
        cli::cli_abort(c(x = "Regulators: If provided as a named list, all names should be targets."))
      }
      regulators <- unique(unlist(regulators))
    }
    if (!is.null(regulators)) {
      if(length(regulators) < 2) cli::cli_abort(c(x = "Provide at least 2 potential regulators."))

      if (!is.vector(regulators)) {
        cli::cli_abort(c(x = "Parameter 'regulators' must a vector of indices or gene names."))
      }

      if (is.numeric(regulators)) {
        if(max(regulators) > nrow(matrix)) {
          cli::cli_abort(c(x = "At least one index in 'regulators' exceeds the number of genes."))
        }
        if(min(regulators) <= 0) cli::cli_abort(c(x = "The indexes in 'regulators' should be >=1."))
      }

      if(any(table(regulators) > 1)) cli::cli_abort(c(x = "Please provide each regulator only once."))

      if (is.character(regulators)) {
        regulatorsInMatrix <- intersect(regulators, colnames(matrix))
        if(length(regulatorsInMatrix) < 2) {
          cli::cli_abort(c(x = "Fewer than 2 regulators in the columns of expression matrix."))
        }

        if(length(regulatorsInMatrix) < length(regulators)) {
          message.warning("Only {length(regulatorsInMatrix)} out of {length(regulators)} candidate regulators are in the expression matrix.")
        }
      }
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    cli::cli_abort(c(x = "Parameter cores should be a stricly positive integer."))
  }

  if (verbose) message.success("All parameters check done.")

  if (verbose) message.success("Using '", penalty, "' penalty.")
  if (verbose & crossValidation) message.success("Using cross validation.")
}

#' Print information
#'
#' @param ... Input information
#'
#' @return NULL
#' @export
#'
message.success <- function(...) {
  cli::cli_alert_success(cli::col_cyan(...))
}

#' Print warning information
#'
#' @param ... Input information
#'
#' @return NULL
#' @export
#'
message.warning <- function(...) {
  cli::cli_alert_warning(cli::col_br_blue(...))
}

#' @title Format weight table
#'
#' @param weightDT The weight data table of network
#' @param regulators Regulators list
#'
#' @return Format weight table
#' @export
#'
net.format <- function(weightDT,
                       regulators = NULL) {
  colnames(weightDT) <- c("regulator", "target", "weight")
  weightDT$weight <- as.numeric(weightDT$weight)
  weightDT <- dplyr::filter(weightDT, weight !=0)
  if (!is.null(regulators)) {
    weightDT <- purrr::map_dfr(regulators, function(x) {
      dplyr::filter(weightDT, regulator == x)
    })
  }
  weightDT$Interaction <- "Activation"
  weightDT$Interaction[weightDT$weight < 0] <- "Repression"
  weightDT$weight <- abs(weightDT$weight)
  return(weightDT)
}

#' @title Extracts a specific solution in the regularization path
#'
#' @param object The output of inferCSN.fit or inferCSN.cvfit
#' @param lambda The value of lambda at which to extract the solution
#' @param gamma The value of gamma at which to extract the solution
#' @param supportSize The number of non-zeros each solution extracted will contain
#' @param ... Other parameters
#'
#' @method coef inferCSN
#'
#' @export
#'
coef.inferCSN <- function(object,
                          lambda = NULL,
                          gamma = NULL,
                          supportSize = NULL, ...) {
  if (!is.null(supportSize) && !is.null(lambda)) {
    cli::cli_abort(c(x = "If 'supportSize' is provided to 'coef' only 'gamma' can also be provided."))
  }

  if (is.null(lambda) && is.null(gamma) && is.null(supportSize)) {
    # If all three are null, return all solutions
    t <- do.call(cbind, object$beta)
    if (object$settings$intercept) {
      intercepts <- unlist(object$a0)
      t <- rbind(intercepts, t)
    }
    return(t)
  }

  if (is.null(gamma)) gamma <- object$gamma[1]

  diffGamma <- abs(object$gamma - gamma)
  gammaindex <- which(diffGamma == min(diffGamma))

  indices <- NULL
  if (!is.null(lambda)) {
    diffLambda <- abs(lambda - object$lambda[[gammaindex]])
    indices <- which(diffLambda == min(diffLambda))
  } else if (!is.null(supportSize)) {
    diffSupportSize <- abs(supportSize - object$suppSize[[gammaindex]])
    indices <- which(diffSupportSize == min(diffSupportSize))
  } else {
    indices <- seq_along(object$lambda[[gammaindex]])
  }

  if (object$settings$intercept) {
    t <- rbind(object$a0[[gammaindex]][indices],
               object$beta[[gammaindex]][, indices, drop = FALSE])
    rownames(t) <- c("Intercept",
                     paste0(rep("V", object$p), 1:object$p))
  } else {
    t <- object$beta[[gammaindex]][, indices, drop = FALSE]
    rownames(t) <- paste0(rep("V", object$p), 1:object$p)
  }
  t
}

#' @rdname coef.inferCSN
#'
#' @method coef inferCSNCV
#'
#' @export
#'
coef.inferCSNCV <- function(object,
                            lambda = NULL,
                            gamma = NULL, ...) {
  coef.inferCSN(object$fit, lambda, gamma, ...)
}

#' @title Prints a summary of inferCSN.fit
#'
#' @param x The output of inferCSN.fit or inferCSN.cvfit
#' @param ... Other parameters
#'
#' @method print inferCSN
#'
#' @export
#'
print.inferCSN <- function(x, ...) {
  gammas <- rep(x$gamma, times = lapply(x$lambda, length))
  data.frame(lambda = unlist(x["lambda"]),
             gamma = gammas,
             suppSize = unlist(x["suppSize"]),
             row.names = NULL)
}

#' @rdname print.inferCSN
#'
#' @method print inferCSNCV
#'
#' @export
#'
print.inferCSNCV <- function(x, ...) {
  print.inferCSN(x$fit)
}

#' @title Predict Response
#'
#' @description Predicts the response for a given sample
#'
#' @param object The output of inferCSN.fit
#' @param newx A matrix on which predictions are made. The matrix should have p columns
#' @param lambda The value of lambda to use for prediction.
#' A summary of the lambdas in the regularization path can be obtained using \code{print(fit)}
#' @param gamma The value of gamma to use for prediction.
#' A summary of the gammas in the regularization path can be obtained using \code{print(fit)}
#' @param ... Other parameters
#'
#' @method predict inferCSN
#'
#' @details
#' If both lambda and gamma are not supplied, then a matrix of predictions for all the solutions in the regularization path is returned.
#' If lambda is supplied but gamma is not, the smallest value of gamma is used.
#' In case of logistic regression, probability values are returned
#'
#' @export
#'
predict.inferCSN <- function(object,
                             newx,
                             lambda=NULL,
                             gamma=NULL, ...) {
  beta = coef.inferCSN(object, lambda, gamma)
  if (object$settings$intercept){
    # add a column of ones for the intercept
    x = cbind(1,newx)
  }	else{
    x = newx
  }
  prediction = x%*%beta
  if (object$loss == "Logistic"){
    prediction = 1 / (1 + exp(-prediction))
  }
  prediction
}

#' @rdname predict.inferCSN
#'
#' @method predict inferCSNCV
#'
#' @export
#'
predict.inferCSNCV <- function(object,
                               newx,
                               lambda=NULL,
                               gamma=NULL, ...) {
  predict.inferCSN(object$fit, newx, lambda, gamma, ...)
}

#' @title is.scalar
#'
#' @param x Value
#'
#' @return Null
#' @export
#'
is.scalar <- function(x) {
  is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
}