utils::globalVariables(c("x",
                         "y",
                         "xend",
                         "yend",
                         "weight",
                         "Interaction",
                         "name",
                         "target",
                         "degree",
                         "edges",
                         "curvetype"))

#' Check input arguments
#'
#' @inheritParams inferCSN
#'
#' @return NULL
#' @export
#'
check.arguments <- function(matrix,
                            penalty,
                            algorithm,
                            crossValidation,
                            nFolds,
                            kFolds,
                            rThreshold,
                            regulators,
                            targets,
                            maxSuppSize,
                            verbose,
                            cores) {
  matrixMessage <- paste0("Parameter matrix must be a two-dimensional matrix,
                          where each column corresponds to a gene and each row corresponds to a sample/cell......")
  if (!is.matrix(matrix) && !is.array(matrix)) {
    stop(matrixMessage)
  }

  if (length(dim(matrix)) != 2) {
    stop(matrixMessage)
  }

  if (is.null(colnames(matrix))) {
    stop("Parameter matrix must contain the names of the genes as colnames......")
  }

  # Check the penalty term of the regression model
  if (!any(c("L0", "L0L2") == penalty)) {
    stop("inferCSN does not support '", penalty, "' penalty regression......\n",
         "Please set penalty item as 'L0' or 'L0L2'......")
  }

  # Check the algorithm of the regression model
  if (!any(c("CD", "CDPSI") == algorithm)) {
    stop("inferCSN does not support '", algorithm, "' algorithm......\n",
         "Please set algorithm as 'CD' or 'CDPSI'......")
  }

  if (!is.null(targets)) {
    if (!is.vector(targets)) {
      stop("Parameter 'targets' must a vector (of indices or gene names)......")
    }

    if (is.numeric(targets)) {
      if(max(targets) > nrow(matrix)) stop("At least one index in 'targets' exceeds the number of genes......")
      if(min(targets) <= 0) stop("The indexes in 'targets' should be >=1......")
    }

    if(any(table(targets) > 1)) stop("Please, provide each target (name/ID) only once......")

    if(is.character(targets)){
      targetsInMatrix <- intersect(targets, colnames(matrix))
      if(length(targetsInMatrix) == 0) {
        stop("The genes must contain at least one target......")
      }

      if(length(targetsInMatrix) < length(targets)) {
        warning("Only", length(targetsInMatrix), "out of", length(targets),
                "target genes are in the expression matrix......")
      }
    }
  }

  if (!is.null(regulators)) {
    if(is.list(regulators)) {
      if(!all(names(regulators) %in% targets)) {
        stop("Regulators: If provided as a named list, all names should be targets......")
      }
      regulators <- unique(unlist(regulators))
    }
    if (!is.null(regulators)) {
      if(length(regulators) < 2) stop("Provide at least 2 potential regulators......")

      if (!is.vector(regulators)) {
        stop("Parameter 'regulators' must a vector of indices or gene names......")
      }

      if (is.numeric(regulators)) {
        if(max(regulators) > nrow(matrix)) {
          stop("At least one index in 'regulators' exceeds the number of genes......")
        }
        if(min(regulators) <= 0) stop("The indexes in 'regulators' should be >=1......")
      }

      if(any(table(regulators) > 1)) stop("Please, provide each regulator only once......")

      if (is.character(regulators)) {
        regulatorsInMatrix <- intersect(regulators, colnames(matrix))
        if(length(regulatorsInMatrix) < 2) {
          stop("Fewer than 2 regulators in the columns of expression matrix......")
        }

        if(length(regulatorsInMatrix) < length(regulators)) {
          warning("Only ", length(regulatorsInMatrix), " out of ", length(regulators),
                  " candidate regulators are in the expression matrix......")
        }
      }
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    stop("Parameter cores should be a stricly positive integer......")
  }

  if (verbose) message("All arguments check done......")
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
  if (!is.null(regulators)) {
    weightDT <- purrr::map_dfr(
      regulators, function(x) {
        weightDT[which(weightDT$regulator == x), ]
      }
    )
  }
  weightDT$weight <- as.numeric(weightDT$weight)
  weightDT$Interaction <- "Activation"
  weightDT$Interaction[weightDT$weight < 0] <- "Repression"
  weightDT$weight <- abs(weightDT$weight)
  return(weightDT)
}

#' @title Extracts a specific solution in the regularization path.
#'
#' @param object The output of inferCSN.fit or inferCSN.cvfit
#' @param lambda The value of lambda at which to extract the solution.
#' @param gamma The value of gamma at which to extract the solution.
#' @param supportSize The number of non-zeros each solution extracted will contain
#' @param ... ignore
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
    stop("If `supportSize` is provided to `coef` only `gamma` can also be provided......")
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
#' @param ... ignore
#' @method print inferCSN
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
#' @param object The output of inferCSN.fit or inferCSN.cvfit
#' @param newx A matrix on which predictions are made. The matrix should have p columns
#' @param lambda The value of lambda to use for prediction. A summary of the lambdas in the regularization
#' path can be obtained using \code{print(fit)}
#' @param gamma The value of gamma to use for prediction. A summary of the gammas in the regularization
#' path can be obtained using \code{print(fit)}
#' @param ... ignore
#'
#' @method predict inferCSN
#'
#' @details
#' If both lambda and gamma are not supplied, then a matrix of predictions
#' for all the solutions in the regularization path is returned. If lambda is
#' supplied but gamma is not, the smallest value of gamma is used. In case of
#' of logistic regression, probability values are returned
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
