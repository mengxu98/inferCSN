# import C++ compiled code
#' @useDynLib inferCSN
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom methods is
#' @import Matrix

#' @title Fit an L0-regularized model
#'
#' @description Computes the regularization path for the specified loss function and
#' penalty function (which can be a combination of the L0, L1, and L2 norms).
#' @param x The data matrix.
#' @param y The response vector. For classification, we only support binary vectors.
#' @param loss The loss function. Currently we support the choices "SquaredError" (for regression), "Logistic" (for logistic regression), and "SquaredHinge" (for smooth SVM).
#' @param penalty The type of regularization. This can take either one of the following choices:
#' "L0", "L0L2", and "L0L1".
#' @param algorithm The type of algorithm used to minimize the objective function. Currently "CD" and "CDPSI" are
#' are supported. "CD" is a variant of cyclic coordinate descent and runs very fast. "CDPSI" performs
#' local combinatorial search on top of CD and typically achieves higher quality solutions (at the expense
#' of increased running time).
#' @param maxSuppSize The maximum support size at which to terminate the regularization path. We recommend setting
#' this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small
#' portion of non-zeros.
#' @param nLambda The number of Lambda values to select (recall that Lambda is the regularization parameter
#' corresponding to the L0 norm). This value is ignored if 'lambdaGrid' is supplied.
#' @param nGamma The number of Gamma values to select (recall that Gamma is the regularization parameter
#' corresponding to L1 or L2, depending on the chosen penalty). This value is ignored if 'lambdaGrid' is supplied
#' and will be set to length(lambdaGrid)
#' @param gammaMax The maximum value of Gamma when using the L0L2 penalty. For the L0L1 penalty this is
#' automatically selected.
#' @param gammaMin The minimum value of Gamma when using the L0L2 penalty. For the L0L1 penalty, the minimum
#' value of gamma in the grid is set to gammaMin * gammaMax. Note that this should be a strictly positive quantity.
#' @param partialSort If TRUE partial sorting will be used for sorting the coordinates to do greedy cycling (see our paper for
#' for details). Otherwise, full sorting is used.
#' @param maxIters The maximum number of iterations (full cycles) for CD per grid point.
#' @param rtol The relative tolerance which decides when to terminate optimization (based on the relative change in the objective between iterations).
#' @param atol The absolute tolerance which decides when to terminate optimization (based on the absolute L2 norm of the residuals).
#' @param activeSet If TRUE, performs active set updates.
#' @param activeSetNum The number of consecutive times a support should appear before declaring support stabilization.
#' @param maxSwaps The maximum number of swaps used by CDPSI for each grid point.
#' @param scaleDownFactor This parameter decides how close the selected Lambda values are. The choice should be
#' strictly between 0 and 1 (i.e., 0 and 1 are not allowed). Larger values lead to closer lambdas and typically to smaller
#' gaps between the support sizes. For details, see our paper - Section 5 on Adaptive Selection of Tuning Parameters).
#' @param screenSize The number of coordinates to cycle over when performing initial correlation screening.
#' @param autoLambda Ignored parameter. Kept for backwards compatibility.
#' @param lambdaGrid A grid of Lambda values to use in computing the regularization path. This is by default an empty list and is ignored.
#' When specified, LambdaGrid should be a list of length 'nGamma', where the ith element (corresponding to the ith gamma) should be a decreasing sequence of lambda values
#' which are used by the algorithm when fitting for the ith value of gamma (see the vignette for details).
#' @param excludeFirstK This parameter takes non-negative integers. The first excludeFirstK features in x will be excluded from variable selection,
#' i.e., the first excludeFirstK variables will not be included in the L0-norm penalty (they will still be included in the L1 or L2 norm penalties.).
#' @param intercept If FALSE, no intercept term is included in the model.
#' @param lows Lower bounds for coefficients. Either a scalar for all coefficients to have the same bound or a vector of size p (number of columns of X) where lows[i] is the lower bound for coefficient i.
#' @param highs Upper bounds for coefficients. Either a scalar for all coefficients to have the same bound or a vector of size p (number of columns of X) where highs[i] is the upper bound for coefficient i.
#' @return An S3 object of type "inferCSN" describing the regularization path. The object has the following members.
#' \item{a0}{a0 is a list of intercept sequences. The ith element of the list (i.e., a0[[i]]) is the sequence of intercepts corresponding to the ith gamma value (i.e., gamma[i]).}
#' \item{beta}{This is a list of coefficient matrices. The ith element of the list is a p x \code{length(lambda)} matrix which
#' corresponds to the ith gamma value. The jth column in each coefficient matrix is the vector of coefficients for the jth lambda value.}
#' \item{lambda}{This is the list of lambda sequences used in fitting the model. The ith element of lambda (i.e., lambda[[i]]) is the sequence
#' of Lambda values corresponding to the ith gamma value.}
#' \item{gamma}{This is the sequence of gamma values used in fitting the model.}
#' \item{suppSize}{This is a list of support size sequences. The ith element of the list is a sequence of support sizes (i.e., number of non-zero coefficients)
#' corresponding to the ith gamma value.}
#' \item{converged}{This is a list of sequences for checking whether the algorithm has converged at every grid point. The ith element of the list is a sequence
#' corresponding to the ith value of gamma, where the jth element in each sequence indicates whether the algorithm has converged at the jth value of lambda.}
#'
#' @export
inferCSN.fit <- function(x, y,
                         loss="SquaredError",
                         penalty="L0",
                         algorithm="CD",
                         maxSuppSize=100,
                         nLambda=100,
                         nGamma=10,
                         gammaMax=10,
                         gammaMin=0.0001,
                         partialSort = TRUE,
                         maxIters=200,
                         rtol=1e-6,
                         atol=1e-9,
                         activeSet=TRUE,
                         activeSetNum=3,
                         maxSwaps=100,
                         scaleDownFactor=0.8,
                         screenSize=1000,
                         autoLambda = NULL,
                         lambdaGrid = list(),
                         excludeFirstK=0,
                         intercept = TRUE,
                         lows=-Inf,
                         highs=Inf) {

  if ((rtol < 0) || (rtol >= 1)){
    stop("The specified rtol parameter must exist in [0, 1)")
  }

  if (atol < 0){
    stop("The specified atol parameter must exist in [0, INF)")
  }

  # Some sanity checks for the inputs
  if ( !(loss %in% c("SquaredError","Logistic","SquaredHinge")) ){
    stop("The specified loss function is not supported.")
  }
  if ( !(penalty %in% c("L0","L0L2","L0L1")) ){
    stop("The specified penalty is not supported.")
  }
  if ( !(algorithm %in% c("CD","CDPSI")) ){
    stop("The specified algorithm is not supported.")
  }
  if (loss=="Logistic" | loss=="SquaredHinge"){
    if (dim(table(y)) != 2){
      stop("Only binary classification is supported. Make sure y has only 2 unique values.")
    }
    y = factor(y,labels=c(-1,1)) # returns a vector of strings
    y = as.numeric(levels(y))[y]

    if (penalty == "L0"){

      if ((length(lambdaGrid) != 0) && (length(lambdaGrid) != 1)){
        stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
    			             Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
      }
      penalty = "L0L2"
      nGamma = 1
      gammaMax = 1e-7
      gammaMin = 1e-7
    }
  }

  # Handle Lambda Grids:
  if (length(lambdaGrid) != 0){
    if (!is.null(autoLambda)){
      warning("autoLambda is ignored and inferred if 'lambdaGrid' is supplied", call.=FALSE)
    }
    autoLambda = FALSE
  } else {
    autoLambda = TRUE
    lambdaGrid = list(0)
  }

  if (penalty == "L0" && !autoLambda){
    bad_lambdaGrid = FALSE
    if (length(lambdaGrid) != 1){
      bad_lambdaGrid = TRUE
    }
    current = Inf
    for (nxt in lambdaGrid[[1]]){
      if (nxt > current){
        bad_lambdaGrid = TRUE
        break
      }
      if (nxt < 0){
        bad_lambdaGrid = TRUE
        break
      }
      current = nxt

    }

    if (bad_lambdaGrid){
      stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
                 Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
    }
  }

  if (penalty != "L0" && !autoLambda){
    # Covers L0L1, L0L2 cases
    bad_lambdaGrid = FALSE
    if (length(lambdaGrid) != nGamma){
      warning("nGamma is ignored and replaced with length(lambdaGrid)", call.=FALSE)
      nGamma = length(lambdaGrid)
    }

    for (i in 1:length(lambdaGrid)){
      current = Inf
      for (nxt in lambdaGrid[[i]]){
        if (nxt > current){
          bad_lambdaGrid = TRUE
          break
        }
        if (nxt < 0){
          bad_lambdaGrid = TRUE
          break
        }
        current = nxt
      }
      if (bad_lambdaGrid){
        break
      }
    }

    if (bad_lambdaGrid){
      stop("L0L1 or L0L2 Penalty requires 'lambdaGrid' to be a list of length 'nGamma'.
                 Where lambdaGrid[[i]] is a list or vector of decreasing positive values.")
    }


  }

  is.scalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0 && !is.nan(x) && !is.na(x)

  p = dim(x)[[2]]

  withBounds = FALSE

  if ((!identical(lows, -Inf)) || (!identical(highs, Inf))){
    withBounds = TRUE

    if (algorithm == "CDPSI"){
      if (any(lows != -Inf) || any(highs != Inf)){
        stop("Bounds are not YET supported for CDPSI algorithm.")
      }
    }

    if (is.scalar(lows)){
      lows = lows*rep(1, p)
    } else if (!all(sapply(lows, is.scalar)) || length(lows) != p) {
      stop('Lows must be a vector of real values of length p')
    }

    if (is.scalar(highs)){
      highs = highs*rep(1, p)
    } else if (!all(sapply(highs, is.scalar)) || length(highs) != p) {
      stop('Highs must be a vector of real values of length p')
    }

    if (any(lows >= highs) || any(lows > 0) || any(highs < 0)){
      stop("Bounds must conform to the following conditions: Lows <= 0, Highs >= 0, Lows < Highs")
    }

  }

  M = list()
  if (is(x, "sparseMatrix")){
    M <- .Call('_inferCSN_inferCSNFit_sparse', PACKAGE = 'inferCSN', x, y, loss, penalty,
               algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
               partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
               scaleDownFactor, screenSize, !autoLambda, lambdaGrid,
               excludeFirstK, intercept, withBounds, lows, highs)
  } else{
    M <- .Call('_inferCSN_inferCSNFit_dense', PACKAGE = 'inferCSN', x, y, loss, penalty,
               algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
               partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
               scaleDownFactor, screenSize, !autoLambda, lambdaGrid,
               excludeFirstK, intercept, withBounds, lows, highs)
  }

  settings = list()
  settings[[1]] = intercept # Settings only contains intercept for now. Might include additional elements later.
  names(settings) <- c("intercept")

  # Find potential support sizes exceeding maxSuppSize and remove them (this is due to
  # the C++ core whose last solution can exceed maxSuppSize
  for (i in 1:length(M$SuppSize)){
    last = length(M$SuppSize[[i]])
    if (M$SuppSize[[i]][last] > maxSuppSize){
      if (last == 1){
        warning("Warning! Only 1 element in path with support size > maxSuppSize. \n
                Try increasing maxSuppSize to resolve the issue.")
      } else{
        M$SuppSize[[i]] = M$SuppSize[[i]][-last]
        M$Converged[[i]] = M$Converged[[i]][-last]
        M$lambda[[i]] = M$lambda[[i]][-last]
        M$a0[[i]] = M$a0[[i]][-last]
        M$beta[[i]] = as(M$beta[[i]][,-last], "sparseMatrix")
      }
    }
  }

  G <- list(beta = M$beta,
            lambda=lapply(M$lambda,signif, digits=6),
            a0=M$a0,
            converged = M$Converged,
            suppSize= M$SuppSize,
            gamma=M$gamma,
            penalty=penalty,
            loss=loss,
            settings = settings)

  if (is.null(colnames(x))){
    varnames <- 1:dim(x)[2]
  } else {
    varnames <- colnames(x)
  }
  G$varnames <- varnames

  class(G) <- "inferCSN"
  G$n <- dim(x)[1]
  G$p <- dim(x)[2]
  G
}

#' @title Cross Validation
#'
#' @inheritParams inferCSN.fit
#' @description Computes a regularization path and performs K-fold cross-validation.
#' @param nFolds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @return An S3 object of type "inferCSNCV" describing the regularization path. The object has the following members.
#' \item{cvMeans}{This is a list, where the ith element is the sequence of cross-validation errors corresponding to the ith gamma value, i.e., the sequence
#' cvMeans[[i]] corresponds to fit$gamma[i]}
#' \item{cvSDs}{This a list, where the ith element is a sequence of standard deviations for the cross-validation errors: cvSDs[[i]] corresponds to cvMeans[[i]].}
#' \item{fit}{The fitted model with type "inferCSN", i.e., this is the same object returned by \code{\link{inferCSN.fit}}.}
#'
#' @export
inferCSN.cvfit <- function(x,y, loss="SquaredError", penalty="L0", algorithm="CD",
                           maxSuppSize=100, nLambda=100, nGamma=10, gammaMax=10,
                           gammaMin=0.0001, partialSort = TRUE, maxIters=200,
                           rtol=1e-6, atol=1e-9, activeSet=TRUE, activeSetNum=3, maxSwaps=100,
                           scaleDownFactor=0.8, screenSize=1000, autoLambda=NULL,
                           lambdaGrid = list(), nFolds=10, seed=1, excludeFirstK=0,
                           intercept=TRUE, lows=-Inf, highs=Inf)
{
  set.seed(seed)

  if ((rtol < 0) || (rtol >= 1)){
    stop("The specified rtol parameter must exist in [0, 1)")
  }

  if (atol < 0){
    stop("The specified atol parameter must exist in [0, INF)")
  }

  # Some sanity checks for the inputs
  if ( !(loss %in% c("SquaredError","Logistic","SquaredHinge")) ){
    stop("The specified loss function is not supported.")
  }
  if ( !(penalty %in% c("L0","L0L2","L0L1")) ){
    stop("The specified penalty is not supported.")
  }
  if ( !(algorithm %in% c("CD","CDPSI")) ){
    stop("The specified algorithm is not supported.")
  }
  if (loss=="Logistic" | loss=="SquaredHinge"){
    if (dim(table(y)) != 2){
      stop("Only binary classification is supported. Make sure y has only 2 unique values.")
    }
    y = factor(y,labels=c(-1,1)) # returns a vector of strings
    y = as.numeric(levels(y))[y]

    if (penalty == "L0"){
      if ((length(lambdaGrid) != 0) && (length(lambdaGrid) != 1)){
        stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
    			             Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
      }
      penalty = "L0L2"
      nGamma = 1
      gammaMax = 1e-7
      gammaMin = 1e-7
    }
  }

  # Handle Lambda Grids:
  if (length(lambdaGrid) != 0){
    if (!is.null(autoLambda) && !autoLambda){
      warning("autoLambda is ignored and inferred if 'lambdaGrid' is supplied", call.=FALSE)
    }
    autoLambda = FALSE
  } else {
    autoLambda = TRUE
    lambdaGrid = list(0)
  }

  if (penalty == "L0" && !autoLambda){
    bad_lambdaGrid = FALSE
    if (length(lambdaGrid) != 1){
      bad_lambdaGrid = TRUE
    }
    current = Inf
    for (nxt in lambdaGrid[[1]]){
      if (nxt > current){
        # This must be > instead of >= to allow first iteration L0L1 lambdas of all 0s to be valid
        bad_lambdaGrid = TRUE
        break
      }
      if (nxt < 0){
        bad_lambdaGrid = TRUE
        break
      }
      current = nxt

    }

    if (bad_lambdaGrid){
      stop("L0 Penalty requires 'lambdaGrid' to be a list of length 1.
                 Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
    }
  }

  if (penalty != "L0" && !autoLambda){
    # Covers L0L1, L0L2 cases
    bad_lambdaGrid = FALSE
    if (length(lambdaGrid) != nGamma){
      warning("nGamma is ignored and replaced with length(lambdaGrid)", call.=FALSE)
      nGamma = length(lambdaGrid)
    }

    for (i in 1:length(lambdaGrid)){
      current = Inf
      for (nxt in lambdaGrid[[i]]){
        if (nxt > current){
          # This must be > instead of >= to allow first iteration L0L1 lambdas of all 0s to be valid
          bad_lambdaGrid = TRUE
          break
        }
        if (nxt < 0){
          bad_lambdaGrid = TRUE
          break
        }
        current = nxt
      }
      if (bad_lambdaGrid){
        break
      }
    }

    if (bad_lambdaGrid){
      stop("L0L1 or L0L2 Penalty requires 'lambdaGrid' to be a list of length 'nGamma'.
                 Where lambdaGrid[[i]] is a list or vector of decreasing positive values.")
    }


  }

  is.scalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0 && !is.nan(x) && !is.na(x)

  p = dim(x)[[2]]

  withBounds = FALSE
  if ((!identical(lows, -Inf)) || (!identical(highs, Inf))){
    withBounds=TRUE

    if (algorithm == "CDPSI"){
      if (any(lows != -Inf) || any(highs != Inf)){
        stop("Bounds are not YET supported for CDPSI algorithm.")
      }
    }

    if (is.scalar(lows)){
      lows = lows*rep(1, p)
    } else if (!all(sapply(lows, is.scalar)) || length(lows) != p) {
      stop('Lows must be a vector of real values of length p')
    }

    if (is.scalar(highs)){
      highs = highs*rep(1, p)
    } else if (!all(sapply(highs, is.scalar)) || length(highs) != p) {
      stop('Highs must be a vector of real values of length p')
    }

    if (any(lows >= highs) || any(lows > 0) || any(highs < 0)){
      stop("Bounds must conform to the following conditions: Lows <= 0, Highs >= 0, Lows < Highs")
    }

  }

  M = list()
  if (is(x, "sparseMatrix")){
    M <- .Call('_inferCSN_inferCSNCV_sparse', PACKAGE = 'inferCSN', x, y, loss, penalty,
               algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
               partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
               scaleDownFactor, screenSize, !autoLambda, lambdaGrid, nFolds,
               seed, excludeFirstK, intercept, withBounds, lows, highs)
  } else {
    M <- .Call('_inferCSN_inferCSNCV_dense', PACKAGE = 'inferCSN', x, y, loss, penalty,
               algorithm, maxSuppSize, nLambda, nGamma, gammaMax, gammaMin,
               partialSort, maxIters, rtol, atol, activeSet, activeSetNum, maxSwaps,
               scaleDownFactor, screenSize, !autoLambda, lambdaGrid, nFolds,
               seed, excludeFirstK, intercept, withBounds, lows, highs)
  }

  settings = list()
  settings[[1]] = intercept # Settings only contains intercept for now. Might include additional elements later.
  names(settings) <- c("intercept")

  # the C++ core whose last solution can exceed maxSuppSize
  for (i in 1:length(M$SuppSize)){
    last = length(M$SuppSize[[i]])
    if (M$SuppSize[[i]][last] > maxSuppSize){
      if (last == 1){
        warning("Warning! Only 1 element in path with support size > maxSuppSize. \n
                Try increasing maxSuppSize to resolve the issue.")
      }
      else{
        M$SuppSize[[i]] = M$SuppSize[[i]][-last]
        M$Converged[[i]] = M$Converged[[i]][-last]
        M$lambda[[i]] = M$lambda[[i]][-last]
        M$a0[[i]] = M$a0[[i]][-last]
        M$beta[[i]] = as(M$beta[[i]][,-last], "sparseMatrix") # conversion to sparseMatrix is necessary to handle the case of a single column
        M$CVMeans[[i]] = M$CVMeans[[i]][-last]
        M$CVSDs[[i]] = M$CVSDs[[i]][-last]
      }
    }
  }

  fit <- list(beta = M$beta,
              lambda=lapply(M$lambda,signif, digits=6),
              a0=M$a0,
              converged = M$Converged,
              suppSize= M$SuppSize,
              gamma=M$gamma,
              penalty=penalty,
              loss=loss,
              settings=settings)

  if (is.null(colnames(x))){
    varnames <- 1:dim(x)[2]
  } else {
    varnames <- colnames(x)
  }
  fit$varnames <- varnames
  class(fit) <- "inferCSN"
  fit$n <- dim(x)[1]
  fit$p <- dim(x)[2]
  G <- list(fit=fit, cvMeans=M$CVMeans,cvSDs=M$CVSDs)
  class(G) <- "inferCSNCV"
  G
}

#' @title Extract Solutions
#'
#' @description Extracts a specific solution in the regularization path.
#'
#' @param object The output of inferCSN.fit or inferCSN.cvfit
#' @param lambda The value of lambda at which to extract the solution.
#' @param gamma The value of gamma at which to extract the solution.
#' @param supportSize The number of non-zeros each solution extracted will
#' contain. If no solutions have `supportSize` non-zeros, solutions with
#' the closest number will be extracted.
#' @param ... ignore
#'
#' @method coef inferCSN
#'
#' @export
coef.inferCSN <- function(object,
                          lambda = NULL,
                          gamma = NULL,
                          supportSize = NULL, ...) {
  if (!is.null(supportSize) && !is.null(lambda)) {
    stop("If `supportSize` is provided to `coef` only `gamma` can also be provided.")
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

  if (is.null(gamma)) {
    # if lambda is present but gamma is not, use smallest value of gamma
    gamma <- object$gamma[1]
  }

  diffGamma <- abs(object$gamma - gamma)
  gammaindex <- which(diffGamma == min(diffGamma))

  indices <- NULL
  if (!is.null(lambda)) {
    diffLambda <- abs(lambda - object$lambda[[gammaindex]])
    indices <- which(diffLambda == min(diffLambda))
  } else if(!is.null(supportSize)) {
    diffSupportSize <- abs(supportSize - object$suppSize[[gammaindex]])
    indices <- which(diffSupportSize == min(diffSupportSize))
  } else {
    indices <- seq_along(object$lambda[[gammaindex]])
  }

  if (object$settings$intercept) {
    t <- rbind(object$a0[[gammaindex]][indices],
               object$beta[[gammaindex]][, indices, drop = FALSE])
    rownames(t) <- c("Intercept",
                     paste(rep("V", object$p),
                           1:object$p,
                           sep = ""))
  } else {
    t <- object$beta[[gammaindex]][, indices, drop = FALSE]
    rownames(t) <- paste(rep("V", object$p),
                         1:object$p,
                         sep = "")
  }
  t
}

#' @rdname coef.inferCSN
#' @method coef inferCSNCV
#' @export
coef.inferCSNCV <- function(object, lambda=NULL, gamma=NULL, ...) {
  coef.inferCSN(object$fit, lambda, gamma, ...)
}

#' @title Print inferCSN.fit object
#'
#' @description Prints a summary of inferCSN.fit
#' @param x The output of inferCSN.fit or inferCSN.cvfit
#' @param ... ignore
#' @method print inferCSN
#' @export
print.inferCSN <- function(x, ...) {
  gammas = rep(x$gamma, times=lapply(x$lambda, length) )
  data.frame(lambda = unlist(x["lambda"]),
             gamma = gammas,
             suppSize = unlist(x["suppSize"]),
             row.names = NULL)
}

#' @rdname print.inferCSN
#' @method print inferCSNCV
#' @export
print.inferCSNCV <- function(x, ...) {
  print.inferCSN(x$fit)
}

