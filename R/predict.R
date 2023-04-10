#' @title Predict Response
#'
#' @description Predicts the response for a given sample.
#' @param object The output of inferCSN.fit or inferCSN.cvfit
#' @param ... ignore
#' @param newx A matrix on which predictions are made. The matrix should have p columns.
#' @param lambda The value of lambda to use for prediction. A summary of the lambdas in the regularization
#' path can be obtained using \code{print(fit)}.
#' @param gamma The value of gamma to use for prediction. A summary of the gammas in the regularization
#' path can be obtained using \code{print(fit)}.
#' @method predict inferCSN
#' @details
#' If both lambda and gamma are not supplied, then a matrix of predictions
#' for all the solutions in the regularization path is returned. If lambda is
#' supplied but gamma is not, the smallest value of gamma is used. In case of
#' of logistic regression, probability values are returned.
#'
#' @export
predict.inferCSN <- function(object,newx,lambda=NULL,gamma=NULL, ...)
{
		beta = coef.inferCSN(object, lambda, gamma)
		if (object$settings$intercept){
				# add a column of ones for the intercept
				x = cbind(1,newx)
		}
		else{
				x = newx
		}
		prediction = x%*%beta
		#if (object$loss == "Logistic" || object$loss == "SquaredHinge"){
		#		prediction = sign(prediction)
		#}
		if (object$loss == "Logistic"){
				prediction = 1/(1+exp(-prediction))
		}
		prediction
}

#' @rdname predict.inferCSN
#' @method predict inferCSNCV
#' @export
predict.inferCSNCV <- function(object,newx,lambda=NULL,gamma=NULL, ...)
{
    predict.inferCSN(object$fit,newx,lambda,gamma, ...)
}
