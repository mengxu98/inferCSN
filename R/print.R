#' @title Print inferCSN.fit object
#'
#' @description Prints a summary of inferCSN.fit
#' @param x The output of inferCSN.fit or inferCSN.cvfit
#' @param ... ignore
#' @method print inferCSN
#' @export
print.inferCSN <- function(x, ...) {
		gammas = rep(x$gamma, times=lapply(x$lambda, length) )
		data.frame(lambda = unlist(x["lambda"]), gamma = gammas, suppSize = unlist(x["suppSize"]), row.names = NULL)
}

#' @rdname print.inferCSN
#' @method print inferCSNCV
#' @export
print.inferCSNCV <- function(x, ...) {
    print.inferCSN(x$fit)
}
