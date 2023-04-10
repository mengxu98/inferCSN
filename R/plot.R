#' @title Plot Regularization Path
#'
#' @description Plots the regularization path for a given gamma.
#' @param gamma The value of gamma at which to plot.
#' @param x The output of inferCSN.fit
#' @param showLines If TRUE, the lines connecting the points in the plot are shown.
#' @param ... ignore
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @method plot inferCSN
#' @export
plot.inferCSN <- function(x, gamma=0, showLines=FALSE, ...){
  j = which(abs(x$gamma-gamma)==min(abs(x$gamma-gamma)))
  p = x$p
  allin = c() # contains all the non-zero variables in the path
  for (i in 1:length(x$lambda[[j]])){
    BetaTemp = x$beta[[j]][,i]
    supp = which(as.matrix(BetaTemp != 0))
    allin = c(allin, supp)
  }
  allin = unique(allin)

  #ggplot needs a dataframe
  yy = t(as.matrix(x$beta[[j]][allin,])) # length(lambda) x length(allin) matrix
  data <- as.data.frame(yy)

  colnames(data)  = x$varnames[allin]

  #id variable for position in matrix
  data$id <- x$suppSize[[j]]

  #reshape to long format
  plot_data <- melt(data,id.var="id")

  #breaks = x$suppSize[[j]]

  #plot
  plotObject = ggplot(plot_data, aes_string(x="id",y="value",group="variable",colour="variable")) + geom_point(size=2.5) +
    labs(x = "Support Size", y = "Coefficient") + theme(axis.title=element_text(size=14)) # + scale_x_continuous(breaks = breaks) + theme(axis.text = element_text(size = 12))

  if (showLines == TRUE){
    plotObject = plotObject + geom_line(aes_string(lty="variable"),alpha=0.3)
  }
  plotObject
}

#' @title Plot Cross-validation Errors
#'
#' @description Plots cross-validation errors for a given gamma.
#' @param x The output of inferCSN.cvfit
#' @inheritParams plot.inferCSN
#'
#' @method plot inferCSNCV
#' @export
plot.inferCSNCV <- function(x, gamma=0, ...){
  j = which(abs(x$fit$gamma-gamma)==min(abs(x$fit$gamma-gamma)))
  data = data.frame(x=x$fit$suppSize[[j]], y=x$cvMeans[[j]], sd=x$cvSDs[[j]])
  ggplot(data, aes_string(x="x",y="y")) + geom_point() + geom_errorbar(aes_string(ymin="y-sd", ymax="y+sd"))+
    labs(x = "Support Size", y = "Cross-validation Error") + theme(axis.title=element_text(size=14)) + theme(axis.text = element_text(size = 12))
}
