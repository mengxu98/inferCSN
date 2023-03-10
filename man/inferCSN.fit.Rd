% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inferCSN.fit.R
\name{inferCSN.fit}
\alias{inferCSN.fit}
\title{inferCSN.fit}
\usage{
inferCSN.fit(
  X,
  y,
  crossValidation = crossValidation,
  penalty = penalty,
  maxSuppSize = maxSuppSize,
  nFolds = nFolds,
  nGamma = nGamma,
  verbose = verbose
)
}
\arguments{
\item{X}{The rows are samples and the columns are genes of the matrix}

\item{y}{Target vector}

\item{crossValidation}{Cross validation}

\item{penalty}{penalty = penalty}

\item{maxSuppSize}{maxSuppSize = maxSuppSize}

\item{nFolds}{N folds cross validation}

\item{nGamma}{nGamma = nGamma}

\item{verbose}{Print detailed information}
}
\value{
A vector
}
\description{
inferCSN.fit
}
