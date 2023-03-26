#' @title inferCSN
#' @description A method for inferring cell-type-specific gene regulatory network
#' from single-cell transcriptome data.
#'
#' @param data [Default = NULL] An object
#' @param normalize [Default = FALSE] Data normalize
#' @param penalty [Default = NULL]
#' @param crossValidation [Default = FALSE] Check whether cross validation is used
#' @param nFolds [Default = 10] N folds cross validation
#' @param regulators [Default = NULL] Regulator genes
#' @param targets [Default = NULL] Target genes
#' @param maxSuppSize [Default = NULL] The number of non-zore coef
#' @param nGamma [Default = 5] nGamma
#' @param verbose [Default = FALSE] Print detailed information
#' @param cores [Default = 1] CPU cores
#'
#' @import magrittr
#' @importFrom utils "methods" "read.table" "setTxtProgressBar" "txtProgressBar"
#'
#' @return A list of gene-gene regulatory relationship
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' weightList <- inferCSN(exampleDataMatrix, verbose = TRUE)
inferCSN <- function(data = NULL,
                     normalize = FALSE,
                     penalty = NULL,
                     crossValidation = FALSE,
                     nFolds = 10,
                     regulators = NULL,
                     targets = NULL,
                     maxSuppSize = NULL,
                     nGamma = 5,
                     verbose = FALSE,
                     cores = 1) {
  # Data processing
  if (!is.null(data)) {
    if (verbose) message("Data processing......")
    matrix <- data.processing(data, normalize = normalize, verbose = verbose)
  } else {
    stop("Please ensure provide an object!")
  }

  # Check the penalty terms of the regression model
  if (!is.null(penalty)) {
    if (!any(c("L0", "L0L2") == penalty)) {
      stop(paste(
        "Note: inferCSN does not support", penalty, "penalty regression......\n",
        "Please set penalty item as 'L0' or 'L0L2'......\n"
      ))
    }
  } else {
    penalty <- "L0"
  }
  if (verbose) message(paste("Using", penalty, "penalty regression......"))

  # Check whether cross validation is used
  if (verbose & crossValidation) {
    if (verbose) message(paste("Using ", penalty, "cross validation......"))
  } else {
    if (verbose) message(paste("Using ", penalty, "fit......"))
  }

  if (!is.null(regulators)) {
    regulatorsMatrix <- matrix[, intersect(colnames(matrix), regulators)]
  } else {
    regulatorsMatrix <- matrix
  }
  regulators <- colnames(regulatorsMatrix)

  if (!is.null(targets)) {
    targetsMatrix <- matrix[, intersect(colnames(matrix), targets)]
  } else {
    targetsMatrix <- matrix
  }
  targets <- colnames(targetsMatrix)

  if (is.null(maxSuppSize)) maxSuppSize <- ncol(targetsMatrix)

  if (cores == 1) {
    # Format print progress information
    pb <- progress::progress_bar$new(
      format = "Running [:bar] :percent, No.:current of :total gene,:elapsed......",
      total = length(targets),
      clear = TRUE,
      width = 100
    )

    weightList <- c()
    for (i in 1:length(targets)) {
      X <- as.matrix(regulatorsMatrix[, setdiff(colnames(regulatorsMatrix), targets[i])])

      y <- targetsMatrix[, targets[i]]

      temp <- inferCSN.fit(
        X, y,
        penalty = penalty,
        crossValidation = crossValidation,
        nFolds = nFolds,
        maxSuppSize = maxSuppSize,
        nGamma = nGamma,
        verbose = verbose
      ) %>% as.vector()
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
      if (length(wghts) != ncol(X)) {
        weightd <- data.frame(regulator = colnames(X), target = targets[i], weight = 0)
      } else {
        weightd <- data.frame(regulator = colnames(X), target = targets[i], weight = wghts)
      }
      weightList <- rbind.data.frame(weightList, weightd)

      # Print progress
      pb$tick()
      Sys.sleep(0.05)
    }
  } else {
    cores <- min(parallel::detectCores(logical = FALSE), cores, length(targets))
    cl <- snow::makeSOCKcluster(cores)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(
      min = 1,
      max = length(targets),
      width = 100,
      style = 3
    )
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    "%dopar%" <- foreach::"%dopar%"
    weightList <- foreach::foreach(
      target = targets,
      .combine = "rbind",
      .export = "inferCSN.fit",
      .packages = c("Kendall"),
      .options.snow = opts
    ) %dopar% {
      X <- as.matrix(regulatorsMatrix[, setdiff(colnames(regulatorsMatrix), target)])

      y <- targetsMatrix[, target]

      temp <- inferCSN.fit(
        X, y,
        penalty = penalty,
        crossValidation = crossValidation,
        nFolds = nFolds,
        maxSuppSize = maxSuppSize,
        nGamma = nGamma,
        verbose = verbose
      ) %>% as.vector()
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
      if (length(wghts) != ncol(X)) {
        weightd <- data.frame(regulator = colnames(X), target = target, weight = 0)
      } else {
        weightd <- data.frame(regulator = colnames(X), target = target, weight = wghts)
      }
    }

    # cores <- min(parallel::detectCores(logical = FALSE), cores, length(targets))
    # cl <- parallel::makeCluster(cores)
    # doParallel::registerDoParallel(cl)
    # # doParallel::registerDoParallel(cores = cores)
    # if (verbose) message(paste("Using", foreach::getDoParWorkers(), "cores......"))
    # "%dopar%" <- foreach::"%dopar%"
    # suppressPackageStartupMessages(
    #   weightList <- doRNG::"%dorng%"(
    #     foreach::foreach(
    #       target = targets,
    #       .combine = "rbind",
    #       .export = "inferCSN.fit",
    #       .packages = c("doParallel", "dplyr")
    #     ), {
    #
    #       X <- as.matrix(regulatorsMatrix[, setdiff(colnames(regulatorsMatrix), target)])
    #
    #       y <- targetsMatrix[, target]
    #
    #       temp <- inferCSN.fit(
    #         X, y,
    #         penalty = penalty,
    #         crossValidation = crossValidation,
    #         nFolds = nFolds,
    #         maxSuppSize = maxSuppSize,
    #         nGamma = nGamma,
    #         verbose = verbose
    #       ) %>% as.vector()
    #       wghts <- temp[-1]
    #       wghts <- abs(wghts)
    #       wghts <- wghts / sum(wghts)
    #       if (length(wghts) != ncol(X)) {
    #         weightd <- data.frame(regulator = colnames(X), target = target, weight = 0)
    #       } else {
    #         weightd <- data.frame(regulator = colnames(X), target = target, weight = wghts)
    #       }
    #   })
    # )
    close(pb)
    parallel::stopCluster(cl)
  }
  attr(weightList, "rng") <- NULL
  attr(weightList, "doRNG_version") <- NULL
  weightList <- weightList[order(weightList$weight, decreasing = TRUE), ]
  return(weightList)
}
