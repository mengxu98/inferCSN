utils::globalVariables(c("x", "y", "xend", "yend", "weight", "Interaction", "name", ".", "target", "curvetype"))

#' @title Sparse regression model
#'
#' @param X The data matrix
#' @param y The response vector
#' @inheritParams inferCSN
#'
#' @importFrom stats coef
#'
#' @return The coefficients
#' @export
#'
sparse.regression <- function(X, y,
                              crossValidation = FALSE,
                              penalty = "L0",
                              algorithm = "CD",
                              maxSuppSize = NULL,
                              nFolds = 10,
                              verbose = FALSE) {
  if (crossValidation) {
    message("Using '", penalty, "' penalty and cross validation......")
    fit <- try(inferCSN.cvfit(X, y,
                              penalty = penalty,
                              algorithm = algorithm,
                              maxSuppSize = maxSuppSize,
                              nFolds = nFolds))
    if (class(fit)[1] == "try-error") {
      if (verbose) message("Cross validation error, used fit instead......")
      fit <- inferCSN.fit(X, y,
                          penalty = penalty,
                          algorithm = algorithm,
                          maxSuppSize = maxSuppSize)

      fitInf <- print(fit)
      lambda <- fitInf$lambda[fitInf$suppSize %>% which.max()]
      gamma <- fitInf$gamma[fitInf$suppSize %>% which.max()]
    } else {
      gamma <- fit$fit$gamma[which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))]
      lambdaList <- print(fit) %>% dplyr::filter(gamma == gamma, )
      if (maxSuppSize %in% lambdaList$maxSuppSize) {
        lambda <- lambdaList$maxSuppSize[which(lambdaList$maxSuppSize == maxSuppSize)]
      } else {
        lambda <- min(lambdaList$lambda)
      }
    }
  } else {
    message("Using '", penalty, "' penalty......")
    fit <- inferCSN.fit(X, y,
                        penalty = penalty,
                        algorithm = algorithm,
                        maxSuppSize = maxSuppSize)

    fitInf <- print(fit)
    lambda <- fitInf$lambda[fitInf$suppSize %>% which.max()]
    gamma <- fitInf$gamma[fitInf$suppSize %>% which.max()]
  }
  return(coef(fit, lambda = lambda, gamma = gamma) %>% as.vector() %>% .[-1])
}

#' @title Sparse regression model for single gene
#'
#' @param regulatorsMatrix regulatorsMatrix
#' @param targetsMatrix targetsMatrix
#' @param target target
#' @inheritParams inferCSN
#'
#' @return The weight data table of sub-network
#' @export
#'
sub.inferCSN <- function(regulatorsMatrix,
                         targetsMatrix,
                         target = NULL,
                         crossValidation = FALSE,
                         penalty = "L0",
                         algorithm = "CD",
                         maxSuppSize = NULL,
                         nFolds = 10,
                         verbose = FALSE) {
  X <- as.matrix(regulatorsMatrix[, setdiff(colnames(regulatorsMatrix), target)])
  y <- targetsMatrix[, target]

  if (is.null(maxSuppSize)) maxSuppSize <- ncol(X)

  coefficients <- sparse.regression(X, y,
                                    crossValidation = crossValidation,
                                    penalty = penalty,
                                    algorithm = algorithm,
                                    maxSuppSize = maxSuppSize,
                                    nFolds = nFolds,
                                    verbose = verbose)

  coefficients <- coefficients / sum(abs(coefficients))
  if (length(coefficients) != ncol(X)) coefficients <- 0
  return(data.frame(regulator = colnames(X), target = target, weight = coefficients))
}

#' @title AUC value calculate
#'
#' @param weightDT The weight data table of network
#' @param groundTruth Ground truth for calculate AUC
#' @param plot If true, draw and print figure of AUC
#' @param fileSave The figure name
#' @param interaction If true, consider the positivity or negativity of interaction
#'
#' @import patchwork
#' @import ggplot2
#'
#' @return AUC values and figure
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' data("exampleGroundTruth")
#' weightDT <- inferCSN(exampleMatrix, cores = 1, verbose = TRUE, algorithm = "CDPSI")
#' auc <- auc.calculate(weightDT, exampleGroundTruth, plot = TRUE)
#' head(auc)
#'
auc.calculate <- function(weightDT,
                          groundTruth,
                          plot = FALSE,
                          fileSave = NULL,
                          interaction = FALSE) {
  # Check input data
  colnames(weightDT) <- c("regulator", "target", "weight")
  if (!interaction) weightDT$weight <- abs(weightDT$weight)

  if (ncol(groundTruth) > 2) groundTruth <- groundTruth[, 1:2]
  names(groundTruth) <- c("regulator", "target")

  groundTruth$gold <- rep(1, nrow(groundTruth))
  gold <- merge(weightDT, groundTruth, by = c("regulator", "target"), all.x = TRUE)
  gold$gold[is.na(gold$gold)] <- 0
  aucCurves <- precrec::evalmod(scores = gold$weight, labels = gold$gold)

  auc <- attr(aucCurves, "auc")
  aucMetric <- data.frame(AUROC = rep(0.000, 1), AUPRC = rep(0.000, 1))
  aucMetric[1, "AUROC"] <- sprintf("%0.3f", auc$aucs[1])
  aucMetric[1, "AUPRC"] <- sprintf("%0.3f", auc$aucs[2])

  if (plot) {
    # Subset data to separate prc and roc
    auprcDf <- subset(fortify(aucCurves), curvetype == "PRC")
    aurocDf <- subset(fortify(aucCurves), curvetype == "ROC")

    # Plot
    auroc <- ggplot(aurocDf, aes(x = x, y = y)) +
      geom_line() +
      geom_abline(slope = 1,
                  color = "gray",
                  linetype = "dotted", linewidth = 1) +
      labs(title = paste("AUROC:", aucMetric[1]),
           x = "False positive rate",
           y = "True positive rate") +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    auprc <- ggplot(auprcDf, aes(x = x, y = y)) +
      geom_line() +
      labs(title = paste("AUPRC:", aucMetric[2]),
           x = "Recall",
           y = "Precision") +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    # Combine two plots by patchwork
    p <- auroc + auprc
    print(p)

    # Save figure
    if (!is.null(fileSave)) {
      if (!grepl(".*\\.(pdf|png|jpe?g)$", fileSave)) fileSave <- paste0(fileSave, ".png")
      cowplot::ggsave2(file = fileSave,
                       p,
                       width = 7,
                       height = 3,
                       dpi = 600)
    }
  }
  return(aucMetric)
}

#' @title The heatmap of network
#'
#' @param weightDT The weight data table of network
#' @param switchMatrix switchMatrix
#' @param heatmapSize heatmapSize
#' @param heatmapTitle heatmapTitle
#' @param heatmapColor heatmapColor
#' @param showNames showNames
#' @param legendName legendName
#'
#' @return heatmap
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' data("exampleGroundTruth")
#' weightDT <- inferCSN(exampleMatrix)
#' p1 <- network.heatmap(exampleGroundTruth,
#'                       heatmapTitle = "Ground truth")
#'
#' p2 <- network.heatmap(weightDT,
#'                       legendName = "Weight2",
#'                       heatmapTitle = "inferCSN")
#'
#' ComplexHeatmap::draw(p1 + p2)
#'
#' p3 <- network.heatmap(weightDT,
#'                       heatmapTitle = "inferCSN",
#'                       heatmapColor = c("#20a485", "#410054", "#fee81f"))
#'
#' p4 <- network.heatmap(weightDT,
#'                       heatmapTitle = "inferCSN",
#'                       legendName = "Weight2",
#'                       heatmapColor = c("#20a485", "white", "#fee81f"))
#'
#' ComplexHeatmap::draw(p3 + p4)
#'
#' p5 <- network.heatmap(weightDT,
#'                       heatmapTitle = "inferCSN",
#'                       showNames = TRUE)
#' p5
#'
network.heatmap <- function(weightDT,
                            switchMatrix = TRUE,
                            heatmapSize = NULL,
                            heatmapTitle = NULL,
                            heatmapColor = NULL,
                            showNames = FALSE,
                            legendName = NULL) {
  if (switchMatrix) {
    colnames(weightDT) <- c("regulator", "target", "weight")
    genes <- c(weightDT$regulator, weightDT$target)
    weightMatrix <- .Call("_inferCSN_DT2Matrix", PACKAGE = "inferCSN", weightDT)
  } else {
    genes <- c(rownames(weightDT), colnames(weightDT))
    weightMatrix <- weightDT
  }
  genes <- gtools::mixedsort(unique(genes))
  weightMatrix <- weightMatrix[genes, genes]

  if (is.null(legendName)) legendName <- "Weight"

  if (is.null(heatmapColor)) heatmapColor <- c("#1966ad", "white", "#bb141a")

  if (showNames) {
    if (is.null(heatmapSize)) heatmapSize <- length(genes) / 2
  } else {
    heatmapSize <- 6
  }

  minWeight <- min(weightMatrix)
  maxWeight <- max(weightMatrix)
  if (minWeight >= 0) {
    colorFun <- circlize::colorRamp2(c(minWeight, maxWeight), heatmapColor[-1])
  } else if (maxWeight <= 0) {
    colorFun <- circlize::colorRamp2(c(minWeight, maxWeight), heatmapColor[-3])
  } else {
    colorFun <- circlize::colorRamp2(c(minWeight, 0, maxWeight), heatmapColor)
  }

  p <- ComplexHeatmap::Heatmap(weightMatrix,
                               name = legendName,
                               col = colorFun,
                               column_title = heatmapTitle,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               show_row_names = showNames,
                               show_column_names = showNames,
                               width = unit(heatmapSize, "cm"),
                               height = unit(heatmapSize, "cm"),
                               border = "black")
  return(p)
}

#' @title Plot of dynamic networks
#'
#' @param weightDT weightDT
#' @param regulators regulators
#' @param legend.position legend.position
#'
#' @import ggplot2
#' @import ggnetwork
#'
#' @return A list of ggplot2 objects
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix)
#' ranks <- compute.gene.rank(weightDT)
#' p <- dynamic.networks(weightDT, ranks[1, 1])
#' p
#'
dynamic.networks <- function(weightDT,
                             regulators = NULL,
                             legend.position = "right") {
  # Format input data
  weightDT <- net.format(weightDT,
                         regulators = regulators)

  net <- igraph::graph_from_data_frame(weightDT[, c("regulator", "target", "weight", "Interaction")],
                                       directed = FALSE)

  layout <- igraph::layout_with_fr(net)
  rownames(layout) <- igraph::V(net)$name
  layout_ordered <- layout[igraph::V(net)$name,]
  regulatorNet <- ggnetwork(net,
                            layout = layout_ordered,
                            cell.jitter = 0)

  regulatorNet$isRegulator <- as.character(regulatorNet$name %in% regulators)
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  # Plot
  g <- ggplot() +
    geom_edges(data = regulatorNet,
               aes(x = x, y = y,
                   xend = xend, yend = yend,
                   size = weight,
                   color = Interaction),
               size = 0.75,
               curvature = 0.1,
               alpha = .6) +
    geom_nodes(data = regulatorNet[regulatorNet$isRegulator == "FALSE", ],
               aes(x = x, y = y),
               color = "darkgray",
               size = 3,
               alpha = .5) +
    geom_nodes(data = regulatorNet[regulatorNet$isRegulator == "TRUE", ],
               aes(x = x, y = y),
               color = "#8C4985",
               size = 6,
               alpha = .8) +
    scale_color_manual(values = cols) +
    geom_nodelabel_repel(data = regulatorNet[regulatorNet$isRegulator == "FALSE", ],
                         aes(x = x, y = y, label = name),
                         size = 2,
                         color = "#5A8BAD") +
    geom_nodelabel_repel(data = regulatorNet[regulatorNet$isRegulator == "TRUE", ],
                         aes(x = x, y = y, label = name),
                         size = 3.5,
                         color = "black") +
    theme_blank()
  g <- g + theme(legend.position = legend.position)
  return(g)
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

#' @title Compute and rank TFs in network
#'
#' @param weightDT The weight data table of network.
#' @param directedGraph If GRN is directed or not
#'
#' @return A data.table with three columns
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix)
#' ranks <- compute.gene.rank(weightDT)
#' head(ranks)
#'
compute.gene.rank <- function(weightDT,
                              directedGraph = FALSE) {
  colnames(weightDT) <- c("regulatory", "target", "weight")
  tfnet <- igraph::graph_from_data_frame(weightDT, directed = directedGraph)
  pageRank <- data.frame(igraph::page_rank(tfnet, directed = directedGraph)$vector)
  colnames(pageRank) <- c("pageRank")
  pageRank$gene <- rownames(pageRank)
  pageRank <- pageRank[, c("gene", "pageRank")]
  pageRank <- pageRank[order(pageRank$pageRank, decreasing = TRUE), ]
  pageRank$isRegulator <- FALSE
  pageRank$isRegulator[pageRank$gene %in% unique(weightDT$regulatory)] <- TRUE
  return(pageRank)
}

# import C++ compiled code
#' @useDynLib inferCSN
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom methods is
#' @import Matrix

#' @title Fit a sparse regression model
#' @description Computes the regularization path for the specified loss function and penalty function
#'
#' @param x The data matrix
#' @param y The response vector
#' @param loss The loss function. Currently support the choices "SquaredError" (for regression), "Logistic" (for logistic regression), and "SquaredHinge" (for smooth SVM).
#' @param penalty The type of regularization.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' @param maxSuppSize The maximum support size at which to terminate the regularization path. Recommend setting
#' this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small
#' portion of non-zeros
#' @param nLambda The number of Lambda values to select
#' @param nGamma The number of Gamma values to select
#' @param gammaMax The maximum value of Gamma when using the L0L2 penalty.
#' @param gammaMin The minimum value of Gamma when using the L0L2 penalty.
#' @param partialSort If TRUE partial sorting will be used for sorting the coordinates to do greedy cycling. Otherwise, full sorting is used.
#' @param maxIters The maximum number of iterations (full cycles) for CD per grid point.
#' @param rtol The relative tolerance which decides when to terminate optimization (based on the relative change in the objective between iterations).
#' @param atol The absolute tolerance which decides when to terminate optimization (based on the absolute L2 norm of the residuals).
#' @param activeSet If TRUE, performs active set updates.
#' @param activeSetNum The number of consecutive times a support should appear before declaring support stabilization.
#' @param maxSwaps The maximum number of swaps used by CDPSI for each grid point.
#' @param scaleDownFactor This parameter decides how close the selected Lambda values are.
#' @param screenSize The number of coordinates to cycle over when performing initial correlation screening.
#' @param autoLambda Ignored parameter. Kept for backwards compatibility.
#' @param lambdaGrid A grid of Lambda values to use in computing the regularization path.
#' @param excludeFirstK This parameter takes non-negative integers.
#' @param intercept If FALSE, no intercept term is included in the model.
#' @param lows Lower bounds for coefficients.
#' @param highs Upper bounds for coefficients.
#'
#' @return An S3 object of type "inferCSN" describing the regularization path
#' @export
#'
inferCSN.fit <- function(x, y,
                         loss = "SquaredError",
                         penalty = "L0",
                         algorithm = "CD",
                         maxSuppSize = 100,
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
  if ((rtol < 0) || (rtol >= 1)) stop("The specified rtol parameter must exist in [0, 1)")
  if (atol < 0) stop("The specified atol parameter must exist in [0, INF)")
  if (!(loss %in% c("SquaredError", "Logistic", "SquaredHinge"))) stop("The specified loss function is not supported.")
  if (!(penalty %in% c("L0", "L0L2", "L0L1"))) stop("The specified penalty is not supported.")
  if (!(algorithm %in% c("CD", "CDPSI"))) stop("The specified algorithm is not supported.")
  if (loss == "Logistic" | loss == "SquaredHinge") {
    if (dim(table(y)) != 2) {
      stop("Only binary classification is supported. Make sure y has only 2 unique values.")
    }
    y <- factor(y, labels = c(-1, 1)) # returns a vector of strings
    y <- as.numeric(levels(y))[y]

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

  # Handle Lambda Grids:
  if (length(lambdaGrid) != 0) {
    if (!is.null(autoLambda)) {
      warning("autoLambda is ignored and inferred if 'lambdaGrid' is supplied", call. = FALSE)
    }
    autoLambda <- FALSE
  } else {
    autoLambda <- TRUE
    lambdaGrid <- list(0)
  }

  if (penalty == "L0" && !autoLambda) {
    bad_lambdaGrid <- FALSE
    if (length(lambdaGrid) != 1) {
      bad_lambdaGrid <- TRUE
    }
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
                 Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
    }
  }

  if (penalty != "L0" && !autoLambda) {
    # Covers L0L1, L0L2 cases
    bad_lambdaGrid <- FALSE
    if (length(lambdaGrid) != nGamma) {
      warning("nGamma is ignored and replaced with length(lambdaGrid)", call. = FALSE)
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
      stop("L0L1 or L0L2 Penalty requires 'lambdaGrid' to be a list of length 'nGamma'.
                 Where lambdaGrid[[i]] is a list or vector of decreasing positive values.")
    }
  }

  is.scalar <- function(x) {
    is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
  }

  p <- dim(x)[[2]]

  withBounds <- FALSE

  if ((!identical(lows, -Inf)) || (!identical(highs, Inf))) {
    withBounds <- TRUE

    if (algorithm == "CDPSI") {
      if (any(lows != -Inf) || any(highs != Inf)) {
        stop("Bounds are not YET supported for CDPSI algorithm.")
      }
    }

    if (is.scalar(lows)) {
      lows <- lows * rep(1, p)
    } else if (!all(sapply(lows, is.scalar)) || length(lows) != p) {
      stop("Lows must be a vector of real values of length p")
    }

    if (is.scalar(highs)) {
      highs <- highs * rep(1, p)
    } else if (!all(sapply(highs, is.scalar)) || length(highs) != p) {
      stop("Highs must be a vector of real values of length p")
    }

    if (any(lows >= highs) || any(lows > 0) || any(highs < 0)) {
      stop("Bounds must conform to the following conditions: Lows <= 0, Highs >= 0, Lows < Highs")
    }
  }

  M <- list()
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

  settings <- list()
  settings[[1]] <- intercept # Settings only contains intercept for now. Might include additional elements later.
  names(settings) <- c("intercept")

  # Find potential support sizes exceeding maxSuppSize and remove them (this is due to
  # the C++ core whose last solution can exceed maxSuppSize
  for (i in 1:length(M$SuppSize)) {
    last <- length(M$SuppSize[[i]])
    if (M$SuppSize[[i]][last] > maxSuppSize) {
      if (last == 1) {
        warning("Warning! Only 1 element in path with support size > maxSuppSize. \n
                Try increasing maxSuppSize to resolve the issue.")
      } else {
        M$SuppSize[[i]] <- M$SuppSize[[i]][-last]
        M$Converged[[i]] <- M$Converged[[i]][-last]
        M$lambda[[i]] <- M$lambda[[i]][-last]
        M$a0[[i]] <- M$a0[[i]][-last]
        M$beta[[i]] <- as(M$beta[[i]][, -last], "sparseMatrix")
      }
    }
  }

  G <- list(beta = M$beta,
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
  G$varnames <- varnames
  class(G) <- "inferCSN"
  G$n <- dim(x)[1]
  G$p <- dim(x)[2]
  G
}

#' @title Computes a regularization path and performs K-fold cross-validation
#'
#' @inheritParams inferCSN.fit
#' @param nFolds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation
#'
#' @return An S3 object of type "inferCSNCV" describing the regularization path
#' @export
#'
inferCSN.cvfit <- function(x, y,
                           loss = "SquaredError",
                           penalty = "L0",
                           algorithm = "CD",
                           maxSuppSize = 100,
                           nLambda = 100,
                           nGamma = 10,
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
                           nFolds = 10,
                           seed = 1,
                           excludeFirstK = 0,
                           intercept = TRUE,
                           lows = -Inf,
                           highs = Inf) {
  set.seed(seed)
  if ((rtol < 0) || (rtol >= 1)) stop("The specified rtol parameter must exist in [0, 1)")
  if (atol < 0) stop("The specified atol parameter must exist in [0, INF)")
  if (!(loss %in% c("SquaredError", "Logistic", "SquaredHinge"))) stop("The specified loss function is not supported.")
  if (!(penalty %in% c("L0", "L0L2", "L0L1"))) stop("The specified penalty is not supported.")
  if (!(algorithm %in% c("CD", "CDPSI"))) stop("The specified algorithm is not supported.")
  if (loss == "Logistic" | loss == "SquaredHinge") {
    if (dim(table(y)) != 2) {
      stop("Only binary classification is supported. Make sure y has only 2 unique values.")
    }
    y <- factor(y, labels = c(-1, 1)) # Returns a vector of strings
    y <- as.numeric(levels(y))[y]

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

  # Handle Lambda Grids:
  if (length(lambdaGrid) != 0) {
    if (!is.null(autoLambda) && !autoLambda) {
      warning("autoLambda is ignored and inferred if 'lambdaGrid' is supplied", call. = FALSE)
    }
    autoLambda <- FALSE
  } else {
    autoLambda <- TRUE
    lambdaGrid <- list(0)
  }

  if (penalty == "L0" && !autoLambda) {
    bad_lambdaGrid <- FALSE
    if (length(lambdaGrid) != 1) {
      bad_lambdaGrid <- TRUE
    }
    current <- Inf
    for (nxt in lambdaGrid[[1]]) {
      if (nxt > current) {
        # This must be > instead of >= to allow first iteration L0L1 lambdas of all 0s to be valid
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
                 Where lambdaGrid[[1]] is a list or vector of decreasing positive values.")
    }
  }

  if (penalty != "L0" && !autoLambda) {
    # Covers L0L1, L0L2 cases
    bad_lambdaGrid <- FALSE
    if (length(lambdaGrid) != nGamma) {
      warning("nGamma is ignored and replaced with length(lambdaGrid)", call. = FALSE)
      nGamma <- length(lambdaGrid)
    }

    for (i in 1:length(lambdaGrid)) {
      current <- Inf
      for (nxt in lambdaGrid[[i]]) {
        if (nxt > current) {
          # This must be > instead of >= to allow first iteration L0L1 lambdas of all 0s to be valid
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
      stop("L0L1 or L0L2 Penalty requires 'lambdaGrid' to be a list of length 'nGamma'.
                 Where lambdaGrid[[i]] is a list or vector of decreasing positive values.")
    }
  }

  is.scalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)

  p <- dim(x)[[2]]

  withBounds <- FALSE
  if ((!identical(lows, -Inf)) || (!identical(highs, Inf))) {
    withBounds <- TRUE

    if (algorithm == "CDPSI") {
      if (any(lows != -Inf) || any(highs != Inf)) {
        stop("Bounds are not YET supported for CDPSI algorithm.")
      }
    }

    if (is.scalar(lows)) {
      lows <- lows * rep(1, p)
    } else if (!all(sapply(lows, is.scalar)) || length(lows) != p) {
      stop("Lows must be a vector of real values of length p")
    }

    if (is.scalar(highs)) {
      highs <- highs * rep(1, p)
    } else if (!all(sapply(highs, is.scalar)) || length(highs) != p) {
      stop("Highs must be a vector of real values of length p")
    }

    if (any(lows >= highs) || any(lows > 0) || any(highs < 0)) {
      stop("Bounds must conform to the following conditions: Lows <= 0, Highs >= 0, Lows < Highs")
    }
  }

  M <- list()
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

  settings <- list()
  settings[[1]] <- intercept
  names(settings) <- c("intercept")

  # The C++ core whose last solution can exceed maxSuppSize
  for (i in 1:length(M$SuppSize)) {
    last <- length(M$SuppSize[[i]])
    if (M$SuppSize[[i]][last] > maxSuppSize) {
      if (last == 1) {
        warning("Warning! Only 1 element in path with support size > maxSuppSize. \n
                Try increasing maxSuppSize to resolve the issue......")
      } else {
        M$SuppSize[[i]] <- M$SuppSize[[i]][-last]
        M$Converged[[i]] <- M$Converged[[i]][-last]
        M$lambda[[i]] <- M$lambda[[i]][-last]
        M$a0[[i]] <- M$a0[[i]][-last]
        M$beta[[i]] <- as(M$beta[[i]][, -last], "sparseMatrix") # conversion to sparseMatrix is necessary to handle the case of a single column
        M$CVMeans[[i]] <- M$CVMeans[[i]][-last]
        M$CVSDs[[i]] <- M$CVSDs[[i]][-last]
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
  G <- list(fit = fit, cvMeans = M$CVMeans, cvSDs = M$CVSDs)
  class(G) <- "inferCSNCV"
  G
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
#' @method coef inferCSNCV
#' @export
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
#' @method print inferCSNCV
#' @export
#'
print.inferCSNCV <- function(x, ...) {
  print.inferCSN(x$fit)
}
