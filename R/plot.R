globalVariables(c("x", "y", "xend", "yend", "weight", "Interaction", "name"))

#' @title dynamic.networks
#' @description Plot of dynamic networks
#'
#' @param weightList weightList
#' @param regulators regulators
#' @param order order = NULL
#' @param thresh thresh = NULL
#' @param onlyRegulators onlyregulators
#' @param legend.position legend.position
#'
#' @return A list of ggplot2 objects
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' weightList <- inferCSN(exampleDataMatrix)
#' p <- dynamic.networks(weightList)
#' p
dynamic.networks <- function(weightList,
                             regulators = NULL,
                             onlyRegulators = TRUE,
                             order = NULL,
                             thresh = NULL,
                             legend.position = "right") {
  # Format input data
  weightList <- net.format(
    weightList,
    regulators = regulators
  )

  net <- igraph::graph_from_data_frame(
    weightList[, c("regulator", "target", "weight", "Interaction")],
    directed = FALSE
  )

  layout <- igraph::layout_with_fr(net)
  rownames(layout) <- igraph::V(net)$name
  layout_ordered <- layout[igraph::V(net)$name,]
  regulatorNet <- ggnetwork::ggnetwork(
    net,
    layout = layout_ordered,
    cell.jitter = 0
  )
  regulatorNet$is_regulator <- as.character(regulatorNet$name %in% regulators)
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  # Plot
  g <- ggplot2::ggplot() +
    ggnetwork::geom_edges(
      data = regulatorNet,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, size = weight, color = Interaction),
      size = 0.75,
      curvature = 0.1,
      alpha = .6
    ) +
    ggnetwork::geom_nodes(
      data = regulatorNet[regulatorNet$is_regulator == "FALSE", ],
      ggplot2::aes(x = x, y = y), # , xend = xend, yend = yend
      color = "darkgray",
      size = 3,
      alpha = .5
    ) +
    ggnetwork::geom_nodes(
      data = regulatorNet[regulatorNet$is_regulator == "TRUE", ],
      ggplot2::aes(x = x, y = y), # , xend = xend, yend = yend
      color = "#8C4985",
      size = 6,
      alpha = .8
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggnetwork::geom_nodelabel_repel(
      data = regulatorNet[regulatorNet$is_regulator == "FALSE", ],
      ggplot2::aes(x = x, y = y, label = name),
      size = 2,
      color = "#5A8BAD"
    ) +
    ggnetwork::geom_nodelabel_repel(
      data = regulatorNet[regulatorNet$is_regulator == "TRUE", ],
      ggplot2::aes(x = x, y = y, label = name),
      size = 3.5,
      color = "black"
    ) +
    ggnetwork::theme_blank()
  # ggtitle(names(weightList)[i])
  g <- g + ggplot2::theme(legend.position = legend.position)
}

#' @title net.format
#' @description Format weight table
#'
#' @param weightList weightList of CSN
#' @param regulators Regulators list
#'
#' @return Format weight table
#' @export
#'
net.format <- function(weightList,
                       regulators = NULL) {
  colnames(weightList) <- c("regulator", "target", "weight")
  if (!is.null(regulators)) {
    weightList <- purrr::map_dfr(
      regulators, function(x) {
        weight <- weightList[which(weightList$regulator == x), ]
      }
    )
  }
  weightList$weight <- as.numeric(weightList$weight)
  weightList$Interaction <- "Activation"
  weightList$Interaction[weightList$weight < 0] <- "Repression"
  weightList$weight <- abs(weightList$weight)
  return(weightList)
}
