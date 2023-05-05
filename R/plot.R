globalVariables(c("x", "y", "xend", "yend", "weight", "Interaction", "name"))

#' @title dynamic.networks
#' @description Plot of dynamic networks
#'
#' @param weightDT weightDT
#' @param regulators regulators
#' @param legend.position legend.position
#'
#' @return A list of ggplot2 objects
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' weightDT <- inferCSN(exampleDataMatrix)
#' p <- dynamic.networks(weightDT)
#' p
dynamic.networks <- function(weightDT,
                             regulators = NULL,
                             legend.position = "right") {
  # Format input data
  weightDT <- net.format(
    weightDT,
    regulators = regulators
  )

  net <- igraph::graph_from_data_frame(
    weightDT[, c("regulator", "target", "weight", "Interaction")],
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
      ggplot2::aes(x = x, y = y),
      color = "darkgray",
      size = 3,
      alpha = .5
    ) +
    ggnetwork::geom_nodes(
      data = regulatorNet[regulatorNet$is_regulator == "TRUE", ],
      ggplot2::aes(x = x, y = y),
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
  # ggtitle(names(weightDT)[i])
  g <- g + ggplot2::theme(legend.position = legend.position)
}

#' @title net.format
#' @description Format weight table
#'
#' @param weightDT The weight data table of network.
#' @param regulators Regulators list.
#'
#' @return Format weight table.
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
