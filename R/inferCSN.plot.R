#' @title Plot of dynamic networks
#' @description Plot
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
#'    data("exampleDataMatrix")
#'    weightList <- inferCSN(exampleDataMatrix)
#'    inferCSN.plot.dynamic.networks(weightList)
inferCSN.plot.dynamic.networks <- function(weightList,
                                           regulators = NULL,
                                           onlyRegulators = TRUE,
                                           order = NULL,
                                           thresh = NULL,
                                           legend.position = "right") {
  # Format input data
  df <- net.format(
    weightList,
    regulators = regulators
  )

  net <- igraph::graph_from_data_frame(
    df[, c("regulator", "target", "Interaction")],
    directed = FALSE
  )
  layout <- igraph::layout_with_fr(net)
  rownames(layout) <- igraph::V(net)$name
  layout_ordered <- layout[igraph::V(net)$name, ]
  regulatornet <- ggnetwork::ggnetwork(net,
                                       layout = layout_ordered,
                                       cell.jitter = 0
  )
  regulatornet$is_regulator <- as.character(regulatornet$name %in% regulators)
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")
  g <- ggplot2::ggplot() +
    ggnetwork::geom_edges(
      data = regulatornet,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = Interaction),
      size = 0.75,
      curvature = 0.1,
      alpha = .6
    ) +
    ggnetwork::geom_nodes(
      data = regulatornet[regulatornet$is_regulator == "FALSE", ],
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "darkgray",
      size = 3,
      alpha = .5
    ) +
    ggnetwork::geom_nodes(
      data = regulatornet[regulatornet$is_regulator == "TRUE", ],
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "#8C4985",
      size = 6,
      alpha = .8
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggnetwork::geom_nodelabel_repel(
      data = regulatornet[regulatornet$is_regulator == "FALSE", ],
      ggplot2::aes(x = x, y = y, label = name),
      size = 2,
      color = "#5A8BAD"
    ) +
    ggnetwork::geom_nodelabel_repel(
      data = regulatornet[regulatornet$is_regulator == "TRUE", ],
      ggplot2::aes(x = x, y = y, label = name),
      size = 3.5,
      color = "black"
    ) +
    ggnetwork::theme_blank()
  # ggtitle(names(weightList)[i])
  g <- g + ggplot2::theme(legend.position = legend.position)
}

#' net.format
#'
#' @param weightList weightList of GRN
#' @param regulators Regulators list
#'
#' @return A formated weight list
#' @export
#'
net.format <- function(weightList,
                       regulators = regulators) {
  colnames(weightList) <- c("regulator", "target", "weight")
  if (!is.null(regulators)) {
    weightListNew <- purrr::map_dfr(
      regulators, function(x) {
        weight <- weightList[which(weightList$regulator == x), ]
      }
    )
    weightList <- weightListNew
  }
  weightList$weight <- as.numeric(weightList$weight)
  weightList$Interaction <- "Activation"
  weightList$Interaction[weightList$weight < 0] <- "Repression"
  return(weightList)
}

#' inferCSN.plot
#' @description Plot
#'
#' @param data A long data table
#' @param plotType boxplot
#'
#' @importFrom graphics "boxplot"
#' @return A ggplot2 object
#' @export
#'
inferCSN.plot <- function(data, plotType = NULL) {
  if (is.null(plotType)) {
    plotType <- boxplot
  }
  if (plotType == boxplot) {
    p <- ggplot(
      data,
      aes(x = Method, y = AUPRC)
    ) +
      # geom_violin(aes(fill = Method),
      #     trim = FALSE
      # ) +
      geom_boxplot(aes(fill = Method),
                   width = 0.8
      ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        comparisons = my_comparisons,
        bracket.size = 0.6,
        sizen = 4,
        color = "#6699cc"
      ) +
      scale_fill_manual(values = mycol) +
      # scale_color_manual(values = mycol) +
      scale_x_discrete(labels = methods) +
      labs(x = "Methods", y = "AUPRC") +
      theme(legend.position = "bottom") +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 10
        )
      )
  }
  p
}
