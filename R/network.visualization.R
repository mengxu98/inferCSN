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
#' @return Return a heatmap of ggplot2 object
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
    if (is.null(heatmapSize)) heatmapSize <- 6
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
#' g <- dynamic.networks(weightDT, regulators = weightDT[1, 1])
#' g
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
    theme_blank() +
    theme(legend.position = legend.position)
  return(g)
}

#' @title contrast.networks
#'
#' @param weightDT weightDT
#' @param degreeValue degreeValue
#' @param weightValue weightValue
#' @param legend.position legend.position
#'
#' @import ggplot2
#' @import ggraph
#'
#' @return Return a ggplot2 object
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix)
#' g <- contrast.networks(weightDT[1:50, ])
#' g
#'
contrast.networks <- function(weightDT,
                              degreeValue = 0,
                              weightValue = 0,
                              legend.position = "bottom") {
  weightDT <- net.format(weightDT)

  graph <- tidygraph::as_tbl_graph(weightDT)
  graph <- dplyr::mutate(graph, degree = tidygraph::centrality_degree(mode = 'out'))
  graph <- dplyr::filter(graph, degree > degreeValue)
  graph <- tidygraph::activate(graph, edges)

  g <- ggraph(graph, layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(colour = Interaction,
                      filter = weight > weightValue,
                      edge_width = weight),
                  arrow = arrow(length = unit(3, 'mm')),
                  start_cap = square(3, 'mm'),
                  end_cap = circle(3, 'mm')) +
    scale_edge_width(range=c(0, 1)) +
    facet_edges(~Interaction) +
    geom_node_point(aes(size = degree), colour = '#A1B7CE') +
    geom_node_text(aes(label = name), repel = TRUE) +
    coord_fixed() +
    theme_graph(base_family = "serif",
                foreground = 'steelblue',
                fg_text_colour = 'white') +
    theme(legend.position = legend.position)

  return(g)
}
