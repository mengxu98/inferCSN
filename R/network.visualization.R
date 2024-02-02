#' @title The heatmap of network
#'
#' @param weight_table The weight data table of network
#' @param switch_watrix switch_watrix
#' @param heatmap_size heatmap_size
#' @param heatmap_title heatmap_title
#' @param heatmap_color heatmap_color
#' @param show_names show_names
#' @param legend_name legend_name
#'
#' @return Return a heatmap of ggplot2 object
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' data("example_ground_truth")
#' weight_table <- inferCSN(example_matrix)
#' p1 <- network.heatmap(example_ground_truth,
#'   heatmap_title = "Ground truth"
#' )
#'
#' p2 <- network.heatmap(weight_table,
#'   legend_name = "Weight2",
#'   heatmap_title = "inferCSN"
#' )
#'
#' ComplexHeatmap::draw(p1 + p2)
#'
#' p3 <- network.heatmap(weight_table,
#'   heatmap_title = "inferCSN",
#'   heatmap_color = c("#20a485", "#410054", "#fee81f")
#' )
#'
#' p4 <- network.heatmap(weight_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "Weight2",
#'   heatmap_color = c("#20a485", "white", "#fee81f")
#' )
#'
#' ComplexHeatmap::draw(p3 + p4)
#'
#' p5 <- network.heatmap(
#'   weight_table,
#'   heatmap_title = "inferCSN",
#'   show_names = TRUE
#' )
#' p5
network.heatmap <- function(
    weight_table,
    switch_watrix = TRUE,
    heatmap_size = NULL,
    heatmap_title = NULL,
    heatmap_color = NULL,
    show_names = FALSE,
    legend_name = NULL) {
  if (switch_watrix) {
    colnames(weight_table) <- c("regulator", "target", "weight")
    genes <- c(weight_table$regulator, weight_table$target)
    weight_matrix <- table.to.matrix(weight_table)
  } else {
    genes <- c(rownames(weight_table), colnames(weight_table))
    weight_matrix <- weight_table
  }
  genes <- gtools::mixedsort(unique(genes))
  weight_matrix <- weight_matrix[genes, genes]

  if (is.null(legend_name)) legend_name <- "Weight"

  if (is.null(heatmap_color)) heatmap_color <- c("#1966ad", "white", "#bb141a")

  if (show_names) {
    if (is.null(heatmap_size)) heatmap_size <- length(genes) / 2
  } else {
    if (is.null(heatmap_size)) heatmap_size <- 6
  }

  min_weight <- min(weight_matrix)
  max_weight <- max(weight_matrix)
  if (min_weight >= 0) {
    color_function <- circlize::colorRamp2(
      c(min_weight, max_weight),
      heatmap_color[-1]
    )
  } else if (max_weight <= 0) {
    color_function <- circlize::colorRamp2(
      c(min_weight, max_weight),
      heatmap_color[-3]
    )
  } else {
    color_function <- circlize::colorRamp2(
      c(min_weight, 0, max_weight),
      heatmap_color
    )
  }

  p <- ComplexHeatmap::Heatmap(
    weight_matrix,
    name = legend_name,
    col = color_function,
    column_title = heatmap_title,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_names,
    show_column_names = show_names,
    width = unit(heatmap_size, "cm"),
    height = unit(heatmap_size, "cm"),
    border = "black"
  )
  return(p)
}

#' @title Plot of dynamic networks
#'
#' @param weight_table weight_table
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
#' \dontrun{
#' # because `ggnetwork` package may be cause error
#' # you can install `igraph 1.6.0` or lower to deal error
#' # or waitting a new version of `ggnetwork` package
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' g <- dynamic.networks(weight_table, regulators = weight_table[1, 1])
#' g
#' }
dynamic.networks <- function(
    weight_table,
    regulators = NULL,
    legend.position = "right") {
  # Format input data
  weight_table <- net.format(
    weight_table,
    regulators = regulators
  )

  net <- igraph::graph_from_data_frame(
    weight_table[, c("regulator", "target", "weight", "Interaction")],
    directed = FALSE
  )

  layout <- igraph::layout_with_fr(net)
  rownames(layout) <- igraph::V(net)$name
  layout_ordered <- layout[igraph::V(net)$name, ]
  regulator_network <- ggnetwork(
    net,
    layout = layout_ordered,
    cell.jitter = 0
  )

  regulator_network$is_regulator <- as.character(
    regulator_network$name %in% regulators
  )
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  # Plot
  g <- ggplot() +
    geom_edges(
      data = regulator_network,
      aes(
        x = x, y = y,
        xend = xend, yend = yend,
        size = weight,
        color = Interaction
      ),
      size = 0.75,
      curvature = 0.1,
      alpha = .6
    ) +
    geom_nodes(
      data = regulator_network[regulator_network$is_regulator == "FALSE", ],
      aes(x = x, y = y),
      color = "darkgray",
      size = 3,
      alpha = .5
    ) +
    geom_nodes(
      data = regulator_network[regulator_network$is_regulator == "TRUE", ],
      aes(x = x, y = y),
      color = "#8C4985",
      size = 6,
      alpha = .8
    ) +
    scale_color_manual(values = cols) +
    geom_nodelabel_repel(
      data = regulator_network[regulator_network$is_regulator == "FALSE", ],
      aes(x = x, y = y, label = name),
      size = 2,
      color = "#5A8BAD"
    ) +
    geom_nodelabel_repel(
      data = regulator_network[regulator_network$is_regulator == "TRUE", ],
      aes(x = x, y = y, label = name),
      size = 3.5,
      color = "black"
    ) +
    theme_blank() +
    theme(legend.position = legend.position)
  return(g)
}

#' @title contrast.networks
#'
#' @param weight_table weight_table
#' @param degree_value degree_value
#' @param weight_value weight_value
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
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' g <- contrast.networks(weight_table[1:50, ])
#' g
contrast.networks <- function(
    weight_table,
    degree_value = 0,
    weight_value = 0,
    legend.position = "bottom") {
  weight_table <- net.format(weight_table)

  graph <- tidygraph::as_tbl_graph(weight_table)
  graph <- dplyr::mutate(
    graph,
    degree = tidygraph::centrality_degree(mode = "out")
  )
  graph <- dplyr::filter(graph, degree > degree_value)
  graph <- tidygraph::activate(graph, edges)

  g <- ggraph(graph, layout = "linear", circular = TRUE) +
    geom_edge_arc(
      aes(
        colour = Interaction,
        filter = weight > weight_value,
        edge_width = weight
      ),
      arrow = arrow(length = unit(3, "mm")),
      start_cap = square(3, "mm"),
      end_cap = circle(3, "mm")
    ) +
    scale_edge_width(range = c(0, 1)) +
    facet_edges(~Interaction) +
    geom_node_point(aes(size = degree), colour = "#A1B7CE") +
    geom_node_text(aes(label = name), repel = TRUE) +
    coord_fixed() +
    theme_graph(
      base_family = "serif",
      foreground = "steelblue",
      fg_text_colour = "white"
    ) +
    theme(legend.position = legend.position)

  return(g)
}
