#' @title The heatmap of network
#'
#' @inheritParams net.format
#' @param switch_matrix Logical value, whether to weight data table to matrix.
#' @param show_names Logical value, whether to show names of row and column.
#' @param heatmap_size The size of heatmap, default set to 5.
#' @param heatmap_height The height of heatmap.
#' @param heatmap_width The width of heatmap.
#' @param heatmap_title The title of heatmap.
#' @param heatmap_color Colors of heatmap.
#' @param legend_name The name of legend.
#' @param row_title The title of row.
#'
#' @return Return a heatmap of ggplot2 object
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' data("example_ground_truth")
#' weight_table <- inferCSN(example_matrix)
#'
#' p1 <- network.heatmap(
#'   example_ground_truth[, 1:3],
#'   heatmap_title = "Ground truth",
#'   legend_name = "Ground truth"
#' )
#' p2 <- network.heatmap(
#'   weight_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "inferCSN"
#' )
#' ComplexHeatmap::draw(p1 + p2)
#'
#' p3 <- network.heatmap(
#'   weight_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "Weight1",
#'   heatmap_color = c("#20a485", "#410054", "#fee81f")
#' )
#' p4 <- network.heatmap(
#'   weight_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "Weight2",
#'   heatmap_color = c("#20a485", "white", "#fee81f")
#' )
#' ComplexHeatmap::draw(p3 + p4)
#'
#' network.heatmap(
#'   weight_table,
#'   heatmap_title = "inferCSN",
#'   show_names = TRUE
#' )
#'
#' network.heatmap(
#'   weight_table,
#'   regulators = c("g1", "g2"),
#'   heatmap_title = "inferCSN",
#'   show_names = TRUE
#' )
#'
#' network.heatmap(
#'   weight_table,
#'   targets = c("g1", "g2"),
#'   heatmap_title = "inferCSN",
#'   show_names = TRUE
#' )
#'
#' network.heatmap(
#'   weight_table,
#'   regulators = c("g1", "g3", "g5"),
#'   targets = c("g3", "g6", "g9"),
#'   heatmap_title = "inferCSN",
#'   show_names = TRUE
#' )
network.heatmap <- function(
    weight_table,
    regulators = NULL,
    targets = NULL,
    switch_matrix = TRUE,
    show_names = FALSE,
    heatmap_size = 5,
    heatmap_height = NULL,
    heatmap_width = NULL,
    heatmap_title = NULL,
    heatmap_color = c("#1966ad", "white", "#bb141a"),
    legend_name = "Weight",
    row_title = "Regulators",
    abs_weight = FALSE) {
  if (switch_matrix) {
    weight_table <- net.format(
      weight_table,
      regulators = regulators,
      targets = targets,
      abs_weight = FALSE
    )[, 1:3]

    regulators <- weight_table$regulator
    targets <- weight_table$target

    weight_matrix <- table.to.matrix(weight_table)
  } else {
    if (is.null(regulators)) {
      regulators <- rownames(weight_matrix)
    }
    if (is.null(targets)) {
      targets <- colnames(weight_matrix)
    }

    weight_matrix <- weight_table
  }
  weight_matrix[is.na(weight_matrix)] <- 0

  unique_regulators <- gtools::mixedsort(unique(regulators))
  unique_targets <- gtools::mixedsort(unique(targets))
  weight_matrix <- weight_matrix[unique_regulators, unique_targets]

  if (show_names) {
    if (is.null(heatmap_height) || is.null(heatmap_width)) {
      heatmap_height <- length(unique_regulators) / 2
      heatmap_width <- length(unique_targets) / 2
    }
  } else {
    if (is.null(heatmap_height) || is.null(heatmap_width)) {
      heatmap_height <- heatmap_size * length(unique_regulators) / length(unique_targets)
      heatmap_width <- heatmap_size
    }
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
    row_title = row_title,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_names,
    show_column_names = show_names,
    column_names_rot = 45,
    width = unit(heatmap_width, "cm"),
    height = unit(heatmap_height, "cm"),
    border = "gray"
  )

  return(p)
}

#' @title Plot of dynamic networks
#'
#' @inheritParams net.format
#' @param legend_position The position of legend.
#'
#' @import ggplot2
#' @import ggnetwork
#'
#' @return A list of ggplot2 objects
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' dynamic.networks(
#'   weight_table,
#'   regulators = weight_table[1, 1]
#' )
#' dynamic.networks(
#'   weight_table,
#'   targets = weight_table[1, 1]
#' )
#' dynamic.networks(
#'   weight_table,
#'   regulators = weight_table[1, 1],
#'   targets = weight_table[1, 2]
#' )
dynamic.networks <- function(
    weight_table,
    regulators = NULL,
    targets = NULL,
    legend_position = "right") {
  # Format input data
  weight_table <- net.format(
    weight_table,
    regulators = regulators,
    targets = targets
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
    theme(legend.position = legend_position)

  return(g)
}

#' @title contrast.networks
#'
#' @inheritParams dynamic.networks
#' @param degree_value degree_value
#' @param weight_value weight_value
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
#' contrast.networks(weight_table[1:50, ])
contrast.networks <- function(
    weight_table,
    degree_value = 0,
    weight_value = 0,
    legend_position = "bottom") {
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
    theme(legend.position = legend_position)

  return(g)
}
