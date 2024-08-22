#' @title Plot network heatmap
#'
#' @inheritParams network_format
#' @param switch_matrix Logical value, default is *`TRUE`*, whether to weight data table to matrix.
#' @param show_names Logical value, default is *`FALSE`*, whether to show names of row and column.
#' @param heatmap_size_lock Lock the size of heatmap.
#' @param heatmap_size Default is *`5`*. The size of heatmap.
#' @param heatmap_height The height of heatmap.
#' @param heatmap_width The width of heatmap.
#' @param heatmap_title The title of heatmap.
#' @param heatmap_color Colors of heatmap.
#' @param border_color Default is *`gray`*. Color of heatmap border.
#' @param rect_color Default is *`NA`*. Color of heatmap rect.
#' @param anno_width Width of annotation.
#' @param anno_height Height of annotation.
#' @param row_anno_type Default is *`NULL`*,
#' could add a annotation plot to row,
#' choose one of *`boxplot`*, *`barplot`*, *`histogram`*, *`density`*, *`lines`*, *`points`*, and *`horizon`*.
#' @param column_anno_type Default is *`NULL`*,
#' could add a annotation plot to column,
#' choose one of *`boxplot`*, *`barplot`*, *`histogram`*, *`density`*, *`lines`*, and *`points`*.
#' @param legend_name The name of legend.
#' @param row_title The title of row.
#'
#' @md
#' @return A heatmap
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#'
#' p1 <- plot_network_heatmap(
#'   example_ground_truth[, 1:3],
#'   heatmap_title = "Ground truth",
#'   legend_name = "Ground truth"
#' )
#' p2 <- plot_network_heatmap(
#'   network_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "inferCSN"
#' )
#' ComplexHeatmap::draw(p1 + p2)
#'
#' p3 <- plot_network_heatmap(
#'   network_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "Weight1",
#'   heatmap_color = c("#20a485", "#410054", "#fee81f")
#' )
#' p4 <- plot_network_heatmap(
#'   network_table,
#'   heatmap_title = "inferCSN",
#'   legend_name = "Weight2",
#'   heatmap_color = c("#20a485", "white", "#fee81f")
#' )
#' ComplexHeatmap::draw(p3 + p4)
#'
#' plot_network_heatmap(
#'   network_table,
#'   show_names = TRUE,
#'   rect_color = "gray90",
#'   row_anno_type = "density",
#'   column_anno_type = "barplot"
#' )
#'
#' plot_network_heatmap(
#'   network_table,
#'   regulators = c("g1", "g2"),
#'   show_names = TRUE
#' )
#'
#' plot_network_heatmap(
#'   network_table,
#'   targets = c("g1", "g2"),
#'   row_anno_type = "boxplot",
#'   column_anno_type = "histogram",
#'   show_names = TRUE
#' )
#'
#' plot_network_heatmap(
#'   network_table,
#'   regulators = c("g1", "g3", "g5"),
#'   targets = c("g3", "g6", "g9"),
#'   show_names = TRUE
#' )
plot_network_heatmap <- function(
    network_table,
    regulators = NULL,
    targets = NULL,
    switch_matrix = TRUE,
    show_names = FALSE,
    heatmap_size_lock = TRUE,
    heatmap_size = 5,
    heatmap_height = NULL,
    heatmap_width = NULL,
    heatmap_title = NULL,
    heatmap_color = c("#1966ad", "white", "#bb141a"),
    border_color = "gray",
    rect_color = NA,
    anno_width = 1,
    anno_height = 1,
    row_anno_type = NULL,
    column_anno_type = NULL,
    legend_name = "Weight",
    row_title = "Regulators") {
  if (switch_matrix) {
    weight_matrix <- table_to_matrix(
      network_table,
      regulators = regulators,
      targets = targets
    )
  } else {
    weight_matrix <- network_table
  }
  weight_matrix <- filter_sort_matrix(
    weight_matrix,
    regulators = regulators,
    targets = targets
  )

  unique_regulators <- rownames(weight_matrix)
  unique_targets <- colnames(weight_matrix)

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

  if (!is.null(row_anno_type)) {
    row_anno_type <- match.arg(
      row_anno_type,
      c("boxplot", "barplot", "histogram", "density", "lines", "points", "horizon")
    )
    row_anno <- switch(row_anno_type,
      "boxplot" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_boxplot(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      ),
      "barplot" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_barplot(
          abs(weight_matrix),
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:ncol(weight_matrix))
        )
      ),
      "histogram" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_histogram(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      ),
      "density" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_density(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      ),
      "lines" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_lines(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      ),
      "points" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_points(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      ),
      "horizon" = ComplexHeatmap::rowAnnotation(
        Anno = ComplexHeatmap::anno_horizon(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      )
    )
  } else {
    row_anno <- NULL
  }

  if (!is.null(column_anno_type)) {
    column_anno_type <- match.arg(
      column_anno_type,
      c("boxplot", "barplot", "histogram", "density", "lines", "points")
    )
    column_anno <- switch(column_anno_type,
      "boxplot" = ComplexHeatmap::columnAnnotation(
        Anno = ComplexHeatmap::anno_boxplot(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = 1:ncol(weight_matrix))
        )
      ),
      "barplot" = ComplexHeatmap::columnAnnotation(
        Anno = ComplexHeatmap::anno_barplot(
          abs(weight_matrix),
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = 1:nrow(weight_matrix))
        )
      ),
      "histogram" = ComplexHeatmap::columnAnnotation(
        Anno = ComplexHeatmap::anno_histogram(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = 1:ncol(weight_matrix))
        )
      ),
      "density" = ComplexHeatmap::columnAnnotation(
        Anno = ComplexHeatmap::anno_density(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = 1:ncol(weight_matrix))
        )
      ),
      "lines" = ComplexHeatmap::columnAnnotation(
        Anno = ComplexHeatmap::anno_lines(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = 1:ncol(weight_matrix))
        )
      ),
      "points" = ComplexHeatmap::columnAnnotation(
        Anno = ComplexHeatmap::anno_points(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = 1:ncol(weight_matrix))
        )
      )
    )
  } else {
    column_anno <- NULL
  }

  if (heatmap_size_lock) {
    width <- grid::unit(heatmap_width, "cm")
    height <- grid::unit(heatmap_height, "cm")
  } else {
    width <- NULL
    height <- NULL
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
    border = border_color,
    rect_gp = grid::gpar(col = rect_color),
    width = width,
    height = height,
    top_annotation = column_anno,
    left_annotation = row_anno
  )

  return(p)
}

#' @title Plot dynamic networks
#'
#' @inheritParams network_format
#' @param legend_position The position of legend.
#'
#' @import ggnetwork
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' plot_static_networks(
#'   network_table,
#'   regulators = network_table[1, 1]
#' )
#' plot_static_networks(
#'   network_table,
#'   targets = network_table[1, 1]
#' )
#' plot_static_networks(
#'   network_table,
#'   regulators = network_table[1, 1],
#'   targets = network_table[1, 2]
#' )
plot_static_networks <- function(
    network_table,
    regulators = NULL,
    targets = NULL,
    legend_position = "right") {
  network_table <- network_format(
    network_table,
    regulators = regulators,
    targets = targets
  )

  net <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "weight", "Interaction")],
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

#' @title Plot contrast networks
#'
#' @inheritParams plot_static_networks
#' @param degree_value Degree value to filter nodes.
#' @param weight_value Weight value to filter edges.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' plot_contrast_networks(network_table[1:50, ])
plot_contrast_networks <- function(
    network_table,
    degree_value = 0,
    weight_value = 0,
    legend_position = "bottom") {
  network_table <- network_format(network_table)

  graph <- tidygraph::as_tbl_graph(network_table)
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

#' @title Plot dynamic networks
#'
#' @inheritParams network_format
#' @param celltypes_order The order of cell types.
#' @param ntop The number of top genes to plot.
#' @param title The title of figure.
#' @param theme_type Default is \code{theme_void}, the theme of figure,
#' could be \code{theme_void}, \code{theme_blank} or \code{theme_facet}.
#' @param plot_type Default is \code{"ggplot"}, the type of figure,
#' could be \code{ggplot}, \code{animate} or \code{ggplotly}.
#' @param layout Default is \code{"fruchtermanreingold"}, the layout of figure,
#' could be \code{fruchtermanreingold} or \code{kamadakawai}.
#' @param nrow The number of rows of figure.
#' @param figure_save Default is \code{FALSE},
#' Logical value, whether to save the figure file.
#' @param figure_name The name of figure file.
#' @param figure_width The width of figure.
#' @param figure_height The height of figure.
#' @param seed Default is \code{1}, the seed random use to plot network.
#'
#' @return A dynamic figure object
#' @export
#'
#' @examples
#' data("example_matrix")
#' network <- inferCSN(example_matrix)[1:100, ]
#' network$celltype <- c(
#'   rep("cluster1", 20),
#'   rep("cluster2", 20),
#'   rep("cluster3", 20),
#'   rep("cluster5", 20),
#'   rep("cluster6", 20)
#' )
#'
#' celltypes_order <- c(
#'   "cluster5", "cluster3",
#'   "cluster2", "cluster1",
#'   "cluster6"
#' )
#'
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order
#' )
#'
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order[1:3]
#' )
#'
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order,
#'   plot_type = "ggplotly"
#' )
#'
#' \dontrun{
#' # If setting `plot_type = "animate"` to plot and save `gif` figure,
#' # please install `gifski` package first.
#' plot_dynamic_networks(
#'   network,
#'   celltypes_order = celltypes_order,
#'   plot_type = "animate"
#' )
#' }
plot_dynamic_networks <- function(
    network_table,
    celltypes_order,
    ntop = 10,
    title = NULL,
    theme_type = "theme_void",
    plot_type = "ggplot",
    layout = "fruchtermanreingold",
    nrow = 2,
    figure_save = FALSE,
    figure_name = NULL,
    figure_width = 6,
    figure_height = 6,
    seed = 1) {
  names(network_table) <- c("regulator", "target", "weight", "celltype")
  network_table$regulator <- as.character(network_table$regulator)
  network_table$target <- as.character(network_table$target)
  celltypes_list <- unique(intersect(celltypes_order, network_table$celltype))

  network_table <- purrr::map_dfr(
    celltypes_list,
    .f = function(x) {
      network_table[which(network_table$celltype == x), ]
    }
  )

  # Get nodes information
  nodes <- unique(c(network_table$regulator, network_table$target))
  dnodes <- data.frame(id = 1:length(nodes), label = nodes)
  edges <- dplyr::left_join(
    network_table,
    dnodes,
    by = c("regulator" = "label")
  )
  edges <- dplyr::rename(edges, from = id)
  edges <- dplyr::left_join(
    edges,
    dnodes,
    by = c("target" = "label")
  )
  edges <- dplyr::rename(edges, to = id)
  edges <- dplyr::select(edges, from, to, weight, celltype)
  edges$Interaction <- ifelse(
    edges$weight > 0, "Activation", "Repression"
  )
  edges$weight <- abs(edges$weight)

  dedges <- unique(edges)
  dnodes$label <- gsub("\\.", "-", dnodes$label)
  network_data <- network::network(
    dedges,
    vertex.attr = dnodes,
    matrix.type = "edgelist",
    ignore.eval = FALSE,
    directed = TRUE,
    multiple = TRUE
  )

  set.seed(seed)
  layout <- match.arg(
    layout,
    c("fruchtermanreingold", "kamadakawai")
  )
  ggnetwork_data <- ggnetwork(
    network_data,
    arrow.size = 0.1,
    arrow.gap = 0.015,
    by = "celltype",
    weights = "weight",
    layout = layout
  )

  # Get out-degree for each regulator node per celltype
  nodes_data <- purrr::map_dfr(
    celltypes_list,
    .f = function(x) {
      nodes_data_celltype <- network_table[which(network_table$celltype == x), ]
      nodes_data_celltype <- dplyr::group_by(
        nodes_data_celltype,
        regulator
      )
      nodes_data_celltype <- dplyr::summarise(
        nodes_data_celltype,
        targets_num = dplyr::n()
      )
      nodes_data_celltype <- dplyr::arrange(
        nodes_data_celltype,
        dplyr::desc(targets_num)
      )
      nodes_data_celltype <- as.data.frame(nodes_data_celltype)
      nodes_data_celltype$label_genes <- as.character(
        nodes_data_celltype$regulator
      )
      if (nrow(nodes_data_celltype) > ntop) {
        cf <- nodes_data_celltype$targets_num[ntop]
        nodes_data_celltype$label_genes[which(nodes_data_celltype$targets_num < cf)] <- ""
      } else if (nrow(nodes_data_celltype) == 0) {
        return()
      }
      nodes_data_celltype$celltype <- x

      return(nodes_data_celltype)
    }
  )

  names(nodes_data)[1] <- "label"
  ggnetwork_data <- merge(
    ggnetwork_data,
    nodes_data,
    by = c("label", "celltype"),
    all.x = T
  )
  ggnetwork_data$targets_num[which(is.na(ggnetwork_data$targets_num))] <- 0
  ggnetwork_data$label_genes[which(is.na(ggnetwork_data$label_genes))] <- ""
  ggnetwork_data$celltype <- factor(
    ggnetwork_data$celltype,
    levels = celltypes_order
  )
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  plot_type <- match.arg(
    plot_type,
    c("ggplot", "animate", "ggplotly")
  )
  p <- ggplot(ggnetwork_data, aes(x, y, xend = xend, yend = yend))
  if (plot_type == "ggplotly") {
    p <- p + geom_edges(
      aes(color = Interaction),
      size = 0.7,
      arrow = arrow(length = unit(3, "pt"), type = "closed")
    )
  } else {
    p <- p + geom_edges(
      aes(color = Interaction, alpha = weight),
      size = 0.7,
      arrow = arrow(length = unit(3, "pt"), type = "closed")
    )
  }
  p <- p +
    geom_nodes(
      aes(size = targets_num),
      color = "darkgray",
      alpha = 0.9
    ) +
    geom_nodetext(
      aes(label = label_genes), # , size = targets_num - 1
      color = "black"
    ) +
    theme(aspect.ratio = 2, legend.position = "bottom") +
    scale_color_manual(values = cols)

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  theme_type <- match.arg(
    theme_type,
    c("theme_void", "theme_facet", "theme_blank")
  )
  p <- switch(
    EXPR = theme_type,
    "theme_void" = p + theme_void(),
    "theme_blank" = p + theme_blank(),
    "theme_facet" = p + theme_facet()
  )
  if (plot_type == "ggplot") {
    p <- p +
      facet_wrap(~celltype, nrow = nrow)
    if (figure_save) {
      if (is.null(figure_name)) {
        figure_name <- "networks.pdf"
      }
      ggsave(
        figure_name,
        p,
        width = figure_width,
        height = figure_height
      )
    }
  }

  if (plot_type == "animate") {
    p <- p + gganimate::transition_states(states = celltype)
    p <- gganimate::animate(
      p,
      render = gganimate::gifski_renderer()
    )
    if (figure_save) {
      if (is.null(figure_name)) {
        figure_name <- "networks.gif"
      }
      gganimate::anim_save(figure_name, animation = p)
    }
  }

  if (plot_type == "ggplotly") {
    p <- p +
      facet_wrap(~celltype, nrow = nrow)
    p <- plotly::ggplotly(p)
  }

  return(p)
}
