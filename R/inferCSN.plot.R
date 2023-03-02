#' Plot of dynamic networks
#'
#' @param weightdf
#' @param tfs
#' @param onlyTFs
#' @param order
#' @param thresh
#'
#' @return
#' @export
#'
#' @examples
inferCSN.plot.dynamic.network <- function(weightdf, tfs = NULL, onlyTFs = TRUE, order = NULL, thresh = NULL){
  library("ggnetwork")
  df <- net.format(weightdf)
  net <- graph_from_data_frame(df[, c("TF", "TG", "interaction")], directed=FALSE)
  #tfnet<-ggnetwork(net,layout="fruchtermanreingold",cell.jitter=0)
  layout <- layout_with_fr(net)
  rownames(layout) <- V(net)$name
  layout_ordered <- layout[V(net)$name,]
  tfnet <- ggnetwork(net, layout = layout_ordered, cell.jitter=0)
  tfnet$is_regulator <- as.character(tfnet$name %in% tfs)
  cols<-c("activation" = "blue", "repression" = "red")
  g<-ggplot()+
    geom_edges(data=tfnet,
               aes(x=x, y=y, xend=xend, yend=yend, color=interaction),
               size=0.75,
               curvature=0.1,
               alpha=.6)+
    geom_nodes(data=tfnet,
               aes(x=x, y=y, xend=xend, yend=yend),
               color="darkgray",
               size=6,
               alpha=.5)+
    geom_nodes(data=tfnet[tfnet$is_regulator=="TRUE",],
               aes(x=x, y=y, xend=xend, yend=yend),
               color="#8C4985",
               size=6,
               alpha=.8)+
    #geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=6, color="#8856a7")+
    scale_color_manual(values=cols)+
    geom_nodelabel_repel(data=tfnet,
                         aes(x=x, y=y, label=name),
                         size=2.5,
                         color="#5A8BAD")+
    theme_blank() #+
  #ggtitle(names(grn)[i])
  g <- g + theme(legend.position="none")
}

#' Title
#'
#' @param data
#' @param plotType
#'
#' @return
#' @export
#'
#' @examples
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

#' net.format
#'
#' @param grn
#'
#' @return
#' @export
#'
#' @examples
net.format <- function(grn){
  colnames(grn) <- c("TF","TG","weight")
  grn$weight <- as.numeric(grn$weight)
  grn$interaction <- "activation"
  grn$interaction[grn$weight < 0] <- "repression"
  return(grn)
}
