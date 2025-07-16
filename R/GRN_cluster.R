#' GRN_clustering_focus function
#'
#' It is a function that represents the interaction of genes in a cluster plot using a clustering algorithm.
#' @param nodes means genes
#' @param edges indicates gene-gene interactions
#' @examples
#' x <- yeast_SNP
#' y <- yeast_exp
#'
#' result <- cor_filter(x, y)
#' result <- cis_eff_filter(x, y, result)
#' MRggi <- MGI(x, y, calc.idx = result,c.size = 5, GRN = FALSE)
#'
#' grn <- GRN_clustering_focus(MRggi$nodes, MRggi$edges)
#'
#' @export
GRN_clustering_focus <- function(nodes, edges){

  # nodes color
  pal.nodes = "Set3"
  pal.idx.nodes = rev(colorspace::qualitative_hcl(n = nlevels(as.factor(nodes$group)), palette = pal.nodes))  # sequential palette
  clu.idx = levels(as.factor(nodes$group))

  color = rep(0,nrow(nodes))
  for (i in 1:nlevels(as.factor(nodes$group))){
    color[which(nodes$group==clu.idx[i])] = pal.idx.nodes[i]
  }
  nodes$color = color

  # edges color
  pal.edges = "Blue-Red"
  pal.idx = colorspace::diverging_hcl(n = 7, palette = pal.edges)
  Bgg.idx = sort(c(Inf, as.vector(summary(abs(edges$Bgg))[c("1st Qu.","Median","3rd Qu.")]),
                   -as.vector(summary(abs(edges$Bgg))[c("1st Qu.","Median","3rd Qu.")]), -Inf))
  temp = edges$Bgg

  edges$color = sapply(1:nrow(edges), function(x){
    for (i in 1:7){
      if (temp[x]>=Bgg.idx[i] & temp[x]<Bgg.idx[i+1]){col = pal.idx[i]}
    } ; col
  })
  lnodes <- data.frame(label = clu.idx,
                       color = pal.idx.nodes)

    visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visNodes(size = 10, shape = "circle") %>%
    visNetwork::visEdges(arrows ="to", width = 2) %>%
    visNetwork::visOptions(manipulation = TRUE) %>%
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visNetwork::visLegend(addNodes = lnodes, useGroups = FALSE)
}
