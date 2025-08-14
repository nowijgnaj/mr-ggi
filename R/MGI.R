#' MGI function
#'
#' This function identifies gene-gene interactions for all gene-gene pairs that genes can have.
#' @param x an genotype pair in list form
#' @param y an expression data (n x k matrix, n = sample number, k = gene number)
#' @param calc.idx the data result value after correlation filtering and cis-effect filtering
#' @param p.adjust.method multiple correlation method
#' @param fdrth FDR threshold
#' @param c.size the number of genes constructing one cluster
#' @param GRN show clusters using clustering algorithm
#' @examples
#'
#' x <- yeast_SNP
#' y <- yeast_exp
#'
#' result <- cor_filter(x, y)
#' result <- cis_eff_filter(x, y, result)
#'
#' MRggi <- MGI(x, y, calc.idx = result, fdrth = 0.03, c.size = 5, GRN = TRUE)
#'
#' head(MRggi$data$edges) ; head(MRggi$data$nodes)
#' MRggi$GRN
#'
#' @export
MGI<- function(x, y, calc.idx, p.adjust.method = "bonferroni", fdrth = 0.05, c.size = 10, GRN = F) {
  # x = genotype/cis-SNP | y = gene expression/gene
  g1 <- c() ; g2 <- c() ; GGcor <- c()
  max_Bs1g1_vec <- c() ; min_Bs1g1_vec <- c()
  max_Bs2g2_vec <- c() ; min_Bs2g2_vec <- c()
  Bg1g2 <- c() ; Bg2g1 <- c()
  pval_Bg1g2 <- c() ; pval_Bg2g1  <- c()

  for (k in 1:nrow(calc.idx)) {
    i <- calc.idx[k, 1] ; j <- calc.idx[k, 2]
    x1 <- x[[i]] ; x2 <- x[[j]]
    y1 <- y[, i] ; y2 <- y[, j]

    Bs1g1 <- Multiple_regression(x1, y1) ; Bs2g2 <- Multiple_regression(x2, y2)

    Bg1g2.data <- Quartet_2SLS(x1, x2, y1, y2)
    Bg2g1.data <- Quartet_2SLS(x2, x1, y2, y1)

    g1 <- append(g1, colnames(y)[i]) ; g2 <- append(g2, colnames(y)[j])
    GGcor <- append(GGcor, stats::cor(y1, y2))
    max_Bs1g1_vec <- append(max_Bs1g1_vec, max(Bs1g1)) ; min_Bs1g1_vec <- append(min_Bs1g1_vec, min(Bs1g1))
    max_Bs2g2_vec <- append(max_Bs2g2_vec, max(Bs2g2)) ; min_Bs2g2_vec <- append(min_Bs2g2_vec, min(Bs2g2))
    Bg1g2 <- append(Bg1g2, Bg1g2.data$Bxy) ; pval_Bg1g2 <- append(pval_Bg1g2, Bg1g2.data$Pvalue)
    Bg2g1 <- append(Bg2g1, Bg2g1.data$Bxy) ; pval_Bg2g1 <- append(pval_Bg2g1, Bg2g1.data$Pvalue)
  }

  FDR_Bg1g2 <- rep(0, length(g1))
  g1.lv <- levels(as.factor(g1))
  for (i in g1.lv) {
    idx <- which(g1 == i)
    pval.idx <- pval_Bg1g2[idx]
    FDR.idx <- stats::p.adjust(pval.idx, method = stats::p.adjust.methods)
    FDR_Bg1g2[idx] <- FDR.idx
  }
  FDR_Bg2g1 <- rep(0, length(g2))
  g2.lv <- levels(as.factor(g2))
  for (i in g2.lv) {
    idx <- which(g2 == i)
    pval.idx <- pval_Bg2g1[idx]
    FDR.idx <- stats::p.adjust(pval.idx, method = stats::p.adjust.methods)
    FDR_Bg2g1[idx] <- FDR.idx
  }

  res <- data.frame(g1 = g1, g2 = g2, GGcor = GGcor, max_Bs1g1 = max_Bs1g1_vec, min_Bs1g1 = min_Bs1g1_vec,
                    max_Bs2g2 = max_Bs2g2_vec, min_Bs2g2 = min_Bs2g2_vec, Bg1g2 = Bg1g2, pval_Bg1g2 = pval_Bg1g2,
                    FDR_Bg1g2 = FDR_Bg1g2, Bg2g1 = Bg2g1, pval_Bg2g1 = pval_Bg2g1, FDR_Bg2g1 = FDR_Bg2g1)

  sig_res <- res[res$FDR_Bg1g2 < fdrth | res$FDR_Bg2g1 < fdrth, ]

  message("filtering < ", fdrth, "...")
  if (nrow(sig_res) == 0){
    stop("There is no network.")
  }

  sig_res$fdrth.g1g2 <- sapply(1:length(sig_res$FDR_Bg1g2), function(x) {
    if (sig_res$FDR_Bg1g2[x] < fdrth) {1} else {0}
  })
  sig_res$fdrth.g2g1 <- sapply(1:length(sig_res$FDR_Bg2g1), function(x) {
    if (sig_res$FDR_Bg2g1[x] < fdrth) {1} else {0}
  })

  temp1 <- sig_res[sig_res$fdrth.g1g2 == 1, c("g1", "g2", "Bg1g2")]
  temp2 <- sig_res[sig_res$fdrth.g2g1 == 1, c("g2", "g1", "Bg2g1")]
  colnames(temp1) <- colnames(temp2) <- c("from", "to", "Bgg")
  edges <- rbind(temp1, temp2)
  edges[, c(1,2)] <- apply(edges[,c(1,2)], 2, as.character)

  nodes <- data.frame(id = sort(unique(c(edges$from, edges$to))),
                      label = sort(unique(c(edges$from, edges$to))))
  message("nodes: ", nrow(nodes), " / edges: ", nrow(edges))

  graph <- influential::graph_from_data_frame(edges, directed = FALSE)
  cluster <- igraph::cluster_louvain(graph)

  cluster_df <- as.data.frame(t(data.frame(as.list(igraph::membership(cluster)))))
  cluster_df$label <- rownames(cluster_df)
  cluster_df <- cluster_df[order(cluster_df$label),]
  temp <- as.data.frame(table(cluster_df$V1))
  temp.idx <- temp$Var1[ifelse(temp$Freq > c.size, TRUE, FALSE)]
  cluster_df$temp <- ifelse(cluster_df$V1 %in% temp.idx, cluster_df$V1, 0)

  nodes <- nodes[order(nodes$id),]
  nodes$temp <- cluster_df$temp
  temp.group.idx <- as.factor(nodes$temp)
  nodes$group <- sapply(1:nrow(nodes), function(x) {
    for (i in 1:nlevels(temp.group.idx)) {
      if (nodes$temp[x] == as.double(levels(temp.group.idx))[i]) {
        a <- paste("cluster", i-1)
      } else if (nodes$temp[x] == 0) {
        a <- "-"
      }
    }
    a
  }) ; nodes <- nodes[, -3]


  # -------------------- temp -------------------- #

  res = list(edges = edges, nodes = nodes)

  if (GRN == TRUE){
    a = GRN_clustering_focus(nodes = res$nodes, edges = res$edges)  # GRNcluster_temp.R
    res = list(data = res, GRN = a)

    return(res)
  } else{
    return(res)
  }
}
