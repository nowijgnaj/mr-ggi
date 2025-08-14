rm(list = ls())


cis_eff_filter <- function(x, y, calc.idx, Bsg.thr = 0) { # x = genotype/cis-SNP | y = gene expression/gene
  y.temp <- scale(y)
  x.temp <- lapply(x, scale)

  Bsg <- sapply(1:ncol(y.temp), function(x){
    lm.fit <- summary(stats::lm(y.temp[, x] ~ x.temp[[x]]))
    stats::coef(lm.fit)[2,1]
  })
  Bsg.idx <- which(abs(Bsg) > Bsg.thr)
  calc.idx <- calc.idx[calc.idx[, 1] %in% Bsg.idx & calc.idx[, 2] %in% Bsg.idx,]
  calc.idx <- as.matrix(calc.idx)
  message("cis-effect : >", Bsg.thr)
  message("Possible network : ", nrow(calc.idx))

  return(calc.idx = calc.idx)
}

cor_filter <- function(x, y, cor.thr = 0) { # x = genotype/cis-SNP | y = gene expression/gene
  x.temp <- lapply(x, scale)
  y.temp <- scale(y)

  cor.y <- stats::cor(y.temp)
  cor.mat <- matrix(0, nrow(cor.y), ncol(cor.y))
  cor.mat[which(abs(cor.y) > cor.thr)] = 1
  cor.mat[lower.tri(cor.mat, diag = TRUE)] = 0
  calc.idx <- which(cor.mat == 1, arr.ind = TRUE)
  message("gene-gene correlation : >", cor.thr)

  return(calc.idx)
}

FineMapping <- function(x, y, pip.thr = 0){ # x = genotype/SNP | y = gene expression/gene
  x.FineMap <- lapply(1:ncol(y), function(i){

    res.chr <- susie(x, y[,i])
    res.cs <- res.chr$sets$cs
    res.pip <- res.chr$pip
    res.pip.idx <- which(res.pip > pip.thr)

    snp.idx <- unlist(lapply(1:length(res.cs), function(x){
      temp1 <- res.cs[[x]][res.cs[[x]] %in% res.pip.idx]

      if(length(temp1) > 0){
        temp1.pip <- res.pip[temp1]
        temp1.max.pip.idx <- which(temp1.pip == max(temp1.pip))
        temp1.snp.idx <- temp1.max.pip.idx[which.max(temp1.max.pip.idx)]
        a <- temp1[temp1.snp.idx]
        a
      }
    }))
    if (is.null(snp.idx)){
      snp.idx <- max(which(res.pip == max(res.pip)))
    }

    if (is.null(res.cs)){
      if (max(res.pip) != 0){
        snp.idx <- max(which(res.pip == max(res.pip)))
      } else if (max(res.pip) == 0){
        stat <- univariate_regression(x, y[,i])
        z <- stat$betahat/stat$sebetahat
        snp.idx <- max(which(z == max(z)))
      }
    }
    as.matrix(x[,snp.idx])
  }) ; names(x.FineMap) <- paste0("S", 1:length(x.FineMap))

  x.FineMap
  return(x.FineMap)
}

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


yeast_SNP <- function(){
  data_path <- system.file("data", "/data/yeast_SNP.txt", package = "MRggi")
  data_load <- as.matrix(data_path)
  return(data_load)
}


#'
#' @name expression
#' @title Example gene expression data set
#' @description This dataset contains data for demonstaration.
#' @format
#' A data frame with 1012 rows and 82 columns
#' @source <https://figshare.com/s/83bddc1ddf3f97108ad4>
yeast_exp <- function(){
  data_path <- system.file("data", "/data/yeast_exp.txt", package = "MRggi")
  data_loaded <- as.matrix(data_path)
  return(data_loaded)
}

MGI<- function(x, y, calc.idx, p.adjust.method = "bonferroni", fdrth = 0.05, c.size = 10, GRN = F) {

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

Multiple_regression <- function(x, y) { # x = genotype/cis-SNP (about gene Sn) | y = gene expression/gene
  Beta_fit <- summary(stats::lm(x ~ y))
  Beta_coef <- stats::coef(Beta_fit)
  Beta = Beta_coef[, 1][-1]

  return(Beta)
}

Quartet_2SLS <- function(x1, x2, y1, y2){ # x = genotype/cis-SNP (about gene Sn) | y = gene expression/gene
  Bzx_fit <- summary(stats::lm(x1 ~ y1))
  Bzx_coef <- stats::coef(Bzx_fit)
  Bzx <- Bzx_coef[,1][-1]

  x2.resid <- stats::resid(stats::lm(x2~y2))
  Bzy_fit <- summary(stats::lm(x2.resid ~ y1))
  Bzy_coef <- stats::coef(Bzy_fit)
  Bzy <- Bzy_coef[,1][-1]
  se_Bzy <- Bzy_coef[,2][-1]
  IV <- 1/(se_Bzy^2)

  Bxy <- sum(Bzx * Bzy * IV)/sum(Bzx * Bzx * IV)
  se_Bxy <- sqrt(1/sum(Bzx * Bzx * IV))

  pval <- 2*stats::pnorm(abs(Bxy/se_Bxy), lower.tail = FALSE)

  return(list(Bxy = Bxy, se_Bxy = se_Bxy, Pvalue = pval))
}

library(dplyr)
library(data.table)
library(susieR)
library(magrittr)

gene = fread("./data/yeast_exp.txt") %>% as.matrix()  # Gene expression
SNP = fread("./data/yeast_SNP.txt") %>% as.matrix()  # SNP genotype
SNP.double <- data.matrix(SNP)
storage.mode(SNP.double) <- "double"
cis = FineMapping(SNP.double,gene) #cis SNP 찾기
# x = genotype/SNP | y = gene expression/gene

result <- cor_filter(cis, gene)
# x = genotype/cis-SNP | y = gene expression/gene
head(result)
nrow(result)

result2 <- cor_filter(cis, gene, cor.thr = 0.8)
# x = genotype/cis-SNP | y = gene expression/gene
head(result2)
nrow(result2)

#table(sapply(cis, ncol))


cis_eff <- cis_eff_filter(cis, gene, result2)
# x = genotype/cis-SNP | y = gene expression/gene
cis_eff
head(cis_eff)

cis_eff2 <- cis_eff_filter(cis,gene, result2, Bsg.thr = 0.1)
head(cis_eff2)


set.seed(1000)

#여기는 랜덤 데이터 생성 같으니 넘기자
beta_coef <- Multiple_regression(cis$S1, gene) #변수 순서가 중요
# x = genotype/cis-SNP (about gene Sn) | y = gene expression/gene

coef<-Quartet_2SLS(cis$S1, cis$S1, gene, gene)
# x = genotype/cis-SNP (about gene Sn) | y = gene expression/gene


head(beta_coef)
print(coef)


mic <- MGI(cis, gene, calc.idx = result2, c.size = 5)
# x = genotype/cis-SNP | y = gene expression/gene

head(mic$edges); head(mic$nodes)

mic2 <- MGI(cis, gene, calc.idx = result2, c.size = 5, GRN = TRUE)

print(mic2$GRN)




# eMatrix.RData와 gList.RData 파일을 로드
load("./data/eMatrix.RData")
load("./data/gList.RData")

# 데이터 확인
ls() # 로드된 객체 리스트를 확인
head(eMatrix) # eMatrix 데이터의 첫 몇 줄 확인
head(gList) # gList 데이터의 첫 몇 줄 확인

result <- cor_filter(gList, eMatrix, cor.thr = 0.8)
head(result)
cis_eff <- cis_eff_filter(gList,eMatrix, result)
cis_eff
head(cis_eff)

