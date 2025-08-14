#' FineMapping function
#'
#' This function identifies causal variants to be used as instrumental variables.
#' @param x an genotype data in double matrix
#' @param y an expression data in duble matrix
#' @param pip.thr an PIP threshold
#' @keywords fine mapping
#' @export
#' @examples
#'
#' library("dplyr")
#' set.seed(1000)
#' x = sapply(1:10, function(i){
#'   rbinom(100, 2, 0.3) %>% as.double()
#'   }) ; colnames(x) = paste0("rs", 1:ncol(x))
#' y = matrix(0.5*scale(x)[,1] + 0.5*scale(x)[,2] + rnorm(100))
#' x_paris = FineMapping(x, y)
#' head(x_paris)
#'
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
