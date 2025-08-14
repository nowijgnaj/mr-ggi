#' filtering gene-gene correlations
#'
#' This function is used to filter out small values of gene-gene correlation.
#' @param x an gentype pair data in list form
#' @param y an expression data in matrix form
#' @param cor.thr an threshold to filter gene-gene correlations that are too small, default is 0
#' @examples
#' x <- yeast_SNP
#' y <- yeast_exp
#'
#' cor_filter(x = x, y = y, cor.thr = 0 )
#'
#' @export
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
