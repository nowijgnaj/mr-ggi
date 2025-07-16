#' filtering cisSNP-effect
#'
#' This function is used to filter out small values of cisSNP-effect.
#' @param x an gentype pair data in list form
#' @param y an expression data in matrix form
#' @param calc.idx the result value obtained by using the cor_filter function
#' @param Bsg.thr an threshold to filter cis-effect that are too small, default is 0
#' @examples
#' x <- yeast_SNP
#' y <- yeast_exp
#'
#' result <- cor_filter(x = x, y = y, cor.thr = 0 )
#'
#' cis_eff_filter(x = x, y = y, calc.idx=result, Bsg.thr = 0)
#'
#' @export
cis_eff_filter <- function(x, y, calc.idx, Bsg.thr = 0) {
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

