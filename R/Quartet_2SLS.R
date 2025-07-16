#' Quartet_2SLS Function
#'
#' This function is used for regression analysis in MR calculation and obtains p-values using the Wald test.
#' @param x1 an genotype pair data 
#' @param x2 an genotype pair data 
#' @param y1 an expression data
#' @param y2 an expression data
#' @examples 
#' 
#' library(dplyr)
#' set.seed(1000)
#' x = sapply(1:10, function(i){
#'   rbinom(100, 2, 0.3) %>% as.double()
#'   }) ; colnames(x) = paste0("rs", 1:ncol(x))
#' y = matrix(0.5*scale(x)[,1] + 0.5*scale(x)[,2] + rnorm(100))
#' x_paris = FineMapping(x, y)
#' 
#' coef <- Quartet_2SLS(x1 = x_paris$S1, x2 = x_paris$S1, y1 = y, y2 = y)
#' 
#' @export

Quartet_2SLS <- function(x1, x2, y1, y2){
  Bzx_fit <- summary(stats::lm(y1 ~ x1))
  Bzx_coef <- stats::coef(Bzx_fit)
  Bzx <- Bzx_coef[,1][-1]
  
  y2.resid <- stats::resid(stats::lm(y2~x2))
  Bzy_fit <- summary(stats::lm(y2.resid ~ x1))
  Bzy_coef <- stats::coef(Bzy_fit)
  Bzy <- Bzy_coef[,1][-1]
  se_Bzy <- Bzy_coef[,2][-1]
  IV <- 1/(se_Bzy^2)
  
  Bxy <- sum(Bzx * Bzy * IV)/sum(Bzx * Bzx * IV)
  se_Bxy <- sqrt(1/sum(Bzx * Bzx * IV))
  
  pval <- 2*stats::pnorm(abs(Bxy/se_Bxy), lower.tail = FALSE)
  
  return(list(Bxy = Bxy, se_Bxy = se_Bxy, Pvalue = pval))
  
}
