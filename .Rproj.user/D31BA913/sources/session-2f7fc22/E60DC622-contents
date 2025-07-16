#' Multiple_regression Function
#'
#' This function is used for regression analysis in MR calculations and requires one or more variants.
#' @param x an genotype pair data 
#' @param y an expression data
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
#' beta_coef = Multiple_regression(x = x_paris$S1, y = y)
#' 
#' @export
Multiple_regression <- function(x, y) {
  Beta_fit <- summary(stats::lm(y ~ x))
  Beta_coef <- stats::coef(Beta_fit)
  Beta = Beta_coef[, 1][-1]
  
  return(Beta) 
}
