#'
#'
#' @name gene_pair
#' @title Example gene-gene pair Dataset
#' @description This dataset contains data for demonstaration.
#' @format
#' A data frame with 1012 rows and 780 columns
#' @source <https://figshare.com/s/83bddc1ddf3f97108ad4>
yeast_SNP <- function(){
  data_path <- system.file("data", "yeast_SNP.txt", package = "MRggi")
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
  data_path <- system.file("data", "yeast_exp.txt", package = "MRggi")
  data_loaded <- as.matrix(data_path)
  return(data_loaded)
}
