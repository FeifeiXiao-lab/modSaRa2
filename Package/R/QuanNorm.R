#' Quantile normalization of the original intentisites
#'
#' This function runs the quantile normalization procedure for each sequence separately as a preprocessing step.
#' @param lrr the matrix of intensities. Each column represents a sequence or subject, each row represents a single marker
#' @return This function generates a vector of signal intensities for a single sequence after quantile normalization
#' @examples
#' # Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform
#' data(example.data.lrr)
#' lrr.qn <- QuanNorm(example.data.lrr) ## quantile normalization of the intensities
#' @export
 QuanNorm <- function(lrr) {
  T <- dim(lrr)[1]
  n <- dim(lrr)[2]
  lrr.new <- matrix(NA,T,n)
  for (i in 1:n) {
  sample <- rnorm(T)
  sample <- sort(sample)
  lrr.rank <- rank(lrr[,i],ties.method="random")
  lrr.new[,i] <- sample[lrr.rank]
  }
  return(lrr = lrr.new)
 }
