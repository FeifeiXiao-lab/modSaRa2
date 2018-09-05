#' Smooth the original intensities to remove outliers
#'
#' This function runs the smoothing procedure in the original intensities to remove outliers.
#' @param lrr the matrix of the signal intensities. Each column presents a sequence or subject, each row represents a single marker
#’ @param R predefined parameter for smoothing region. For position i in the sequence, the smoothing region was defined as {i-R,...,i,...,i+R}, defaults to 10
#' @param t the tuning parameter for smoothing region, defaults to 2
#' @examples
#' # Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform.
#' data(example.data.lrr)
#' data(example.data.baf)
#' # Use LRR and BAF information of ten samples to calculate eCN
#' eCN.cal <- eCN(lrr=example.data.lrr,baf=example.data.baf)
#' e_CN <- eCN.cal$e_CN
#' # This returns a matrix of smoothed signal intensites of copy number estimatese
#' e_CN.smo <- smooth(e_CN, R=10, t=2)
#' @export
# f <- function(v, first, last) {
#   .Call('f', PACKAGE = 'modSaRa2', v, first, last)
# }

smooth <- function(lrr, R = 10, t = 2) {
  .Call('_modSaRa2_smooth', PACKAGE = 'modSaRa2', lrr, R, t)
}