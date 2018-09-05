#' Screening procedure processing single sequence to find local maximizers of the local dignostic statistic
#'
#' This function runs the screening step under multiple bandwidths processing a single sequence.
#' @param Y the vector of the intensities of markers
#' @param h1 the bandwidth 1 for the screening procedure, defaults to 5
#' @param h2 the bandwidth 2 for the screening procedure, defaults to 10
#' @param h3 the bandwidth 3 for the screening procedure, defaults to 15
#' @param FINV the inverse CDF of the local minimum p-values, approximated by function fInverse()
#' @param precise the precision of the inverse CDF of local min p-values. This will be used only if FINV is not specified. Defaults to 10000
#' @param sigma the standard deviation for the intensities between two adjacent change-points, defaults to NULL
#' @return The return is a list of index with local minimum p values at each bandwidth.
#' @seealso \link{SARA} for processing the screening and ranking steps using single bandwidth
#' @export
multiSaRa <-function(Y, h1 = 3*round(log10(length(Y))), h2 = 2*round(log10(length(Y))),
           h3 = round(log10(length(Y))),FINV = NULL, precise =10000, sigma = NULL) {
    cp = vector("list", 3)
    index = vector("list", 3)
    pv = vector("list", 3)
    h = c(h1, h2, h3)
    for (i in 1 : 3) {
      local.test = SARA(Y = Y, h = h[i], FINV = FINV[[i]], precise = precise, sigma = sigma)
      index[[i]] = local.test$index
      pv[[i]]    = local.test$pV
      cp[[i]]    = matrix(NA, length(index[[i]]), 2)
      cp[[i]][,1] = index[[i]]
      cp[[i]][,2] = pv[[i]]
    }
    return (cp = cp)
  }

