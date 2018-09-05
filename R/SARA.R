#' Screening procedure processing single sequence and p value calculation of local maximizers
#'
#' This function runs the screening and ranking algorithm under single bandwidth processing a single sequence. Local min p-values are corrected by approximation emprically by the distribution of the local minimizers in a long standard normal sequence.
#' @param Y the numeric vector of the intensities of markers
#' @param h the bandwidth for the screening procedure, defaults to 10
#' @param hh the bandwidth for the local minimum procedure
#' @param FINV the inverse CDF of the local min p-values, approximated by function fInverse
#' @param sigma the standard deviation for the intensities between two adjacent change points, defaults to NULL
#' @param precise the precision of the inverse CDF of local minimum p-values. This will be used only if FINV is not specified. Defaults to 10000
#' @return The return is a vector of corrected local p-value minimum and the marker index of these minimum.
#' @return index index of markers for the corrected local p-value minimizers
#' @return pV corrected local min p-values
#' @seealso \link{SARAp} for processing the screening step using single bandwidth to find local maximizers of the diagnostic statistic
#' @export
SARA <-function(Y, h = 10, hh = 2 * h, FINV = NULL, sigma = NULL, precise = 10000){
   object = SARAp(Y = Y, h = h, hh = hh, sigma = sigma)
   index  = object$index
   pV     = object$pV
   if (is.null(FINV)) FINV = fInverse(n =  length(Y), h = h, hh = hh, precise = precise, simT=100)
   pVcorr = pV
   for (i in 1:length(pV)){
         #pVcorr[i] = (max(length(which(FINV<pV[i])),1)-1)/precise
         pVcorr[i] = (length(which(FINV < pV[i])))/ (precise + 1)
   }
   return (list(index = index,pV = pVcorr))
   }

