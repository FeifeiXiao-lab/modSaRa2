#' Screening procedure processing single sequence to find local maximizers of the dignostic statistic
#'
#' This function runs the screening procedure under single bandwidth processing a single sequence. Local maximizers of the diagnostic statistic or the corresponding minimum p-values within bandwidth h are identified.
#' @param Y the numeric vector of the intensities of markers
#' @param h the bandwidth for the screening procedure, defaults to 10
#' @param hh the bandwidth for the local minimum procedure
#' @param sigma the standard deviation for the intensities between two adjacent change points, defaults to NULL
#' @return The return is a vector of local p-value minimum and the marker index of these minimum.
#' @return index numeric vector a position index for all the local p-value minimum
#' @return pV local minimum p-values
#' @seealso \link{estimateSigma} for estimation of the standard deviation between two adjacent change-points. \link{localMax} for calculation of the local maximizers of local diagnostic function
#' @export
SARAp <-function(Y, h, hh = 2 * h, sigma = NULL){
   n        = length(Y)
   LDF      =  localDiagnostic(Y, h)
   LDF.pos      =  LDF
   LDF.neg      =  -LDF

   if (is.null(sigma)) sigma = estimateSigma(Y, h = max(3, 2 * floor(log(n))))
   pV.pos       = 1 - 2 * pnorm(LDF.pos / (sqrt(2 / h) * sigma))
   LocalMax     = localMax(LDF.pos, span = hh)                            # aviod the case that several P-values are zero
   LocalMax        = clean(LocalMax, h)
   LocalMaxValue = pV.pos[LocalMax]                                       # all local min value

   pV.neg       = 1 - 2 * pnorm(LDF.neg / (sqrt(2 / h)*sigma))
   LocalMin = localMax(LDF.neg, span = hh)                               # aviod the case that several P-values are zero
   LocalMin = clean(LocalMin, h)

   LocalMinValue = pV.neg[LocalMin]                                      # all local min value

   LocalExt <- c(LocalMax, LocalMin)
   LocalExtValue <- c(LocalMaxValue, LocalMinValue)

   LocalExtValue <- LocalExtValue[order(LocalExt)]
   LocalExt <- sort(LocalExt)

  return (list(index = LocalExt, pV = LocalExtValue))

}

clean <- function(LocalM, h) {
   len <- length(LocalM)
   rm.list <- NULL
   for (i in 1 : (len - 1)) {
    if(LocalM[i] >= LocalM[i + 1] - h){
       rm.list <- c(rm.list, i)
      }
   }
   if(length(rm.list) > 0) {
   LocalM <- LocalM[-as.vector(rm.list)]
   }
   return (LocalM = LocalM)
   }
