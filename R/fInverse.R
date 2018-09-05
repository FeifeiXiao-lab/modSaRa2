#' Get the inverse cumulative distribution function of local min p-values
#'
#' This function approximates the inverse cumulative distribution function (CDF) of the local min p-values via simulation.
#' @param n the length of the data to simulate
#' @param h the bandwidth for the screening procedure. Defaults to 10
#' @param hh the bandwidth for the local minimum procedure
#' @param precise the precision of the approximation (the number of quantiles to use)
#' @param simT number of simulations
#' @return the quantiles of the approximated inverse CDF
#' @export
fInverse <-
function(n = 10000, h= 10, hh = 2 * h, precise = 10000, simT = 100){        #need to be faster
   empirical = NULL
   for (i in 1 : simT){
     Y   = rnorm(n)
     LDF  =  localDiagnostic(Y, h)
     LDF.pos      =  LDF
     LDF.neg      =  -LDF
     sigma = 1
     index.pos = localMax(y = LDF.pos, span = hh)
     pV.pos   = 1 - 2 * pnorm(LDF.pos[index.pos] / (sqrt(2 / h) * sigma))
     index.neg = localMax(y = LDF.neg, span = hh)
     pV.neg   = 1 - 2 * pnorm(LDF.neg[index.neg] / (sqrt(2 / h) * sigma))
     #empirical = c(empirical, pV.pos, pV.neg)
     index <- c(index.pos, index.neg)
     pv <- c(pV.pos, pV.neg)
     pv <- pv[order(index)]
     index <- sort(index)
     len <- length(index)
     rm.list <- NULL
     for (j in 1 : (len - 1)) {
       if(index[j] >= index[j + 1] - h){
         rm.list <- c(rm.list, j)
        }
     }
     if (length(rm.list) > 0) {
      pv <- pv[-rm.list]
     }
     empirical <- c(empirical, pv)
     if (length(empirical) > 10 * precise) break
   }
   return(quantile(empirical, probs = c(0 : precise) / precise))
}
