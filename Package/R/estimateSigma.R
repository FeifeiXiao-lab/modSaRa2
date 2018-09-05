#' Estimate the standard deviation of the intensities between two adjacent change-points
#'
#' This function estimates the standard deviation between two adjacent change-points using a local smoother.
#' @param Y the numeric vector of the intensities of markers
#' @param h the bandwidth of the local smoother
#' @return This function estimates the standard deviation of the intensities between two adjacent change points
#' @export
estimateSigma <-
function(Y,h=10){                   #constant case
  n     = length(Y)                                 # can we make it faster?
  YBar  = rep(0,n)
  for (i in 1:n) {
     a       = min(n,i+h)
     b       = max(1,i-h)
     YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y-YBar)*(2*h+1)/(2*h)))
}
