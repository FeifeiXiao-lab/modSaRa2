#' CNVcluster
#'
#' This function uses different priors settings to achieve convergence of the Expectation-Maximization algorithm, and then determine the best clustering results by applying the modified BIC.
#' @param Y the numeric vector of the intensities of markers
#' @param cp the numeric vector of the position index for the identified change-points
#' @param L repeat times in the EM algorithm, defaults to 100
#' @return The return is the clustered CNV segments by presenting the start position and end position, length of the CNV and the copy number states (duplication or deletion). It also returns a vector of final candidates of change-points.
#' @return  newcp the final list of change-points
#' @return  h the bandwidth used for the identification of change-points
#' @return  cnv.state copy number state for each CNV
#' @return  cnv.start start position of each CNV
#' @return  cnv.end end position of each CNV
#' @seealso \link{gausianMixture} for clustering of CNVs using Expectation-Maximization algorithm.
#' @export
CNVcluster <-function(Y, cp, L) {
  bic.v <- vector()
  EM <- vector("list", 3)
  st1 = 2
  p1 = rep(1/st1, st1)
  sigma1 = rep(0.5, st1)
  mu1 = c(0.5, 2)
  priors1 = list(p = p1, mu = mu1, sigma = sigma1)
  EM[[1]] = gausianMixture(Y, cp, priors1, L, st1)
  st2 = 2
  p2 = rep(1/st2, st2)
  sigma2 = rep(0.5, st2)
  mu2 = c(2, 4)
  priors2 = list(p = p2, mu = mu2, sigma = sigma2)
  EM[[2]] = gausianMixture(Y, cp, priors2, L, st2)
  EM[[2]]$state.new[which(EM[[2]]$state.new == 2)] = 0
  EM[[2]]$state.new[which(EM[[2]]$state.new == 1)] = 2
  EM[[2]]$state.new[which(EM[[2]]$state.new == 0)] = 1
  st3 = 3
  mu3 = c(0.5, 2, 5)
  p3 = rep(1/st3, st3)
  sigma3 = rep(0.5, st3)
  priors3 <- list(p = p3, mu = mu3, sigma = sigma3)
  EM[[3]] = gausianMixture(Y, cp, priors3, L, st3)
  for (i in 1:3) {
    bic.v[i] <- getOneBIC(Y, EM[[i]]$cp.final)$bic
  }
  if (bic.v[3] == min(bic.v)) {
    mins = 3
  } else {
    mins = which.min(bic.v)
  }
  EM.f = EM[[mins]]
  newcp = EM.f$cp.final
  h = EM.f$index.final
  cnv.state <- getState(EM = EM.f, mins)
  return(list(newcp = newcp, h = h, cnv.state = cnv.state$cnv.state, 
              cnv.start = cnv.state$cnv.start, cnv.end = cnv.state$cnv.end))
}

getState <-function (EM = EM, mins) {
  state      = EM$state.new
  cp.f       = EM$cp.final
  start.index  = which(state!=2)  
  cnv.start     = cp.f[start.index]
  cnv.start = cnv.start[!is.na(cnv.start)]
  cnv.end = cp.f[start.index+1]
  cnv.end = cnv.end[!is.na(cnv.end)]
  cnv.state = state[start.index]

  if (TRUE %in% (mins == 3)) {
    if(EM$mu.final[1] < 2 & EM$mu.final[3] > 2) {
      cnv.state[which(cnv.state == 1)] = "del"
      cnv.state[which(cnv.state == 3)] = "dup"
    }else if(EM$mu.final[1] < 2 & EM$mu.final[3] < 2){
      cnv.state[which(cnv.state != 2)] = "del"
    }else if(EM$mu.final[1] > 2 & EM$mu.final[3] > 2){
      cnv.state[which(cnv.state != 2)] = "dup"    
    }else if(EM$mu.final[1] > 2 & EM$mu.final[3] < 2){
      cnv.state[which(cnv.state == 1)] = "dup"
      cnv.state[which(cnv.state == 3)] = "del"
    }
  }
  if (TRUE %in% (mins == 1)) {
    cnv.state[which(cnv.state == 1)] = "del"
  }
  else if (TRUE %in%(mins == 2)) {
    cnv.state[which(cnv.state == 1)] = "dup"
  }
  return (list(cnv.state = cnv.state, cnv.start = cnv.start, cnv.end = cnv.end))
}
