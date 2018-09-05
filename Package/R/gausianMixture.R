#' Clustering of CNVs using Expectation-Maximization algorithm
#'
#' This function clusters the identified change-points to make final CNV calling.The potential CNV segments between two neighbor candidate change-points are assigned to different copy number states according to the average intensities in the segment intervals. We use three clusters including duplication, normal state and deletion. Each cluster is presented by Gaussian distribution with unknown mean and variance. Expectation-Maximization (EM) algorithm is applied for the mixture of Gaussians to assign each segment to the most probable cluster/state. Two physically linked candidate CNV segments in the same group are merged to one unique CNV segment.
#' @param x the vector of the intensities of markers
#' @param cp the vector of the marker index of the identified change-points
#' @param priors given initial parameters for the EM algorithm
#' @param L repeat times in the EM algorithm. Defaults to 100
#' @param st number of assumed states in the EM algorithm
#' @return The return is the clustered CNV segments by presenting the start position and end position using SNP or CNV marker index, and the copy number states. It also returns a vector of final candidates of change-points.
#' @return p.final probability of falling into each state for each CNV segment after convergence
#' @return mu.final segment means of each state after convergence
#' @return cp.final list of change-points after EM algorithm
#' @return index.final the bandwidth of change-points
#' @return state.new assigned copy number state for each CNV
#' @export
gausianMixture <- function(x, cp, priors, L,st) {
  T      <- length(x)
  cp.new <- cp[which(cp!=T)]
  index  <- names(cp.new)
  start  <- cp.new+1
  start <- start[-(length(start))]
  end    <- cp.new[-1]
  len    <- end - start + 1
  N      <- length(len)
  means  <- vector()
  sum.x.sq <- vector()
  state    <- vector()
  for (i in 1:N) {
    ##Find segment means
    means[i] <- mean(x[start[i]:end[i]])
    ##Find segment sum.x.sq
    sum.x.sq[i] <- sum(x[start[i]:end[i]]^2)
  }
  para.new <- updateEM(priors$p, priors$mu, priors$sigma,means, sum.x.sq, N, len,st)
  for (iter in 1:L){
    para.new <- updateEM(para.new$p.new, para.new$mu.new, para.new$sigma.new, means, sum.x.sq, N, len,st)
  }
  for (pt in 1:N) {
    state[pt] <- which.max(para.new$p[pt,])
  }
  index.change <- which((state[1:(N-1)] == state[2:N])== "TRUE")
  if (length(index.change) == 0) {
    cp.final     <- cp.new
    index.final  <- index
    state.new    <- state
  }else{
    cp.final     <- cp.new[-(index.change+1)]
    index.final  <- index[-(index.change+1)]
    state.new    <- state[-index.change]
  }
  return(list(p.final = para.new$p.new, mu.final = para.new$mu.new, sigma.final = para.new$sigma.new, cp.final = cp.final, index.final = index.final, state.new = state.new))
}

updateEM <-
  function (p.in, mu.in, sigma.in, means, sum.x.sq, N,len,st) {
    ##Calcute the prob of each segment belong to each state
    p <- dens <- matrix(NA, N, st)
    for (i in 1:N) {
      a <- rep(NA, st)
      for (j in 1:st) {
        dens[i,j] <- dnorm(means[i], mu.in[j], sqrt(sigma.in[j]/len[i]), log=TRUE)
        a[j]      <- log(p.in[j]) + dens[i,j]
      }
      max.a     <- max(a)
      for (k in 1:st) {
        p[i,k] <- (exp(a[k]-max.a))/sum(exp(a-max.a))
      }
    }
    ##Update p, mu, sigma
    p <- na.omit(p)
    p.new <- mu.new <- sigma.new <- rep(NA, st)
    for (k in 1:st) {
      p.new[k]  <- sum(p[,k])/N
      mu.new[k] <- (sum(means*len*p[,k]))/(sum(len*p[,k]))
      sigma.new[k] <- sum(p[,k]*(sum.x.sq-2*len*mu.new[k]*means+len*mu.new[k]^2))/(sum(p[,k]*len))
    }
    return(list(p.new = p.new, mu.new = mu.new, sigma.new = sigma.new, p = p))
  }
