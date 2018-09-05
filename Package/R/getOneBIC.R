#' The modified Bayesian Information Criterion step to remove false positives in change-points
#'
#' @param x the vector of the intensities of markers
#' @param cp the vector of the marker index of the identified change-points
#' @return a list of change-points after filtering false positives
#' @export
getOneBIC <- function(x, cp, mod = FALSE) {
  T <- length(x)
  cp.new <- unique(cp[which(cp != T)])
  index  <- vector()
  for (i in 1 : length(cp.new)) {
    index[i] <- max(as.numeric(names(cp[which(cp == cp.new[i])])))
  }
  J <- length(cp.new)

  start <- c(1, cp.new + 1)
  end <- c(cp.new, T)
  len <- end - start + 1
  means <- rep(NA, J + 1)
  for (i in 1:(J + 1)) {
    means[i] <- mean(x[start[i]:end[i]])
  }
  SST <- sum(x^2)
  sigma2 <- (SST - sum(len * means^2))/T

  if (mod) {
    bic <- T/2 * log(sigma2) + 3/2 * J * log(T) + 1/2 * sum(log(len/T))
  } else {
    bic <- T/2 * log(sigma2) + J * log(T)
  }
  return(list(bic = bic, means = means, len = len, SST = SST, cp.new = cp.new, index = index, mod = mod))
}

getAllBICs <- function(x, cp.new, mod = FALSE) {
  N <- nrow(x)
  if (is.list(cp.new) && length(cp.new) != N) {
    stop("length(cp.new) != nrow(x)\n")
  }
  bic <- list()
  if(is.list(cp.new)) {
    for(i in 1:N) {
      bic[[i]] <- getOneBIC(x[i, ], cp.new[[i]], mod)
    }
  } else {
    for(i in 1:N) {
      bic[[i]] <- getOneBIC(x[i, ], cp.new, mod)
    }
  }
  return(bic)
}

removeOneCP <- function(BIC) {
  T <- sum(BIC$len)
  J <- length(BIC$len) - 1
  if(J < 1) stop("No change point left!\n")
  if (BIC$mod) {
    delta.term3 <- rep(NA, J)
    for (i in 1:J) {
      delta.term3[i] <- 1/2 * (log((BIC$len[i] + BIC$len[i + 1])/T) - sum(log(BIC$len[i:(i + 1)]/T)))
    }
    delta0 <- -3/2 * log(T) + delta.term3[i]
  } else {
    delta0 <- -log(T)
  }
  s2 <- (BIC$SST - sum(BIC$len * BIC$means^2))/T
  delta1 <- rep(NA, J)
  for(i in 1:J) {
    means.temp <- BIC$means
    means.temp[i] <- means.temp[i + 1] <- (BIC$means[i] * BIC$len[i] + BIC$means[i + 1] * BIC$len[i + 1])/sum(BIC$len[i:(i + 1)])
    s2.new <- (BIC$SST - sum(BIC$len * means.temp^2))/T
    delta1[i] <- T/2 * (log(s2.new) - log(s2))
  }
  delta <- delta0 + delta1
  bic <- BIC$bic + delta
  bic.new <- min(bic)
  cp.new.remove <- which.min(bic)
  index.remove <- cp.new.remove
  len.new <- BIC$len[-cp.new.remove]
  len.new[cp.new.remove] <- BIC$len[cp.new.remove] + BIC$len[cp.new.remove + 1]
  means.new <- BIC$means[-cp.new.remove]
  means.new[cp.new.remove] <- (BIC$means[cp.new.remove] * BIC$len[cp.new.remove] + BIC$means[cp.new.remove + 1] * BIC$len[cp.new.remove + 1])/sum(BIC$len[cp.new.remove:(cp.new.remove + 1)])
  BIC.new <- list(bic = bic.new, means = means.new, len = len.new, SST = BIC$SST, mod = BIC$mod, cp.new = BIC$cp.new[-cp.new.remove], index = BIC$index[-index.remove])
  return(list(BIC.new = BIC.new, cp.new.remove = cp.new.remove, delta = delta))
}

iterRemove <- function(BIC) {
  try <- TRUE
  while (try && length(BIC$cp.new) > 0) {
    BIC.next <- removeOneCP(BIC)
    if (min(BIC.next$delta) < 0) {
      BIC <- BIC.next$BIC.new
    } else {
      try <- FALSE
    }
  }
  return(BIC)
}

thresholding <- function(x, cp.news, lambda, s) {
  T <- ncol(x)
  N <- nrow(x)

  cp.new.all <- sort(unique(unlist(cp.news)))
  keep <- list()
  delta <- list()
  for (k in 1:length(cp.news)) {
    cp.new <- cp.news[[k]]
    J <- length(cp.new)
    start <- c(1, cp.new + 1)
    end <- c(cp.new, T)
    means <- matrix(NA, N, J + 1)
    for (i in 1:(J + 1)) {
      means[, i] <- rowMeans(x[, start[i]:end[i], drop = FALSE])
    }

    kp <- matrix(FALSE, N, J)
    de <- matrix(0, N, J)
    for(i in 1:N) {
      kk <- 1:J
      mm <- means[i, ]
      ll <- end - start + 1
      dd <- mm[2:(J + 1)] - mm[1:J]
      while(length(kk) > 0 && any(abs(dd) < lambda[i, k])) {
        j <- which.min(abs(dd))
        ll[j] <- ll[j] + ll[j + 1]
        mm[j] <- (mm[j] * ll[j] + mm[j + 1] * ll[j + 1])/(ll[j] + ll[j + 1])
        ll <- ll[-(j + 1)]
        mm <- mm[-(j + 1)]
        kk <- kk[-j]
        dd <- mm[2:length(mm)] - mm[1:(length(mm) - 1)]
      }
      kp[i, kk] <- TRUE
      de[i, kk] <- dd
    }

    keep[[k]] <- kp
    delta[[k]] <- de
    #mu <- sign(means) * pmax(abs(means) - lambda, 0)
    #keep[[k]] <- (mu[, 1:J, drop = FALSE] != 0) | (mu[, 2:(J + 1), drop = FALSE] != 0)
    #delta[[k]] <- mu[, 2:(J + 1), drop = FALSE] - mu[, 1:J, drop = FALSE]
  }

  res <- list()
  for (i in 1:N) {
    cp.new <- list()
    del <- list()
    for(k in 1:length(cp.news)) {
      cp.new[[k]] <- cp.news[[k]][keep[[k]][i, ]]
      del[[k]] <- delta[[k]][i, keep[[k]][i, ]]
    }
    temp <- cp.new[[1]]
    temp.d <- del[[1]]
    if(length(cp.news) > 1) for(k in 2:length(cp.news)) {
      if(length(cp.new[[k]]) > 0) for(j in 1:length(cp.new[[k]])) {
        cc <- cp.new[[k]][j]
        dd <- del[[k]][j]
        ii <- which((temp %in% (cc + ((-s):s))) & (temp.d * dd > 0))
        if(length(ii) > 0) {
          tog <- c(cc, temp[ii])
          tog.d <- c(dd, temp.d[ii])
          temp <- temp[-ii]
          temp.d <- temp.d[-ii]
          jj <- which.max(abs(tog.d))
          temp <- c(temp, tog[jj])
          or <- order(temp)
          temp <- temp[or]
          temp.d <- c(temp.d, tog.d[[jj]])[or]
        } else {
          temp <- c(temp, cc)
          temp.d <- c(temp.d, dd)
          or <- order(temp)
          temp.d <- temp.d[or]
        }
      }
    }
    res[[i]] <- temp
  }
  return(res)
}
