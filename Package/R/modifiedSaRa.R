#' CNV detection processing multiple sequences using the modified SaRa algorithm
#'
#' This function runs the modified SaRa algorithm and cluster the change-points to CNVs processing multiple sequences.
#' @param Y the numeric vector of the intensities of markers
#' @param alpha the significance levels for the test to accept change-points
#' @param L number of iterations in the EM algorithm for CNV clustering
#' @param h1 the bandwidth 1 for the screening procedure, defaults to 5
#' @param h2 the bandwidth 2 for the screening procedure, defaults to 10
#' @param h3 the bandwidth 3 for the screening procedure, defaults to 15
#' @param precise the precision of the inverse CDF of local min p-values. This will be used only if FINV is not specified. Defaults to 10000
#' @param simT number of simulations in getting the inverse CDF of the local minimum p values
#' @param sigma the standard deviation for the intensities between two adjacent change-points, defaults to NULL
#' @return This function generates a list of detected change-points and clustered CNVs for all samples.
#' @return newcp a list of vectors presenting detected change-points, which is in marker index units. Length of the list is the number of samples or sequences
#' @return h a list of vectors presenting the bandwidth used for this detected change-points. Length of the list is the number of samples
#' @return cnv.state state of detected CNV segments, duplication or deletion
#' @return cnv.start a list of vectors presenting the start position of CNV segments
#' @return cnv.end a list of vectors presenting the end position of CNV segments
#' @seealso \link{multiSaRa} for processing the screening and ranking steps for single sequence
#' @export
modifiedSaRa <-function(Y, alpha = 0.01, h1 = 5, h2 = 10, h3 = 15, L =100, sigma = NULL, precise=10000, FINV=NULL) {
  newcp  <-  h  <-  cnv.state <- cnv.start <-  cnv.end   <- vector() 
  allcp    = multiSaRa(Y, h1 = h1, h2 = h2, h3 = h3, FINV = FINV, precise=precise, sigma = sigma)
  changeP1 = CutCp(allcp[[1]], cutoff = alpha, h = h1)
  changeP2 = CutCp(allcp[[2]], cutoff = alpha, h = h2)
  changeP3 = CutCp(allcp[[3]], cutoff = alpha, h = h3)
  cp.s     = sort(c(changeP1, changeP2, changeP3))
  if (TRUE%in% (length(cp.s)==0)) {
    newcp <- NULL
    h <- NULL
    cnv.state <- NULL
    cnv.start <- NULL
    cnv.end <- NULL
  }else{
    BIC      = getOneBIC(x = Y, cp.s)
    if (BIC$bic != "-Inf") {
      BIC   = iterRemove(BIC)
    }
    cp.BIC   = BIC$cp.new
    names(cp.BIC) = BIC$index
    Gaus.cluster = CNVcluster(Y, cp = cp.BIC, L = L)
    newcp <- Gaus.cluster$newcp
    h <- Gaus.cluster$h
    cnv.state <- Gaus.cluster$cnv.state
    cnv.start <- Gaus.cluster$cnv.start
    cnv.end <- Gaus.cluster$cnv.end
  }
  return (list(newcp =newcp, h = h, cnv.state = cnv.state, cnv.start = cnv.start, cnv.end = cnv.end))
}


CutCp <- function(cp.all, cutoff, h) {
  cp.pv   <- cp.all[order(cp.all[,2]),]
  cut.len <- round(cutoff*dim(cp.pv)[1])
  cut.len <- length(cp.pv[which(cp.pv[, 2]<=cutoff), 2])
  cut.cp  <- cp.pv[c(1 : cut.len), 1]
  names(cut.cp) <- rep(h,length(cut.cp))
  return (cut.cp = cut.cp)
}

