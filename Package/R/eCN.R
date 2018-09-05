#' eCN
#'
#' This function calculates the copy number estimate eCN using the LRR and BAF intensities.
#' @param lrr the matrix of the log R ratio intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @param baf the matrix of the B Allele Frequency intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @return laf the matrix of lesser allele frequency intensities
#' @return e_CN.smo the matrix of copy number estimate intensities after smoothing
#' @return e_CN the matrix of copy number estimate intensities
#' @return cn.est.L the list of likelihood of being each copy number state for each position in the sequence
#' @seealso \link{Likeli.single} for calculating the likelihood of being each copy number state
#' @seealso \link{smooth} for smoothing the intensities 
#' @examples
#' # Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform.
#' # The map file displays annotation of the markers including the chromosome and location
#' # information of each SNP or CNV marker.
#' data(example.data.lrr)
#' data(example.data.baf)
#' # Use LRR and BAF information of ten samples to calculate eCN
#' eCN.cal <- eCN(lrr=example.data.lrr,baf=example.data.baf)
#' e_CN <- eCN.cal$e_CN
#' # This returns a matrix of new signal intensites of copy number estimates
#' @export
eCN <-function(lrr, baf){
  laf <- baf
  laf[laf > 0.5 & !is.na(laf)] <- 1-laf[laf > 0.5 & !is.na(laf)]
  #cp.f <- vector("list", dim(lrr)[2])
  e_CN <-  matrix(NA, dim(lrr)[1], dim(lrr)[2])
  cn.est.L <- list()
  for (i in 1:dim(lrr)[2]) {
    cn.est.L[[i]] <- list()  
    for (k in 1:length(lrr[,i])) {
            cn.est <-  Likeli.single(lrr=lrr[k,i],laf=laf[k,i])
  	    e_CN[k,i] <- cn.est$X
	    cn.est.L[[i]][[k]] <- NaN*5
	    cn.est.L[[i]][[k]]<- cn.est$L
    }
  }
    e_CN.smo <- smooth(e_CN, R=10, t=2)   
   return (list(laf = laf, e_CN.smo = e_CN.smo,e_CN= e_CN, cn.est.L = cn.est.L))
}
