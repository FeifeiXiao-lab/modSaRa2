#' Likeli.single
#'
#' This function integrates the LRR and BAF intensities to calculate the likelihood of being each copy number state for a position in a sequence
#' @param lrr the matrix of the log R ratio intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @param baf the matrix of the B Allele Frequency intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @return X the calculated copy number estimate
#' @return L the vector of likelihood of being each copy number state
#' @return L_geno the vector of likelihood of being each genotype
#' @export
Likeli.single <- function(lrr,laf) {
  mu_lrr <- c(-0.45,0,0,0.3,0.3,0.75,0.75)
  sd_lrr <- rep(0.18,7)
  mu_laf <- c(0,0,0.5,0,1/3,0,0.25)
  sd_laf <- rep(0.03,7)
  L <- L_geno <- L_lrr <- L_laf <- vector()
  for (i in 1:7) {
	L_lrr[i] <- (1/(sd_lrr[i]*sqrt(2*pi)))*exp(-(lrr-mu_lrr[i])^2/(2*(sd_lrr[i])^2))
        L_laf[i] <- (1/(sd_laf[i]*sqrt(2*pi)))*exp(-(laf-mu_laf[i])^2/(2*(sd_laf[i])^2))
	L_geno[i] <- L_lrr[i]*L_laf[i] 
  }
  L["st0"] <-  (1/(2*sqrt(2*pi)))*exp(-(lrr+5)^2/8) ## being double deletions 
  L["st1"] <- L_geno[1] ## being single deletions
  L["st2"] <- L_geno[2]+L_geno[3] ## being normal
  L["st3"] <- L_geno[4]+L_geno[5] ## being single duplication
  L["st4"] <- L_geno[6]+L_geno[7] ##being double duplications
  X <- (L[2]+2*L[3]+3*L[4]+4*L[5])/sum(L)
  return(list(X=X,L=L,L_geno=L_geno))
}
