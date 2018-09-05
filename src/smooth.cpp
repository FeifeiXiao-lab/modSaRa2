#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::vec f(const arma::vec & v, int first, int last) {
  arma::vec out = v.subvec(first, last);
  return out;
}


//[[Rcpp::export]]
arma::mat smooth(arma::mat lrr, int R=10, int t=2){
  
  int T = lrr.n_rows;
  int n=lrr.n_cols;
  arma::mat lrr_new(T,n);
 
  for (int j = 0; j <= (n-1); j++){
    arma::vec x=lrr.col(j);
   arma::vec sigma(T);
   arma::vec mean(T);
   arma::vec median(T);
   for (int i = 0; i <= (R-1); i++){
     arma::vec tem;
     tem=f(x,0,(i+R));
        sigma(i) = arma::stddev(tem);
      
       mean[i] = arma::mean(tem);
       median[i] = arma::median(tem);
   }
   for (int i = (R); i <= (T-R-1); i++){
     arma::vec tem;
     tem=f(x,(i-R),(i+R));
       sigma[i] = arma::stddev(tem);
     mean[i] = arma::mean(tem);
     median[i] = arma::median(tem);
       }
   for (int i = (T-R); i <= (T-1); i++){
     arma::vec tem;
     tem=f(x,(i-R),(T-1));
       sigma[i] = arma::stddev(tem);
     mean[i] = arma::mean(tem);
     median[i] = arma::median(tem);
    }
   for(int i=0; i<=(T-1);i++){
       if (x(i)>mean(i) + t * sigma(i) || x(i) < mean(i) - t * sigma(i)) {
           x(i) = median(i);
       }
   }
   lrr_new.col(j)=x;
  }
   return(lrr_new);
}







