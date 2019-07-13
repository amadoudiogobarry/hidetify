#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]


// Compute a single influential point detection using the asymHIM measure
// [[Rcpp::export]]
arma::mat rcpp_asymHIM_sdetect(const arma::mat& x, const arma::colvec& y, const arma::mat& xquant, const arma::colvec& yquant, const arma::colvec& inv_rob_sdx, const double rob_sdy, 
                               const arma::colvec& asymvec, arma::uvec& inf_set, arma::uvec& non_inf_set)
{
  int nrow = x.n_rows, ncol = x.n_cols, nasym = asymvec.n_elem, size_inf_set = inf_set.n_elem, size_non_inf_set = non_inf_set.n_elem;
  //int nrow = x.n_rows, ncol = x.n_cols, nasym = asymvec.n_elem;

  inf_set = inf_set - 1;
  non_inf_set = non_inf_set - 1;

  arma::mat x1(nrow, ncol);
  arma::mat xnew1((size_non_inf_set+1), ncol);  
  arma::mat xnew2(size_non_inf_set, ncol);
  arma::mat TT(size_inf_set, nasym);
  arma::vec yquant_vec(nrow);
  arma::vec y1(nrow);
  arma::uvec newSet1((size_non_inf_set+1));
  arma::vec ynew1((size_non_inf_set+1));
  arma::vec ynew2(size_non_inf_set);
  arma::vec rhat1(ncol);
  arma::vec rhat2(ncol);
  arma::uvec new_obs(1);

  for (int i = 0; i < size_inf_set; i++)
  {
    new_obs = inf_set[i];
    for (int j = 0; j < nasym; j++) {
      x1 = x - arma::ones(nrow,1)*xquant.row(j);
      x1 = x1*arma::diagmat(inv_rob_sdx);
          yquant_vec.fill(yquant[j]);
      y1 = (y-yquant_vec)/rob_sdy;

      newSet1 = arma::join_cols(new_obs, non_inf_set);
      xnew1 = x1.rows(newSet1);
      ynew1 = y1.elem(newSet1);
      rhat1 = (xnew1.t()*ynew1)/(size_non_inf_set+1);

      xnew2 = x1.rows(non_inf_set);
      ynew2 = y1.elem(non_inf_set);
      rhat2 = xnew2.t()*ynew2/size_non_inf_set;
      TT(i,j) = pow((size_non_inf_set+1),2)*(sum(pow((rhat1-rhat2),2))/ncol);
    }
  }

  return(TT);
}
