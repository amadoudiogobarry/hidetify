#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]

// Package codes.
// Identify single high dimensional influential observation using the asymHIM measure
// Name of the function: shidetify
// Package name: hidetify: High(Hi) dimensional(d) (e)Influential Observations - Identify(tify)

// [[Rcpp::export]]
arma::mat rcpp_shidetify(const arma::mat x, const arma::colvec y, const arma::mat xquant, const arma::colvec yquant, const arma::colvec inv_rob_sdx,
                         const double rob_sdy, const arma::colvec asymvec, const arma::uvec row_indice)
{
  int nrow = x.n_rows, ncol = x.n_cols, nasym = asymvec.n_elem;

  arma::vec yquant_vec(nrow);
  arma::vec y1(nrow);
  arma::vec y2(nrow-1);
  arma::vec rhat1(ncol);
  arma::vec rhat2(ncol);
  arma::mat x1(nrow, ncol);
  arma::mat x2(nrow-1, ncol);
  arma::mat asymHIM(nrow, nasym);


  for (int j = 0; j < nasym; j++)
  {
    x1 = x - arma::ones(nrow,1)*xquant.row(j);
    x1 = x1*arma::diagmat(inv_rob_sdx);
    yquant_vec.fill(yquant[j]);
    y1 = (y-yquant_vec)/rob_sdy;
    rhat1 = (x1.t()*y1)/nrow;
    for (int i = 0; i < nrow; i++) {
      arma::uvec iout = arma::find(row_indice != i);
      x2 = x1.rows(iout);
      y2 = y1.elem(iout);
      rhat2 = x2.t()*y2/(nrow-1);
      asymHIM(i,j) = pow(nrow,2)*(sum(pow((rhat1-rhat2),2))/ncol);
    }
  }

  return(asymHIM);
}

// Il faut l'inclure dans une fonction qui genere les p-valeurs et l'indice des outliers

