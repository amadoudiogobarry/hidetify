#include "RcppArmadillo.h"

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::uvec arma_sample(arma::uvec& x, int size, bool replace) {
  int nx = x.n_elem;
  arma::vec alpha_prob(nx); 
  alpha_prob.fill(0.5);
  arma::uvec out = sample(x, size, replace, alpha_prob);
  return out;
}


// setdiff function
// https://github.com/oldregan/RcppNotes
// [[Rcpp::export]]
arma::uvec rcpp_setdiff(arma::uvec& x, arma::uvec& y){
  arma::uvec x1 = unique(x);
  arma::uvec y1 = unique(y);
  int ny1 = y1.n_elem;
  //Rcout<<"x1:"<<x1<<endl;
  //Rcout<<"y1"<<y1<<endl;
  arma::uvec x_ret = x1;
  arma::uvec q1;
  if(x.is_empty() && y.is_empty()){
    return(x);
  }else{
    for (int j = 0; j < ny1; j++) {
      if(x_ret.is_empty()){
        break;
      }else{
        q1 = arma::find(x_ret == y1(j));
        if(q1.is_empty()){
          break;
        }
      }
      x_ret.shed_row(q1(0));
    }
    return(x_ret);
  }
}


// This function compute the min of the min and the max of the sum of the asymmetric influence measure

// [[Rcpp::export]]
arma::mat rcpp_mask_swamp_stat(const arma::mat& x, const arma::colvec& y, const arma::mat& xquant, const arma::colvec& yquant, const arma::colvec& inv_rob_sdx, const double rob_sdy,
                               const int number_subset, const int size_subset, arma::uvec& est_clean_set, const arma::colvec& asymvec )
{
  int nrow = x.n_rows, ncol = x.n_cols, nasym = asymvec.n_elem, size_est_clean_set = est_clean_set.n_elem;
  //int nrow = x.n_rows, ncol = x.n_cols, nasym = asymvec.n_elem;
  
  est_clean_set = est_clean_set - 1;
  arma::mat x1(nrow, ncol);
  arma::mat xnew1((size_subset+1), ncol);
  arma::mat xnew2(size_subset, ncol);
  arma::mat TT(number_subset, nasym);
  arma::vec max_sum_Him(size_est_clean_set);
  arma::vec min_min_Him(size_est_clean_set);
  arma::vec min_TT(number_subset);
  arma::vec sum_TT(number_subset);
  arma::vec yquant_vec(nrow);
  arma::vec y1(nrow);
  arma::uvec newSet(size_subset);
  arma::uvec newSet1((size_subset+1));
  arma::vec ynew1((size_subset+1));
  arma::vec ynew2(size_subset);
  arma::vec rhat1(ncol);
  arma::vec rhat2(ncol);
  arma::uvec new_obs(1);
  
  
  for (int i = 0; i < size_est_clean_set; i++)
  {
    //double i = 0;
    new_obs = est_clean_set[i];
    arma::uvec index_sample = rcpp_setdiff(est_clean_set, new_obs);
    
    for (int k = 0; k < number_subset; k++)
    {
      newSet = arma_sample(index_sample, size_subset, FALSE);
      
      for (int j = 0; j < nasym; j++) {
        x1 = x - arma::ones(nrow,1)*xquant.row(j);
        x1 = x1*arma::diagmat(inv_rob_sdx);
        yquant_vec.fill(yquant[j]);
        y1 = (y-yquant_vec)/rob_sdy;
        
        
        newSet1 = arma::join_cols(new_obs, newSet);
        
        xnew1 = x1.rows(newSet1);
        ynew1 = y1.elem(newSet1);
        rhat1 = (xnew1.t()*ynew1)/(size_subset+1);
        
        xnew2 = x1.rows(newSet);
        ynew2 = y1.elem(newSet);
        rhat2 = xnew2.t()*ynew2/size_subset;
        TT(k,j) = pow((size_subset+1),2)*(sum(pow((rhat1-rhat2),2))/ncol);
      }
      
      min_TT[k] = min(TT.row(k)); // Prendre le min de chaque ligne
      sum_TT[k] = sum(TT.row(k)); // Prendre la somme de chaque ligne
    }
    
    min_min_Him[i] = min(min_TT);
    max_sum_Him[i] = max(sum_TT);
  }
  
  return(arma::join_rows(min_min_Him,max_sum_Him));
  //return(TT);
}
