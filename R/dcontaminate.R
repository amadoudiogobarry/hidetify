### This function estimate the degree of contamination of the data using the very conservative statistics minOfmin

dcontaminate = function(x, y, number_subset, size_subset, asymvec, ep=0.1, alpha_dcon){
  
  xquant <- apply(x,2,quantile,asymvec)
  yquant <- quantile(y,asymvec)
  inv_rob_sdx <- 1/apply(x,2,mad)
  rob_sdy = mad(y);
  xnrow <- nrow(x)
  est_clean_set <- c(1:xnrow)
  clean_set = ease_swamping(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset, size_subset, est_clean_set, 
                            asymvec, ep, alpha_dcon)$est_clean_set
  non_inf_set <- clean_set
  inf_set <- setdiff(c(1:xnrow), non_inf_set)
  dataframe <- data.frame(ind = 1:xnrow, outlier_ind = rep(0,xnrow))
  dataframe$outlier_ind[which(dataframe$ind%in%inf_set)] = 1
  dataframe
  
}