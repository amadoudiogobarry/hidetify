### Nous avons remplace la fonction asymMIP par hidetify


mhidetify <- function(x, y, number_subset, size_subset, asymvec, ep=0.1, alpha_swamp, alpha_mask, alpha_validate)
{
  df <- length(asymvec)
  xquant <- apply(x,2,quantile,asymvec); yquant <- quantile(y,asymvec); inv_rob_sdx <- 1/apply(x,2,mad); rob_sdy = mad(y);
  xnrow <- nrow(x)
  est_clean_set <- c(1:xnrow)
  est_clean_set <- ease_swamping(x,y,xquant, yquant, inv_rob_sdx, rob_sdy,number_subset,size_subset,est_clean_set,asymvec,ep=0.1,alpha_swamp)$est_clean_set
  clean_set <- ease_masking(x,y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset,size_subset,est_clean_set,asymvec, alpha_mask)$clean_set

  while (length(clean_set) < xnrow/2)
  {
    est_clean_set <- ease_swamping(x,y,xquant, yquant, inv_rob_sdx, rob_sdy,number_subset,size_subset,est_clean_set,asymvec,ep=0.1,alpha_swamp)$est_clean_set
    clean_set <- ease_masking(x,y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset,size_subset,est_clean_set,asymvec, alpha_mask)$clean_set
  }

  non_inf_set <- clean_set
  inf_set <- setdiff(c(1:xnrow), non_inf_set)
  inf_setfinal <- vhidetify(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, asymvec, inf_set, non_inf_set, alpha_validate)$inf_setfinal
  dataframe <- data.frame(ind = 1:xnrow, outlier_ind = rep(0,xnrow))
  dataframe$outlier_ind[which(dataframe$ind%in%inf_setfinal)] = 1
  dataframe

}

