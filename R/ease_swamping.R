### Nous avons remplace la fonction shrink_swamping_effect par ease_swamping

ease_swamping = function(x,y,xquant, yquant, inv_rob_sdx, rob_sdy,number_subset,size_subset,est_clean_set,asymvec,ep=0.1,alpha)
{
  n = nrow(x)
  min_sum_Him = rcpp_mask_swamp_stat(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset, size_subset, est_clean_set, asymvec) # rcpp function
  min_sum_Him = c(min_sum_Him[,1])
  pvv = pchisq(min_sum_Him, df=1, lower.tail=FALSE)
  Spvv = sort.int(pvv,index.return=TRUE)
  Sii = Spvv$ix
  dpv = Spvv$x-alpha/length(min_sum_Him)

  In = which(dpv <= 0)
  if (length(In) == 0)
  {clean_set=est_clean_set}
  else
  {
    rin = max(In)
    inf_setv = est_clean_set[Sii[1:min(floor(ep*n),rin)]]
    clean_set=setdiff(est_clean_set,inf_setv)
  }

  est_clean_set = clean_set
  list(est_clean_set = est_clean_set)
}
