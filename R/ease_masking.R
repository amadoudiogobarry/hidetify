### Nous avons remplace la fonction shrink_masking_effect par ease_masking


ease_masking = function(x,y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset,size_subset,est_clean_set,asymvec,alpha)
{
  #max_sum_Him = (rcpp_mask_swamp_stat(x,y,number_subset,size_subset,est_clean_set,asymvec))$max_sum_Him
  max_sum_Him = rcpp_mask_swamp_stat(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset, size_subset, est_clean_set, asymvec) # rcpp function
  max_sum_Him = c(max_sum_Him[,2])
  df =  length(asymvec)
  pv = pchisq(max_sum_Him, df=df, lower.tail=FALSE)
  Spv=sort.int(pv,index.return=TRUE)  # sorted p value
  Si=Spv$ix
  dp=Spv$x-alpha/length(max_sum_Him)

  In = which(dp<=0)    # BH procedure to control the error rate
  if (length(In) == 0)
  {clean_set = est_clean_set}
  else
  {
    rin = max(In)
    inf_set = est_clean_set[Si[1:rin]]
    clean_set = setdiff(est_clean_set,inf_set)
  }

  list(clean_set=clean_set)
}
