### Nous avons remplace la fonction final_shrinking_swamping_effect par vhidetify. v for validate.
### Validate the influential observation set.


vhidetify <- function(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, asymvec, inf_set, non_inf_set, alpha){
  n = nrow(x)
  n_inf = length(inf_set)
  asymHimMat = rcpp_asymHIM_sdetect(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, asymvec, inf_set, non_inf_set) # rcpp function
  vectHim_new = apply(asymHimMat,1,sum)
  df = length(asymvec)
  pv_inf = pchisq(vectHim_new, df=df, lower.tail=FALSE)
  Spv_inf=sort.int(pv_inf,index.return=TRUE)
  Si=Spv_inf$ix
  dp=Spv_inf$x-alpha/n_inf

  In = which(dp<=0)
  if (length(In)==0){
    clean_setfinal = c(1:n)
    inf_setfinal = setdiff(c(1:n),clean_setfinal)
  } else{
    rin=max(In)
    inf_setfinal = inf_set[Si[1:rin]]
    clean_setfinal = setdiff(c(1:n),inf_setfinal)
  }

  list(inf_setfinal = inf_setfinal)
}
