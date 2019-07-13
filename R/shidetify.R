shidetify <- function(x, y, asymvec, alpha)
{

  xquant <- apply(x,2,quantile,asymvec); yquant <- quantile(y,asymvec); inv_rob_sdx <- 1/apply(x,2,mad); rob_sdy = mad(y);
  df <- length(asymvec)
  n <- nrow(x)
  asymHimMat <- rcpp_shidetify(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, asymvec, row_indice)
  vectHim_new = apply(asymHimMat,1,sum)
  pv_inf = pchisq(vectHim_new, df=df, lower.tail=FALSE)
  Spv_inf=sort.int(pv_inf,index.return=TRUE)
  Si=Spv_inf$ix
  dp=Spv_inf$x-alpha/n

  In = which(dp<=0)
  rin=max(In)
  ind_index = 1:n
  inf_set = ind_index[Si[1:rin]]

  dataframe <- data.frame(ind = ind_index, outlier_ind = rep(0,n))
  dataframe$outlier_ind[which(dataframe$ind%in%inf_set)] = 1
  dataframe

}







