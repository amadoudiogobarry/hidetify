#### main function
hidetify = function(predictors, response, nsample=5, ssize=floor(nrow(x)/2), vtau=c(0.25,0.5,0.75), alpha_shide = 0.05, 
                    alpha_swamp = 0.1, alpha_mask = 0.01, alpha_validate = 0.01, method = c("single", "multiple"))
{

  method = match.arg(method)

  if (!is.matrix(predictors)){
    warning("predictors has to be a matrix object")
  } else {
    x <- predictors
  }

  if (!is.vector(response)){
    warning("response has to be a vector object")
  } else {
    y <- response
  }

  if (any(is.na(vtau)) || !is.vector(vtau) || any(vtau > 1) || any(vtau < 0)) {
    asymvec <- c(0.25,0.5,0.75)
  } else {
    asymvec <- vtau
  }



    
  if (method == "single") {
    
    dfout <- shidetify(x, y, asymvec, alpha_shide)
    
  } else if (method == "multiple") {
    
    if (any(nsample %% 1 != 0 || ssize %% 1 != 0)){
      warning("nsample and ssize are integers")
    } else {
      number_subset <- nsample
      size_subset <- ssize
    }
    dfout <- mhidetify(x, y, number_subset, size_subset, asymvec, ep=0.1, alpha_swamp, alpha_mask, alpha_validate)
    
  }
    



  dfout

}
