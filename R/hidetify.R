#### main function
hidetify = function(predictors, response, nsample=5, ssize, vtau=c(0.25,0.5,0.75), valpha = 0.05, method = c("single", "multiple"))
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

  if (is.na(valpha) || valpha >= 1 || valpha <= 0 ) {
    warning("valpha has to be a value between 0 and 1")
  } else {
    alpha <- valpha
  }

  if (method == "single") {

    dfout <- shidetify(x, y, asymvec, alpha)

  } else if (method == "multiple") {

    if (any(nsample %% 1 != 0 || ssize %% 1 != 0)){
      warning("nsample and ssize are integers")
    } else {
      number_subset <- nsample
      size_subset <- ssize
    }
    dfout <- mhidetify(x, y, number_subset, size_subset, asymvec, ep=0.1, alpha)
  }

  dfout

}
