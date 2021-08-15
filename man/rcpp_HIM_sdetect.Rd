\name{rcpp_HIM_sdetect}
\alias{rcpp_HIM_sdetect}
\docType{package}
\title{Arma Single Detection Statistic}
\description{
\code{rcpp_HIM_sdetect} computes the influence measure
statistic for the single detection technique and for \eqn{\tau=0.5}.
}
\usage{
rcpp_HIM_sdetect(
  x, 
  y, 
  xquant, 
  yquant, 
  inv_rob_sdx, 
  rob_sdy,
  inf_set, 
  non_inf_set)

}

\arguments{

  \item{x}{a matrix of elements.}
  
  \item{y}{a vector of elements.}
  
  \item{xquant}{quantiles of the columns of \code{x} stacked in the matrix \code{xquant}.}
  
  \item{yquant}{quantiles vector of the vector \code{y}.}
  
  \item{inv_rob_sdx}{inverse of the median absolute deviation of the matrix \code{x}.}
  
  \item{rob_sdy}{median absolute deviation of the vector \code{y}.}

  \item{inf_set}{influential set.}
  
  \item{non_inf_set}{non-influential set.}

}

\value{
A vector of the influence measure statistic.  
}

\examples{

\dontrun{
## Simulate a dataset where the first 10 observations are influentials
require("MASS")

# the parameter of interest
beta_param <- c(3,1.5,0,0,2,rep(0,1000-5))

# the contamination parameter 
gama_param <- c(0,0,1,1,0,rep(1,1000-5))

# Covariance matrice for the predictors distribution 
sigmain <- diag(rep(1,1000))
for (i in 1:1000)
{
  for (j in i:1000) 
  {
    sigmain[i,j] <- 0.5^(abs(j-i))
    sigmain[j,i] <- sigmain[i,j]
  }
}

# set the seed
set.seed(13)

# the predictor matrix
x  <- mvrnorm(100, rep(0, 1000), sigmain)

# the error variable
error_var <- rnorm(100)

# the response variable
y  <- x \%*\% beta_param + error_var
y = as.numeric(y)

### Generate influential observations

# the contaminated response variable
youtlier = y
youtlier[1:10] <- x[1:10,] \%*\% (beta_param +  1.2*gama_param)  + error_var[1:10]
youtlier = as.numeric(youtlier)

# the quantile of the predictors
xquant <- apply(x,2,quantile,0.5)

# the quantile of contaminated response variable
yquant <- quantile(youtlier,0,5)

# the inverse of the mad predictors
inv_rob_sdx <- 1/apply(x,2,mad)

# the mad contaminated response variable
rob_sdy = mad(youtlier)

# influential set
inf_set = 1:20

# non-influential set
non_inf_set = 21:100

rcpp_HIM_sdetect(x, youtlier, xquant, yquant, inv_rob_sdx, rob_sdy, inf_set, non_inf_set)

}
}
