\name{ease_swanping}
\alias{ease_swanping}
\title{ Compute the min of the min of a sequence of asymmetric influence measure}
\description{
This function is part of the algorithm which identify multiple influential observations in high dimension linear regression.It computes the min of the min of the asymmetric influence measure to ease the swanping effect}
\usage{
ease_swamping(x,y,xquant, yquant, inv_rob_sdx, rob_sdy,number_subset,size_subset,est_clean_set,asymvec,ep=0.1,alpha)
}

\arguments{
  \item{x}{ Matrix with the values of the predictors.}
  \item{y}{ Numertic vector of the response variable.   }

  \item{xquant }{ Matrix with the quantiles of the predictors. }

  \item{yquant }{ Numertic vector of the quantiles of the response variable.}

  \item{inv_rob_sdx}{ Numertic vector of the inverse of the median absolute deviation of the predictors. }

  \item{rob_sdy}{ Median absolute deviation of the response variable. }

  \item{number_subset}{ Number of random subsets.}

  \item{size_subset}{ Size of the random subsets. The default is half of the initial sample size.}

  \item{est_clean_set}{ The subjct id of the estimated clean subset. The default is the initial sample.}

  \item{asymvec}{ Numertic vector of the asymmetric values.}

  \item{ep}{ Threshold value to ensure that the estimated clean set is not empty. The default value is 0.1.}

  \item{alpha}{Significance level. }
}

\value{
A index vector identifying the estimated non-influential observations using a conservative approach
}
\author{
Amadou Barry \email{barryhafia@gmail.com}}

\examples{
## Simulate a dataset where the first 10 observations are influentials
require("MASS")
asymvec  <- c(0.25,0.5,0.75)
beta_param <- c(3,1.5,0,0,2,rep(0,1000-5))
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

set.seed(13)
x  <- mvrnorm(100, rep(0, 1000), sigmain)
error_var <- rnorm(100)
y  <- x %*% beta_param + error_var
# Generate influential observations
youtlier = y
youtlier[1:10] <- x[1:10,]%*%(beta_param +  1.2*gama_param)  + error_var[1:10]
xquant <- apply(x,2,quantile,asymvec)
yquant <- quantile(y,asymvec)
inv_rob_sdx <- 1/apply(x,2,mad)
rob_sdy = mad(y)
number_subset = 5
size_subset = 100/2
est_clean_set = 1:100
alpha = 0.05
est_clean_set_ease_swamping = ease_swamping(x,youtlier,xquant, yquant, inv_rob_sdx, rob_sdy,number_subset,size_subset,est_clean_set,asymvec,ep=0.1,alpha)
}

\references{
Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2020). \emph{Asymmetric
influence measure for high dimensional regression}. Communications in Statistics - Theory and
Methods.

Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2021). \emph{An algorithm-based multiple detection influence measure for high dimensional regression using expectile}. 
arXiv: 2105.12286 [stat]. arXiv: 2105.12286.
}
