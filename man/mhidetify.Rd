\name{mhidetify}
\alias{mhidetify}
\title{
Compute the asymmetric influence measure for high dimensional regression.
}
\description{
Compute the asymmetric influential measure to identify influential observations in high dimensional regression.
}
\usage{
mhidetify(x, y, number_subset, size_subset, asymvec, ep=0.1, alpha_swamp, alpha_mask, alpha_validate)
}
\arguments{
  \item{x}{ Matrix with the values of the predictors.}
  \item{y}{ Numertic vector of the response variable.}
  \item{number_subset}{Number of random subsets.}
  \item{size_subset}{ Size of the random subsets. The default is half of the initial sample size.}
  \item{asymvec}{ Numertic vector of the asymmetric values.}
  \item{ep}{ Threshold value to ensure that the estimated clean set is not empty. The default value is 0.1.}
  \item{alpha_swamp}{Significance level for the swamping stage.}
  \item{alpha_mask}{Significance level for the masking stage.}
  \item{alpha_validate}{Significance level for the validation stage.}
}

\value{
A dataframe with two variables.
\item{ind }{Index of the subjects of the sample}
\item{outlier_ind }{Influential observations indicator: 1 is influential and 0 otherwise}
}
\references{
Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2020). \emph{Asymmetric
influence measure for high dimensional regression}. Communications in Statistics - Theory and
Methods.

Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2021). \emph{An algorithm-based multiple detection influence measure for high dimensional regression using expectile}. 
arXiv: 2105.12286 [stat]. arXiv: 2105.12286.

}
\author{Amadou Barry \email{barryhafia@gmail.com}}

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
number_subset = 5
size_subset = 100/2
asymvec=c(0.25,0.5,0.75)
out = mhidetify(x, youtlier, number_subset, size_subset, asymvec, ep=0.1, alpha_swamp = 0.1, 
alpha_mask = 0.01, alpha_validate = 0.01)
}



