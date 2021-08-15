\name{mhidetify}
\alias{mhidetify}
\author{Amadou Barry \email{barryhafia@gmail.com}}
\title{Multiple detection asymmetric influential measure 
for high dimensional linear regression.}
\description{
The function computes the asymmetric influential measure to identify influential observations 
in high dimensional linear regression using the multiple detection approach.
}
\usage{
mhidetify(
  x, 
  y, 
  number_subset, 
  size_subset, 
  asymvec, 
  ep=0.1, 
  alpha_swamp, 
  alpha_mask, 
  alpha_validate
  )
}

\arguments{
  \item{x}{Matrix of the predictors.}
  
  \item{y}{Numertic vector of the response variable.}
  
  \item{number_subset}{Number of random subsets, default is 5.}
  
  \item{size_subset}{Size of the random subsets. 
  The default is half of the initial sample size.}
  
  \item{asymvec}{Numeric vector of the asymmetric values.
    It is suggested to choose 3 asymetric points within the 
  quartile.}
  
  \item{ep}{Threshold value to ensure that the estimated clean set is not empty. 
  The default value is 0.1.}
  
  \item{alpha_swamp}{Significance level for the swamping stage.}
  
  \item{alpha_mask}{Significance level for the masking stage.}
  
  \item{alpha_validate}{Significance level for the validation stage.}
}

\value{
A dataframe with two variables.
\item{ind }{Index of the subjects of the sample}
\item{outlier_ind }{Influential observations indicator: 1 is influential and 0 otherwise}
}

\examples{
## Simulate a dataset where the first 10 observations are influentials
require("MASS")
# the vector of asymetric point
asymvec  <- c(0.25,0.5,0.75)

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


# number of random subsets
number_subset = 5

# the size of the random subset
size_subset = 100/2

# the significance level for the swamping stage
alpha_swamp = 0.1

# the significance level for the masking stage
alpha_mask = 0.01

# the significance level for the validation stage
alpha_validate = 0.01

# Threshold value to ensure that the estimated clean set is not empty. 
ep = 0.1

out = 
  mhidetify(
    x, 
    youtlier, 
    number_subset, 
    size_subset, 
    asymvec, 
    ep, 
    alpha_swamp, 
    alpha_mask,
    alpha_validate)
}


\references{
Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2020). \emph{Asymmetric
influence measure for high dimensional regression}. Communications in Statistics - Theory and
Methods.

Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2021). \emph{An algorithm-based multiple detection influence measure for high dimensional regression using expectile}. 
arXiv: 2105.12286 [stat]. arXiv: 2105.12286.

}


