\name{shidetify}
\alias{shidetify}
\author{Amadou Barry \email{barryhafia@gmail.com}}
\title{Single detection asymmetric influential measure 
for high dimensional linear regression.}
\description{
The function computes the asymmetric influential measure to identify influential observations 
in high dimensional linear regression using the single detection approach.
}
\usage{
shidetify(
  x, 
  y, 
  asymvec, 
  alpha
  )
}

\arguments{
  \item{x}{ Matrix of the predictors.}
  
  \item{y}{ Numertic vector of the response variable.}
  
  \item{asymvec}{Numeric vector of the asymmetric values.
    It is suggested to choose 3 asymetric points within the 
  quartile.}
  
  \item{alpha}{Significance level. }
}

\value{
A dataframe with two variables.
\item{ind }{Index of the observations}
\item{outlier_ind }{Influential observations indicator: 1 if influential and 0 otherwise}
}

\examples{
\dontrun{

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

# the significance level 
alpha = 0.05

df_single_influential = 
  shidetify(
    x, 
    youtlier, 
    asymvec, 
    alpha)
}
}

\references{
Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2020). \emph{Asymmetric
influence measure for high dimensional regression}. Communications in Statistics - Theory and
Methods.

Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2021). \emph{An algorithm-based multiple detection influence measure for high dimensional regression using expectile}. 
arXiv: 2105.12286 [stat]. arXiv: 2105.12286.

}

