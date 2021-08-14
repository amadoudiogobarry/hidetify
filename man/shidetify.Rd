\name{shidetify}
\alias{shidetify}
\title{
Compute the single influence measure for high dimensional regression.
}
\description{
Compute the single influence measure to identify influential observations in high dimensional regression
}
\usage{
shidetify(x, y, asymvec, alpha)
}

\arguments{
  \item{x}{ Matrix with the values of the predictors.}
  \item{y}{ Numertic vector of the response variable.}
  \item{asymvec}{ Numertic vector of the asymmetric values.}
  \item{alpha}{Significance level.}
}

\value{
A dataframe with two variables.
\item{ind }{Index of the observations}
\item{outlier_ind }{Influential observations indicator: 1 if influential and 0 otherwise}
}
\references{
Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2020). \emph{Asymmetric
influence measure for high dimensional regression}. Communications in Statistics - Theory and
Methods.

Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2021). \emph{An algorithm-based multiple detection influence measure for high dimensional regression using expectile}. 
arXiv: 2105.12286 [stat]. arXiv: 2105.12286.

}
\author{
Amadou Barry \email{barryhafia@gmail.com} }

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
alpha = 0.05
df_single_influential = shidetify(x, youtlier, asymvec, alpha)
}

