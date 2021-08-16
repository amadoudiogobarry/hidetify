\name{vhidetify}
\alias{vhidetify}
\author{Amadou Barry \email{barryhafia@gmail.com}}
\title{Compute the single influence measure to validate the estimated influential set.}
\description{
This function is part of the algorithm which indentify multiple influential 
observations in high dimension. 
It applies a single detection technique to validate 
the  estimated influential set.
}

\usage{
  vhidetify(
    x, 
    y, 
    xquant, 
    yquant, 
    inv_rob_sdx, 
    rob_sdy, 
    asymvec, 
    inf_set, 
    non_inf_set, 
    alpha
    )
}

\arguments{
  \item{x}{Matrix of the predictors.}
  
  \item{y}{Numertic vector of the response variable.}
  
  \item{xquant}{Matrix of the quantiles of the predictors.}
  
  \item{yquant}{Numertic vector of the quantiles of the response variable.}
  
  \item{inv_rob_sdx}{Numertic vector of the inverse of 
  the median absolute deviation of the predictors.}
  
  \item{rob_sdy}{Median absolute deviation of the response variable.}
  
  \item{asymvec}{Numeric vector of the asymmetric values or percentiles.
    It is suggested to choose 3 asymetric points within the 
  quartile.}
  
  \item{inf_set}{Estimated set of influential observations}
  
  \item{non_inf_set}{Estimated set of non-influential observations}
  
  \item{alpha}{Significance level. }
}

\value{
Vector of index values containing the influential observations id.
\item{inf_setfinal}{The final set of influential observations.}
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

# the quantile of the predictors
xquant <- apply(x,2,quantile,asymvec)

# the quantile of contaminated response variable
yquant <- quantile(youtlier,asymvec)

# the inverse of the mad predictors
inv_rob_sdx <- 1/apply(x,2,mad)

# the mad contaminated response variable
rob_sdy = mad(youtlier)

# the number of random subsets
number_subset = 5

# the size of random subsets
size_subset = 100/2

# the initial clean set
est_clean_set = 1:100

# influential set
inf_set = 1:20

# non-influential set
non_inf_set = 21:100

# the significance level
alpha = 0.05

final_inf_set = 
  vhidetify(
    x, 
    youtlier, 
    xquant, 
    yquant, 
    inv_rob_sdx, 
    rob_sdy, 
    asymvec, 
    inf_set, 
    non_inf_set, 
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
