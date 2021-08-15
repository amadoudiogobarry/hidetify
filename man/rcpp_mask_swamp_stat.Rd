\name{rcpp_mask_swamp_stat}
\alias{rcpp_mask_swamp_stat}
\docType{package}
\title{Arma Masking and Swamping Statistics}
\description{
\code{rcpp_mask_swamp_stat} computes the \code{min} of the \code{min} 
and the \code{max} of the \code{sum} statistics of the asymmetric influence measure
which represent the influential score of the observations in the swamping and 
masking steps, respectively.
}
\usage{
rcpp_mask_swamp_stat(
  x, 
  y, 
  xquant, 
  yquant, 
  inv_rob_sdx, 
  rob_sdy,
  number_subset, 
  size_subset, 
  est_clean_set, 
  asymvec )
}
\arguments{

  \item{x}{a matrix of elements.}
  
  \item{y}{a vector of elements.}
  
  \item{xquant}{quantiles of the columns of \code{x} stacked in the matrix \code{xquant}.}
  
  \item{yquant}{quantiles vector of the vector \code{y}.}
  
  \item{inv_rob_sdx}{inverse of the median absolute deviation of the matrix \code{x}.}
  
  \item{rob_sdy}{median absolute deviation of the vector \code{y}.}
  
  \item{number_subset}{number of random subsets.}
  
  \item{size_subset}{size of random subsets.}
  
  \item{est_clean_set}{estimated cleaned set.}
  
  \item{asymvec}{vector of asymetric points or percentiles.}
}

\value{
A matrix with the \code{min} of the \code{min} 
and the \code{max} of the \code{sum} statistics
in the first and second column for each observation, respectively. 
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

rcpp_mask_swamp_stat(
  x, 
  y, 
  xquant, 
  yquant, 
  inv_rob_sdx, 
  rob_sdy,
  number_subset, 
  size_subset, 
  est_clean_set, 
  asymvec )
}
}
