\name{hidetify}
\alias{hidetify}
\author{Amadou Barry \email{barryhafia@gmail.com}}
\title{ Identify the influential observations in high dimensional regression}
\description{
This function proposes two detection methods to identify influential
observations in high dimensional regression setting:
a single detection technique and a multiple detection technique.
}
\usage{
hidetify(
  predictors, 
  response, 
  nsample=5, 
  ssize=floor(length(response)/2), 
  vtau=c(0.25,0.5,0.75), 
  alpha_shide = 0.05, 
  alpha_swamp = 0.1, 
  alpha_mask = 0.01, 
  alpha_validate = 0.01,
  method = c("single", "multiple")
  )
}
\arguments{
  \item{predictors}{Matrix with the values of the predictors.}
  
  \item{response}{Numertic vector of the response variable.}
  
  \item{nsample}{Number of random subsets, default is 5.}
  
  \item{ssize}{Size of the random subsets. The default is 
  half of the initial sample size.}
  
  \item{vtau}{Numertic vector of the asymmetric values.
  It is suggested to choose 3 asymetric points within the 
  quartile.}
  
  \item{alpha_shide}{Significance level for the single detection method.
  The default is set to 0.05.}
  
  \item{alpha_swamp}{Significance level for the swamping stage.
  The default is set to 0.1.}
  
  \item{alpha_mask}{Significance level for the masking stage.
  The default is set to 0.01.}
  
  \item{alpha_validate}{Significance level for the validation stage.
  The default is set to 0.01.}
  
  \item{method}{The parameter option for the detection method. 
  There is two options: single ou multiple.}
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
vtau  <- c(0.25,0.5,0.75)

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
nsample = 5

# the size of the random subset
ssize = 100/2

# initial clean set
est_clean_set = 1:100

# the significance level for the single detection method
alpha_shide = 0.05

# the significance level for the swamping stage
alpha_swamp = 0.1

# the significance level for the masking stage
alpha_mask = 0.01

# the significance level for the validation stage
alpha_validate = 0.01

# the method of detection
method = "single"

out = 
  hidetify(
    x, 
    youtlier, 
    nsample, 
    ssize, 
    vtau, 
    alpha_shide, 
    alpha_swamp, 
    alpha_mask, 
    alpha_validate, 
    method = "single")
}

\references{
Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2020). \emph{Asymmetric
influence measure for high dimensional regression}. Communications in Statistics - Theory and
Methods.

Barry, A., Bhagwat, N., Misic, B., Poline, J.-B., and Greenwood, C. M. T. (2021). \emph{An algorithm-based multiple detection influence measure for high dimensional regression using expectile}. 
arXiv: 2105.12286 [stat]. arXiv: 2105.12286.
}
