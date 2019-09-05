\name{ease_masking}
\alias{ease_masking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the max of the sum of a sequence of asymmetric influence measure
%%  ~~function to do ... ~~
}
\description{
This function is part of the algorithm which identify multiple influential observations in high dimension.
It computes the max statistic of the asymmetric influence measure to ease the masking effect
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ease_masking(x,y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset, size_subset, est_clean_set, asymvec,alpha=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Matrix with the values of the predictors.
%%     ~~Describe \code{x} here~~
}
  \item{y}{ Numertic vector of the response variable.   }

  \item{xquant }{ Matrix with the quantiles of the predictors. }

  \item{yquant }{ Numertic vector of the quantiles of the response variable.}

  \item{inv_rob_sdx}{ Numertic vector of the inverse of the median absolute deviation of the predictors. }

  \item{rob_sdy}{ Median absolute deviation of the response variable. }

  \item{number_subset}{ Number of random subsets.}

  \item{size_subset}{ Size of the random subsets. The default is half of the initial sample size.}

  \item{est_clean_set}{ The subjct id of the estimated clean subset. The default is the initial sample}

  \item{asymvec}{ Numertic vector of the asymmetric values. The three points of the quartile (0.25,0.5,0.75) give best results}

  \item{alpha}{Scalar value for the Type I error level. The default value is 0.05}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A index vector identifying the estimated non-influential observations using an agressive approach
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Amadou Barry \email{barryhafia@gmail.com}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
youtlier[1:10] <- x[1:10,]%*%(beta_param +  1.2*gama_param)  + error_var[1:n_tilde]
xquant <- apply(x,2,quantile,asymvec)
yquant <- quantile(y,asymvec)
inv_rob_sdx <- 1/apply(x,2,mad)
rob_sdy = mad(y)
number_subset = 5
size_subset = 100/2
est_clean_set = 1:100
alpha = 0.05
est_clean_set_ease_masking = ease_masking(x,y, xquant, yquant, inv_rob_sdx, rob_sdy, number_subset,size_subset,est_clean_set,asymvec,alpha)

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
