\name{vhidetify}
\alias{vhidetify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute the single influence measure to validate the estimated influential set 
}
\description{
This function is part of the algorithm which indentify multiple influential observations in high dimension.
It applies a single detection technique to validate the  estimated influential set.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  vhidetify(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, asymvec, inf_set, non_inf_set, alpha)
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

  \item{asymvec}{ Numertic vector of the asymmetric values. The three points of the quartile (0.25,0.5,0.75) give best results.}

  \item{inf_set}{Estimated of influential set}

  \item{non_inf_set}{Estimated of non-influential set}

  \item{alpha}{Scalar value for the Type I error level. The default value is 0.05}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
Vector of index values containing the influential observations id.
\item{inf_setfinal }{The final set of influential observations }
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
inf_set = 1:20
non_inf_set = 21:100
final_inf_set = vhidetify(x, y, xquant, yquant, inv_rob_sdx, rob_sdy, asymvec, inf_set, non_inf_set, alpha)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
