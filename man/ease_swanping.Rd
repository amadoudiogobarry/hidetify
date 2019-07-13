\name{ease_swanping}
\alias{ease_swanping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the min of the min of a sequence of asymmetric influence measure
%%  ~~function to do ... ~~
}
\description{
This function is part of the algorithm which identify multiple influential observations in high dimension linear regression.
It computes the min of the min of the asymmetric influence measure to ease the swanping effect
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ease_swamping(x,y,xquant, yquant, inv_rob_sdx, rob_sdy,number_subset,size_subset,est_clean_set,asymvec,ep,alpha)
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

  \item{est_clean_set}{ The subjct id of the estimated clean subset. The default is the initial sample.}

  \item{asymvec}{ Numertic vector of the asymmetric values. The three points of the quartile (0.25,0.5,0.75) give best results.}

  \item{ep}{ Threshold value to ensure that the estimated clean set is not empty. The default value is 0.1.}

  \item{alpha}{Scalar value for the Type I error level. The default value is 0.05}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of index, identifying the estimated non-influential subjects
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
Amadou Barry
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

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
