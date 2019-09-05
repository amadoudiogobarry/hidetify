\name{shidetify}
\alias{shidetify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the single influence measure for high dimensional regression.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the single influence measure to identify influential observations in high dimensional regression
}
\usage{
shidetify(x, y, asymvec, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Matrix with the values of the predictors.
%%     ~~Describe \code{x} here~~
}
  \item{y}{ Numertic vector of the response variable.   }

  \item{asymvec}{ Numertic vector of the asymmetric values. The three points of the quartile (0.25,0.5,0.75) give best results.}

  \item{alpha}{Scalar value for the Type I error level. The default value is 0.05}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe with two variables.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ind }{Index of the observations}
\item{outlier_ind }{Influential observations indicator: 1 if influential and 0 otherwise}
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
alpha = 0.05
df_single_influential = shidetify(x, y, asymvec, alpha)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
