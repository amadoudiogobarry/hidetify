\name{shidetify}
\alias{shidetify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the single influential points identifier measure
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the single influential points measure to identify influential observations in high dimensional regression
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
\item{ind }{Index of the subjects of the sample}
\item{outlier_ind }{Influential observations indicator: 1 is influential and 0 otherwise}
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
