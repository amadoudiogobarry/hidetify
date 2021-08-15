\name{rcpp_setdiff}
\alias{rcpp_setdiff}
\docType{package}
\title{Arma Set Difference of Subsets}
\description{
\code{rcpp_setdiff} calculates the (nonsymmetric) set difference of 
subsets of a probability space.
}
\usage{
rcpp_setdiff(x, y)
}
\arguments{

  \item{x}{a vector of integers.}
  
  \item{y}{a vector of integers.}
  
}

\value{
A subset of a probability space of the same type as its arguments.
}

\examples{
\dontrun{
x = 1:20
y = 1:12
rcpp_setdiff(x, y)
}
}
