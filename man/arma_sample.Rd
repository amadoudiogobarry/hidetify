\name{arma_sample}
\alias{arma_sample}
\docType{package}
\title{Arma Samples and Permutations}
\description{
\code{arma_sample} takes a sample of the specified size from 
the elements of \code{x} using either with or without replacement.
}
\usage{
arma_sample(x, size, replace)
}
\arguments{

  \item{x}{a vector of elements from which to choose.}
  
  \item{size}{a non-negative integer giving the number of items to choose.}
  
  \item{replace}{should sampling be with replacement.}
}

\value{
A vector of length \code{size} with elements drawn from \code{x}. 
}

\examples{
\dontrun{
x = 1:12
size = 5
arma_sample(x, size, replace = TRUE)
}
}
