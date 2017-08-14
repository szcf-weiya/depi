\name{eigenK}
\alias{eigenK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the eigen values and eigen vectors of matrix K.
}
\description{
K is a matrix constructed from a dataset. This function returns the eigen values and eigen vectors of matrix K.
}
\usage{
eigenK(inputFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{
  the file name of .ped format
%%     ~~Describe \code{inputFile} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{eigen values}{eigen values of matrix K}
\item{eigen vectors}{eigen vectors of matrix K}
}

\author{
Lijun Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

## eigenK("X.PED")
}

\keyword{ eigenK }
