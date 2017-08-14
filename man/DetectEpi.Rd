\name{DetectEpi}
\alias{DetectEpi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Detect Epi
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Detect Epi
}
\usage{
DetectEpi(inputfile_X, inputfile_Y, inputfile_COV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfile_X}{
  X's filename
%%     ~~Describe \code{x} here~~
}
\item{inputfile_Y}{Y's filename}
\item{inputfile_COV}{COV's filename. It can be NULL when does not consider covariates, otherwise, it would consider 4 covariates default.}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The p-value of iteraction.
}
\author{
%%  ~~who you are~~
Lijun Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Consider covariates

## DetectEpi("data/Data.txt", "data/Y.FA.txt", NULL)

## DO NOT consider covariates.
## DetectEpi("data/Data.txt", "data/Y.FA.txt", "data/cov.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ DetectEpi }% __ONLY ONE__ keyword per line
