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
DetectEpi(inputfile_X, inputfile_Y, inputfile_COV, MAX_ROW, MAX_COL, IS_EXIST_COV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfile_X}{
  X's filename
%%     ~~Describe \code{x} here~~
}
\item{inputfile_Y}{Y's filename}
\item{inputfile_COV}{COV's filename}
\item{MAX_ROW}{The row number of X}
\item{MAX_COL}{The col number of X}
\item{IS_EXIST_COV}{Whether consider cov's effect}
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

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ DetectEpi }% __ONLY ONE__ keyword per line
