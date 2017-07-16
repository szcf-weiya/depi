# Instructions

## install Rcpp and RcppGSL packages

```
install.packages('Rcpp')
install.packages('RcppGSL')
```

## install gsl

[https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)

## usage

```
library(Rcpp)
sourceCpp('main_r.cpp')
DetectEpi(XFILE, YFILE, COVFILE, MAXROW, MAXCOL, FALSE)
## if consider covariate, the last parameter is TRUE, otherwise FALSE
```
