# Instructions

## install Rcpp packages

```
install.packages('Rcpp')
```

## install gsl

[https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)

### ubuntu


```
cd gsl-****
./configure
make
make install
```

### windows 
todo



## usage

```
library(Rcpp)
sourceCpp('main_r.cpp')
DetectEpi(XFILE, YFILE, COVFILE, MAXROW, MAXCOL, FALSE)
## if consider covariate, the last parameter is TRUE, otherwise FALSE
```
