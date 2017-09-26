# epi

Use gsl in cpp to speed the computing, and call the cpp function in R by the packages Rcpp. In Linux, it is very easy to implement this procedure, while on windows we have to do more things.

## Linux

```
Rcpp::sourceCpp('./v_rcpp/main_r.cpp')
```

## Windows

### Compile GSL for Windows

1. [Download MinGW](www.mingw.org) and install.
2. In MinGW Installation Manager, choose all packages in basic setup and install.
3. In cmd, type 'msys.bat MSYS' to open a unix-like terminal.
4. [Download gsl](ftp://ftp.gnu.org/gnu/gsl/) to the file folder of MSYS.
5. In terminal, type

```
./configure
make
make install
```
6. copy 'local/lib' to 'R-3.4.1/bin/i386'; copy 'local/include' to 'R-3.4.1/include'.  


### Install Rtools

[Rtools](https://cran.r-project.org/bin/windows/Rtools/)

Remember to choose edit the path.

Now, you can call .cpp in R using the following code:
```
Rcpp::sourceCpp('./v_rcpp/main_r.cpp')
```

### Others

DO NOT use stderr in .cpp on Windows, I do not know why.
