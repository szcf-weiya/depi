# depi

A R package for detecting interaction of epi.

Use the following command in R to install the latest development version from Github.

```
devtools::install_github('szcf-weiya/depi')
```

## build log @ win

### 2017.09.08 
1. git clone repo to win
2. open restudio and try to rebuild, alert to install rtools.
3. setup lib and includes of gsl in Makevars.win
4. rename stderr
5. check cannot pass, because it will compile with arch i386 and x64. Check on arch i386 can succeed, while it failed on arch x64.
```
* installing *source* package 'depi' ...
** libs

*** arch - i386
C:/RBuildTools/3.4/mingw_32/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c DetectEpi.cpp -o DetectEpi.o
DetectEpi.cpp:205:0: warning: ignoring #pragma omp sections [-Wunknown-pragmas]
 # pragma omp sections
 ^
DetectEpi.cpp:207:0: warning: ignoring #pragma omp section [-Wunknown-pragmas]
 # pragma omp section
 ^
DetectEpi.cpp:209:0: warning: ignoring #pragma omp section [-Wunknown-pragmas]
 # pragma omp section
 ^
DetectEpi.cpp:211:0: warning: ignoring #pragma omp section [-Wunknown-pragmas]
 # pragma omp section
 ^
DetectEpi.cpp:217:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
     # pragma omp parallel for schedule(dynamic) // faster!!!
 ^
DetectEpi.cpp: In function 'void readData(std::string, gsl_matrix*, int, int, std::vector<std::basic_string<char> >&, std::vector<std::basic_string<char> >&)':
DetectEpi.cpp:89:24: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
   for (size_t i = 0; i < col; i++)
                        ^
DetectEpi.cpp:96:24: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
   for (size_t i = 0; i < row; i++)
                        ^
DetectEpi.cpp:102:28: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
       for (size_t j = 0; j < col; j++)
                            ^
DetectEpi.cpp: In function 'Rcpp::List DetectEpi(SEXP, SEXP, SEXP)':
DetectEpi.cpp:203:21: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
   for (int i = 0; i < ncol; i++)
                     ^
DetectEpi.cpp:218:25: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
     for (int j = i+1; j < ncol; j++)
                         ^
DetectEpi.cpp:189:10: warning: unused variable 'min_pair_geno' [-Wunused-variable]
   double min_pair_geno;
          ^
C:/RBuildTools/3.4/mingw_32/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c RcppExports.cpp -o RcppExports.o
C:/RBuildTools/3.4/mingw_32/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c rcpp_hello.cpp -o rcpp_hello.o
C:/RBuildTools/3.4/mingw_32/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c readPED.cpp -o readPED.o
C:/RBuildTools/3.4/mingw_32/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c safeGetline.cpp -o safeGetline.o
C:/RBuildTools/3.4/mingw_32/bin/g++ -shared -s -static-libgcc -o depi.dll tmp.def DetectEpi.o RcppExports.o rcpp_hello.o readPED.o safeGetline.o -LC:\Users\weiya\Documents\GitHub\GSLwin\lib -lgsl -lgslcblas -Ld:/Compiler/gcc-4.9.3/local330/lib/i386 -Ld:/Compiler/gcc-4.9.3/local330/lib -LC:/PROGRA~1/R/R-34~1.1/bin/i386 -lR
Warning messages:
1: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
2: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
3: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
4: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
5: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
installing to C:/Users/weiya/Documents/GitHub/win_depi/depi.Rcheck/depi/libs/i386

*** arch - x64
C:/RBuildTools/3.4/mingw_64/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c DetectEpi.cpp -o DetectEpi.o
DetectEpi.cpp:205:0: warning: ignoring #pragma omp sections [-Wunknown-pragmas]
 # pragma omp sections
 ^
DetectEpi.cpp:207:0: warning: ignoring #pragma omp section [-Wunknown-pragmas]
 # pragma omp section
 ^
DetectEpi.cpp:209:0: warning: ignoring #pragma omp section [-Wunknown-pragmas]
 # pragma omp section
 ^
DetectEpi.cpp:211:0: warning: ignoring #pragma omp section [-Wunknown-pragmas]
 # pragma omp section
 ^
DetectEpi.cpp:217:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
     # pragma omp parallel for schedule(dynamic) // faster!!!
 ^
DetectEpi.cpp: In function 'void readData(std::string, gsl_matrix*, int, int, std::vector<std::basic_string<char> >&, std::vector<std::basic_string<char> >&)':
DetectEpi.cpp:89:24: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
   for (size_t i = 0; i < col; i++)
                        ^
DetectEpi.cpp:96:24: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
   for (size_t i = 0; i < row; i++)
                        ^
DetectEpi.cpp:102:28: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
       for (size_t j = 0; j < col; j++)
                            ^
DetectEpi.cpp: In function 'Rcpp::List DetectEpi(SEXP, SEXP, SEXP)':
DetectEpi.cpp:203:21: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
   for (int i = 0; i < ncol; i++)
                     ^
DetectEpi.cpp:218:25: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
     for (int j = i+1; j < ncol; j++)
                         ^
DetectEpi.cpp:189:10: warning: unused variable 'min_pair_geno' [-Wunused-variable]
   double min_pair_geno;
          ^
C:/RBuildTools/3.4/mingw_64/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c RcppExports.cpp -o RcppExports.o
C:/RBuildTools/3.4/mingw_64/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c rcpp_hello.cpp -o rcpp_hello.o
C:/RBuildTools/3.4/mingw_64/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c readPED.cpp -o readPED.o
C:/RBuildTools/3.4/mingw_64/bin/g++  -I"C:/PROGRA~1/R/R-34~1.1/include" -DNDEBUG -I"C:\Users\weiya\Documents\GitHub\GSLwin\include" -I../inst/include -I"C:/Users/weiya/R/win-library/3.4/Rcpp/include"   -I"d:/Compiler/gcc-4.9.3/local330/include"     -O2 -Wall  -mtune=core2 -c safeGetline.cpp -o safeGetline.o
C:/RBuildTools/3.4/mingw_64/bin/g++ -shared -s -static-libgcc -o depi.dll tmp.def DetectEpi.o RcppExports.o rcpp_hello.o readPED.o safeGetline.o -LC:\Users\weiya\Documents\GitHub\GSLwin\lib -lgsl -lgslcblas -Ld:/Compiler/gcc-4.9.3/local330/lib/x64 -Ld:/Compiler/gcc-4.9.3/local330/lib -LC:/PROGRA~1/R/R-34~1.1/bin/x64 -lR
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: skipping incompatible C:\Users\weiya\Documents\GitHub\GSLwin\lib/libgsl.a when searching for -lgsl
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: skipping incompatible C:\Users\weiya\Documents\GitHub\GSLwin\lib\libgsl.a when searching for -lgsl
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: skipping incompatible C:\Users\weiya\Documents\GitHub\GSLwin\lib/libgsl.a when searching for -lgsl
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: cannot find -lgsl
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: skipping incompatible C:\Users\weiya\Documents\GitHub\GSLwin\lib/libgslcblas.a when searching for -lgslcblas
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: skipping incompatible C:\Users\weiya\Documents\GitHub\GSLwin\lib\libgslcblas.a when searching for -lgslcblas
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: skipping incompatible C:\Users\weiya\Documents\GitHub\GSLwin\lib/libgslcblas.a when searching for -lgslcblas
C:/RBuildTools/3.4/mingw_64/bin/../lib/gcc/x86_64-w64-mingw32/4.9.3/../../../../x86_64-w64-mingw32/bin/ld.exe: cannot find -lgslcblas
collect2.exe: error: ld returned 1 exit status
Warning messages:
1: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
2: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
3: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
4: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
5: In FUN(X[[i]], ...) : this requires 'nm' to be on the PATH
no DLL was created
ERROR: compilation failed for package 'depi'
* removing 'C:/Users/weiya/Documents/GitHub/win_depi/depi.Rcheck/depi'
```