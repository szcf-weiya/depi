// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<stdlib.h>
#include<stdio.h>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_combination.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_cdf.h>
#include<omp.h>

using namespace std;

//#define MAX_ROW 305
//#define MAX_COL 5000
#define COV_COL 4
#define MAX_PAIR_CORR 0.98
#define MAX_EPS 1e+12
#define MIN_P_VALUE 0.001


#define coef(i) (gsl_vector_get(coef, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))

// readData when unknown nrow and ncol
void readData(string FILE, vector<vector<double> > &mv, vector<string> &rowname, vector<string> &colname, size_t *nrow, size_t *ncol)
{
  ifstream input(FILE.c_str());
  string line;
  string tmp, tmpp;
  double val;
  *nrow = 0;
  *ncol = 0;
  // the first line
  getline(input, line);
  stringstream ss(line);
  while(!ss.eof())
  {
    ss >> tmp;
    colname.push_back(tmp);
    (*ncol)++;
  }
  //cout << colname.size() << endl;

  while (getline(input, line)) {
    mv[*nrow].resize(*ncol);
    if(mv.size() == *nrow+1)
      mv.resize(mv.size()*2);
    ss.clear();
    stringstream ss(line);
    ss >> tmp; // rowname
    rowname.push_back(tmp);
    for (size_t j = 0; j < *ncol; j++)
    {
      ss >> tmp;
      val = strtod(tmp.c_str(), NULL); // transfer str to double
      //mv[*nrow].push_back(val);
      mv[*nrow][j] = val;
    }
    (*nrow)++;
  }
  input.close();
  mv.resize(*nrow);
}

// read data with known nrow and known ncol
void readData(string FILE, gsl_matrix *m, int row, int col, vector<string> &rowname, vector<string> &colname)
{
  ifstream input(FILE.c_str());
  string tmp, tmpp;
  double val;
  // colnames
  string line;
  // the first line
  getline(input, line);
  stringstream ss(line);
  //ss >> tmp; // null no need!!!
  for (size_t i = 0; i < col; i++)
    {
      ss >> tmp;
      tmpp = tmp.substr(1, tmp.size() - 2);
      colname.push_back(tmpp);
    }

  for (size_t i = 0; i < row; i++)
    {
      getline(input, line);
      stringstream ss(line);
      ss >> tmp; // rowname
      rowname.push_back(tmp);
      for (size_t j = 0; j < col; j++)
      	{
      	  ss >> tmp;
      	  val = strtod(tmp.c_str(), NULL); // transfer str to double
      	  gsl_matrix_set(m, i, j, val);
      	}
    }
  input.close();
}

//[[Rcpp::export]]

Rcpp::List DetectEpi(SEXP inputfile_X, SEXP inputfile_Y, SEXP inputfile_COV, bool IS_EXIST_COV)
{
  // convert file name
  Rcpp::CharacterVector r_inputfile_X(inputfile_X);
  Rcpp::CharacterVector r_inputfile_Y(inputfile_Y);
  Rcpp::CharacterVector r_inputfile_COV(inputfile_COV);
  string file_X = Rcpp::as<string>(r_inputfile_X);
  string file_Y = Rcpp::as<string>(r_inputfile_Y);
  string file_COV = Rcpp::as<string>(r_inputfile_COV);

  // read data X
  //gsl_matrix *G = gsl_matrix_alloc(MAX_ROW, MAX_COL);
  vector<string> G_colname;
  vector<string> G_rowname;
  //readData(file_X.c_str(), G, MAX_ROW, MAX_COL, G_rowname, G_colname);
  vector<vector<double> > mv(1, vector<double>(1));
  size_t nrow, ncol;
  readData(file_X.c_str(), mv, G_rowname, G_colname, &nrow, &ncol);
  cout << ncol << endl;
  gsl_matrix *G = gsl_matrix_calloc(nrow, ncol);
  cout << mv.size() << " " << mv[0].size() << endl;
  cout << nrow << " " << ncol << endl;
  cout << G->size1 << " " << G->size2 << endl;
  cout << mv[0][0] << endl;
  cout << mv[0][10] << endl;
  for (size_t i = 0; i < nrow; i++)
    for (size_t j = 0; j < ncol; j++)
      gsl_matrix_set(G, i, j, mv[i][j]);

  int p;
  if (IS_EXIST_COV)
    p = 4 + COV_COL;
  else
    p = 4;
  gsl_matrix *X = gsl_matrix_alloc(nrow, p);
  if (IS_EXIST_COV)
  {
    gsl_matrix *coveriate = gsl_matrix_alloc(nrow, COV_COL);
    // read data cov
    vector<string> cov_colname;
    vector<string> cov_rowname;
    readData(file_COV.c_str(), coveriate, nrow, COV_COL, cov_rowname, cov_colname);
    gsl_vector *coveriate_vec = gsl_vector_alloc(nrow);
    for (int ip = 0; ip < COV_COL; ip++)
    {
      gsl_matrix_get_col(coveriate_vec, coveriate, ip);
      gsl_matrix_set_col(X, 4+ip, coveriate_vec);
    }
    gsl_vector_free(coveriate_vec);
    gsl_matrix_free(coveriate);
  }

  // read data Y
  ifstream input2(file_Y.c_str());
  double val;
  gsl_vector *Y = gsl_vector_alloc(nrow);
  for (size_t i = 0; i < nrow; i++)
  {
    input2 >> val;
    gsl_vector_set(Y, i, val);
  }
  input2.close();


  // fit
  gsl_vector *x0 = gsl_vector_alloc(nrow);
  gsl_vector_set_all(x0, 1);
  //  gsl_vector *x1_copy = gsl_vector_alloc(MAX_ROW);


  int max_val, min_val;

  double min_pair_geno;

  // cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  // cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";
  //RcppGSL::vector<double> r1, r2, r1_p, r2_p, r1r2_p;
  vector<string> r1, r2;
  vector<double> r1_p, r2_p, r1r2_p;


//  fprintf(output, "r1\tr2\tr1.p.value\tr2.p.value\tr1*r2.p.value\n");
  //# pragma omp parallel
  //# pragma omp for schedule(dynamic) // relatively slow
  gsl_matrix_set_col(X, 0, x0);
  gsl_vector *x1 = gsl_vector_alloc(nrow);
  cout << "hello "<< endl;
  for (int i = 0; i < ncol/100; i++)
  {
# pragma omp sections
{
# pragma omp section
  gsl_matrix_get_col(x1, G, i);
# pragma omp section
  max_val = gsl_vector_max(x1);
# pragma omp section
  min_val = gsl_vector_min(x1);
}
    if (max_val == min_val)
      continue;
    gsl_matrix_set_col(X, 1, x1);
    # pragma omp parallel for schedule(dynamic) // faster!!!
    for (int j = i+1; j < ncol; j++)
    {
      gsl_vector *x2 = gsl_vector_alloc(nrow);

      size_t flag, df;

      double stderr, t;
      double pvalue[3];
      df = nrow - p;

      gsl_matrix_get_col(x2, G, j);
      // if same
      max_val = gsl_vector_max(x2);
      min_val = gsl_vector_min(x2);
      if (max_val == min_val)
      {
        gsl_vector_free(x2);
        continue;
      }
      gsl_matrix_set_col(X, 2, x2);

      double res_corr = (float)gsl_stats_correlation(x1->data, 1, x2->data, 1, nrow);
      if (fabs(res_corr) >= MAX_PAIR_CORR)
        continue;
      gsl_vector *x3 = gsl_vector_alloc(nrow);
      gsl_matrix_get_col(x3, G, i);
      gsl_vector_mul(x3, x2);
      gsl_vector_free(x2);

      max_val = gsl_vector_max(x3);
      //	  min_val = gsl_vector_min(x3);
      if (max_val == 0)
      {
        gsl_vector_free(x3);
        continue;
      }

      gsl_matrix_set_col(X, 3, x3);

      gsl_vector_free(x3);

      gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nrow, p);
      gsl_matrix *cov;
      gsl_vector *coef;
      double chisq;
      cov = gsl_matrix_alloc(p, p);
      coef = gsl_vector_alloc(p);

      gsl_multifit_linear(X, Y, coef, cov, &chisq, work);
      gsl_multifit_linear_free(work);
      flag = 1;

      for (int k = 1; k < p; k++){
        if (coef(k) > MAX_EPS)
          flag = 0;
      }
      if (flag == 0)
        continue;
      // compute p-value for interaction first
      stderr = sqrt(COV(3, 3));
      t = coef(3) / stderr;
      pvalue[2] = 2*(t < 0 ? (1 - gsl_cdf_tdist_P(-t, df)) : (1 - gsl_cdf_tdist_P(t, df)));
      if (pvalue[2] > MIN_P_VALUE)
      {
        gsl_matrix_free(cov);
        gsl_vector_free(coef);
        continue;
      }
      stderr = sqrt(COV(2, 2));
      t = coef(2) / stderr;
      pvalue[1] = 2*(t < 0 ? (1 - gsl_cdf_tdist_P(-t, df)) : (1 - gsl_cdf_tdist_P(t, df)));
      stderr = sqrt(COV(1, 1));
      t = coef(1) / stderr;
      pvalue[0] = 2*(t < 0 ? (1 - gsl_cdf_tdist_P(-t, df)) : (1 - gsl_cdf_tdist_P(t, df)));
      gsl_matrix_free(cov);
      gsl_vector_free(coef);

      r1.push_back(G_colname[i]);
      r2.push_back(G_colname[j]);
      r1_p.push_back(pvalue[0]);
      r2_p.push_back(pvalue[1]);
      r1r2_p.push_back(pvalue[2]);

    }
  }
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("r1") = r1,
                                                   Rcpp::Named("r2") = r2,
                                                   Rcpp::Named("r1.p.value") = r1_p,
                                                   Rcpp::Named("r2.p.value") = r2_p,
                                                   Rcpp::Named("r1*r2.p.value") = r1r2_p);

  gsl_vector_free(x0);
  gsl_vector_free(x1);
  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_matrix_free(G);

  return output;
}
