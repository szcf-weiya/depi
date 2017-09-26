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
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_blas.h>
#include<omp.h>

#include "safeGetline.h"

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
  safeGetline(input, line);
  stringstream ss(line);
  while(!ss.eof())
  {
    ss >> tmp;
    colname.push_back(tmp);
    (*ncol)++;
  }
  //cout << colname.size() << endl;

  while (!safeGetline(input, line).eof()) {
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
  safeGetline(input, line);
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
      safeGetline(input, line);
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

Rcpp::List DetectEpi(SEXP inputfile_X, SEXP inputfile_Y, SEXP inputfile_COV)
{
  // convert file name
  Rcpp::CharacterVector r_inputfile_X(inputfile_X);
  Rcpp::CharacterVector r_inputfile_Y(inputfile_Y);

  string file_X = Rcpp::as<string>(r_inputfile_X);
  string file_Y = Rcpp::as<string>(r_inputfile_Y);

  bool IS_EXIST_COV = FALSE;
  string file_COV;
  // check cov file is null or not
  if(!Rf_isNull(inputfile_COV))
  {
    Rcpp::CharacterVector r_inputfile_COV(inputfile_COV);
    file_COV = Rcpp::as<string>(r_inputfile_COV);
    IS_EXIST_COV = TRUE;
  }
  // read data X
  //gsl_matrix *G = gsl_matrix_alloc(MAX_ROW, MAX_COL);
  vector<string> G_colname;
  vector<string> G_rowname;
  //readData(file_X.c_str(), G, MAX_ROW, MAX_COL, G_rowname, G_colname);
  vector<vector<double> > mv(1, vector<double>(1));
  size_t nrow, ncol;
  readData(file_X.c_str(), mv, G_rowname, G_colname, &nrow, &ncol);
  gsl_matrix *G = gsl_matrix_calloc(nrow, ncol);

  for (size_t i = 0; i < nrow; i++)
    for (size_t j = 0; j < ncol; j++)
      gsl_matrix_set(G, i, j, mv[i][j]);

  int p;
  if (IS_EXIST_COV)
    p = 4 + COV_COL;
  else
    p = 4;
  
  gsl_matrix *b = gsl_matrix_alloc(nrow, COV_COL+2);
  gsl_matrix *B = gsl_matrix_alloc(COV_COL+2, COV_COL+2);
  gsl_matrix *invB = gsl_matrix_alloc(COV_COL+2, COV_COL+2);
  
  gsl_matrix *X = gsl_matrix_alloc(nrow, p);
  if (IS_EXIST_COV)
  {
    gsl_matrix *covariate = gsl_matrix_alloc(nrow, COV_COL);
    // read data cov
    vector<string> cov_colname;
    vector<string> cov_rowname;
    readData(file_COV.c_str(), covariate, nrow, COV_COL, cov_rowname, cov_colname);
    gsl_vector *covariate_vec = gsl_vector_alloc(nrow);
    for (int ip = 0; ip < COV_COL; ip++)
    {
      gsl_matrix_get_col(covariate_vec, covariate, ip);
      gsl_matrix_set_col(b, 2+ip, covariate_vec);
    }
    gsl_vector_free(covariate_vec);
    gsl_matrix_free(covariate);
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

  // cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  // cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";
  //RcppGSL::vector<double> r1, r2, r1_p, r2_p, r1r2_p;
  vector<string> r1, r2;
  vector<double> r1_p, r2_p, r1r2_p;

  gsl_permutation *permutation_B = gsl_permutation_alloc(B->size1);
  int status;
  
//  fprintf(output, "r1\tr2\tr1.p.value\tr2.p.value\tr1*r2.p.value\n");
  //# pragma omp parallel
  //# pragma omp for schedule(dynamic) // relatively slow
  gsl_matrix_set_col(X, 0, x0);
  gsl_vector *x1 = gsl_vector_alloc(nrow);
  for (int i = 0; i < ncol; i++)
  {
    gsl_matrix_get_col(x1, G, i);
    max_val = gsl_vector_max(x1);
    min_val = gsl_vector_min(x1);

    if (max_val == min_val)
      continue;
    
    // construct x1
    gsl_matrix_set_col(X, 1, x1);
    
    // B = b'b
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, b, b, 0.0, B);
    
    // inv(B) by LU decomposition
    gsl_linalg_LU_decomp(B, permutation_B, &status);
    gsl_linalg_LU_invert(B, permutation_B, invB);
    
    # pragma omp parallel for schedule(dynamic) // faster!!!
    for (int j = i+1; j < ncol/100; j++)
    {
      gsl_vector *x2 = gsl_vector_alloc(nrow);

      size_t flag, df;

      double stderr, t;
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
      double A_1i_11, A_1i_12, A_1i_22;
      double a11, a12, a21, a22, a_det;
      gsl_matrix *invA_1i = gsl_matrix_alloc(2, 2);
      gsl_matrix *a_1i = gsl_matrix_alloc(nrow, 2);
      gsl_matrix *V_1i = gsl_matrix_alloc(COV_COL+2, 2); 
      gsl_matrix *invB_mul_V_1i = gsl_matrix_alloc(COV_COL+2, 2);
      gsl_matrix *m_tmp = gsl_matrix_alloc(2, 2);
      gsl_matrix *B_1 = gsl_matrix_alloc(2, 2);
      gsl_matrix *invD = gsl_matrix_alloc(2, 2);
      gsl_matrix *m_tmp2 = gsl_matrix_alloc(COV_COL+2, 2);
      gsl_matrix *m_tmp3 = gsl_matrix_alloc(2, 2);
      gsl_matrix *m_tmp4 = gsl_matrix_alloc(COV_COL+2, 2);
      gsl_matrix *invXX_11 = gsl_matrix_alloc(2, 2);
      gsl_matrix *invXX_22 = gsl_matrix_alloc(COV_COL+2, COV_COL+2);
      gsl_matrix *invXX_21 = gsl_matrix_alloc(COV_COL+2, 2);
      gsl_vector *XY_1 = gsl_vector_alloc(2);
      gsl_vector *XY_2 = gsl_vector_alloc(COV_COL+2);
      gsl_vector *beta_1 = gsl_vector_alloc(2);
      gsl_vector *beta_2 = gsl_vector_alloc(COV_COL+2);
      gsl_vector *Yhat = gsl_vector_alloc(nrow);
      double rss, zscore, pvalue;
      
      // A_1i
      gsl_blas_ddot(x2, x2, &A_1i_11);
      gsl_blas_ddot(x2, x3, &A_1i_12);
      gsl_blas_ddot(x3, x3, &A_1i_22);
      // invA_1i
      a_det = A_1i_11*A_1i_22-A_1i_12*A_1i_12;
      
      a11 = A_1i_22/a_det;
      a12 = -A_1i_12/a_det;
      a22 = A_1i_11/a_det;
      gsl_matrix_set(invA_1i, 0, 0, a11);
      gsl_matrix_set(invA_1i, 1, 1, a22);
      gsl_matrix_set(invA_1i, 0, 1, a12);
      gsl_matrix_set(invA_1i, 1, 0, a12);
      
      // construct a_1i
      gsl_matrix_set_col(a_1i, 0, x2);
      gsl_matrix_set_col(a_1i, 1, x3); 
      // V_1i = b'a_1i
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, b, a_1i, 0.0, V_1i);
      gsl_vector_free(x2);
      gsl_vector_free(x3);
      // invB_mul_V_1i
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invB, V_1i, 0.0, invB_mul_V_1i);
      // B_1 = V_1i' mul invB mul V_1i
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, invB_mul_V_1i, V_1i, 0.0, B_1);
      
      // D = I - B_1 * invA_1i
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, B_1, invA_1i, 0.0, m_tmp);
      
      a11 = gsl_matrix_get(m_tmp, 0, 0);
      a22 = gsl_matrix_get(m_tmp, 1, 1);
      a12 = gsl_matrix_get(m_tmp, 0, 1);
      
      // D is noy symmetric !!
      a21 = gsl_matrix_get(m_tmp, 1, 0);
      a11 += 1;
      a22 += 1;
      a_det = a11*a22-a12*a21;
      
      gsl_matrix_set(invD, 0, 0, a22/a_det);
      gsl_matrix_set(invD, 1, 1, a11/a_det);
      gsl_matrix_set(invD, 1, 0, -a21/a_det);
      gsl_matrix_set(invD, 0, 1, -a12/a_det);
      
      // invXX_11
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invA_1i, B_1, 0.0, m_tmp3);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, m_tmp3, invD, 0.0, m_tmp); // Do not replace m_tmp with m_tmp3 !! 
      gsl_matrix_memcpy(invXX_11, invA_1i);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m_tmp, invA_1i, 1.0, invXX_11);
      
      // invXX_22
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invB_mul_V_1i, invA_1i, 0.0, m_tmp2);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, m_tmp2, invD, 0.0, m_tmp4);
      gsl_matrix_memcpy(invXX_22, invB);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, m_tmp4, invB_mul_V_1i, 1.0, invXX_22);
      
      // invXX_21
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invB, V_1i, 0.0, m_tmp2);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, m_tmp2, invD, 0.0, m_tmp4);
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, m_tmp4, invA_1i, 0.0, invXX_21);
      
      // X'Y
      gsl_blas_dgemv(CblasTrans, 1.0, a_1i, Y, 0.0, XY_1);
      gsl_blas_dgemv(CblasTrans, 1.0, b, Y, 0.0, XY_2);
      
      // beta
      gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_11, XY_1, 0.0, beta_1);
      gsl_blas_dgemv(CblasTrans, 1.0, invXX_21, XY_2, 1.0, beta_1);
      gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_21, XY_1, 0.0, beta_2);
      gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_22, XY_2, 1.0, beta_2);
      
      // RSS
      gsl_blas_dgemv(CblasNoTrans, 1.0, a_1i, beta_1, 0.0, Yhat);
      gsl_blas_dgemv(CblasNoTrans, 1.0, b, beta_2, 1.0, Yhat);
      
      gsl_vector_sub(Yhat, Y);
      
      gsl_blas_ddot(Yhat, Yhat, &rss);
      
      // zscore
      zscore = gsl_vector_get(beta_1, 1)/(sqrt(rss/df*gsl_matrix_get(invXX_11, 1, 1)));
      pvalue = 2*(zscore < 0 ? (1 - gsl_cdf_tdist_P(-zscore, df)) : (1 - gsl_cdf_tdist_P(zscore, df)));
      
      //if (pvalue > MIN_P_VALUE)
      //		continue;
      gsl_vector_free(XY_1);
      gsl_vector_free(XY_2);
      gsl_vector_free(beta_1);
      gsl_vector_free(beta_2);
      gsl_vector_free(Yhat);
      gsl_matrix_free(invXX_11);
      gsl_matrix_free(invXX_22);
      gsl_matrix_free(invXX_21);
      gsl_matrix_free(invA_1i);
      gsl_matrix_free(a_1i);
      gsl_matrix_free(V_1i);
      gsl_matrix_free(invB_mul_V_1i);
      gsl_matrix_free(m_tmp);
      gsl_matrix_free(B_1);
      gsl_matrix_free(invD);
      gsl_matrix_free(m_tmp2);
      gsl_matrix_free(m_tmp3);
      gsl_matrix_free(m_tmp4);
      
      r1.push_back(G_colname[i]);
      r2.push_back(G_colname[j]);
      r1r2_p.push_back(pvalue);

    }
  }
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("r1") = r1,
                                                   Rcpp::Named("r2") = r2,
      //                                             Rcpp::Named("r1.p.value") = r1_p,
        //                                           Rcpp::Named("r2.p.value") = r2_p,
                                                   Rcpp::Named("r1*r2.p.value") = r1r2_p);

  gsl_vector_free(x0);
  gsl_vector_free(x1);
  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_matrix_free(G);

  return output;
}
