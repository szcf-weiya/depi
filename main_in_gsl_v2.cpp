#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
clock_t start, stop;
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_combination.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_cdf.h>
#include<omp.h>

using namespace std;

#define MAX_ROW 305
#define MAX_COL 5000
#define COV_COL 4
#define MAX_PAIR_CORR 0.98
#define MAX_EPS 1e+12
#define IS_EXIST_COV 1

#define coef(i) (gsl_vector_get(coef, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))

void readData(char* FILE, gsl_matrix *m, int row, int col)
{
  ifstream input(FILE);
  string tmp;
  vector<string> rowname(row);
  vector<string> colname(col);
  
  double val;
  // colnames
  string line;
  // the first line
  getline(input, line);
  stringstream ss(line);
  ss >> tmp; // null
  for (size_t i = 0; i < row; i++)
    {
      ss >> tmp;
      colname.push_back(tmp);
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

int main()
{
  start = clock();
  /////////////////////////////////////
  //
  //  read data G
  //
  /////////////////////////////////////
  //ofstream output;
  //output.open("res_gsl.txt");
  FILE *output;
  output = fopen("res_gsl.txt", "w");
  char FILE[] = "Data.txt";
  gsl_matrix *G = gsl_matrix_alloc(MAX_ROW, MAX_COL);
  readData(FILE, G, MAX_ROW, MAX_COL);
  //  ifstream input3("cov.txt");
  char FILE3[] = "cov.txt";
  gsl_matrix *coveriate = gsl_matrix_alloc(MAX_ROW, COV_COL);
  readData(FILE3, coveriate, MAX_ROW, COV_COL);
  /////////////////////////////////
  //
  //   read data Y.FA
  //
  /////////////////////////////////
  
  ifstream input2("Y.FA.txt");
  double val;
  gsl_vector *Y = gsl_vector_alloc(MAX_ROW);
  for (size_t i = 0; i < MAX_ROW; i++)
    {
      input2 >> val;
      gsl_vector_set(Y, i, val);
    }
  input2.close();
  ///////////////////////////////////////
  //
  // combination
  //
  //////////////////////////////////////
  /*
  gsl_combination *c;
  vector<size_t> r1, r2;
  //  size_t k = 2;
  size_t *res;
  res = (size_t*)malloc(sizeof(size_t)*2);
  c = gsl_combination_calloc(MAX_COL, 2);
  do{
    res = gsl_combination_data(c);
    r1.push_back(res[0]);
    r2.push_back(res[1]);
  }while(gsl_combination_next(c) == GSL_SUCCESS);
  cout << r1.size() << " " << r2.size() << endl;
  */
#ifndef _OPENMP
  fprintf(stderr, "OpenMP not supported");
#endif
  ////////////////////////////////////////
  //
  //  main
  //
  ////////////////////////////////////////
  gsl_vector *x0 = gsl_vector_alloc(MAX_ROW);
  gsl_vector_set_all(x0, 1);
  gsl_vector *x1 = gsl_vector_alloc(MAX_ROW);
  gsl_vector *x1_copy = gsl_vector_alloc(MAX_ROW);
  gsl_vector *x2 = gsl_vector_alloc(MAX_ROW);
  gsl_vector *x3 = gsl_vector_alloc(MAX_ROW);
  gsl_matrix *X;
  if (!IS_EXIST_COV)
    X = gsl_matrix_alloc(MAX_ROW, 4);
  else
    X = gsl_matrix_alloc(MAX_ROW, 4 + COV_COL); 
  //gsl_matrix *X = gsl_matrix_alloc(MAX_ROW, 4); // TO DO WITH COV
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";
  //output << fixed << setprecision(6); // keep two the decimal places
  //output << "r1\tr2\tr1.p.value\tr2.p.value\tr1*r2.p.value"<<endl;
  fprintf(output, "r1\tr2\tr1.p.value\tr2.p.value\tr1*r2.p.value\n");
  
  # pragma omp parallel for
  
  for (size_t i = 0; i < 10; i++)
    {
      gsl_matrix_get_col(x1, G, i);
      //cout << "i" << i << endl;
      # pragma omp parallel for
      for (size_t j = i+1; j < MAX_COL; j++)
	{
	  //      gsl_matrix_get_col(x1, G, r1.at(i));
	    int max_val, min_val;
	    double res_corr;
	    double min_pair_geno;
	    gsl_matrix *cov;
	    gsl_vector *coef, *coveriate_vec;
	    coveriate_vec = gsl_vector_alloc(MAX_ROW);
	    double chisq;
	    size_t flag, df;
	    int p;
	    double stderr, t;

      gsl_matrix_get_col(x2, G, j);
      gsl_vector_memcpy(x3, x1);
      gsl_vector_mul(x3, x2);
      // if same
      max_val = gsl_vector_max(x1);
      min_val = gsl_vector_min(x1);
      if (max_val == min_val)
	continue;
      max_val = gsl_vector_max(x2);
      min_val = gsl_vector_min(x2);
      if (max_val == min_val)
	continue;
      max_val = gsl_vector_max(x3);
      if (max_val == 0)
	continue;
      res_corr = (float)gsl_stats_correlation(x1->data, 1, x2->data, 1, MAX_ROW);
      //      cout << "cor " << res_corr << endl;
      if (abs(res_corr) < MAX_PAIR_CORR)
	{
	  gsl_matrix_set_col(X, 0, x0);
	  gsl_matrix_set_col(X, 1, x1);
	  gsl_matrix_set_col(X, 2, x2);
	  gsl_matrix_set_col(X, 3, x3);
	  if (!IS_EXIST_COV)
	    p = 4;
	  else
	    {
	      p = 4 + COV_COL;
	      
	      for (size_t ip = 0; ip < COV_COL; ip++)
		{
		  gsl_matrix_get_col(coveriate_vec, coveriate, ip);
		  gsl_matrix_set_col(X, 4+ip, coveriate_vec);
		}
		
	    }
	  df = MAX_ROW - p;
	  cov = gsl_matrix_alloc(p, p);
	  coef = gsl_vector_alloc(p);
	  double pvalue[3]; 
	  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(MAX_ROW, p);
	  gsl_multifit_linear(X, Y, coef, cov, &chisq, work);
	  flag = 1;
	  
	  # pragma omp parallel \
	    shared (coef, flag)	
	  # pragma omp for
	  for (size_t k = 1; k < p; k++)
	    {
	      if (coef(k) > MAX_EPS)
		flag = 0;
	    }
	  if (flag == 0)
	    continue;
	  for (size_t ii = 1; ii < 4; ii++)
	    {
	      stderr = sqrt(COV(ii, ii));
	      t = coef(ii) / stderr;
	      pvalue[ii-1] = 2*(t < 0 ? (1 - gsl_cdf_tdist_P(-t, df)) : (1 - gsl_cdf_tdist_P(t, df)));
	    }
	  // fprintf(output, "%lu\t%lu\t%.6f\t%.6f\t%.6f\n", i, j, pvalue[0], pvalue[1], pvalue[2]);
	  //output << r1.at(i) << "\t" << r2.at(i) << "\t" << pvalue[0] << "\t" << pvalue[1] << "\t" << pvalue[2] <<endl;
	  /*
	  output << "#######################" <<endl
		 << "r1 = " << r1.at(i) << "; r2 = " << r2.at(i) <<  endl
		 << "# best fit: Y = " << coef(0) << " + " << coef(1) << "x1 + " << coef(2) << "x2 + " << coef(3) << "x3" << endl
		 << "# chisq = " << chisq << endl;
	  */
	  gsl_multifit_linear_free(work);
	}
	}
    }
  
  stop = clock();
  
  cout << "Done!" << (float)(stop-start)/CLOCKS_PER_SEC << endl;
  gsl_vector_free(x1);
  gsl_vector_free(x2);
  gsl_vector_free(Y);
  gsl_matrix_free(G);
  fclose(output);
  // output.close();
  return 0;
}

