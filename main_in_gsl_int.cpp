#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_combination.h>
#include<gsl/gsl_statistics.h>

using namespace std;

#define MAX_ROW 305
#define MAX_COL 5000
#define MAX_PAIR_CORR 0.98

int main()
{
  /////////////////////////////////////
  //
  //  read data G
  //
  /////////////////////////////////////

  gsl_matrix_int *G = gsl_matrix_int_alloc(MAX_ROW, MAX_COL);
  vector<string> rowname(MAX_ROW);
  vector<string> colname(MAX_COL);
  
  ifstream input("Data.csv");
  string tmp;
  string comma;
  int val;
  input >> comma >> comma;
  // colnames
  for (int i = 0; i < MAX_COL; i++)
    {
      if (i == MAX_COL - 1)
	input >> tmp;
      else
	input >> tmp >> comma;
      rowname.push_back(tmp);
    }

  for (int i = 0; i < MAX_ROW; i++)
    {
      input >> tmp;
      rowname.push_back(tmp);
      for (int j = 0; j < MAX_COL; j++)
	{
	  input >> comma >> val;
	  gsl_matrix_int_set(G, i, j, val);
	}
    }
  input.close();
  /////////////////////////////////
  //
  //   read data Y.FA
  //
  /////////////////////////////////

  ifstream input2("Y.FA.txt");
  gsl_vector_int *Y = gsl_vector_int_alloc(MAX_ROW);
  for (int i = 0; i < MAX_ROW; i++)
    {
      input2 >> val;
      gsl_vector_int_set(Y, i, val);
    }
  input2.close();
  ///////////////////////////////////////
  //
  // combination
  //
  //////////////////////////////////////
  gsl_combination *c;
  vector<size_t> r1, r2;
  //  size_t k = 2;
  size_t *res;
  res = (size_t*)malloc(sizeof(size_t)*2);
  c = gsl_combination_calloc(MAX_COL, 2);
  do{
    res = gsl_combination_data(c);
    //    cout << res[0] << " "  << res[1] << endl;
    r1.push_back(res[0]);
    r2.push_back(res[1]);
  }while(gsl_combination_next(c) == GSL_SUCCESS);
  cout << r1.size() << " " << r2.size() << endl;

  ////////////////////////////////////////
  //
  //  main
  //
  ////////////////////////////////////////
  gsl_vector_int *x1 = gsl_vector_int_alloc(MAX_ROW);
  gsl_vector_int *x1_copy = gsl_vector_int_alloc(MAX_ROW);
  gsl_vector_int *x2 = gsl_vector_int_alloc(MAX_ROW);
  gsl_vector_int *x3 = gsl_vector_int_alloc(MAX_ROW);
  int max_val, min_val;
  double res_corr;
  for (size_t i = 0; i < r1.size(); i++)
    {
      gsl_matrix_int_get_col(x1, G, r1.at(i));
      gsl_matrix_int_get_col(x2, G, r2.at(i));
      gsl_vector_int_memcpy(x1_copy, x1);
      gsl_vector_int_mul(x1_copy, x2);
      gsl_vector_int_memcpy(x3, x1_copy);
      // if same
      max_val = gsl_vector_int_max(x1);
      min_val = gsl_vector_int_min(x1);
      if (max_val == min_val)
	continue;
      max_val = gsl_vector_int_max(x2);
      min_val = gsl_vector_int_min(x2);
      if (max_val == min_val)
	continue;
      res_corr = gsl_stats_correlation((gsl_vector) x1, 1, (gsl_vector) x2, 1, MAX_ROW);
      
      
      
    }

  
  cout << "Done!" << endl;
  gsl_vector_int_free(x1);
  gsl_vector_int_free(x2);
  gsl_vector_int_free(Y);
  gsl_matrix_int_free(G);
  return 0;
}
