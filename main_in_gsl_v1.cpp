#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_combination.h>

using namespace std;

#define MAX_ROW 305
#define MAX_COL 5000

int main()
{
  /////////////////////////////////////
  //
  //  read data G
  //
  /////////////////////////////////////
  /* 
  double** G;
  G = (double**)malloc(sizeof(double)*MAX_ROW);
  for (int i = 0; i < MAX_ROW; i++)
    G[i] = (double*)malloc(sizeof(double)*MAX_COL);
  */

  // gsl implementation
  gsl_matrix *G = gsl_matrix_alloc(MAX_ROW, MAX_COL);
  vector<string> rowname(MAX_ROW);
  vector<string> colname(MAX_COL);
  
  ifstream input("Data.csv");
  string tmp;
  string comma;
  double val;
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
	  gsl_matrix_set(G, i, j, val);
	}
    }
  input.close();
  /////////////////////////////////
  //
  //   read data Y.FA
  //
  /////////////////////////////////

  ifstream input2("Y.FA.txt");
  double Y[MAX_ROW];
  for (int i = 0; i < MAX_ROW; i++)
    input2 >> Y[i];
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

  


  
  cout << "Done!" << endl;
  gsl_matrix_free(G);
  return 0;
}
