#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
using namespace std;

void readPED(string FILE, vector<vector<double> > &mv, size_t &n0, size_t &nh, size_t &n1)
{
  ifstream input(FILE.c_str());
  string line;
  string tmp, tmp1, tmp2;
  stringstream ss;
  double val;
  //vector<vector<double> > mv;
  size_t ncol = 0, nrow = 0;
  n0 = 0;
  nh = 0;
  n1 = 0;
  while(getline(input, line))
  {
    nrow++;
    if (nrow != 1)
      mv[nrow-1].resize(ncol);
    if (mv.size() == nrow)
      mv.resize(mv.size()*2);
    ss.clear();
    ss.str(line);
    // discard the first six elements
    for (size_t i = 0; i < 6; i++)
      ss >> tmp;
    while (!ss.eof()) {
      ss >> tmp1 >> tmp2;
      // must judge eof status!!!
      if(ss.eof())
        break;
      tmp = tmp1 + tmp2;
      if (strcmp(tmp.c_str(), "11"))
      {
        val = 0;
        n0++;
      }
      else if (strcmp(tmp.c_str(), "12"))
      {
        val = 0.5;
        nh++;
      }
      else if (strcmp(tmp.c_str(), "22"))
      {
        val = 1;
        n1++;
      }

      //val = strtod(tmp.c_str(), NULL);
      if (nrow == 1)
      {
        if (mv[nrow-1].size() == ncol)
          mv[nrow-1].resize(ncol*2);
        ncol++;
      }
      mv[nrow-1].push_back(val);
    }
    if (nrow == 1)
      mv[nrow-1].resize(ncol);
  }
  mv.resize(nrow);
  cout << "ncol = " << ncol << endl
       << "nrow = " << nrow << endl;
}

void kinship(gsl_matrix* m, char* method, char* use, size_t &nh)
{
  if (strcmp("dominant",  method))
  {

  }
  else if (strcmp("recessive", method))
  {

  }
  else if (strcmp("additive", method))
  {
    //TODO
  }
}

void rowMeans(gsl_matrix *m, gsl_vector *rmu)
{
  gsl_vector *tmpRow = gsl_vector_calloc(m->size2);
  double mu;
  for (size_t i = 0; i < m->size1; i++)
  {
    gsl_matrix_get_row(tmpRow, m, i);
    mu = gsl_stats_mean(tmpRow->data, 1, tmpRow->size);
    if (mu > 0.5)
      gsl_vector_set(rmu, i, 1);
    else
      gsl_vector_set(rmu, i, 0);
  }
}


int main()
{
  string filename = "../data/X.PED";
  vector<vector<double> > mv(1, vector<double>(1));
  size_t n0, nh, n1;
  readPED(filename, mv, n0, nh, n1);
  gsl_matrix *m;
  size_t nrow = mv.size();
  size_t ncol = mv[0].size();
  m = gsl_matrix_calloc(nrow, ncol);
  for (size_t i = 0; i < nrow; i++)
  {
    for (size_t j = 0; j < ncol; j++)
      gsl_matrix_set(m, i, j, mv[i][j]);
  }
  cout << m->size1 <<"*" << m->size2 << endl;
  gsl_matrix_free(m);
  return 0;
}
