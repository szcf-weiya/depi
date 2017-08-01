#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include <gsl/gsl_matrix.h>
using namespace std;

void readPED(string FILE, vector<vector<double> > &mv)
{
  ifstream input(FILE.c_str());
  string line;
  string tmp, tmp1, tmp2;
  stringstream ss;
  double val;
  //vector<vector<double> > mv;
  size_t ncol = 0, nrow = 0;
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
        val = 0;
      else if (strcmp(tmp.c_str(), "12"))
        val = 1;
      else if (strcmp(tmp.c_str(), "22"))
        val = 2;
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

int main()
{
  string filename = "../data/X.PED";
  vector<vector<double> > mv(1, vector<double>(1));
  readPED(filename, mv);
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
