#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

#define ERR 1e-10

void readPED(string FILE, vector<vector<double> > &mv, size_t &n0, size_t &nh, size_t &n1, size_t &nNA)
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
  nNA = 0;
  while(getline(input, line))
  {
    nrow++;
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
      if (!strcmp(tmp.c_str(), "11"))
      {
        val = 0;
        n0++;
      }
      else if (!strcmp(tmp.c_str(), "12"))
      {
        val = 0.5;
        nh++;
      }
      else if (!strcmp(tmp.c_str(), "22"))
      {
        val = 1;
        n1++;
      }
      else // NA values
      {
        val = 0.0/0.0;
        nNA++;
      }

      //val = strtod(tmp.c_str(), NULL);
      if (nrow == 1)
        ncol++;
      mv[nrow-1].push_back(val);
    }
    if (nrow == 1)
      mv[nrow-1].resize(ncol);
  }
  input.close();
  mv.resize(nrow);
  cout << "ncol = " << ncol << endl
       << "nrow = " << nrow << endl;
}


void subNanValue(gsl_matrix *m, gsl_vector *mu)
{
  for (size_t j = 0; j < m->size2; j++)
    for (size_t i = 0; i < m->size1; i++)
      if (gsl_isnan(gsl_matrix_get(m, i, j)))
        gsl_matrix_set(m, i, j, gsl_vector_get(mu, i));
}

void subValue(gsl_matrix *m, gsl_vector *mu)
{
  double tmp;
  for (size_t j = 0; j < m->size2; j++)
    for (size_t i = 0; i < m->size1; i++)
    {
      tmp = gsl_matrix_get(m, i, j);
      if (gsl_isnan(tmp))
        continue;
      else
      {
        if (!gsl_fcmp(tmp, 0.5, ERR))
          gsl_matrix_set(m, i, j, gsl_vector_get(mu, i));
      }
    }
}

double noNanMean(const gsl_vector *v)
{
  double res = 0, tmp;
  size_t num = 0;
  for (size_t i = 0; i < v->size; i++)
  {
    tmp = gsl_vector_get(v, i);
    if (gsl_isnan(tmp))
      continue;
    num++;
    res += tmp;
  }
  return res/num*1.0;
}

void rowMeans(gsl_matrix *m, gsl_vector *rmu)
{
  gsl_vector *tmpRow = gsl_vector_calloc(m->size2);
  double mu;
  for (size_t i = 0; i < m->size1; i++)
  {
    gsl_matrix_get_row(tmpRow, m, i);
    //mu = gsl_stats_mean(tmpRow->data, 1, tmpRow->size);
    mu = noNanMean(tmpRow);
    if (mu > 0.5)
      gsl_vector_set(rmu, i, 1);
    else
      gsl_vector_set(rmu, i, 0);
  }
  gsl_vector_free(tmpRow);
}

bool isNanInVector(gsl_vector *v)
{
  for (size_t i = 0; i < v->size; i++)
  {
    if (gsl_isnan(gsl_vector_get(v, i)))
      return true;
  }
  return false;
}

void rbind(const gsl_matrix *m1, const gsl_matrix *m2, gsl_matrix *m)
{
  gsl_vector *tmp = gsl_vector_alloc(m->size2);
  for (size_t i = 0; i < m1->size1; i++)
  {
    gsl_matrix_get_row(tmp, m1, i);
    gsl_matrix_set_row(m, i, tmp);
  }
  for (size_t i = 0; i < m2->size1; i++)
  {
    gsl_matrix_get_row(tmp, m2, i);
    gsl_matrix_set_row(m, i+m1->size1, tmp);
  }
  gsl_vector_free(tmp);
}

void gini(gsl_vector *xi, const gsl_vector *xj)
{
  gsl_vector *tmp = gsl_vector_calloc(xi->size);
  gsl_vector_memcpy(tmp, xi);
  gsl_vector_mul(tmp, xj);
  gsl_vector_scale(tmp, 2);
  gsl_vector_add_constant(tmp, 1);
  gsl_vector_sub(tmp, xi);
  gsl_vector_sub(tmp, xj);
  gsl_vector_memcpy(xi, tmp);
  gsl_vector_free(tmp);
}

void kinship(gsl_matrix* m, const char* method, const char* use, size_t &nh, gsl_matrix* K)
{
  gsl_vector *snp = gsl_vector_calloc(m->size1);
  gsl_vector *snp2 = gsl_vector_calloc(m->size1);
  rowMeans(m, snp);
  gsl_vector_sub(snp2, snp); // snp2 = 1-snp
  gsl_vector_add_constant(snp2, 1);

  gsl_matrix_set(K, 0, 0, 1);

  if (!strcmp("additive", method) && (nh > 0))
  {
    gsl_matrix *m2 = gsl_matrix_calloc(m->size1, m->size2); // DO NOT forget
    gsl_matrix_memcpy(m2, m);
    subValue(m, snp);
    subValue(m2, snp2);
    gsl_matrix *G = gsl_matrix_calloc(m->size1+m2->size1, m->size2);
    gsl_vector *rmu = gsl_vector_calloc(m->size1+m2->size1);
    rbind(m, m2, G);
    if (!strcmp(use, "all"))
    {
      rowMeans(G, rmu);
  //    subNanValue(G, rmu);
    }
    else if (!strcmp(use, "complete.obs"))
    {
      ; // no need to remove nan, in the following step, the nan value will be removed
    }
    gsl_vector *xi = gsl_vector_calloc(G->size1);
    gsl_vector *xj = gsl_vector_calloc(G->size1);
    for (size_t i = 1; i < G->size2; i++)
    {
      for (size_t j = 0; j < i; j++)
      {
        gsl_matrix_get_col(xi, G, i);
        gsl_matrix_get_col(xj, G, j);
        gini(xi, xj);
        gsl_matrix_set(K, i, j, noNanMean(xi));
        gsl_matrix_set(K, j, i, noNanMean(xi));
      }
      gsl_matrix_set(K, i, i, 1);
    }
    gsl_vector_free(xi);
    gsl_vector_free(xj);
    gsl_vector_free(rmu);
    gsl_matrix_free(m2);
    gsl_matrix_free(G);
  }
  else
  {
    if (!strcmp("dominant",  method))
      subValue(m, snp);
    else if (!strcmp("recessive", method))
      subValue(m, snp2);
    gsl_vector *rmu = gsl_vector_calloc(m->size1);
    gsl_vector *xi = gsl_vector_calloc(m->size1);
    gsl_vector *xj = gsl_vector_calloc(m->size1);
    if (!strcmp(use, "all"))
    {
      rowMeans(m, rmu);
      subNanValue(m, rmu);
    }
    else if (!strcmp(use, "complete.obs"))
    {
      ; // no need to remove nan, in the following step, the nan value will be removed
    }
    for (size_t i = 1; i < m->size2; i++)
    {
      for (size_t j = 0; j < i-1; j++)
      {
        gsl_matrix_get_col(xi, m, i);
        gsl_matrix_get_col(xj, m, j);
        gini(xi, xj);
        gsl_matrix_set(K, i, j, noNanMean(xi));
        gsl_matrix_set(K, j, i, noNanMean(xi));
      }
      gsl_matrix_set(K, i, i, 1);
    }
    gsl_vector_free(xi);
    gsl_vector_free(xj);
    gsl_vector_free(rmu);
  }
  gsl_vector_free(snp);
  gsl_vector_free(snp2);
}

void eigen(const gsl_matrix *K, gsl_vector *eval, gsl_matrix *evec)
{
  // matrix K is symmetric
  gsl_matrix *KK = gsl_matrix_alloc(K->size1, K->size2);
  gsl_matrix_memcpy(KK, K);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(K->size1);
  gsl_eigen_symmv(KK, eval, evec, w);
  gsl_eigen_symmv_free(w);
  gsl_matrix_free(KK);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
}

int main()
{
  string filename = "../data/X.PED";
  vector<vector<double> > mv(1, vector<double>(1));
  size_t n0, nh, n1, nNA;
  readPED(filename, mv, n0, nh, n1, nNA);
  if (n0+nh+n1+nNA != mv.size()*mv[0].size())
  {
    cout << "stop!" << endl;
    return 0;
  }
  gsl_matrix *m;
  size_t nrow = mv.size();
  size_t ncol = mv[0].size();
  m = gsl_matrix_calloc(nrow, ncol);
  for (size_t i = 0; i < nrow; i++)
  {
    for (size_t j = 0; j < ncol; j++)
    {
      gsl_matrix_set(m, i, j, mv[i][j]);
    }
  }
  cout << m->size1 <<"*" << m->size2 << endl;
  cout << n0 << endl
       << nh << endl
       << n1 << endl
       << nNA << endl;

  gsl_matrix *K = gsl_matrix_calloc(m->size2, m->size2);

  string method = "additive";
  string use = "all";
  kinship(m, method.c_str(), use.c_str(), nh, K);
  gsl_matrix *evec = gsl_matrix_alloc(K->size1, K->size2);
  gsl_vector *eval = gsl_vector_alloc(K->size1);
  eigen(K, eval, evec);
  for (size_t i = 0; i < K->size1; i++)
    cout << gsl_vector_get(eval, i) << endl;
  /*
  ofstream output("../data/res_K.dat");
  for (size_t i = 0; i < K->size1; i++)
    for (size_t j = 0; j < K->size2; j++)
    {
      if (j!=K->size2-1)
        output << gsl_matrix_get(K, i, j) << ",";
      else
        output << gsl_matrix_get(K, i, j) << endl;
    }
  output.close();
  */
  //FILE *output = fopen("../data/res_K.dat", "w");
  //gsl_matrix_fprintf(output, K, "%f");
  //gsl_matrix_fwrite(output, K);
  //fclose(output);

  gsl_matrix_free(m);
  gsl_matrix_free(K);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  return 0;
}
