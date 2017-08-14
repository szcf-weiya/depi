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
#define MIN_P_VALUE 0.001


#define coef(i) (gsl_vector_get(coef, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))


// https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

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
  //getline(input, line);
  safeGetline(input, line);
  stringstream ss(line);
  while(!ss.eof())
  {
    ss >> tmp;
    colname.push_back(tmp);
    (*ncol)++;
  }
  cout << colname.size() << endl;

  while (safeGetline(input, line)) {
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
      mv[*nrow].push_back(val);
    }
    (*nrow)++;
  }
  input.close();
  mv.resize(*nrow);
}

void readData(char* FILE, gsl_matrix *m, int row, int col, vector<string> &rowname, vector<string> &colname)
{
  ifstream input(FILE);
  string tmp, tmpp;
  //  vector<string> rowname;
  //vector<string> colname;

  double val;
  // colnames
  string line;
  // the first line
  safeGetline(input, line);
  stringstream ss(line);
  // ss >> tmp; // null no need!!
  //  cout << tmp << endl;
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

int main()
{
  start = clock();
  /////////////////////////////////////
  //
  //  read data G
  //
  /////////////////////////////////////
  FILE *output;
  output = fopen("res_gsl.txt", "w");
  char FILE[] = "../data/Data.txt";
  gsl_matrix *G = gsl_matrix_alloc(MAX_ROW, MAX_COL);
  vector<string> G_colname;
  vector<string> G_rowname;
  size_t nrow, ncol;
  vector<vector<double> > mv(1, vector<double>(1));
  readData(FILE, mv, G_rowname, G_colname, &nrow, &ncol);
  //readData(FILE, G, MAX_ROW, MAX_COL, G_rowname, G_colname);
  cout << G_colname[0] << endl
       << G_colname[1] << endl;

/*
  //  ifstream input3("cov.txt");
  char FILE3[] = "cov.txt";
  gsl_matrix *coveriate = gsl_matrix_alloc(MAX_ROW, COV_COL);
  vector<string> colname;
  vector<string> rowname;
  readData(FILE3, coveriate, MAX_ROW, COV_COL, rowname, colname);
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
  ////////////////////////////////////////
  //
  //  main
  //
  ////////////////////////////////////////
  gsl_vector *x0 = gsl_vector_alloc(MAX_ROW);
  gsl_vector_set_all(x0, 1);
  //  gsl_vector *x1_copy = gsl_vector_alloc(MAX_ROW);

  int p;
  int max_val, min_val;

  double min_pair_geno;
  gsl_vector *coveriate_vec;
  coveriate_vec = gsl_vector_alloc(MAX_ROW);
  if (!IS_EXIST_COV)
      p = 4;
  else
      p = 4 + COV_COL;
  gsl_matrix *X;
  X = gsl_matrix_alloc(MAX_ROW, p);
  if (IS_EXIST_COV)
    {
        for (int ip = 0; ip < COV_COL; ip++)
	  {
	    gsl_matrix_get_col(coveriate_vec, coveriate, ip);
	    gsl_matrix_set_col(X, 4+ip, coveriate_vec);
	  }
    }
  gsl_vector_free(coveriate_vec);
  gsl_matrix_free(coveriate);

  // cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  // cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  fprintf(output, "r1\tr2\tr1.p.value\tr2.p.value\tr1*r2.p.value\n");
  //# pragma omp parallel
  //# pragma omp for schedule(dynamic) // relatively slow
  gsl_matrix_set_col(X, 0, x0);
  gsl_vector *x1 = gsl_vector_alloc(MAX_ROW);
  for (int i = 0; i < MAX_COL/100; i++)
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
      for (int j = i+1; j < MAX_COL; j++)
	{
	  gsl_vector *x2 = gsl_vector_alloc(MAX_ROW);

	  size_t flag, df;

	  double stderr, t;
	  double pvalue[3];
	  df = MAX_ROW - p;

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

	  double res_corr = (float)gsl_stats_correlation(x1->data, 1, x2->data, 1, MAX_ROW);
	  if (fabs(res_corr) >= MAX_PAIR_CORR)
	    continue;
	  gsl_vector *x3 = gsl_vector_alloc(MAX_ROW);
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

	  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(MAX_ROW, p);
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

	  fprintf(output, "%d\t%d\t%.6f\t%.6f\t%.10f\n", i, j, pvalue[0], pvalue[1], pvalue[2]);
	}
    }
  //for (int i = 0; i < MAX_ROW; i++)
  // cout << x1055->data[i] << ", " << x1056->data[i] << ", " << x1057->data[i] << endl;

  //  stop = clock();
  //  cout << "Done!" << (float)(stop-start)/CLOCKS_PER_SEC << endl;
  gsl_vector_free(x0);
  gsl_vector_free(x1);
  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_matrix_free(G);
  fclose(output);
  //  output.close();
  */
  return 0;

}
