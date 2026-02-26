/* Tests the performance of double double matrix multiplication
 * on regular cores of the GPU. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "ddmm_host.h"
#include "ddmm_kernels.h"

using namespace std;

double max_error
 ( int nrows, int ncols,
   double **Ahi, double **Alo, double **Bhi, double **Blo );
/*
 * Given two nrows-by-ncols double double matrices (Ahi, Alo)
 * and (Bhi, Blo), returns the maximum difference between them. */

void seconds_flops ( int nrows, int ncols, int dim, double seconds );
/*
 * Reports the number of flops for the matrix multiplication
 * if seconds > 0.0. */

void millisec_flops ( int nrows, int ncols, int dim, double millisec );
/*
 * Reports the number of flops for the matrix multiplication
 * if millisec > 0.0. */

int main ( int argc, char **argv )
{
   int seed = time(NULL);
   srand(seed);

   cout << "Give #rows of the product A*B : ";
   int m; cin >> m;
   cout << "Give #columns of the product A*B : ";
   int n; cin >> n;
   cout << "Give #columns of A, #rows of B : ";
   int k; cin >> k;

   double **Ahi = new double*[m];
   double **Alo = new double*[m];
   double **Bhi = new double*[k];
   double **Blo = new double*[k];
   double **Chi = new double*[m];
   double **Clo = new double*[m];

   random_dd_matrices(m, n, k, Ahi, Alo, Bhi, Blo, Chi, Clo, 1);

   cout << "computing double double matrix product plainly ..." << endl;
   clock_t start = clock();
   double_double_matmatmul(m, n, k, Ahi, Alo, Bhi, Blo, Chi, Clo);
   clock_t end = clock();
   double timelapsec = double(end - start)/CLOCKS_PER_SEC;
   cout << fixed << setprecision(3);
   cout << "elapsed CPU time " << timelapsec << " seconds" << endl;
   seconds_flops(m, n, k, timelapsec);
  
   double **Thi = new double*[n];
   double **Tlo = new double*[n];
   for(int i=0; i<n; i++)
   {
      Thi[i] = new double[k];
      Tlo[i] = new double[k];
   }
   transpose_dd_matrix(k, n, Bhi, Blo, Thi, Tlo);

   cout << "computing dd product using transposed B  ..." << endl;
   start = clock();
   double_double_transposed_mm(m, n, k, Ahi, Alo, Thi, Tlo, Chi, Clo);
   end = clock();
   timelapsec = double(end - start)/CLOCKS_PER_SEC;
   cout << fixed << setprecision(3);
   cout << "elapsed CPU time " << timelapsec << " seconds" << endl;
   seconds_flops(m, n, k, timelapsec);

   double **Dhi = new double*[m];
   double **Dlo = new double*[m];
   for(int i=0; i<m; i++)
   {
      Dhi[i] = new double[n];
      Dlo[i] = new double[n];
   } 
   float lapsedms;

   cout << "computing product on GPU ..." << endl;
   start = clock();
   GPU_dd_matmatmul(m, n, k, Ahi, Alo, Thi, Tlo, Dhi, Dlo, &lapsedms);
   end = clock();
   timelapsec = double(end - start)/CLOCKS_PER_SEC;
   cout << fixed << setprecision(3);
   cout << "elapsed CPU time " << timelapsec << " seconds" << endl;
   cout << "elapsed kernel time " << lapsedms << " milliseconds" << endl;
   // seconds_flops(m, n, k, timelapsec);
   millisec_flops(m, n, k, ((double) lapsedms));

   double errmax = max_error(m, n, Chi, Clo, Dhi, Dlo);

   cout << scientific << setprecision(3)
        << "error : " << errmax << endl;

   cout << "seed used : " << seed << endl;

   return 0;
}

double max_error
 ( int nrows, int ncols,
   double **Ahi, double **Alo, double **Bhi, double **Blo )
{
   double err = 0.0;

   cout << scientific << setprecision(16);

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         // cout << "Ahi[" << i << "][" << j << "] : " << Ahi[i][j] << endl;
         // cout << "Bhi[" << i << "][" << j << "] : " << Bhi[i][j] << endl;
         // cout << "Alo[" << i << "][" << j << "] : " << Alo[i][j] << endl;
         // cout << "Blo[" << i << "][" << j << "] : " << Blo[i][j] << endl;
         double d = abs(Ahi[i][j] - Bhi[i][j])
                  + abs(Alo[i][j] - Blo[i][j]);
         if(d > err) err = d;
      }

   return err;
}

void seconds_flops ( int nrows, int ncols, int dim, double seconds )
{
   if(seconds > 0.0)
   {
      long long int add,mul;
      flopcount_dd_matmatmul(nrows, ncols, dim, &add, &mul);
      const long long int totflop = add + mul;
      cout << "number of FP64 ops : " << totflop;

      double wallflops = ((double) totflop)/seconds;
      const int gigacnt = pow(2.0,30);
      cout << ", Flops : "
            << scientific << setprecision(3) << wallflops;
      cout << fixed << setprecision(3)
           << " = " << wallflops/gigacnt << " Gigaflops" << endl;
   }
}

void millisec_flops ( int nrows, int ncols, int dim, double millisec )
{
   if(millisec > 0.0)
   {
      long long int add,mul;
      flopcount_dd_matmatmul(nrows, ncols, dim, &add, &mul);
      const long long int totflop = add + mul;
      cout << "number of FP64 ops : " << totflop;

      double flops = 1000.0*((double) totflop)/millisec;
      const int gigacnt = pow(2.0,30);
      cout << ", Flops : "
            << scientific << setprecision(3) << flops;
      cout << fixed << setprecision(3)
           << " = " << flops/gigacnt << " Gigaflops" << endl;
   }
}
