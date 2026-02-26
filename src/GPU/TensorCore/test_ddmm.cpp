/* Tests the performance of double double matrix multiplication
 * on regular cores of the GPU. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "ddmm_host.h"

using namespace std;

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

   if(timelapsec != 0.0)
   {
      long long int add,mul;
      flopcount_dd_matmatmul(m, n, k, &add, &mul);
      const long long int totflop = add + mul;
      cout << "number of FP64 ops : " << totflop << endl;

      double wallflops = ((double) totflop)/timelapsec;
      const int gigacnt = pow(2.0,30);
      cout << "Wall Time Flops : "
            << scientific << setprecision(3) << wallflops;
      cout << fixed << setprecision(3)
           << " = " << wallflops/gigacnt << " Gigaflops" << endl;
   }
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

   if(timelapsec != 0.0)
   {
      long long int add,mul;
      flopcount_dd_matmatmul(m, n, k, &add, &mul);
      const long long int totflop = add + mul;
      cout << "number of FP64 ops : " << totflop << endl;

      double wallflops = ((double) totflop)/timelapsec;
      const int gigacnt = pow(2.0,30);
      cout << "Wall Time Flops : "
            << scientific << setprecision(3) << wallflops;
      cout << fixed << setprecision(3)
           << " = " << wallflops/gigacnt << " Gigaflops" << endl;
   }

   cout << "seed used : " << seed << endl;

   return 0;

}
