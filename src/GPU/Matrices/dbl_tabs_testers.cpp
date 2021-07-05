// The file dbl_tabs_testers.cpp defines the function with prototypes
// in dbl_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random_matrices.h"
#include "dbl_factorizations.h"
#include "dbl_tabs_host.h"
#include "dbl_tabs_kernels.h"
#include "dbl_tabs_testers.h"

using namespace std;

double dbl_Difference_Sum ( int n, double *x, double *y )
{
   double result = 0.0;

   for(int i=0; i<n; i++) result = result + abs(x[i] - y[i]);

   return result;
}

double dbl_Column_Sum ( int dim, int col, double **A )
{
   double result = 0.0;

   for(int i=0; i<dim; i++) result = result + abs(A[i][col]);

   return result;
}

double dbl_Max_Column_Sum ( int dim, double **A )
{
   double result = dbl_Column_Sum(dim,0,A);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = dbl_Column_Sum(dim,j,A);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double dbl_condition ( int dim, double **A, double **invA )
{
   double Amaxcolsum = dbl_Max_Column_Sum(dim,A);
   double invAmaxcolsum = dbl_Max_Column_Sum(dim,invA);

   return Amaxcolsum*invAmaxcolsum;
}

double dbl_Matrix_Difference_Sum ( int n, double **A, double **B )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(A[i][j] - B[i][j]);

   return result;
}

void dbl_random_upper_factor ( int dim, double **A )
{
   random_dbl_matrix(dim,dim,A);

   int *pivots = new int[dim];

   CPU_dbl_factors_lufac(dim,A,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++) A[i][j] = 0.0;

   free(pivots);
}

void test_real_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   // random_dbl_upper_matrix(dim,dim,A);
   dbl_random_upper_factor(dim,A);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhs[i] << endl;
   }
   double **invA_h = new double*[dim];
   for(int i=0; i<dim; i++) invA_h[i] = new double[dim];

   CPU_dbl_upper_inverse(dim,A,invA_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invA_h[i][j] << endl;
   }
   double **invA_d = new double*[dim];
   for(int i=0; i<dim; i++) invA_d[i] = new double[dim];

   GPU_dbl_upper_inverse(dim,A,invA_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invA_d[i][j] << endl;
   }
   cout << "   Sum of errors : "
        << dbl_Matrix_Difference_Sum(dim,invA_h,invA_d) << endl;

   double *x = new double[dim];
   for(int i=0; i<dim; i++)
   {
      x[i] = 0.0;
      for(int j=0; j<dim; j++) x[i] = x[i] + invA_h[i][j]*rhs[j];
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << x[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : " << dbl_Difference_Sum(dim,sol,x) << endl;
   cout << "Condition number : "
        << dbl_condition(dim,A,invA_h) << endl;
}

void test_real_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   // random_dbl_upper_matrix(dim,dim,A);
   dbl_random_upper_factor(dim,A);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhs[i] << endl;
   }
   double *x = new double[dim];

   CPU_dbl_upper_tiled_solver(dim,sizetile,numtiles,A,rhs,x);

   if(verbose > 0)
   {
      cout << "The solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << x[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : " << dbl_Difference_Sum(dim,sol,x) << endl;
}
