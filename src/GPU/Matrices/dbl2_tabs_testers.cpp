// The file dbl2_tabs_testers.cpp defines the functions specified in
// the file dbl2_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "double_double_functions.h"
#include "random2_matrices.h"
#include "dbl2_factorizations.h"
#include "dbl2_tabs_host.h"
#include "dbl2_tabs_kernels.h"

using namespace std;

double dbl2_Difference_Sum
 ( int n, double *xhi, double *xlo, double *yhi, double *ylo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xhi[i] - yhi[i]) + abs(xlo[i] - ylo[i]);

   return result;
}

double dbl2_Column_Sum ( int dim, int col, double **Ahi, double **Alo )
{
   double resulthi = 0.0;
   double resultlo = 0.0;

   for(int i=0; i<dim; i++)
      ddf_inc(&resulthi,&resultlo,abs(Ahi[i][col]),abs(Alo[i][col]));

   return resulthi;
}

double dbl2_Max_Column_Sum ( int dim, double **Ahi, double **Alo )
{
   double result = dbl2_Column_Sum(dim,0,Ahi,Alo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = dbl2_Column_Sum(dim,j,Ahi,Alo);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double dbl2_condition
 ( int dim, double **Ahi, double **Alo, double **invAhi, double **invAlo )
{
   double Amaxcolsum = dbl2_Max_Column_Sum(dim,Ahi,Alo);
   double invAmaxcolsum = dbl2_Max_Column_Sum(dim,invAhi,invAlo);

   return Amaxcolsum*invAmaxcolsum;
}

double dbl2_Matrix_Difference_Sum
 ( int n, double **Ahi, double **Alo, double **Bhi, double **Blo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Ahi[i][j] - Bhi[i][j])
                         + abs(Alo[i][j] - Blo[i][j]);

   return result;
}

void dbl2_random_upper_factor ( int dim, double **Ahi, double **Alo )
{
   random_dbl2_matrix(dim,dim,Ahi,Alo);

   int *pivots = new int[dim];

   CPU_dbl2_factors_lufac(dim,Ahi,Alo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Ahi[i][j] = 0.0;
         Alo[i][j] = 0.0;
      }

   free(pivots);
}

void test_real2_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }
   // random_dbl2_upper_matrix(dim,dim,Ahi,Alo);
   dbl2_random_upper_factor(dim,Ahi,Alo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   double *solhi = new double[dim];
   double *sollo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solhi[i] = 1.0;
      sollo[i] = 0.0;
   }
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];
   double acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
      for(int j=0; j<dim; j++)  // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhshi[i] << "  " << rhslo[i] << endl;
   }
   double **invAhi_h = new double*[dim];
   double **invAlo_h = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAhi_h[i] = new double[dim];
      invAlo_h[i] = new double[dim];
   }
   CPU_dbl2_upper_inverse(dim,Ahi,Alo,invAhi_h,invAlo_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invAhi_h[i][j] << "  " << invAlo_h[i][j] << endl;
   }
   double **invAhi_d = new double*[dim];
   double **invAlo_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAhi_d[i] = new double[dim];
      invAlo_d[i] = new double[dim];
   }
   GPU_dbl2_upper_inverse(dim,Ahi,Alo,invAhi_d,invAlo_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAhi_d[i][j] << "  " << invAlo_d[i][j] << endl;
   }
   cout << "   Sum of errors : "
        << dbl2_Matrix_Difference_Sum
              (dim,invAhi_h,invAlo_h,invAhi_d,invAlo_d)
        << endl;

   double *xhi = new double[dim];
   double *xlo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      xhi[i] = 0.0;
      xlo[i] = 0.0;
      for(int j=0; j<dim; j++)   // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         ddf_mul(invAhi_h[i][j],invAlo_h[i][j],rhshi[j],rhslo[j],
                 &acchi,&acclo);
         ddf_inc(&xhi[i],&xlo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhi[i] << "  " << xlo[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi,xlo) << endl;
   cout << "Condition number : "
        << dbl2_condition(dim,Ahi,Alo,invAhi_h,invAlo_h) << endl;
}

void test_real2_upper_tiling ( void )
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

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++) 
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }

   // random_dbl_upper_matrix(dim,dim,A);
   dbl2_random_upper_factor(dim,Ahi,Alo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   double *solhi = new double[dim];
   double *sollo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solhi[i] = 1.0;
      sollo[i] = 0.0;
   }
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];
   double acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshi[i] << "  " << rhslo[i] << endl;
   }
   double *xhi = new double[dim];
   double *xlo = new double[dim];

   CPU_dbl2_upper_tiled_solver
      (dim,sizetile,numtiles,Ahi,Alo,rhshi,rhslo,xhi,xlo);

   if(verbose > 0)
   {
      cout << "The solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhi[i] << "  " << xlo[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi,xlo) << endl;
}
