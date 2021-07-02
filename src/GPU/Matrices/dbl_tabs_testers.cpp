// The file dbl_tabs_testers.cpp defines the function with prototypes
// in dbl_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random_matrices.h"
#include "dbl_tabs_host.h"
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

void test_real_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   random_dbl_upper_matrix(dim,dim,A);

   cout << scientific << setprecision(16);

   cout << "A random upper triangular matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;

   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : " << rhs[i] << endl;

   double **invA = new double*[dim];
   for(int i=0; i<dim; i++) invA[i] = new double[dim];

   CPU_dbl_upper_inverse(dim,A,invA);

   cout << "The inverse of the upper triangular matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "invA[" << i << "][" << j << "] : " << invA[i][j] << endl;

   double *x = new double[dim];
   for(int i=0; i<dim; i++)
   {
      x[i] = 0.0;
      for(int j=0; j<dim; j++) x[i] = x[i] + invA[i][j]*rhs[j];
   }
   cout << "The solution computed with the inverse :" << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : " << x[i] << endl;

   cout << scientific << setprecision(2);
   cout << "   Sum of errors : " << dbl_Difference_Sum(dim,sol,x) << endl;
   cout << "Condition number : " << dbl_condition(dim,A,invA) << endl;
}

void test_real_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   random_dbl_upper_matrix(dim,dim,A);

   cout << scientific << setprecision(16);

   cout << "A random upper triangular matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;

   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : " << rhs[i] << endl;

   double *x = new double[dim];

   CPU_dbl_upper_tiled_solver(dim,sizetile,numtiles,A,rhs,x);

   cout << "The solution computed with tiling :" << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : " << x[i] << endl;

   cout << scientific << setprecision(2);
   cout << "   Sum of errors : " << dbl_Difference_Sum(dim,sol,x) << endl;
}
