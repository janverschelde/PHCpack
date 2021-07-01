// The file dbl_tabs_testers.cpp defines the function with prototypes
// in dbl_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "random_matrices.h"
#include "dbl_tabs_host.h"
#include "dbl_tabs_testers.h"

using namespace std;

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
}
