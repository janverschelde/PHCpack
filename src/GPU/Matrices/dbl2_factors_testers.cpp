/* The file dbl2_factors_testers.cpp define the functions specified in
   the file dbl2_factors_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "double_double_functions.h"
#include "random2_matrices.h"
#include "dbl2_factorizations.h"

using namespace std;

void test_factors_real2_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }
   random_dbl2_matrix(dim,dim,Ahi,Alo);

   cout << scientific << setprecision(16);

   cout << "A random triangular matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "A[" << i << "][" << j << "] : "
              << Ahi[i][j] << "  " << Alo[i][j] << endl;

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
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : "
           << rhshi[i] << "  " << rhslo[i] << endl;

   double *xhi = new double[dim];
   double *xlo = new double[dim];
   int *pivots = new int[dim];

   CPU_dbl2_factors_lusolve(dim,Ahi,Alo,pivots,rhshi,rhslo,xhi,xlo);

   cout << "The computed solution :" << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : "
           << xhi[i] << "  " << xlo[i] << endl;
}
