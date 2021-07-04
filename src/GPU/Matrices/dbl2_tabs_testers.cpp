// The file dbl2_tabs_testers.cpp defines the functions specified in
// the file dbl2_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "double_double_functions.h"
#include "random2_matrices.h"
#include "dbl2_tabs_host.h"

using namespace std;

void test_real2_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }
   random_dbl2_upper_matrix(dim,dim,Ahi,Alo);

   cout << scientific << setprecision(16);

   cout << "A random upper triangular matrix :" << endl;
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
      for(int j=0; j<dim; j++)  // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : " << rhshi[i] << "  " << rhslo[i] << endl;

   double **invAhi_h = new double*[dim];
   double **invAlo_h = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAhi_h[i] = new double[dim];
      invAlo_h[i] = new double[dim];
   }
   CPU_dbl2_upper_inverse(dim,Ahi,Alo,invAhi_h,invAlo_h);

   cout << "The CPU inverse of the upper triangular matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "invA_h[" << i << "][" << j << "] : "
              << invAhi_h[i][j] << "  " << invAlo_h[i][j] << endl;

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
   cout << "The solution computed with the CPU inverse :" << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : "
           << xhi[i] << "  " << xlo[i] << endl;
}
