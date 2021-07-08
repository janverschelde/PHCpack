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

   cout << "A random matrix :" << endl;
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

void test_factors_cmplx2_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Arehi = new double*[dim];
   double **Arelo = new double*[dim];
   double **Aimhi = new double*[dim];
   double **Aimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Arehi[i] = new double[dim];
      Arelo[i] = new double[dim];
      Aimhi[i] = new double[dim];
      Aimlo[i] = new double[dim];
   }
   random_cmplx2_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);

   cout << scientific << setprecision(16);

   cout << "A random matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         cout << "A[" << i << "][" << j << "]re : "
              << Arehi[i][j] << "  " << Arelo[i][j] << endl;
         cout << "A[" << i << "][" << j << "]im : "
              << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
      }

   double *solrehi = new double[dim];
   double *solrelo = new double[dim];
   double *solimhi = new double[dim];
   double *solimlo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solrehi[i] = 1.0; solrelo[i] = 0.0;
      solimhi[i] = 0.0; solimlo[i] = 0.0;
   }
   double *rhsrehi = new double[dim];
   double *rhsrelo = new double[dim];
   double *rhsimhi = new double[dim];
   double *rhsimlo = new double[dim];
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;

   for(int i=0; i<dim; i++)
   {
      rhsrehi[i] = 0.0; rhsrelo[i] = 0.0;
      rhsimhi[i] = 0.0; rhsimlo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Arehi[i][j],Arelo[i][j],solrehi[j],solrelo[j],
                 &acc1hi,&acc1lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solimhi[j],solimlo[j],
                 &acc2hi,&acc2lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solrehi[j],solrelo[j],
                 &acc3hi,&acc3lo);
         ddf_mul(Arehi[i][j],Arelo[i][j],solimhi[j],solimlo[j],
                 &acc4hi,&acc4lo);
         ddf_inc(&rhsrehi[i],&rhsrelo[i],acc1hi,acc1lo);
         ddf_dec(&rhsrehi[i],&rhsrelo[i],acc2hi,acc2lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc3hi,acc3lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc4hi,acc4lo);
      }
   }
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "b[" << i << "]re : "
           << rhsrehi[i] << "  " << rhsrelo[i] << endl;
      cout << "b[" << i << "]im : "
           << rhsimhi[i] << "  " << rhsimlo[i] << endl;
   }
   double *xrehi = new double[dim];
   double *xrelo = new double[dim];
   double *ximhi = new double[dim];
   double *ximlo = new double[dim];
   int *pivots = new int[dim];

   CPU_cmplx2_factors_lusolve
      (dim,Arehi,Arelo,Aimhi,Aimlo,pivots,
       rhsrehi,rhsrelo,rhsimhi,rhsimlo,xrehi,xrelo,ximhi,ximlo);

   cout << "The computed solution :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "x[" << i << "]re : "
           << xrehi[i] << "  " << xrelo[i] << endl;
      cout << "x[" << i << "]im : "
           << ximhi[i] << "  " << ximlo[i] << endl;
   }
}
