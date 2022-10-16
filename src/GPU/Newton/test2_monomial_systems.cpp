/* Tests the making of a monomial systems in double double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "random2_series.h"
#include "unimodular_matrices.h"
#include "dbl2_monomial_systems.h"

using namespace std;

int main ( void )
{
   cout << "testing the making of a monomial system ..." << endl;

   const int vrblvl = 2;
   int seed,dim,deg,size,nbritr;

   cout << "-> give the seed (0 for time) : "; cin >> seed;
   cout << "-> give the number of series : "; cin >> dim;
   while(true)
   {
      cout << "-> give the truncation degree : "; cin >> deg;
      if(deg > 0) break;
      cout << "The degree must be one or larger.  Retry" << endl;
   }
   cout << "-> give the size of the numbers : "; cin >> size;
   cout << "-> give the number of iterations : "; cin >> nbritr;

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

// make the unimodular matrix of the system

   int **rowsA = new int*[dim];  // exponents in the rows
   int *nvr = new int[dim];      // number of variables in each monomial
   int **idx = new int*[dim];    // indexes of variables in each monomial
   int **exp = new int*[dim];    // exponents of the variables
   int *nbrfac = new int[dim];   // number of exponents > 1 in each monomial
   int **expfac = new int*[dim]; // exponents of the common factors

   make_monomial_system
      (dim,size,1,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);

   int *expsol = new int[dim];
   int sing = exponents_check(dim,rowsA,expsol,vrblvl);

// generate the solution series

   double **solrehi = new double*[dim];
   double **solrelo = new double*[dim];
   double **solimhi = new double*[dim];
   double **solimlo = new double*[dim];

   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      solrehi[i] = new double[degp1];
      solrelo[i] = new double[degp1];
      solimhi[i] = new double[degp1];
      solimlo[i] = new double[degp1];
   }
   make_complex2_exponentials(dim,deg,solrehi,solrelo,solimhi,solimlo);

   cout << scientific << setprecision(16);

   for(int j=0; j<degp1; j++)
   {
      cout << "coefficients of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "c[" << i << "][" << j << "] : ";     
         cout << solrehi[i][j] << "  " << solrelo[i][j] << endl
              << "  "
              << solimhi[i][j] << "  " << solimlo[i][j] << endl;
      }
   }

// compute the right hand sides via evaluation

   double **rhsrehi = new double*[dim];
   double **rhsrelo = new double*[dim];
   double **rhsimhi = new double*[dim];
   double **rhsimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhsrehi[i] = new double[degp1];
      rhsrelo[i] = new double[degp1];
      rhsimhi[i] = new double[degp1];
      rhsimlo[i] = new double[degp1];

      rhsrehi[i][0] = 1.0;     // initialize product to one
      rhsrelo[i][0] = 0.0;
      rhsimhi[i][0] = 0.0;
      rhsimlo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         rhsrehi[i][k] = 0.0;
         rhsrelo[i][k] = 0.0;
         rhsimhi[i][k] = 0.0;
         rhsimlo[i][k] = 0.0;
      }
   }
   evaluate_complex2_monomials
      (dim,deg,rowsA,solrehi,solrelo,solimhi,solimlo,
                     rhsrehi,rhsrelo,rhsimhi,rhsimlo);

   cout << "the evaluated right hand sides :" << endl;

   for(int j=0; j<degp1; j++)
   {
      cout << "coefficients of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "r[" << i << "][" << j << "] : ";     
         cout << rhsrehi[i][j] << "  " << rhsrelo[i][j] << endl
              << "  "
              << rhsimhi[i][j] << "  " << rhsimlo[i][j] << endl;
      }
   }
   return 0;
}
