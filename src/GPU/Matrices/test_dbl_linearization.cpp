/* Tests operations on series vectors in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "random_matrices.h"
#include "dbl_linearization.h"

using namespace std;

int test_dbl_linearization ( int dim, int deg );
/*
 * Tests the linearization of a random vectors of dimension dim,
 * with series truncated at degree deg, for real data. */

int test_cmplx_linearization ( int dim, int deg );
/*
 * Tests the linearization of a random vectors of dimension dim,
 * with series truncated at degree deg, for complex data. */

int main ( void )
{
   srand(time(NULL));

   cout << "Give the dimension of the vector : ";
   int dim; cin >> dim;

   cout << "Give the truncation degree of the series : ";
   int deg; cin >> deg;

   cout << "Testing on real data ..." << endl;
   test_dbl_linearization(dim,deg);

   cout << "Testing on complex data ..." << endl;
   test_cmplx_linearization(dim,deg);

   return 0;
}

int test_dbl_linearization ( int dim, int deg )
{
   double *x = new double[dim];
   double **v = new double*[dim];

   for(int i=0; i<dim; i++) v[i] = new double[deg+1];

   random_dbl_series_vector(dim,deg,x,v,false);

   cout << scientific << setprecision(16);

   cout << "A vector of random series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "v[" << i << "] :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << v[i][j];
         if(j == 0)
            cout << endl;
         else if(j == 1)
            cout << "*t" << endl;
         else
            cout << "*t^" << j << endl;
      }
   }
   double **w = new double*[deg+1];

   for(int i=0; i<=deg; i++) w[i] = new double[dim];

   dbl_linear_series_vector(dim,deg,v,w);

   cout << "The coefficient vectors of the series :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << "w[" << i << "] :" << endl;
      for(int j=0; j<dim; j++)
         cout << w[i][j] << endl;
   }
   return 0;
}

int test_cmplx_linearization ( int dim, int deg )
{
   double *xre = new double[dim];
   double *xim = new double[dim];
   double **vre = new double*[dim];
   double **vim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      vre[i] = new double[deg+1];
      vim[i] = new double[deg+1];
   }
   random_cmplx_series_vector(dim,deg,xre,xim,vre,vim,false);

   cout << scientific << setprecision(16);

   cout << "A vector of random series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "v[" << i << "] :" << endl;
      for(int j=0; j<=deg; j++)
      {
         if(j > 0) cout << "+ ( ";
         cout << vre[i][j] << "  " << vim[i][j];
         if(j == 0)
            cout << endl;
         else if(j == 1)
            cout << " )*t" << endl;
         else
            cout << " )*t^" << j << endl;
      }
   }
   double **wre = new double*[deg+1];
   double **wim = new double*[deg+1];

   for(int i=0; i<=deg; i++)
   {
      wre[i] = new double[dim];
      wim[i] = new double[dim];
   }
   cmplx_linear_series_vector(dim,deg,vre,vim,wre,wim);

   cout << "The coefficient vectors of the series :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << "w[" << i << "] :" << endl;
      for(int j=0; j<dim; j++)
         cout << wre[i][j] << "  " << wim[i][j] << endl;
   }
   return 0;
}
