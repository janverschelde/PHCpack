/* Tests operations on series vectors in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "random_matrices.h"
#include "dbl_convolutions_host.h"

using namespace std;

void inner_product
 ( int dim, int deg, double **x, double **y, double *z );
/*
 * DESCRIPTION :
 *   Computes the product of two vectors x and y of power series
 *   and assigns the result to z.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors x and y;
 *   deg      truncation degree of the series;
 *   x        dim series truncated to degree deg;
 *   y        dim series truncated to degree deg;
 *   z        space for deg+1 doubles.
 *
 * ON RETURN :
 *   z        the sum of all x[k]*y[k] for k from 0 to dim-1,
 *            as a power series truncated to the degree deg. */

int main ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   double *x = new double[dim];
   double **px = new double*[dim];
   double **mx = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      px[i] = new double[deg+1];
      mx[i] = new double[deg+1];
   }
   random_dbl_series_vectors(dim,deg,x,px,mx);

   cout << scientific << setprecision(16);

   for(int k=0; k<dim; k++)
   {
      cout << "a random x : " << x[k] << endl;
      cout << "series for exp(+x) :" << endl; 
      for(int i=0; i<=deg; i++) cout << px[k][i] << endl;
      cout << "series for exp(-x) :" << endl; 
      for(int i=0; i<=deg; i++) cout << mx[k][i] << endl;
   }

   double *ip = new double[deg+1];

   inner_product(dim,deg,px,mx,ip);

   cout << "the inner product :" << endl;
   for(int i=0; i<=deg; i++) cout << ip[i] << endl;
  
   return 0;
}

void inner_product
 ( int dim, int deg, double **x, double **y, double *z )
{
   double *prod = new double[deg+1];

   for(int i=0; i<=deg; i++) z[i] = 0.0;

   for(int k=0; k<dim; k++)
   {
      CPU_dbl_product(deg,x[k],y[k],prod);
      for(int i=0; i<=deg; i++)
         z[i] = z[i] + prod[i];
   }
   free(prod);
}
