// The file dbl_monomial_systems.cpp defines the functions specified in
// the file dbl_monomial_systems.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random_series.h"
#include "dbl_convolutions_host.h"
#include "dbl_monomial_systems.h"

using namespace std;

void make_complex_exponentials
 ( int dim, int deg, double *angles, double **sre, double **sim )
{
   double xre,xim;

   for(int i=0; i<dim; i++)
   {
      angles[i] = random_angle();
      xre = cos(angles[i]);
      xim = sin(angles[i]);

      cmplx_exponential(deg,xre,xim,sre[i],sim[i]);
   }
}

void evaluate_complex_monomials
 ( int dim, int deg, int **rowsA,
   double **xre, double **xim, double **rhsre, double **rhsim )
{
   const int degp1 = deg+1;

   double *accre = new double[degp1]; // accumulates product
   double *accim = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx_product
                  (deg,xre[j],xim[j],rhsre[i],rhsim[i],accre,accim);

               for(int L=0; L<degp1; L++)
               {
                  rhsre[i][L] = accre[L];
                  rhsim[i][L] = accim[L];
               }
            }
         }
      }
   }
   free(accre); free(accim);
}
