// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "random2_vectors.h"
#include "random2_series.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomial_systems.h"

using namespace std;

void make_complex2_exponentials
 ( int dim, int deg,
   double **srehi, double **srelo, double **simhi, double **simlo )
{
   double xrehi,xrelo,ximhi,ximlo,yhi,ylo;

   cout << scientific << setprecision(16);

   for(int i=0; i<dim; i++)
   {
      random_double_double(&xrehi,&xrelo); // cosine of some angle
 
      cout << "  xre[" << i << "] : " << xrehi << "  " << xrelo << endl;

      ddf_sqr(xrehi,xrelo,&yhi,&ylo);        // y = cos^2

      cout << "xre^2[" << i << "] : " << yhi << "  " << ylo << endl;

      ddf_minus(&yhi,&ylo);                  // y = -cos^2
      ddf_inc_d(&yhi,&ylo,1.0);              // y = 1 - cos^2

      cout << "    y[" << i << "] : " << yhi << "  " << ylo << endl;

      ddf_sqrt(yhi,ylo,&ximhi,&ximlo);       // sin is sqrt(1-cos^2)

      cout << "  xim[" << i << "] : " << ximhi << "  " << ximlo << endl;

      cmplx2_exponential
         (deg,xrehi,xrelo,ximhi,ximlo,srehi[i],srelo[i],simhi[i],simlo[i]);
   }
}

void evaluate_complex2_monomials
 ( int dim, int deg, int **rowsA,
   double **srehi, double **srelo, double **simhi, double **simlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo )
{
   const int degp1 = deg+1;

   double *accrehi = new double[degp1]; // accumulates product
   double *accrelo = new double[degp1];
   double *accimhi = new double[degp1];
   double *accimlo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx2_product
                  (deg,srehi[j],srelo[j],simhi[j],simlo[j],
                   rhsrehi[i],rhsrelo[i],rhsimhi[i],rhsimlo[i],
                   accrehi,accrelo,accimhi,accimlo);

               for(int L=0; L<degp1; L++)
               {
                  rhsrehi[i][L] = accrehi[L];
                  rhsrelo[i][L] = accrelo[L];
                  rhsimhi[i][L] = accimhi[L];
                  rhsimlo[i][L] = accimlo[L];
               }
            }
         }
      }
   }
   free(accrehi); free(accimhi);
   free(accrelo); free(accimlo);
}
