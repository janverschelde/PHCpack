// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "random2_vectors.h"
#include "random2_series.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomial_systems.h"

using namespace std;

void make_real2_exponentials
 ( int dim, int  deg, double **shi, double **slo )
{
   double rndhi,rndlo,acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      random_double_double(&rndhi,&rndlo);        // rnd is in [-1, +1]
      ddf_div(rndhi,rndlo,2.0,0.0,&acchi,&acclo); // acc is in [-0.5, +0.5]
      
      if(rndhi < 0)
      {
         // rnd = rnd - 1.5; if -0.5 <= rnd < 0, rnd - 1.5 is in [-2, -1.5]
         ddf_sub(acchi,acclo,1.5,0.0,&rndhi,&rndlo);
      }
      else
      {
         // rnd = rnd + 1.5; if  0 < rnd <= 0.5, rnd + 1.5 is in [+1.5, +2]
         ddf_add(acchi,acclo,1.5,0.0,&rndhi,&rndlo);
      }
      dbl2_exponential(deg,rndhi,rndlo,shi[i],slo[i]);
   }
}

void make_complex2_exponentials
 ( int dim, int deg,
   double **srehi, double **srelo, double **simhi, double **simlo )
{
   double xrehi,xrelo,ximhi,ximlo,yhi,ylo;

   for(int i=0; i<dim; i++)
   {
      random_double_double(&xrehi,&xrelo); // cosine of some angle
 
      ddf_sqr(xrehi,xrelo,&yhi,&ylo);        // y = cos^2
      ddf_minus(&yhi,&ylo);                  // y = -cos^2
      ddf_inc_d(&yhi,&ylo,1.0);              // y = 1 - cos^2
      ddf_sqrt(yhi,ylo,&ximhi,&ximlo);       // sin is sqrt(1-cos^2)

      cmplx2_exponential
         (deg,xrehi,xrelo,ximhi,ximlo,srehi[i],srelo[i],simhi[i],simlo[i]);
   }
}

void evaluate_real2_monomials
 ( int dim, int deg, int **rowsA,
   double **shi, double **slo, double **rhshi, double **rhslo )
{
   const int degp1 = deg+1;

   double *acchi = new double[degp1]; // accumulates product
   double *acclo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_dbl2_product
                  (deg,shi[j],slo[j],rhshi[i],rhslo[i],acchi,acclo);

               for(int L=0; L<degp1; L++)
               {
                  rhshi[i][L] = acchi[L];
                  rhslo[i][L] = acclo[L];
               }
            }
         }
      }
   }
   free(acchi); free(acclo);
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
