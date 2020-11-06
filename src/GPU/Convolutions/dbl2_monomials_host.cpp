/* The file dbl2_monomials_host.cpp defines functions specified
 * in dbl2_monomials_host.h. */

/* The algorithm to compute forward, backward, and cross products
 * (denoted respectively by arrays f, b, and c)
 * for a monomial cff*x[0]*x[1]* .. *x[n-1] goes as follows:
 *
 * f[0] := cff*x[0]
 * for i from 1 to n-1 do f[i] := f[i-1]*x[i]
 * if n > 2 then
 *    b[0] := x[n-1]*x[n-2]
 *    for i from 1 to n-3 do b[i] := b[i-1]*x[n-2-i]
 *    b[n-3] := b[n-3]*cff
 *    if n = 3 then
 *       c[0] = f[0]*x[2]
 *    else
 *       for i from 0 to n-4 do c[i] := f[i]*b[n-4-i]
 *       c[n-3] := f[n-3]*x[n-1]
 *
 * Compared to the evaluation and differentiation of a product of variables,
 * (without coefficient cff), two extra multiplications must be done,
 * but this is better than n+1 multiplications with cff afterwards. */

#include <cstdlib>
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"

void CPU_dbl2_speel
 ( int nvr, int deg, int *idx, double *cffhi, double *cfflo,
   double **inputhi, double **inputlo, double **forwardhi,
   double **forwardlo, double **backwardhi, double **backwardlo,
   double **crosshi, double **crosslo )
{
   int ix1 = idx[0];
   int ix2;

   CPU_dbl2_product(deg,cffhi,cfflo,inputhi[ix1],inputlo[ix1],
                   forwardhi[0],forwardlo[0]);           // f[0] = cff*x[0] 

   for(int i=1; i<nvr; i++)
   {                                               // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_dbl2_product(deg,forwardhi[i-1],forwardlo[i-1],
                       inputhi[ix2],inputlo[ix2],forwardhi[i],forwardlo[i]);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1];                            // b[0] = x[n-1]*x[n-2]
      ix2 = idx[nvr-2];
      CPU_dbl2_product(deg,inputhi[ix1],inputlo[ix1],inputhi[ix2],
                       inputlo[ix2],backwardhi[0],backwardlo[0]);
      for(int i=1; i<nvr-2; i++)
      {                                            // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];
         CPU_dbl2_product(deg,backwardhi[i-1],backwardlo[i-1],inputhi[ix2],
                          inputlo[ix2],backwardhi[i],backwardlo[i]);
      }
                                                   // b[n-3] = cff*b[n-3]
      CPU_dbl2_product(deg,backwardhi[nvr-3],backwardlo[nvr-3],cffhi,cfflo,
                       crosshi[0],crosslo[0]);
      // cross[0] is work space, cannot write into backward[nvr-3]
      for(int i=0; i<=deg; i++)
      {
         backwardhi[nvr-3][i] = crosshi[0][i];
         backwardlo[nvr-3][i] = crosslo[0][i];
      }
      if(nvr == 3)
      {                                            // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_dbl2_product(deg,forwardhi[0],forwardlo[0],inputhi[ix2],
                          inputlo[ix2],crosshi[0],crosslo[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                         // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_dbl2_product(deg,forwardhi[i],forwardlo[i],backwardhi[ix2],
                             backwardlo[ix2],crosshi[i],crosslo[i]);
         }
                                                   // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_dbl2_product(deg,forwardhi[nvr-3],forwardlo[nvr-3],inputhi[ix2],
                          inputlo[ix2],crosshi[nvr-3],crosslo[nvr-3]);
      }
   }
}

void CPU_cmplx2_speel
 ( int nvr, int deg, int *idx,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo,
   double **inputrehi, double **inputrelo, double **inputimhi,
   double **inputimlo, double **forwardrehi, double **forwardrelo,
   double **forwardimhi, double **forwardimlo, double **backwardrehi,
   double **backwardrelo, double **backwardimhi, double **backwardimlo,
   double **crossrehi, double **crossrelo, double **crossimhi,
   double **crossimlo )
{
}

void CPU_dbl2_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cffhi, double *cfflo,
   double **inputhi, double **inputlo, double **outputhi, double **outputlo )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_dbl2_product(deg,inputhi[ix],inputlo[ix],cffhi,cfflo,
                       outputhi[dim],outputlo[dim]);

      for(int i=0; i<=deg; i++)
      {
         outputhi[ix][i] = cffhi[i];
         outputlo[ix][i] = cfflo[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_dbl2_product(deg,inputhi[ix1],inputlo[ix1],inputhi[ix2],
                       inputlo[ix2],outputhi[dim],outputlo[dim]);
      CPU_dbl2_product(deg,outputhi[dim],outputlo[dim],cffhi,cfflo,
                       outputhi[dim],outputlo[dim]);

      CPU_dbl2_product(deg,cffhi,cfflo,inputhi[ix1],inputlo[ix1],
                       outputhi[ix2],outputlo[ix2]);
      CPU_dbl2_product(deg,cffhi,cfflo,inputhi[ix2],inputlo[ix2],
                       outputhi[ix1],outputlo[ix1]);
   }
   else
   {
      double **forwardhi = new double*[nvr];
      double **forwardlo = new double*[nvr];
      double **backwardhi = new double*[nvr-2];
      double **backwardlo = new double*[nvr-2];
      double **crosshi = new double*[nvr-2];
      double **crosslo = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardhi[i] = new double[deg+1];
         forwardlo[i] = new double[deg+1];
         backwardhi[i] = new double[deg+1];
         backwardlo[i] = new double[deg+1];
         crosshi[i] = new double[deg+1];
         crosslo[i] = new double[deg+1];
      }
      forwardhi[nvr-2] = new double[deg+1];
      forwardlo[nvr-2] = new double[deg+1];
      forwardhi[nvr-1] = new double[deg+1];
      forwardlo[nvr-1] = new double[deg+1];

      CPU_dbl2_speel(nvr,deg,idx,cffhi,cfflo,inputhi,inputlo,forwardhi,
                     forwardlo,backwardhi,backwardlo,crosshi,crosslo);

      // assign the value of the monomial

      for(int i=0; i<deg+1; i++)
      {
         outputhi[dim][i] = forwardhi[nvr-1][i];
         outputlo[dim][i] = forwardlo[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];       // derivative with respect to x[n-1]
         for(int i=0; i<deg+1; i++)
         {
            outputhi[ix][i] = forwardhi[nvr-2][i];
            outputlo[ix][i] = forwardlo[nvr-2][i];
         }
         ix = idx[0];               // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputhi[ix][i] = backwardhi[nvr-3][i];
            outputlo[ix][i] = backwardlo[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {                          // derivative with respect to x[k]
            ix = idx[k];
            for(int i=0; i<deg+1; i++)
            {
               outputhi[ix][i] = crosshi[k-1][i];
               outputlo[ix][i] = crosslo[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardhi[i]); free(backwardhi[i]); free(crosshi[i]);
         free(forwardlo[i]); free(backwardlo[i]); free(crosslo[i]);
      }
      free(forwardhi[nvr-2]); free(forwardhi[nvr-1]);
      free(forwardlo[nvr-2]); free(forwardlo[nvr-1]);
      free(forwardhi); free(backwardhi); free(crosshi);
      free(forwardlo); free(backwardlo); free(crosslo);
   }
}

void CPU_cmplx2_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo,
   double **inputrehi, double **inputrelo, double **inputimhi,
   double **inputimlo, double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo )
{
}
