/* The file dbl4_monomials_host.cpp defines the funnctions specified
 * in dbl4_monomials_host.h. */

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
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"

void CPU_dbl4_speel
 ( int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo )
{
   int ix1 = idx[0];
   int ix2;                                            // f[0] = cff*x[0] 

   CPU_dbl4_product(deg,cffhihi,cfflohi,cffhilo,cfflolo,
      inputhihi[ix1],inputlohi[ix1],inputhilo[ix1],inputlolo[ix1],
      forwardhihi[0],forwardlohi[0],forwardhilo[0],forwardlolo[0]); 

   for(int i=1; i<nvr; i++)
   {                                                // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_dbl4_product(deg,
         forwardhihi[i-1],forwardlohi[i-1],forwardhilo[i-1],forwardlolo[i-1],
         inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
         forwardhihi[i],forwardlohi[i],forwardhilo[i],forwardlolo[i]);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1];                            // b[0] = x[n-1]*x[n-2]
      ix2 = idx[nvr-2];
      CPU_dbl4_product(deg,
         inputhihi[ix1],inputlohi[ix1],inputhilo[ix1],inputlolo[ix1],
         inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
         backwardhihi[0],backwardlohi[0],backwardhilo[0],backwardlolo[0]);

      for(int i=1; i<nvr-2; i++)
      {                                            // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];
         CPU_dbl4_product(deg,backwardhihi[i-1],backwardlohi[i-1],
                              backwardhilo[i-1],backwardlolo[i-1],
            inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
            backwardhihi[i],backwardlohi[i],backwardhilo[i],backwardlolo[i]);
      }
                                                   // b[n-3] = cff*b[n-3]
      CPU_dbl4_product(deg,backwardhihi[nvr-3],backwardlohi[nvr-3],
                           backwardhilo[nvr-3],backwardlolo[nvr-3],
         cffhihi,cfflohi,cffhilo,cfflolo,
         crosshihi[0],crosslohi[0],crosshilo[0],crosslolo[0]);
      // cross[0] is work space, cannot write into backward[nvr-3]
      for(int i=0; i<=deg; i++)
      {
         backwardhihi[nvr-3][i] = crosshihi[0][i];
         backwardlohi[nvr-3][i] = crosslohi[0][i];
         backwardhilo[nvr-3][i] = crosshilo[0][i];
         backwardlolo[nvr-3][i] = crosslolo[0][i];
      }
      if(nvr == 3)
      {                                            // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_dbl4_product(deg,
            forwardhihi[0],forwardlohi[0],forwardhilo[0],forwardlolo[0],
            inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
            crosshihi[0],crosslohi[0],crosshilo[0],crosslolo[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                         // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_dbl4_product(deg,
               forwardhihi[i],forwardlohi[i],forwardhilo[i],forwardlolo[i],
               backwardhihi[ix2],backwardlohi[ix2],
               backwardhilo[ix2],backwardlolo[ix2],
               crosshihi[i],crosslohi[i],crosshilo[i],crosslolo[i]);
         }
                                                   // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_dbl4_product(deg,
            forwardhihi[nvr-3],forwardlohi[nvr-3],
            forwardhilo[nvr-3],forwardlolo[nvr-3],
            inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
            crosshihi[nvr-3],crosslohi[nvr-3],
            crosshilo[nvr-3],crosslolo[nvr-3]);
      }
   }
}

void CPU_cmplx4_speel
 ( int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **backwardrehihi, double **backwardrelohi,
   double **backwardrehilo, double **backwardrelolo,
   double **backwardimhihi, double **backwardimlohi,
   double **backwardimhilo, double **backwardimlolo,
   double **crossrehihi, double **crossrelohi,
   double **crossrehilo, double **crossrelolo,
   double **crossimhihi, double **crossimlohi,
   double **crossimhilo, double **crossimlolo )
{
}

void CPU_dbl4_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_dbl4_product(deg,
         inputhihi[ix],inputlohi[ix],inputhilo[ix],inputlolo[ix],
         cffhihi,cfflohi,cffhilo,cfflolo,
         outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim]);

      for(int i=0; i<=deg; i++)
      {
         outputhihi[ix][i] = cffhihi[i];
         outputlohi[ix][i] = cfflohi[i];
         outputhilo[ix][i] = cffhilo[i];
         outputlolo[ix][i] = cfflolo[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_dbl4_product(deg,
         inputhihi[ix1],inputlohi[ix1],inputhilo[ix1],inputlolo[ix1],
         inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
         outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim]);
      CPU_dbl4_product(deg,
         outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim],
         cffhihi,cfflohi,cffhilo,cfflolo,
         outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim]);

      CPU_dbl4_product(deg,cffhihi,cfflohi,cffhilo,cfflolo,
         inputhihi[ix1],inputlohi[ix1],inputhilo[ix1],inputlolo[ix1],
         outputhihi[ix2],outputlohi[ix2],outputhilo[ix2],outputlolo[ix2]);
      CPU_dbl4_product(deg,cffhihi,cfflohi,cffhilo,cfflolo,
         inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
         outputhihi[ix1],outputlohi[ix1],outputhilo[ix1],outputlolo[ix1]);
   }
   else
   {
      double **forwardhihi = new double*[nvr];
      double **forwardlohi = new double*[nvr];
      double **forwardhilo = new double*[nvr];
      double **forwardlolo = new double*[nvr];
      double **backwardhihi = new double*[nvr-2];
      double **backwardlohi = new double*[nvr-2];
      double **backwardhilo = new double*[nvr-2];
      double **backwardlolo = new double*[nvr-2];
      double **crosshihi = new double*[nvr-2];
      double **crosslohi = new double*[nvr-2];
      double **crosshilo = new double*[nvr-2];
      double **crosslolo = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardhihi[i] = new double[deg+1];
         forwardlohi[i] = new double[deg+1];
         forwardhilo[i] = new double[deg+1];
         forwardlolo[i] = new double[deg+1];
         backwardhihi[i] = new double[deg+1];
         backwardlohi[i] = new double[deg+1];
         backwardhilo[i] = new double[deg+1];
         backwardlolo[i] = new double[deg+1];
         crosshihi[i] = new double[deg+1];
         crosslohi[i] = new double[deg+1];
         crosshilo[i] = new double[deg+1];
         crosslolo[i] = new double[deg+1];
      }
      forwardhihi[nvr-2] = new double[deg+1];
      forwardlohi[nvr-2] = new double[deg+1];
      forwardhilo[nvr-2] = new double[deg+1];
      forwardlolo[nvr-2] = new double[deg+1];
      forwardhihi[nvr-1] = new double[deg+1];
      forwardlohi[nvr-1] = new double[deg+1];
      forwardhilo[nvr-1] = new double[deg+1];
      forwardlolo[nvr-1] = new double[deg+1];

      CPU_dbl4_speel(nvr,deg,idx,cffhihi,cfflohi,cffhilo,cfflolo,
         inputhihi,inputlohi,inputhilo,inputlolo,
         forwardhihi,forwardlohi,forwardhilo,forwardlolo,
         backwardhihi,backwardlohi,backwardhilo,backwardlolo,
         crosshihi,crosslohi,crosshilo,crosslolo);

      // assign the value of the monomial

      for(int i=0; i<deg+1; i++)
      {
         outputhihi[dim][i] = forwardhihi[nvr-1][i];
         outputlohi[dim][i] = forwardlohi[nvr-1][i];
         outputhilo[dim][i] = forwardhilo[nvr-1][i];
         outputlolo[dim][i] = forwardlolo[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];       // derivative with respect to x[n-1]
         for(int i=0; i<deg+1; i++)
         {
            outputhihi[ix][i] = forwardhihi[nvr-2][i];
            outputlohi[ix][i] = forwardlohi[nvr-2][i];
            outputhilo[ix][i] = forwardhilo[nvr-2][i];
            outputlolo[ix][i] = forwardlolo[nvr-2][i];
         }
         ix = idx[0];               // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputhihi[ix][i] = backwardhihi[nvr-3][i];
            outputlohi[ix][i] = backwardlohi[nvr-3][i];
            outputhilo[ix][i] = backwardhilo[nvr-3][i];
            outputlolo[ix][i] = backwardlolo[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {                          // derivative with respect to x[k]
            ix = idx[k];
            for(int i=0; i<deg+1; i++)
            {
               outputhihi[ix][i] = crosshihi[k-1][i];
               outputlohi[ix][i] = crosslohi[k-1][i];
               outputhilo[ix][i] = crosshilo[k-1][i];
               outputlolo[ix][i] = crosslolo[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardhihi[i]); free(forwardlohi[i]);
         free(forwardhilo[i]); free(forwardlolo[i]);
         free(backwardhihi[i]); free(backwardlohi[i]);
         free(backwardhilo[i]); free(backwardlolo[i]);
         free(crosshihi[i]); free(crosslohi[i]);
         free(crosshilo[i]); free(crosslolo[i]);
      }
      free(forwardhihi[nvr-2]); free(forwardlohi[nvr-2]);
      free(forwardhilo[nvr-2]); free(forwardlolo[nvr-2]);
      free(forwardhihi[nvr-1]); free(forwardlohi[nvr-1]);
      free(forwardhilo[nvr-1]); free(forwardlolo[nvr-1]);
      free(forwardhihi); free(forwardlohi);
      free(forwardhilo); free(forwardlolo);
      free(backwardhihi); free(backwardlohi);
      free(backwardhilo); free(backwardlolo);
      free(crosshihi); free(crosslohi);
      free(crosshilo); free(crosslolo);
   }
}

void CPU_cmplx4_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo,
   double **inputrehi, double **inputrelo, double **inputimhi,
   double **inputimlo, double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo )
{
}

