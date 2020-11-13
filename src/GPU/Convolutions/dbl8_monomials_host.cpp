/* The file dbl8_monomials_host.cpp defines the functions specified
 * in dbl8_monomials_host.h. */

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
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"

void CPU_dbl8_speel
 ( int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi, double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo, double *cffhilolo, double *cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **forwardhihihi, double **forwardlohihi,
   double **forwardhilohi, double **forwardlolohi,
   double **forwardhihilo, double **forwardlohilo,
   double **forwardhilolo, double **forwardlololo,
   double **backwardhihihi, double **backwardlohihi,
   double **backwardhilohi, double **backwardlolohi,
   double **backwardhihilo, double **backwardlohilo,
   double **backwardhilolo, double **backwardlololo,
   double **crosshihihi, double **crosslohihi,
   double **crosshilohi, double **crosslolohi,
   double **crosshihilo, double **crosslohilo,
   double **crosshilolo, double **crosslololo )
{
   int ix1 = idx[0];
   int ix2;                                            // f[0] = cff*x[0] 

   CPU_dbl8_product(deg,
      cffhihihi,cfflohihi,cffhilohi,cfflolohi,
      cffhihilo,cfflohilo,cffhilolo,cfflololo,
      inputhihihi[ix1],inputlohihi[ix1],
      inputhilohi[ix1],inputlolohi[ix1],
      inputhihilo[ix1],inputlohilo[ix1],
      inputhilolo[ix1],inputlololo[ix1],
      forwardhihihi[0],forwardlohihi[0],
      forwardhilohi[0],forwardlolohi[0],
      forwardhihilo[0],forwardlohilo[0],
      forwardhilolo[0],forwardlololo[0]); 

   for(int i=1; i<nvr; i++)
   {                                                // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_dbl8_product(deg,
         forwardhihihi[i-1],forwardlohihi[i-1],
         forwardhilohi[i-1],forwardlolohi[i-1],
         forwardhihilo[i-1],forwardlohilo[i-1],
         forwardhilolo[i-1],forwardlololo[i-1],
         inputhihihi[ix2],inputlohihi[ix2],
         inputhilohi[ix2],inputlolohi[ix2],
         inputhihilo[ix2],inputlohilo[ix2],
         inputhilolo[ix2],inputlololo[ix2],
         forwardhihihi[i],forwardlohihi[i],
         forwardhilohi[i],forwardlolohi[i],
         forwardhihilo[i],forwardlohilo[i],
         forwardhilolo[i],forwardlololo[i]);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1];                            // b[0] = x[n-1]*x[n-2]
      ix2 = idx[nvr-2];
      CPU_dbl8_product(deg,
         inputhihihi[ix1],inputlohihi[ix1],
         inputhilohi[ix1],inputlolohi[ix1],
         inputhihilo[ix1],inputlohilo[ix1],
         inputhilolo[ix1],inputlololo[ix1],
         inputhihihi[ix2],inputlohihi[ix2],
         inputhilohi[ix2],inputlolohi[ix2],
         inputhihilo[ix2],inputlohilo[ix2],
         inputhilolo[ix2],inputlololo[ix2],
         backwardhihihi[0],backwardlohihi[0],
         backwardhilohi[0],backwardlolohi[0],
         backwardhihilo[0],backwardlohilo[0],
         backwardhilolo[0],backwardlololo[0]);

      for(int i=1; i<nvr-2; i++)
      {                                            // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];
         CPU_dbl8_product(deg,
            backwardhihihi[i-1],backwardlohihi[i-1],
            backwardhilohi[i-1],backwardlolohi[i-1],
            backwardhihilo[i-1],backwardlohilo[i-1],
            backwardhilolo[i-1],backwardlololo[i-1],
            inputhihihi[ix2],inputlohihi[ix2],
            inputhilohi[ix2],inputlolohi[ix2],
            inputhihilo[ix2],inputlohilo[ix2],
            inputhilolo[ix2],inputlololo[ix2],
            backwardhihihi[i],backwardlohihi[i],
            backwardhilohi[i],backwardlolohi[i],
            backwardhihilo[i],backwardlohilo[i],
            backwardhilolo[i],backwardlololo[i]);
      }
                                                   // b[n-3] = cff*b[n-3]
      CPU_dbl8_product(deg,
         backwardhihihi[nvr-3],backwardlohihi[nvr-3],
         backwardhilohi[nvr-3],backwardlolohi[nvr-3],
         backwardhihilo[nvr-3],backwardlohilo[nvr-3],
         backwardhilolo[nvr-3],backwardlololo[nvr-3],
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         crosshihihi[0],crosslohihi[0],crosshilohi[0],crosslolohi[0],
         crosshihilo[0],crosslohilo[0],crosshilolo[0],crosslololo[0]);
      // cross[0] is work space, cannot write into backward[nvr-3]
      for(int i=0; i<=deg; i++)
      {
         backwardhihihi[nvr-3][i] = crosshihihi[0][i];
         backwardlohihi[nvr-3][i] = crosslohihi[0][i];
         backwardhilohi[nvr-3][i] = crosshilohi[0][i];
         backwardlolohi[nvr-3][i] = crosslolohi[0][i];
         backwardhihilo[nvr-3][i] = crosshihilo[0][i];
         backwardlohilo[nvr-3][i] = crosslohilo[0][i];
         backwardhilolo[nvr-3][i] = crosshilolo[0][i];
         backwardlololo[nvr-3][i] = crosslololo[0][i];
      }
      if(nvr == 3)
      {                                            // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_dbl8_product(deg,
            forwardhihihi[0],forwardlohihi[0],
            forwardhilohi[0],forwardlolohi[0],
            forwardhihilo[0],forwardlohilo[0],
            forwardhilolo[0],forwardlololo[0],
            inputhihihi[ix2],inputlohihi[ix2],
            inputhilohi[ix2],inputlolohi[ix2],
            inputhihilo[ix2],inputlohilo[ix2],
            inputhilolo[ix2],inputlololo[ix2],
            crosshihihi[0],crosslohihi[0],crosshilohi[0],crosslolohi[0],
            crosshihilo[0],crosslohilo[0],crosshilolo[0],crosslololo[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                         // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_dbl8_product(deg,
               forwardhihihi[i],forwardlohihi[i],
               forwardhilohi[i],forwardlolohi[i],
               forwardhihilo[i],forwardlohilo[i],
               forwardhilolo[i],forwardlololo[i],
               backwardhihihi[ix2],backwardlohihi[ix2],
               backwardhilohi[ix2],backwardlolohi[ix2],
               backwardhihilo[ix2],backwardlohilo[ix2],
               backwardhilolo[ix2],backwardlololo[ix2],
               crosshihihi[i],crosslohihi[i],
               crosshilohi[i],crosslolohi[i],
               crosshihilo[i],crosslohilo[i],
               crosshilolo[i],crosslololo[i]);
         }
                                                   // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_dbl8_product(deg,
            forwardhihihi[nvr-3],forwardlohihi[nvr-3],
            forwardhilohi[nvr-3],forwardlolohi[nvr-3],
            forwardhihilo[nvr-3],forwardlohilo[nvr-3],
            forwardhilolo[nvr-3],forwardlololo[nvr-3],
            inputhihihi[ix2],inputlohihi[ix2],
            inputhilohi[ix2],inputlolohi[ix2],
            inputhihilo[ix2],inputlohilo[ix2],
            inputhilolo[ix2],inputlololo[ix2],
            crosshihihi[nvr-3],crosslohihi[nvr-3],
            crosshilohi[nvr-3],crosslolohi[nvr-3],
            crosshihilo[nvr-3],crosslohilo[nvr-3],
            crosshilolo[nvr-3],crosslololo[nvr-3]);
      }
   }
}

void CPU_cmplx8_speel
 ( int nvr, int deg, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **backwardrehihihi, double **backwardrelohihi,
   double **backwardrehilohi, double **backwardrelolohi,
   double **backwardrehihilo, double **backwardrelohilo,
   double **backwardrehilolo, double **backwardrelololo,
   double **backwardimhihihi, double **backwardimlohihi,
   double **backwardimhilohi, double **backwardimlolohi,
   double **backwardimhihilo, double **backwardimlohilo,
   double **backwardimhilolo, double **backwardimlololo,
   double **crossrehihihi, double **crossrelohihi,
   double **crossrehilohi, double **crossrelolohi,
   double **crossrehihilo, double **crossrelohilo,
   double **crossrehilolo, double **crossrelololo,
   double **crossimhihihi, double **crossimlohihi,
   double **crossimhilohi, double **crossimlolohi,
   double **crossimhihilo, double **crossimlohilo,
   double **crossimhilolo, double **crossimlololo )
{
}

void CPU_dbl8_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi, double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo, double *cffhilolo, double *cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_dbl8_product(deg,
         inputhihihi[ix],inputlohihi[ix],inputhilohi[ix],inputlolohi[ix],
         inputhihilo[ix],inputlohilo[ix],inputhilolo[ix],inputlololo[ix],
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         outputhihihi[dim],outputlohihi[dim],
         outputhilohi[dim],outputlolohi[dim],
         outputhihilo[dim],outputlohilo[dim],
         outputhilolo[dim],outputlololo[dim]);

      for(int i=0; i<=deg; i++)
      {
         outputhihihi[ix][i] = cffhihihi[i];
         outputlohihi[ix][i] = cfflohihi[i];
         outputhilohi[ix][i] = cffhilohi[i];
         outputlolohi[ix][i] = cfflolohi[i];
         outputhihilo[ix][i] = cffhihilo[i];
         outputlohilo[ix][i] = cfflohilo[i];
         outputhilolo[ix][i] = cffhilolo[i];
         outputlololo[ix][i] = cfflololo[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_dbl8_product(deg,
         inputhihihi[ix1],inputlohihi[ix1],inputhilohi[ix1],inputlolohi[ix1],
         inputhihilo[ix1],inputlohilo[ix1],inputhilolo[ix1],inputlololo[ix1],
         inputhihihi[ix2],inputlohihi[ix2],inputhilohi[ix2],inputlolohi[ix2],
         inputhihilo[ix2],inputlohilo[ix2],inputhilolo[ix2],inputlololo[ix2],
         outputhihihi[dim],outputlohihi[dim],
         outputhilohi[dim],outputlolohi[dim],
         outputhihilo[dim],outputlohilo[dim],
         outputhilolo[dim],outputlololo[dim]);
      CPU_dbl8_product(deg,
         outputhihihi[dim],outputlohihi[dim],
         outputhilohi[dim],outputlolohi[dim],
         outputhihilo[dim],outputlohilo[dim],
         outputhilolo[dim],outputlololo[dim],
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         outputhihihi[dim],outputlohihi[dim],
         outputhilohi[dim],outputlolohi[dim],
         outputhihilo[dim],outputlohilo[dim],
         outputhilolo[dim],outputlololo[dim]);

      CPU_dbl8_product(deg,
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         inputhihihi[ix1],inputlohihi[ix1],
         inputhilohi[ix1],inputlolohi[ix1],
         inputhihilo[ix1],inputlohilo[ix1],
         inputhilolo[ix1],inputlololo[ix1],
         outputhihihi[ix2],outputlohihi[ix2],
         outputhilohi[ix2],outputlolohi[ix2],
         outputhihilo[ix2],outputlohilo[ix2],
         outputhilolo[ix2],outputlololo[ix2]);
      CPU_dbl8_product(deg,
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         inputhihihi[ix2],inputlohihi[ix2],
         inputhilohi[ix2],inputlolohi[ix2],
         inputhihilo[ix2],inputlohilo[ix2],
         inputhilolo[ix2],inputlololo[ix2],
         outputhihihi[ix1],outputlohihi[ix1],
         outputhilohi[ix1],outputlolohi[ix1],
         outputhihilo[ix1],outputlohilo[ix1],
         outputhilolo[ix1],outputlololo[ix1]);
   }
   else
   {
      double **forwardhihihi = new double*[nvr];
      double **forwardlohihi = new double*[nvr];
      double **forwardhilohi = new double*[nvr];
      double **forwardlolohi = new double*[nvr];
      double **forwardhihilo = new double*[nvr];
      double **forwardlohilo = new double*[nvr];
      double **forwardhilolo = new double*[nvr];
      double **forwardlololo = new double*[nvr];
      double **backwardhihihi = new double*[nvr-2];
      double **backwardlohihi = new double*[nvr-2];
      double **backwardhilohi = new double*[nvr-2];
      double **backwardlolohi = new double*[nvr-2];
      double **backwardhihilo = new double*[nvr-2];
      double **backwardlohilo = new double*[nvr-2];
      double **backwardhilolo = new double*[nvr-2];
      double **backwardlololo = new double*[nvr-2];
      double **crosshihihi = new double*[nvr-2];
      double **crosslohihi = new double*[nvr-2];
      double **crosshilohi = new double*[nvr-2];
      double **crosslolohi = new double*[nvr-2];
      double **crosshihilo = new double*[nvr-2];
      double **crosslohilo = new double*[nvr-2];
      double **crosshilolo = new double*[nvr-2];
      double **crosslololo = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardhihihi[i] = new double[deg+1];
         forwardlohihi[i] = new double[deg+1];
         forwardhilohi[i] = new double[deg+1];
         forwardlolohi[i] = new double[deg+1];
         forwardhihilo[i] = new double[deg+1];
         forwardlohilo[i] = new double[deg+1];
         forwardhilolo[i] = new double[deg+1];
         forwardlololo[i] = new double[deg+1];
         backwardhihihi[i] = new double[deg+1];
         backwardlohihi[i] = new double[deg+1];
         backwardhilohi[i] = new double[deg+1];
         backwardlolohi[i] = new double[deg+1];
         backwardhihilo[i] = new double[deg+1];
         backwardlohilo[i] = new double[deg+1];
         backwardhilolo[i] = new double[deg+1];
         backwardlololo[i] = new double[deg+1];
         crosshihihi[i] = new double[deg+1];
         crosslohihi[i] = new double[deg+1];
         crosshilohi[i] = new double[deg+1];
         crosslolohi[i] = new double[deg+1];
         crosshihilo[i] = new double[deg+1];
         crosslohilo[i] = new double[deg+1];
         crosshilolo[i] = new double[deg+1];
         crosslololo[i] = new double[deg+1];
      }
      forwardhihihi[nvr-2] = new double[deg+1];
      forwardlohihi[nvr-2] = new double[deg+1];
      forwardhilohi[nvr-2] = new double[deg+1];
      forwardlolohi[nvr-2] = new double[deg+1];
      forwardhihilo[nvr-2] = new double[deg+1];
      forwardlohilo[nvr-2] = new double[deg+1];
      forwardhilolo[nvr-2] = new double[deg+1];
      forwardlololo[nvr-2] = new double[deg+1];
      forwardhihihi[nvr-1] = new double[deg+1];
      forwardlohihi[nvr-1] = new double[deg+1];
      forwardhilohi[nvr-1] = new double[deg+1];
      forwardlolohi[nvr-1] = new double[deg+1];
      forwardhihilo[nvr-1] = new double[deg+1];
      forwardlohilo[nvr-1] = new double[deg+1];
      forwardhilolo[nvr-1] = new double[deg+1];
      forwardlololo[nvr-1] = new double[deg+1];

      CPU_dbl8_speel(nvr,deg,idx,
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         inputhihihi,inputlohihi,inputhilohi,inputlolohi,
         inputhihilo,inputlohilo,inputhilolo,inputlololo,
         forwardhihihi,forwardlohihi,forwardhilohi,forwardlolohi,
         forwardhihilo,forwardlohilo,forwardhilolo,forwardlololo,
         backwardhihihi,backwardlohihi,backwardhilohi,backwardlolohi,
         backwardhihilo,backwardlohilo,backwardhilolo,backwardlololo,
         crosshihihi,crosslohihi,crosshilohi,crosslolohi,
         crosshihilo,crosslohilo,crosshilolo,crosslololo);

      // assign the value of the monomial

      for(int i=0; i<deg+1; i++)
      {
         outputhihihi[dim][i] = forwardhihihi[nvr-1][i];
         outputlohihi[dim][i] = forwardlohihi[nvr-1][i];
         outputhilohi[dim][i] = forwardhilohi[nvr-1][i];
         outputlolohi[dim][i] = forwardlolohi[nvr-1][i];
         outputhihilo[dim][i] = forwardhihilo[nvr-1][i];
         outputlohilo[dim][i] = forwardlohilo[nvr-1][i];
         outputhilolo[dim][i] = forwardhilolo[nvr-1][i];
         outputlololo[dim][i] = forwardlololo[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];       // derivative with respect to x[n-1]
         for(int i=0; i<deg+1; i++)
         {
            outputhihihi[ix][i] = forwardhihihi[nvr-2][i];
            outputlohihi[ix][i] = forwardlohihi[nvr-2][i];
            outputhilohi[ix][i] = forwardhilohi[nvr-2][i];
            outputlolohi[ix][i] = forwardlolohi[nvr-2][i];
            outputhihilo[ix][i] = forwardhihilo[nvr-2][i];
            outputlohilo[ix][i] = forwardlohilo[nvr-2][i];
            outputhilolo[ix][i] = forwardhilolo[nvr-2][i];
            outputlololo[ix][i] = forwardlololo[nvr-2][i];
         }
         ix = idx[0];               // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputhihihi[ix][i] = backwardhihihi[nvr-3][i];
            outputlohihi[ix][i] = backwardlohihi[nvr-3][i];
            outputhilohi[ix][i] = backwardhilohi[nvr-3][i];
            outputlolohi[ix][i] = backwardlolohi[nvr-3][i];
            outputhihilo[ix][i] = backwardhihilo[nvr-3][i];
            outputlohilo[ix][i] = backwardlohilo[nvr-3][i];
            outputhilolo[ix][i] = backwardhilolo[nvr-3][i];
            outputlololo[ix][i] = backwardlololo[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {                          // derivative with respect to x[k]
            ix = idx[k];
            for(int i=0; i<deg+1; i++)
            {
               outputhihihi[ix][i] = crosshihihi[k-1][i];
               outputlohihi[ix][i] = crosslohihi[k-1][i];
               outputhilohi[ix][i] = crosshilohi[k-1][i];
               outputlolohi[ix][i] = crosslolohi[k-1][i];
               outputhihilo[ix][i] = crosshihilo[k-1][i];
               outputlohilo[ix][i] = crosslohilo[k-1][i];
               outputhilolo[ix][i] = crosshilolo[k-1][i];
               outputlololo[ix][i] = crosslololo[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardhihihi[i]); free(forwardlohihi[i]);
         free(forwardhilohi[i]); free(forwardlolohi[i]);
         free(forwardhihilo[i]); free(forwardlohilo[i]);
         free(forwardhilolo[i]); free(forwardlololo[i]);
         free(backwardhihihi[i]); free(backwardlohihi[i]);
         free(backwardhilohi[i]); free(backwardlolohi[i]);
         free(backwardhihilo[i]); free(backwardlohilo[i]);
         free(backwardhilolo[i]); free(backwardlololo[i]);
         free(crosshihihi[i]); free(crosslohihi[i]);
         free(crosshilohi[i]); free(crosslolohi[i]);
         free(crosshihilo[i]); free(crosslohilo[i]);
         free(crosshilolo[i]); free(crosslololo[i]);
      }
      free(forwardhihihi[nvr-2]); free(forwardlohihi[nvr-2]);
      free(forwardhilohi[nvr-2]); free(forwardlolohi[nvr-2]);
      free(forwardhihilo[nvr-2]); free(forwardlohilo[nvr-2]);
      free(forwardhilolo[nvr-2]); free(forwardlololo[nvr-2]);
      free(forwardhihihi[nvr-1]); free(forwardlohihi[nvr-1]);
      free(forwardhilohi[nvr-1]); free(forwardlolohi[nvr-1]);
      free(forwardhihilo[nvr-1]); free(forwardlohilo[nvr-1]);
      free(forwardhilolo[nvr-1]); free(forwardlololo[nvr-1]);
      free(forwardhihihi); free(forwardlohihi);
      free(forwardhilohi); free(forwardlolohi);
      free(forwardhihilo); free(forwardlohilo);
      free(forwardhilolo); free(forwardlololo);
      free(backwardhihihi); free(backwardlohihi);
      free(backwardhilohi); free(backwardlolohi);
      free(backwardhihilo); free(backwardlohilo);
      free(backwardhilolo); free(backwardlololo);
      free(crosshihihi); free(crosslohihi);
      free(crosshilohi); free(crosslolohi);
      free(crosshihilo); free(crosslohilo);
      free(crosshilolo); free(crosslololo);
   }
}

void CPU_cmplx8_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo )
{
}
