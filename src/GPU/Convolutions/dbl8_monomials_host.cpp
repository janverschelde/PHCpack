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
   double **forwardrehihihi, double **forwardrelohihi,
   double **forwardrehilohi, double **forwardrelolohi,
   double **forwardrehihilo, double **forwardrelohilo,
   double **forwardrehilolo, double **forwardrelololo,
   double **forwardimhihihi, double **forwardimlohihi,
   double **forwardimhilohi, double **forwardimlolohi,
   double **forwardimhihilo, double **forwardimlohilo,
   double **forwardimhilolo, double **forwardimlololo,
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
   int ix1 = idx[0];
   int ix2;
                                                           // f[0] = cff*x[0]
   CPU_cmplx8_product(deg,
      cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
      cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
      cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
      cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
      inputrehihihi[ix1],inputrelohihi[ix1],
      inputrehilohi[ix1],inputrelolohi[ix1],
      inputrehihilo[ix1],inputrelohilo[ix1],
      inputrehilolo[ix1],inputrelololo[ix1],
      inputimhihihi[ix1],inputimlohihi[ix1],
      inputimhilohi[ix1],inputimlolohi[ix1],
      inputimhihilo[ix1],inputimlohilo[ix1],
      inputimhilolo[ix1],inputimlololo[ix1],
      forwardrehihihi[0],forwardrelohihi[0],
      forwardrehilohi[0],forwardrelolohi[0],
      forwardrehihilo[0],forwardrelohilo[0],
      forwardrehilolo[0],forwardrelololo[0],
      forwardimhihihi[0],forwardimlohihi[0],
      forwardimhilohi[0],forwardimlolohi[0],
      forwardimhihilo[0],forwardimlohilo[0],
      forwardimhilolo[0],forwardimlololo[0]);

   for(int i=1; i<nvr; i++)
   {                                                    // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_cmplx8_product(deg,
         forwardrehihihi[i-1],forwardrelohihi[i-1],
         forwardrehilohi[i-1],forwardrelolohi[i-1],
         forwardrehihilo[i-1],forwardrelohilo[i-1],
         forwardrehilolo[i-1],forwardrelololo[i-1],
         forwardimhihihi[i-1],forwardimlohihi[i-1],
         forwardimhilohi[i-1],forwardimlolohi[i-1],
         forwardimhihilo[i-1],forwardimlohilo[i-1],
         forwardimhilolo[i-1],forwardimlololo[i-1],
         inputrehihihi[ix2],inputrelohihi[ix2],
         inputrehilohi[ix2],inputrelolohi[ix2],
         inputrehihilo[ix2],inputrelohilo[ix2],
         inputrehilolo[ix2],inputrelololo[ix2],
         inputimhihihi[ix2],inputimlohihi[ix2],
         inputimhilohi[ix2],inputimlolohi[ix2],
         inputimhihilo[ix2],inputimlohilo[ix2],
         inputimhilolo[ix2],inputimlololo[ix2],
         forwardrehihihi[i],forwardrelohihi[i],
         forwardrehilohi[i],forwardrelolohi[i],
         forwardrehihilo[i],forwardrelohilo[i],
         forwardrehilolo[i],forwardrelololo[i],
         forwardimhihihi[i],forwardimlohihi[i],
         forwardimhilohi[i],forwardimlolohi[i],
         forwardimhihilo[i],forwardimlohilo[i],
         forwardimhilolo[i],forwardimlololo[i]);
   }
   if(nvr > 2)
   {                                                  // b[0] = x[n-1]*x[n-2]
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];
      CPU_cmplx8_product(deg,
         inputrehihihi[ix1],inputrelohihi[ix1],
         inputrehilohi[ix1],inputrelolohi[ix1],
         inputrehihilo[ix1],inputrelohilo[ix1],
         inputrehilolo[ix1],inputrelololo[ix1],
         inputimhihihi[ix1],inputimlohihi[ix1],
         inputimhilohi[ix1],inputimlolohi[ix1],
         inputimhihilo[ix1],inputimlohilo[ix1],
         inputimhilolo[ix1],inputimlololo[ix1],
         inputrehihihi[ix2],inputrelohihi[ix2],
         inputrehilohi[ix2],inputrelolohi[ix2],
         inputrehihilo[ix2],inputrelohilo[ix2],
         inputrehilolo[ix2],inputrelololo[ix2],
         inputimhihihi[ix2],inputimlohihi[ix2],
         inputimhilohi[ix2],inputimlolohi[ix2],
         inputimhihilo[ix2],inputimlohilo[ix2],
         inputimhilolo[ix2],inputimlololo[ix2],
         backwardrehihihi[0],backwardrelohihi[0],
         backwardrehilohi[0],backwardrelolohi[0],
         backwardrehihilo[0],backwardrelohilo[0],
         backwardrehilolo[0],backwardrelololo[0],
         backwardimhihihi[0],backwardimlohihi[0],
         backwardimhilohi[0],backwardimlolohi[0],
         backwardimhihilo[0],backwardimlohilo[0],
         backwardimhilolo[0],backwardimlololo[0]);

      for(int i=1; i<nvr-2; i++)
      {                                             // b[i] = b[i-1]*x[x-2-i]
         ix2 = idx[nvr-2-i];
         CPU_cmplx8_product(deg,
            backwardrehihihi[i-1],backwardrelohihi[i-1],
            backwardrehilohi[i-1],backwardrelolohi[i-1],
            backwardrehihilo[i-1],backwardrelohilo[i-1],
            backwardrehilolo[i-1],backwardrelololo[i-1],
            backwardimhihihi[i-1],backwardimlohihi[i-1],
            backwardimhilohi[i-1],backwardimlolohi[i-1],
            backwardimhihilo[i-1],backwardimlohilo[i-1],
            backwardimhilolo[i-1],backwardimlololo[i-1],
            inputrehihihi[ix2],inputrelohihi[ix2],
            inputrehilohi[ix2],inputrelolohi[ix2],
            inputrehihilo[ix2],inputrelohilo[ix2],
            inputrehilolo[ix2],inputrelololo[ix2],
            inputimhihihi[ix2],inputimlohihi[ix2],
            inputimhilohi[ix2],inputimlolohi[ix2],
            inputimhihilo[ix2],inputimlohilo[ix2],
            inputimhilolo[ix2],inputimlololo[ix2],
            backwardrehihihi[i],backwardrelohihi[i],
            backwardrehilohi[i],backwardrelolohi[i],
            backwardrehihilo[i],backwardrelohilo[i],
            backwardrehilolo[i],backwardrelololo[i],
            backwardimhihihi[i],backwardimlohihi[i],
            backwardimhilohi[i],backwardimlolohi[i],
            backwardimhihilo[i],backwardimlohilo[i],
            backwardimhilolo[i],backwardimlololo[i]);
      }
                                                       // b[n-3] = b[n-3]*cff
      CPU_cmplx8_product(deg,
         backwardrehihihi[nvr-3],backwardrelohihi[nvr-3],
         backwardrehilohi[nvr-3],backwardrelolohi[nvr-3],
         backwardrehihilo[nvr-3],backwardrelohilo[nvr-3],
         backwardrehilolo[nvr-3],backwardrelololo[nvr-3],
         backwardimhihihi[nvr-3],backwardimlohihi[nvr-3],
         backwardimhilohi[nvr-3],backwardimlolohi[nvr-3],
         backwardimhihilo[nvr-3],backwardimlohilo[nvr-3],
         backwardimhilolo[nvr-3],backwardimlololo[nvr-3],
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         crossrehihihi[0],crossrelohihi[0],
         crossrehilohi[0],crossrelolohi[0],
         crossrehihilo[0],crossrelohilo[0],
         crossrehilolo[0],crossrelololo[0],
         crossimhihihi[0],crossimlohihi[0],
         crossimhilohi[0],crossimlolohi[0],
         crossimhihilo[0],crossimlohilo[0],
         crossimhilolo[0],crossimlololo[0]); 
                                                       // cross is work space
      for(int i=0; i<=deg; i++)
      {
         backwardrehihihi[nvr-3][i] = crossrehihihi[0][i];
         backwardrelohihi[nvr-3][i] = crossrelohihi[0][i];
         backwardrehilohi[nvr-3][i] = crossrehilohi[0][i];
         backwardrelolohi[nvr-3][i] = crossrelolohi[0][i];
         backwardrehihilo[nvr-3][i] = crossrehihilo[0][i];
         backwardrelohilo[nvr-3][i] = crossrelohilo[0][i];
         backwardrehilolo[nvr-3][i] = crossrehilolo[0][i];
         backwardrelololo[nvr-3][i] = crossrelololo[0][i];
         backwardimhihihi[nvr-3][i] = crossimhihihi[0][i];
         backwardimlohihi[nvr-3][i] = crossimlohihi[0][i];
         backwardimhilohi[nvr-3][i] = crossimhilohi[0][i];
         backwardimlolohi[nvr-3][i] = crossimlolohi[0][i];
         backwardimhihilo[nvr-3][i] = crossimhihilo[0][i];
         backwardimlohilo[nvr-3][i] = crossimlohilo[0][i];
         backwardimhilolo[nvr-3][i] = crossimhilolo[0][i];
         backwardimlololo[nvr-3][i] = crossimlololo[0][i];
      }
      if(nvr == 3)
      {                                                   // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_cmplx8_product(deg,
            forwardrehihihi[0],forwardrelohihi[0],
            forwardrehilohi[0],forwardrelolohi[0],
            forwardrehihilo[0],forwardrelohilo[0],
            forwardrehilolo[0],forwardrelololo[0],
            forwardimhihihi[0],forwardimlohihi[0],
            forwardimhilohi[0],forwardimlolohi[0],
            forwardimhihilo[0],forwardimlohilo[0],
            forwardimhilolo[0],forwardimlololo[0],
            inputrehihihi[ix2],inputrelohihi[ix2],
            inputrehilohi[ix2],inputrelolohi[ix2],
            inputrehihilo[ix2],inputrelohilo[ix2],
            inputrehilolo[ix2],inputrelololo[ix2],
            inputimhihihi[ix2],inputimlohihi[ix2],
            inputimhilohi[ix2],inputimlolohi[ix2],
            inputimhihilo[ix2],inputimlohilo[ix2],
            inputimhilolo[ix2],inputimlololo[ix2],
            crossrehihihi[0],crossrelohihi[0],
            crossrehilohi[0],crossrelolohi[0],
            crossrehihilo[0],crossrelohilo[0],
            crossrehilolo[0],crossrelololo[0],
            crossimhihihi[0],crossimlohihi[0],
            crossimhilohi[0],crossimlolohi[0],
            crossimhihilo[0],crossimlohilo[0],
            crossimhilolo[0],crossimlololo[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                            // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_cmplx8_product(deg,
               forwardrehihihi[i],forwardrelohihi[i],
               forwardrehilohi[i],forwardrelolohi[i],
               forwardrehihilo[i],forwardrelohilo[i],
               forwardrehilolo[i],forwardrelololo[i],
               forwardimhihihi[i],forwardimlohihi[i],
               forwardimhilohi[i],forwardimlolohi[i],
               forwardimhihilo[i],forwardimlohilo[i],
               forwardimhilolo[i],forwardimlololo[i],
               backwardrehihihi[ix2],backwardrelohihi[ix2],
               backwardrehilohi[ix2],backwardrelolohi[ix2],
               backwardrehihilo[ix2],backwardrelohilo[ix2],
               backwardrehilolo[ix2],backwardrelololo[ix2],
               backwardimhihihi[ix2],backwardimlohihi[ix2],
               backwardimhilohi[ix2],backwardimlolohi[ix2],
               backwardimhihilo[ix2],backwardimlohilo[ix2],
               backwardimhilolo[ix2],backwardimlololo[ix2],
               crossrehihihi[i],crossrelohihi[i],
               crossrehilohi[i],crossrelolohi[i],
               crossrehihilo[i],crossrelohilo[i],
               crossrehilolo[i],crossrelololo[i],
               crossimhihihi[i],crossimlohihi[i],
               crossimhilohi[i],crossimlolohi[i],
               crossimhihilo[i],crossimlohilo[i],
               crossimhilolo[i],crossimlololo[i]);
         }
                                                    // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_cmplx8_product(deg,
            forwardrehihihi[nvr-3],forwardrelohihi[nvr-3],
            forwardrehilohi[nvr-3],forwardrelolohi[nvr-3],
            forwardrehihilo[nvr-3],forwardrelohilo[nvr-3],
            forwardrehilolo[nvr-3],forwardrelololo[nvr-3],
            forwardimhihihi[nvr-3],forwardimlohihi[nvr-3],
            forwardimhilohi[nvr-3],forwardimlolohi[nvr-3],
            forwardimhihilo[nvr-3],forwardimlohilo[nvr-3],
            forwardimhilolo[nvr-3],forwardimlololo[nvr-3],
            inputrehihihi[ix2],inputrelohihi[ix2],
            inputrehilohi[ix2],inputrelolohi[ix2],
            inputrehihilo[ix2],inputrelohilo[ix2],
            inputrehilolo[ix2],inputrelololo[ix2],
            inputimhihihi[ix2],inputimlohihi[ix2],
            inputimhilohi[ix2],inputimlolohi[ix2],
            inputimhihilo[ix2],inputimlohilo[ix2],
            inputimhilolo[ix2],inputimlololo[ix2],
            crossrehihihi[nvr-3],crossrelohihi[nvr-3],
            crossrehilohi[nvr-3],crossrelolohi[nvr-3],
            crossrehihilo[nvr-3],crossrelohilo[nvr-3],
            crossrehilolo[nvr-3],crossrelololo[nvr-3],
            crossimhihihi[nvr-3],crossimlohihi[nvr-3],
            crossimhilohi[nvr-3],crossimlolohi[nvr-3],
            crossimhihilo[nvr-3],crossimlohilo[nvr-3],
            crossimhilolo[nvr-3],crossimlololo[nvr-3]);
      }
   }
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
      // CPU_dbl8_product(deg,
      //    outputhihihi[dim],outputlohihi[dim],
      //    outputhilohi[dim],outputlolohi[dim],
      //    outputhihilo[dim],outputlohilo[dim],
      //    outputhilolo[dim],outputlololo[dim],
      //    cffhihihi,cfflohihi,cffhilohi,cfflolohi,
      //    cffhihilo,cfflohilo,cffhilolo,cfflololo,
      //    outputhihihi[dim],outputlohihi[dim],
      //    outputhilohi[dim],outputlolohi[dim],
      //    outputhihilo[dim],outputlohilo[dim],
      //    outputhilolo[dim],outputlololo[dim]); // wrong!

      double *acchihihi = new double[deg+1];
      double *acclohihi = new double[deg+1];
      double *acchilohi = new double[deg+1];
      double *acclolohi = new double[deg+1];
      double *acchihilo = new double[deg+1];
      double *acclohilo = new double[deg+1];
      double *acchilolo = new double[deg+1];
      double *acclololo = new double[deg+1];

      CPU_dbl8_product(deg,
         outputhihihi[dim],outputlohihi[dim],
         outputhilohi[dim],outputlolohi[dim],
         outputhihilo[dim],outputlohilo[dim],
         outputhilolo[dim],outputlololo[dim],
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         acchihihi,acclohihi,acchilohi,acclolohi,
         acchihilo,acclohilo,acchilolo,acclololo);

      for(int i=0; i<=deg; i++)
      {
         outputhihihi[dim][i] = acchihihi[i];
         outputlohihi[dim][i] = acclohihi[i];
         outputhilohi[dim][i] = acchilohi[i];
         outputlolohi[dim][i] = acclolohi[i];
         outputhihilo[dim][i] = acchihilo[i];
         outputlohilo[dim][i] = acclohilo[i];
         outputhilolo[dim][i] = acchilolo[i];
         outputlololo[dim][i] = acclololo[i];
      }
      free(acchihihi); free(acclohihi);
      free(acchilohi); free(acclolohi);
      free(acchihilo); free(acclohilo);
      free(acchilolo); free(acclololo);

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
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_cmplx8_product(deg,
         inputrehihihi[ix],inputrelohihi[ix],
         inputrehilohi[ix],inputrelolohi[ix],
         inputrehihilo[ix],inputrelohilo[ix],
         inputrehilolo[ix],inputrelololo[ix],
         inputimhihihi[ix],inputimlohihi[ix],
         inputimhilohi[ix],inputimlolohi[ix],
         inputimhihilo[ix],inputimlohilo[ix],
         inputimhilolo[ix],inputimlololo[ix],
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         outputrehihihi[dim],outputrelohihi[dim],
         outputrehilohi[dim],outputrelolohi[dim],
         outputrehihilo[dim],outputrelohilo[dim],
         outputrehilolo[dim],outputrelololo[dim],
         outputimhihihi[dim],outputimlohihi[dim],
         outputimhilohi[dim],outputimlolohi[dim],
         outputimhihilo[dim],outputimlohilo[dim],
         outputimhilolo[dim],outputimlololo[dim]);

      for(int i=0; i<=deg; i++) 
      {
         outputrehihihi[ix][i] = cffrehihihi[i];
         outputrelohihi[ix][i] = cffrelohihi[i];
         outputrehilolo[ix][i] = cffrehilolo[i];
         outputrelololo[ix][i] = cffrelololo[i];
         outputrehihihi[ix][i] = cffrehihihi[i];
         outputrelohihi[ix][i] = cffrelohihi[i];
         outputrehilolo[ix][i] = cffrehilolo[i];
         outputrelololo[ix][i] = cffrelololo[i];
         outputimhihihi[ix][i] = cffimhihihi[i];
         outputimlohihi[ix][i] = cffimlohihi[i];
         outputimhilolo[ix][i] = cffimhilolo[i];
         outputimlololo[ix][i] = cffimlololo[i];
         outputimhihihi[ix][i] = cffimhihihi[i];
         outputimlohihi[ix][i] = cffimlohihi[i];
         outputimhilolo[ix][i] = cffimhilolo[i];
         outputimlololo[ix][i] = cffimlololo[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0];
      int ix2 = idx[1];

      CPU_cmplx8_product(deg,
         inputrehihihi[ix1],inputrelohihi[ix1],
         inputrehilohi[ix1],inputrelolohi[ix1],
         inputrehihilo[ix1],inputrelohilo[ix1],
         inputrehilolo[ix1],inputrelololo[ix1],
         inputimhihihi[ix1],inputimlohihi[ix1],
         inputimhilohi[ix1],inputimlolohi[ix1],
         inputimhihilo[ix1],inputimlohilo[ix1],
         inputimhilolo[ix1],inputimlololo[ix1],
         inputrehihihi[ix2],inputrelohihi[ix2],
         inputrehilohi[ix2],inputrelolohi[ix2],
         inputrehihilo[ix2],inputrelohilo[ix2],
         inputrehilolo[ix2],inputrelololo[ix2],
         inputimhihihi[ix2],inputimlohihi[ix2],
         inputimhilohi[ix2],inputimlolohi[ix2],
         inputimhihilo[ix2],inputimlohilo[ix2],
         inputimhilolo[ix2],inputimlololo[ix2],
         outputrehihihi[dim],outputrelohihi[dim],
         outputrehilohi[dim],outputrelolohi[dim],
         outputrehihilo[dim],outputrelohilo[dim],
         outputrehilolo[dim],outputrelololo[dim],
         outputimhihihi[dim],outputimlohihi[dim],
         outputimhilohi[dim],outputimlolohi[dim],
         outputimhihilo[dim],outputimlohilo[dim],
         outputimhilolo[dim],outputimlololo[dim]);
 /*
      CPU_cmplx8_product(deg,
         outputrehihihi[dim],outputrelohihi[dim],
         outputrehilohi[dim],outputrelolohi[dim],
         outputrehihilo[dim],outputrelohilo[dim],
         outputrehilolo[dim],outputrelololo[dim],
         outputimhihihi[dim],outputimlohihi[dim],
         outputimhilohi[dim],outputimlolohi[dim],
         outputimhihilo[dim],outputimlohilo[dim],
         outputimhilolo[dim],outputimlololo[dim],
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         outputrehihihi[dim],outputrelohihi[dim],
         outputrehilohi[dim],outputrelolohi[dim],
         outputrehihilo[dim],outputrelohilo[dim],
         outputrehilolo[dim],outputrelololo[dim],
         outputimhihihi[dim],outputimlohihi[dim],
         outputimhilohi[dim],outputimlolohi[dim],
         outputimhihilo[dim],outputimlohilo[dim],
         outputimhilolo[dim],outputimlololo[dim]);
  */

      double *accrehihihi = new double[deg+1];
      double *accrelohihi = new double[deg+1];
      double *accrehilohi = new double[deg+1];
      double *accrelolohi = new double[deg+1];
      double *accrehihilo = new double[deg+1];
      double *accrelohilo = new double[deg+1];
      double *accrehilolo = new double[deg+1];
      double *accrelololo = new double[deg+1];
      double *accimhihihi = new double[deg+1];
      double *accimlohihi = new double[deg+1];
      double *accimhilohi = new double[deg+1];
      double *accimlolohi = new double[deg+1];
      double *accimhihilo = new double[deg+1];
      double *accimlohilo = new double[deg+1];
      double *accimhilolo = new double[deg+1];
      double *accimlololo = new double[deg+1];

      CPU_cmplx8_product(deg,
         outputrehihihi[dim],outputrelohihi[dim],
         outputrehilohi[dim],outputrelolohi[dim],
         outputrehihilo[dim],outputrelohilo[dim],
         outputrehilolo[dim],outputrelololo[dim],
         outputimhihihi[dim],outputimlohihi[dim],
         outputimhilohi[dim],outputimlolohi[dim],
         outputimhihilo[dim],outputimlohilo[dim],
         outputimhilolo[dim],outputimlololo[dim],
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         accrehihihi,accrelohihi,accrehilohi,accrelolohi,
         accrehihilo,accrelohilo,accrehilolo,accrelololo,
         accimhihihi,accimlohihi,accimhilohi,accimlolohi,
         accimhihilo,accimlohilo,accimhilolo,accimlololo);

      for(int i=0; i<=deg; i++)
      {
         outputrehihihi[dim][i] = accrehihihi[i];
         outputrelohihi[dim][i] = accrelohihi[i];
         outputrehilohi[dim][i] = accrehilohi[i];
         outputrelolohi[dim][i] = accrelolohi[i];
         outputrehihilo[dim][i] = accrehihilo[i];
         outputrelohilo[dim][i] = accrelohilo[i];
         outputrehilolo[dim][i] = accrehilolo[i];
         outputrelololo[dim][i] = accrelololo[i];
         outputimhihihi[dim][i] = accimhihihi[i];
         outputimlohihi[dim][i] = accimlohihi[i];
         outputimhilohi[dim][i] = accimhilohi[i];
         outputimlolohi[dim][i] = accimlolohi[i];
         outputimhihilo[dim][i] = accimhihilo[i];
         outputimlohilo[dim][i] = accimlohilo[i];
         outputimhilolo[dim][i] = accimhilolo[i];
         outputimlololo[dim][i] = accimlololo[i];
      }
      free(accrehihihi); free(accrelohihi); 
      free(accrehilohi); free(accrelolohi);
      free(accrehihilo); free(accrelohilo); 
      free(accrehilolo); free(accrelololo);
      free(accimhihihi); free(accimlohihi);
      free(accimhilohi); free(accimlolohi);
      free(accimhihilo); free(accimlohilo);
      free(accimhilolo); free(accimlololo);

      CPU_cmplx8_product(deg,
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         inputrehihihi[ix1],inputrelohihi[ix1],
         inputrehilohi[ix1],inputrelolohi[ix1],
         inputrehihilo[ix1],inputrelohilo[ix1],
         inputrehilolo[ix1],inputrelololo[ix1],
         inputimhihihi[ix1],inputimlohihi[ix1],
         inputimhilohi[ix1],inputimlolohi[ix1],
         inputimhihilo[ix1],inputimlohilo[ix1],
         inputimhilolo[ix1],inputimlololo[ix1],
         outputrehihihi[ix2],outputrelohihi[ix2],
         outputrehilohi[ix2],outputrelolohi[ix2],
         outputrehihilo[ix2],outputrelohilo[ix2],
         outputrehilolo[ix2],outputrelololo[ix2],
         outputimhihihi[ix2],outputimlohihi[ix2],
         outputimhilohi[ix2],outputimlolohi[ix2],
         outputimhihilo[ix2],outputimlohilo[ix2],
         outputimhilolo[ix2],outputimlololo[ix2]);
      CPU_cmplx8_product(deg,
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         inputrehihihi[ix2],inputrelohihi[ix2],
         inputrehilohi[ix2],inputrelolohi[ix2],
         inputrehihilo[ix2],inputrelohilo[ix2],
         inputrehilolo[ix2],inputrelololo[ix2],
         inputimhihihi[ix2],inputimlohihi[ix2],
         inputimhilohi[ix2],inputimlolohi[ix2],
         inputimhihilo[ix2],inputimlohilo[ix2],
         inputimhilolo[ix2],inputimlololo[ix2],
         outputrehihihi[ix1],outputrelohihi[ix1],
         outputrehilohi[ix1],outputrelolohi[ix1],
         outputrehihilo[ix1],outputrelohilo[ix1],
         outputrehilolo[ix1],outputrelololo[ix1],
         outputimhihihi[ix1],outputimlohihi[ix1],
         outputimhilohi[ix1],outputimlolohi[ix1],
         outputimhihilo[ix1],outputimlohilo[ix1],
         outputimhilolo[ix1],outputimlololo[ix1]);
   }
   else
   {
      double **forwardrehihihi = new double*[nvr];
      double **forwardrelohihi = new double*[nvr];
      double **forwardrehilohi = new double*[nvr];
      double **forwardrelolohi = new double*[nvr];
      double **forwardrehihilo = new double*[nvr];
      double **forwardrelohilo = new double*[nvr];
      double **forwardrehilolo = new double*[nvr];
      double **forwardrelololo = new double*[nvr];
      double **forwardimhihihi = new double*[nvr];
      double **forwardimlohihi = new double*[nvr];
      double **forwardimhilohi = new double*[nvr];
      double **forwardimlolohi = new double*[nvr];
      double **forwardimhihilo = new double*[nvr];
      double **forwardimlohilo = new double*[nvr];
      double **forwardimhilolo = new double*[nvr];
      double **forwardimlololo = new double*[nvr];
      double **backwardrehihihi = new double*[nvr-2];
      double **backwardrelohihi = new double*[nvr-2];
      double **backwardrehilohi = new double*[nvr-2];
      double **backwardrelolohi = new double*[nvr-2];
      double **backwardrehihilo = new double*[nvr-2];
      double **backwardrelohilo = new double*[nvr-2];
      double **backwardrehilolo = new double*[nvr-2];
      double **backwardrelololo = new double*[nvr-2];
      double **backwardimhihihi = new double*[nvr-2];
      double **backwardimlohihi = new double*[nvr-2];
      double **backwardimhilohi = new double*[nvr-2];
      double **backwardimlolohi = new double*[nvr-2];
      double **backwardimhihilo = new double*[nvr-2];
      double **backwardimlohilo = new double*[nvr-2];
      double **backwardimhilolo = new double*[nvr-2];
      double **backwardimlololo = new double*[nvr-2];
      double **crossrehihihi = new double*[nvr-2];
      double **crossrelohihi = new double*[nvr-2];
      double **crossrehilohi = new double*[nvr-2];
      double **crossrelolohi = new double*[nvr-2];
      double **crossrehihilo = new double*[nvr-2];
      double **crossrelohilo = new double*[nvr-2];
      double **crossrehilolo = new double*[nvr-2];
      double **crossrelololo = new double*[nvr-2];
      double **crossimhihihi = new double*[nvr-2];
      double **crossimlohihi = new double*[nvr-2];
      double **crossimhilohi = new double*[nvr-2];
      double **crossimlolohi = new double*[nvr-2];
      double **crossimhihilo = new double*[nvr-2];
      double **crossimlohilo = new double*[nvr-2];
      double **crossimhilolo = new double*[nvr-2];
      double **crossimlololo = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardrehihihi[i] = new double[deg+1];
         forwardrelohihi[i] = new double[deg+1];
         forwardrehilohi[i] = new double[deg+1];
         forwardrelolohi[i] = new double[deg+1];
         forwardrehihilo[i] = new double[deg+1];
         forwardrelohilo[i] = new double[deg+1];
         forwardrehilolo[i] = new double[deg+1];
         forwardrelololo[i] = new double[deg+1];
         forwardimhihihi[i] = new double[deg+1];
         forwardimlohihi[i] = new double[deg+1];
         forwardimhilohi[i] = new double[deg+1];
         forwardimlolohi[i] = new double[deg+1];
         forwardimhihilo[i] = new double[deg+1];
         forwardimlohilo[i] = new double[deg+1];
         forwardimhilolo[i] = new double[deg+1];
         forwardimlololo[i] = new double[deg+1];
         backwardrehihihi[i] = new double[deg+1];
         backwardrelohihi[i] = new double[deg+1];
         backwardrehilohi[i] = new double[deg+1];
         backwardrelolohi[i] = new double[deg+1];
         backwardrehihilo[i] = new double[deg+1];
         backwardrelohilo[i] = new double[deg+1];
         backwardrehilolo[i] = new double[deg+1];
         backwardrelololo[i] = new double[deg+1];
         backwardimhihihi[i] = new double[deg+1];
         backwardimlohihi[i] = new double[deg+1];
         backwardimhilohi[i] = new double[deg+1];
         backwardimlolohi[i] = new double[deg+1];
         backwardimhihilo[i] = new double[deg+1];
         backwardimlohilo[i] = new double[deg+1];
         backwardimhilolo[i] = new double[deg+1];
         backwardimlololo[i] = new double[deg+1];
         crossrehihihi[i] = new double[deg+1];
         crossrelohihi[i] = new double[deg+1];
         crossrehilohi[i] = new double[deg+1];
         crossrelolohi[i] = new double[deg+1];
         crossrehihilo[i] = new double[deg+1];
         crossrelohilo[i] = new double[deg+1];
         crossrehilolo[i] = new double[deg+1];
         crossrelololo[i] = new double[deg+1];
         crossimhihihi[i] = new double[deg+1];
         crossimlohihi[i] = new double[deg+1];
         crossimhilohi[i] = new double[deg+1];
         crossimlolohi[i] = new double[deg+1];
         crossimhihilo[i] = new double[deg+1];
         crossimlohilo[i] = new double[deg+1];
         crossimhilolo[i] = new double[deg+1];
         crossimlololo[i] = new double[deg+1];
      }
      forwardrehihihi[nvr-2] = new double[deg+1];
      forwardrelohihi[nvr-2] = new double[deg+1];
      forwardrehilohi[nvr-2] = new double[deg+1];
      forwardrelolohi[nvr-2] = new double[deg+1];
      forwardrehihilo[nvr-2] = new double[deg+1];
      forwardrelohilo[nvr-2] = new double[deg+1];
      forwardrehilolo[nvr-2] = new double[deg+1];
      forwardrelololo[nvr-2] = new double[deg+1];
      forwardimhihihi[nvr-2] = new double[deg+1];
      forwardimlohihi[nvr-2] = new double[deg+1];
      forwardimhilohi[nvr-2] = new double[deg+1];
      forwardimlolohi[nvr-2] = new double[deg+1];
      forwardimhihilo[nvr-2] = new double[deg+1];
      forwardimlohilo[nvr-2] = new double[deg+1];
      forwardimhilolo[nvr-2] = new double[deg+1];
      forwardimlololo[nvr-2] = new double[deg+1];
      forwardrehihihi[nvr-1] = new double[deg+1];
      forwardrelohihi[nvr-1] = new double[deg+1];
      forwardrehilohi[nvr-1] = new double[deg+1];
      forwardrelolohi[nvr-1] = new double[deg+1];
      forwardrehihilo[nvr-1] = new double[deg+1];
      forwardrelohilo[nvr-1] = new double[deg+1];
      forwardrehilolo[nvr-1] = new double[deg+1];
      forwardrelololo[nvr-1] = new double[deg+1];
      forwardimhihihi[nvr-1] = new double[deg+1];
      forwardimlohihi[nvr-1] = new double[deg+1];
      forwardimhilohi[nvr-1] = new double[deg+1];
      forwardimlolohi[nvr-1] = new double[deg+1];
      forwardimhihilo[nvr-1] = new double[deg+1];
      forwardimlohilo[nvr-1] = new double[deg+1];
      forwardimhilolo[nvr-1] = new double[deg+1];
      forwardimlololo[nvr-1] = new double[deg+1];

      CPU_cmplx8_speel(nvr,deg,idx,
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
         inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
         inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
         inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
         forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
         forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
         forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
         forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
         backwardrehihihi,backwardrelohihi,backwardrehilohi,backwardrelolohi,
         backwardrehihilo,backwardrelohilo,backwardrehilolo,backwardrelololo,
         backwardimhihihi,backwardimlohihi,backwardimhilohi,backwardimlolohi,
         backwardimhihilo,backwardimlohilo,backwardimhilolo,backwardimlololo,
         crossrehihihi,crossrelohihi,crossrehilohi,crossrelolohi,
         crossrehihilo,crossrelohilo,crossrehilolo,crossrelololo,
         crossimhihihi,crossimlohihi,crossimhilohi,crossimlolohi,
         crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo);

      for(int i=0; i<deg+1; i++)          // assign value of the monomial
      {
         outputrehihihi[dim][i] = forwardrehihihi[nvr-1][i];
         outputrelohihi[dim][i] = forwardrelohihi[nvr-1][i];
         outputrehilohi[dim][i] = forwardrehilohi[nvr-1][i];
         outputrelolohi[dim][i] = forwardrelolohi[nvr-1][i];
         outputrehihilo[dim][i] = forwardrehihilo[nvr-1][i];
         outputrelohilo[dim][i] = forwardrelohilo[nvr-1][i];
         outputrehilolo[dim][i] = forwardrehilolo[nvr-1][i];
         outputrelololo[dim][i] = forwardrelololo[nvr-1][i];
         outputimhihihi[dim][i] = forwardimhihihi[nvr-1][i];
         outputimlohihi[dim][i] = forwardimlohihi[nvr-1][i];
         outputimhilohi[dim][i] = forwardimhilohi[nvr-1][i];
         outputimlolohi[dim][i] = forwardimlolohi[nvr-1][i];
         outputimhihilo[dim][i] = forwardimhihilo[nvr-1][i];
         outputimlohilo[dim][i] = forwardimlohilo[nvr-1][i];
         outputimhilolo[dim][i] = forwardimhilolo[nvr-1][i];
         outputimlololo[dim][i] = forwardimlololo[nvr-1][i];
      }

      if(nvr > 2)
      {
         int ix = idx[nvr-1];        // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++)
         {
            outputrehihihi[ix][i] = forwardrehihihi[nvr-2][i];
            outputrelohihi[ix][i] = forwardrelohihi[nvr-2][i];
            outputrehilohi[ix][i] = forwardrehilohi[nvr-2][i];
            outputrelolohi[ix][i] = forwardrelolohi[nvr-2][i];
            outputrehihilo[ix][i] = forwardrehihilo[nvr-2][i];
            outputrelohilo[ix][i] = forwardrelohilo[nvr-2][i];
            outputrehilolo[ix][i] = forwardrehilolo[nvr-2][i];
            outputrelololo[ix][i] = forwardrelololo[nvr-2][i];
            outputimhihihi[ix][i] = forwardimhihihi[nvr-2][i];
            outputimlohihi[ix][i] = forwardimlohihi[nvr-2][i];
            outputimhilohi[ix][i] = forwardimhilohi[nvr-2][i];
            outputimlolohi[ix][i] = forwardimlolohi[nvr-2][i];
            outputimhihilo[ix][i] = forwardimhihilo[nvr-2][i];
            outputimlohilo[ix][i] = forwardimlohilo[nvr-2][i];
            outputimhilolo[ix][i] = forwardimhilolo[nvr-2][i];
            outputimlololo[ix][i] = forwardimlololo[nvr-2][i];
         }

         ix = idx[0];                  // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputrehihihi[ix][i] = backwardrehihihi[nvr-3][i];
            outputrelohihi[ix][i] = backwardrelohihi[nvr-3][i];
            outputrehilohi[ix][i] = backwardrehilohi[nvr-3][i];
            outputrelolohi[ix][i] = backwardrelolohi[nvr-3][i];
            outputrehihilo[ix][i] = backwardrehihilo[nvr-3][i];
            outputrelohilo[ix][i] = backwardrelohilo[nvr-3][i];
            outputrehilolo[ix][i] = backwardrehilolo[nvr-3][i];
            outputrelololo[ix][i] = backwardrelololo[nvr-3][i];
            outputimhihihi[ix][i] = backwardimhihihi[nvr-3][i];
            outputimlohihi[ix][i] = backwardimlohihi[nvr-3][i];
            outputimhilohi[ix][i] = backwardimhilohi[nvr-3][i];
            outputimlolohi[ix][i] = backwardimlolohi[nvr-3][i];
            outputimhihilo[ix][i] = backwardimhihilo[nvr-3][i];
            outputimlohilo[ix][i] = backwardimlohilo[nvr-3][i];
            outputimhilolo[ix][i] = backwardimhilolo[nvr-3][i];
            outputimlololo[ix][i] = backwardimlololo[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];               // derivative with respect to x[k]
            for(int i=0; i<deg+1; i++)
            {
               outputrehihihi[ix][i] = crossrehihihi[k-1][i];
               outputrelohihi[ix][i] = crossrelohihi[k-1][i];
               outputrehilohi[ix][i] = crossrehilohi[k-1][i];
               outputrelolohi[ix][i] = crossrelolohi[k-1][i];
               outputrehihilo[ix][i] = crossrehihilo[k-1][i];
               outputrelohilo[ix][i] = crossrelohilo[k-1][i];
               outputrehilolo[ix][i] = crossrehilolo[k-1][i];
               outputrelololo[ix][i] = crossrelololo[k-1][i];
               outputimhihihi[ix][i] = crossimhihihi[k-1][i];
               outputimlohihi[ix][i] = crossimlohihi[k-1][i];
               outputimhilohi[ix][i] = crossimhilohi[k-1][i];
               outputimlolohi[ix][i] = crossimlolohi[k-1][i];
               outputimhihilo[ix][i] = crossimhihilo[k-1][i];
               outputimlohilo[ix][i] = crossimlohilo[k-1][i];
               outputimhilolo[ix][i] = crossimhilolo[k-1][i];
               outputimlololo[ix][i] = crossimlololo[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardrehihihi[i]);  free(forwardrelohihi[i]);
         free(forwardrehilohi[i]);  free(forwardrelolohi[i]);
         free(forwardrehihilo[i]);  free(forwardrelohilo[i]);
         free(forwardrehilolo[i]);  free(forwardrelololo[i]);
         free(forwardimhihihi[i]);  free(forwardimlohihi[i]);
         free(forwardimhilohi[i]);  free(forwardimlolohi[i]);
         free(forwardimhihilo[i]);  free(forwardimlohilo[i]);
         free(forwardimhilolo[i]);  free(forwardimlololo[i]);
         free(backwardrehihihi[i]); free(backwardrelohihi[i]);
         free(backwardrehilohi[i]); free(backwardrelolohi[i]);
         free(backwardrehihilo[i]); free(backwardrelohilo[i]);
         free(backwardrehilolo[i]); free(backwardrelololo[i]);
         free(backwardimhihihi[i]); free(backwardimlohihi[i]);
         free(backwardimhilohi[i]); free(backwardimlolohi[i]);
         free(backwardimhihilo[i]); free(backwardimlohilo[i]);
         free(backwardimhilolo[i]); free(backwardimlololo[i]);
         free(crossrehihihi[i]);    free(crossrelohihi[i]);
         free(crossrehilohi[i]);    free(crossrelolohi[i]);
         free(crossrehihilo[i]);    free(crossrelohilo[i]);
         free(crossrehilolo[i]);    free(crossrelololo[i]);
         free(crossimhihihi[i]);    free(crossimlohihi[i]);
         free(crossimhilohi[i]);    free(crossimlolohi[i]);
         free(crossimhihilo[i]);    free(crossimlohilo[i]);
         free(crossimhilolo[i]);    free(crossimlololo[i]);
      }
      free(forwardrehihihi[nvr-2]); free(forwardrelohihi[nvr-2]);
      free(forwardrehilohi[nvr-2]); free(forwardrelolohi[nvr-2]);
      free(forwardrehihilo[nvr-2]); free(forwardrelohilo[nvr-2]);
      free(forwardrehilolo[nvr-2]); free(forwardrelololo[nvr-2]);
      free(forwardrehihihi[nvr-1]); free(forwardrelohihi[nvr-1]);
      free(forwardrehilohi[nvr-1]); free(forwardrelolohi[nvr-1]);
      free(forwardrehihilo[nvr-1]); free(forwardrelohilo[nvr-1]);
      free(forwardrehilolo[nvr-1]); free(forwardrelololo[nvr-1]);
      free(forwardrehihihi);        free(forwardrelohihi);
      free(forwardrehilohi);        free(forwardrelolohi);
      free(forwardrehihilo);        free(forwardrelohilo);
      free(forwardrehilolo);        free(forwardrelololo);
      free(forwardimhihihi[nvr-2]); free(forwardimlohihi[nvr-2]);
      free(forwardimhilohi[nvr-2]); free(forwardimlolohi[nvr-2]);
      free(forwardimhihilo[nvr-2]); free(forwardimlohilo[nvr-2]);
      free(forwardimhilolo[nvr-2]); free(forwardimlololo[nvr-2]);
      free(forwardimhihihi[nvr-1]); free(forwardimlohihi[nvr-1]);
      free(forwardimhilohi[nvr-1]); free(forwardimlolohi[nvr-1]);
      free(forwardimhihilo[nvr-1]); free(forwardimlohilo[nvr-1]);
      free(forwardimhilolo[nvr-1]); free(forwardimlololo[nvr-1]);
      free(forwardimhihihi);        free(forwardimlohihi);
      free(forwardimhilohi);        free(forwardimlolohi);
      free(forwardimhihilo);        free(forwardimlohilo);
      free(forwardimhilolo);        free(forwardimlololo);
      free(backwardrehihihi);       free(backwardrelohihi);
      free(backwardrehilohi);       free(backwardrelolohi);
      free(backwardrehihilo);       free(backwardrelohilo);
      free(backwardrehilolo);       free(backwardrelololo);
      free(backwardimhihihi);       free(backwardimlohihi);
      free(backwardimhilohi);       free(backwardimlolohi);
      free(backwardimhihilo);       free(backwardimlohilo);
      free(backwardimhilolo);       free(backwardimlololo);
      free(crossrehihihi);          free(crossrelohihi);
      free(crossrehilohi);          free(crossrelolohi);
      free(crossrehihilo);          free(crossrelohilo);
      free(crossrehilolo);          free(crossrelololo);
      free(crossimhihihi);          free(crossimlohihi);
      free(crossimhilohi);          free(crossimlolohi);
      free(crossimhihilo);          free(crossimlohilo);
      free(crossimhilolo);          free(crossimlololo);
   }
}
