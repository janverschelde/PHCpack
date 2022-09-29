/* The file dbl4_monomials_host.cpp defines the functions specified
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
   int ix1 = idx[0];
   int ix2;
                                                           // f[0] = cff*x[0]
   CPU_cmplx4_product(deg,
      cffrehihi,cffrelohi,cffrehilo,cffrelolo,
      cffimhihi,cffimlohi,cffimhilo,cffimlolo,
      inputrehihi[ix1],inputrelohi[ix1],inputrehilo[ix1],inputrelolo[ix1],
      inputimhihi[ix1],inputimlohi[ix1],inputimhilo[ix1],inputimlolo[ix1],
      forwardrehihi[0],forwardrelohi[0],forwardrehilo[0],forwardrelolo[0],
      forwardimhihi[0],forwardimlohi[0],forwardimhilo[0], forwardimlolo[0]);

   for(int i=1; i<nvr; i++)
   {                                                    // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_cmplx4_product(deg,
         forwardrehihi[i-1],forwardrelohi[i-1],
         forwardrehilo[i-1],forwardrelolo[i-1],
         forwardimhihi[i-1],forwardimlohi[i-1],
         forwardimhilo[i-1],forwardimlolo[i-1],
         inputrehihi[ix2],inputrelohi[ix2],inputrehilo[ix2],inputrelolo[ix2],
         inputimhihi[ix2],inputimlohi[ix2],inputimhilo[ix2],inputimlolo[ix2],
         forwardrehihi[i],forwardrelohi[i],forwardrehilo[i],forwardrelolo[i],
         forwardimhihi[i],forwardimlohi[i],forwardimhilo[i],forwardimlolo[i]);
   }
   if(nvr > 2)
   {                                                  // b[0] = x[n-1]*x[n-2]
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];
      CPU_cmplx4_product(deg,
         inputrehihi[ix1],inputrelohi[ix1],inputrehilo[ix1],inputrelolo[ix1],
         inputimhihi[ix1],inputimlohi[ix1],inputimhilo[ix1],inputimlolo[ix1],
         inputrehihi[ix2],inputrelohi[ix2],inputrehilo[ix2],inputrelolo[ix2],
         inputimhihi[ix2],inputimlohi[ix2],inputimhilo[ix2],inputimlolo[ix2],
         backwardrehihi[0],backwardrelohi[0],
         backwardrehilo[0],backwardrelolo[0],
         backwardimhihi[0],backwardimlohi[0],
         backwardimhilo[0],backwardimlolo[0]);

      for(int i=1; i<nvr-2; i++)
      {                                             // b[i] = b[i-1]*x[x-2-i]
         ix2 = idx[nvr-2-i];
         CPU_cmplx4_product(deg,
            backwardrehihi[i-1],backwardrelohi[i-1],
            backwardrehilo[i-1],backwardrelolo[i-1],
            backwardimhihi[i-1],backwardimlohi[i-1],
            backwardimhilo[i-1],backwardimlolo[i-1],
            inputrehihi[ix2],inputrelohi[ix2],
            inputrehilo[ix2],inputrelolo[ix2],
            inputimhihi[ix2],inputimlohi[ix2],
            inputimhilo[ix2],inputimlolo[ix2],
            backwardrehihi[i],backwardrelohi[i],
            backwardrehilo[i],backwardrelolo[i],
            backwardimhihi[i],backwardimlohi[i],
            backwardimhilo[i],backwardimlolo[i]);
      }
                                                       // b[n-3] = b[n-3]*cff
      CPU_cmplx4_product(deg,
         backwardrehihi[nvr-3],backwardrelohi[nvr-3],
         backwardrehilo[nvr-3],backwardrelolo[nvr-3],
         backwardimhihi[nvr-3],backwardimlohi[nvr-3],
         backwardimhilo[nvr-3],backwardimlolo[nvr-3],
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         crossrehihi[0],crossrelohi[0],crossrehilo[0],crossrelolo[0],
         crossimhihi[0],crossimlohi[0],crossimhilo[0],crossimlolo[0]); 
                                                       // cross is work space
      for(int i=0; i<=deg; i++)
      {
         backwardrehihi[nvr-3][i] = crossrehihi[0][i];
         backwardrelohi[nvr-3][i] = crossrelohi[0][i];
         backwardrehilo[nvr-3][i] = crossrehilo[0][i];
         backwardrelolo[nvr-3][i] = crossrelolo[0][i];
         backwardimhihi[nvr-3][i] = crossimhihi[0][i];
         backwardimlohi[nvr-3][i] = crossimlohi[0][i];
         backwardimhilo[nvr-3][i] = crossimhilo[0][i];
         backwardimlolo[nvr-3][i] = crossimlolo[0][i];
      }
      if(nvr == 3)
      {                                                   // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_cmplx4_product(deg,
            forwardrehihi[0],forwardrelohi[0],
            forwardrehilo[0],forwardrelolo[0],
            forwardimhihi[0],forwardimlohi[0],
            forwardimhilo[0],forwardimlolo[0],
            inputrehihi[ix2],inputrelohi[ix2],
            inputrehilo[ix2],inputrelolo[ix2],
            inputimhihi[ix2],inputimlohi[ix2],
            inputimhilo[ix2],inputimlolo[ix2],
            crossrehihi[0],crossrelohi[0],crossrehilo[0],crossrelolo[0],
            crossimhihi[0],crossimlohi[0],crossimhilo[0],crossimlolo[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                            // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_cmplx4_product(deg,
               forwardrehihi[i],forwardrelohi[i],
               forwardrehilo[i],forwardrelolo[i],
               forwardimhihi[i],forwardimlohi[i],
               forwardimhilo[i],forwardimlolo[i],
               backwardrehihi[ix2],backwardrelohi[ix2],
               backwardrehilo[ix2],backwardrelolo[ix2],
               backwardimhihi[ix2],backwardimlohi[ix2],
               backwardimhilo[ix2],backwardimlolo[ix2],
               crossrehihi[i],crossrelohi[i],crossrehilo[i],crossrelolo[i],
               crossimhihi[i],crossimlohi[i],crossimhilo[i],crossimlolo[i]);
         }
                                                    // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_cmplx4_product(deg,
            forwardrehihi[nvr-3],forwardrelohi[nvr-3],
            forwardrehilo[nvr-3],forwardrelolo[nvr-3],
            forwardimhihi[nvr-3],forwardimlohi[nvr-3],
            forwardimhilo[nvr-3],forwardimlolo[nvr-3],
            inputrehihi[ix2],inputrelohi[ix2],
            inputrehilo[ix2],inputrelolo[ix2],
            inputimhihi[ix2],inputimlohi[ix2],
            inputimhilo[ix2],inputimlolo[ix2],
            crossrehihi[nvr-3],crossrelohi[nvr-3],
            crossrehilo[nvr-3],crossrelolo[nvr-3],
            crossimhihi[nvr-3],crossimlohi[nvr-3],
            crossimhilo[nvr-3],crossimlolo[nvr-3]);
      }
   }
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
      // CPU_dbl4_product(deg,
      //    outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim],
      //    cffhihi,cfflohi,cffhilo,cfflolo,
      //    outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim]);
      // wrong, cannot use input and output!
      double *acchihi = new double[deg+1];
      double *acclohi = new double[deg+1];
      double *acchilo = new double[deg+1];
      double *acclolo = new double[deg+1];

      CPU_dbl4_product(deg,
         outputhihi[dim],outputlohi[dim],outputhilo[dim],outputlolo[dim],
         cffhihi,cfflohi,cffhilo,cfflolo,
         acchihi,acclohi,acchilo,acclolo);

      for(int i=0; i<=deg; i++)
      {
         outputhihi[dim][i] = acchihi[i];
         outputlohi[dim][i] = acclohi[i];
         outputhilo[dim][i] = acchilo[i];
         outputlolo[dim][i] = acclolo[i];
      }
      free(acchihi); free(acclohi);
      free(acchilo); free(acclolo);

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
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_cmplx4_product(deg,
         inputrehihi[ix],inputrelohi[ix],inputrehilo[ix],inputrelolo[ix],
         inputimhihi[ix],inputimlohi[ix],inputimhilo[ix],inputimlolo[ix],
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         outputrehihi[dim],outputrelohi[dim],
         outputrehilo[dim],outputrelolo[dim],
         outputimhihi[dim],outputimlohi[dim],
         outputimhilo[dim],outputimlolo[dim]);

      for(int i=0; i<=deg; i++) 
      {
         outputrehihi[ix][i] = cffrehihi[i];
         outputrelohi[ix][i] = cffrelohi[i];
         outputrehilo[ix][i] = cffrehilo[i];
         outputrelolo[ix][i] = cffrelolo[i];
         outputimhihi[ix][i] = cffimhihi[i];
         outputimlohi[ix][i] = cffimlohi[i];
         outputimhilo[ix][i] = cffimhilo[i];
         outputimlolo[ix][i] = cffimlolo[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0];
      int ix2 = idx[1];

      CPU_cmplx4_product(deg,
         inputrehihi[ix1],inputrelohi[ix1],inputrehilo[ix1],inputrelolo[ix1],
         inputimhihi[ix1],inputimlohi[ix1],inputimhilo[ix1],inputimlolo[ix1],
         inputrehihi[ix2],inputrelohi[ix2],inputrehilo[ix2],inputrelolo[ix2],
         inputimhihi[ix2],inputimlohi[ix2],inputimhilo[ix2],inputimlolo[ix2],
         outputrehihi[dim],outputrelohi[dim],
         outputrehilo[dim],outputrelolo[dim],
         outputimhihi[dim],outputimlohi[dim],
         outputimhilo[dim],outputimlolo[dim]);
      // CPU_cmplx4_product(deg,
      //    outputrehihi[dim],outputrelohi[dim],
      //    outputrehilo[dim],outputrelolo[dim],
      //    outputimhihi[dim],outputimlohi[dim],
      //    outputimhilo[dim],outputimlolo[dim],
      //    cffrehihi,cffrelohi,cffrehilo,cffrelolo,
      //    cffimhihi,cffimlohi,cffimhilo,cffimlolo,
      //    outputrehihi[dim],outputrelohi[dim],
      //    outputrehilo[dim],outputrelolo[dim],
      //    outputimhihi[dim],outputimlohi[dim],
      //    outputimhilo[dim],outputimlolo[dim]); // wrong!

      double *accrehihi = new double[deg+1];
      double *accrelohi = new double[deg+1];
      double *accrehilo = new double[deg+1];
      double *accrelolo = new double[deg+1];
      double *accimhihi = new double[deg+1];
      double *accimlohi = new double[deg+1];
      double *accimhilo = new double[deg+1];
      double *accimlolo = new double[deg+1];

      CPU_cmplx4_product(deg,
         outputrehihi[dim],outputrelohi[dim],
         outputrehilo[dim],outputrelolo[dim],
         outputimhihi[dim],outputimlohi[dim],
         outputimhilo[dim],outputimlolo[dim],
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         accrehihi,accrelohi,accrehilo,accrelolo,
         accimhihi,accimlohi,accimhilo,accimlolo);

      for(int i=0; i<=deg; i++)
      {
         outputrehihi[dim][i] = accrehihi[i];
         outputrelohi[dim][i] = accrelohi[i];
         outputrehilo[dim][i] = accrehilo[i];
         outputrelolo[dim][i] = accrelolo[i];
         outputimhihi[dim][i] = accimhihi[i];
         outputimlohi[dim][i] = accimlohi[i];
         outputimhilo[dim][i] = accimhilo[i];
         outputimlolo[dim][i] = accimlolo[i];
      }
      free(accrehihi); free(accrelohi); free(accimhihi); free(accimlohi);
      free(accrehilo); free(accrelolo); free(accimhilo); free(accimlolo);

      CPU_cmplx4_product(deg,
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         inputrehihi[ix1],inputrelohi[ix1],inputrehilo[ix1],inputrelolo[ix1],
         inputimhihi[ix1],inputimlohi[ix1],inputimhilo[ix1],inputimlolo[ix1],
         outputrehihi[ix2],outputrelohi[ix2],
         outputrehilo[ix2],outputrelolo[ix2],
         outputimhihi[ix2],outputimlohi[ix2],
         outputimhilo[ix2],outputimlolo[ix2]);
      CPU_cmplx4_product(deg,
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         inputrehihi[ix2],inputrelohi[ix2],inputrehilo[ix2],inputrelolo[ix2],
         inputimhihi[ix2],inputimlohi[ix2],inputimhilo[ix2],inputimlolo[ix2],
         outputrehihi[ix1],outputrelohi[ix1],
         outputrehilo[ix1],outputrelolo[ix1],
         outputimhihi[ix1],outputimlohi[ix1],
         outputimhilo[ix1],outputimlolo[ix1]);
   }
   else
   {
      double **forwardrehihi = new double*[nvr];
      double **forwardrelohi = new double*[nvr];
      double **forwardrehilo = new double*[nvr];
      double **forwardrelolo = new double*[nvr];
      double **forwardimhihi = new double*[nvr];
      double **forwardimlohi = new double*[nvr];
      double **forwardimhilo = new double*[nvr];
      double **forwardimlolo = new double*[nvr];
      double **backwardrehihi = new double*[nvr-2];
      double **backwardrelohi = new double*[nvr-2];
      double **backwardrehilo = new double*[nvr-2];
      double **backwardrelolo = new double*[nvr-2];
      double **backwardimhihi = new double*[nvr-2];
      double **backwardimlohi = new double*[nvr-2];
      double **backwardimhilo = new double*[nvr-2];
      double **backwardimlolo = new double*[nvr-2];
      double **crossrehihi = new double*[nvr-2];
      double **crossrelohi = new double*[nvr-2];
      double **crossrehilo = new double*[nvr-2];
      double **crossrelolo = new double*[nvr-2];
      double **crossimhihi = new double*[nvr-2];
      double **crossimlohi = new double*[nvr-2];
      double **crossimhilo = new double*[nvr-2];
      double **crossimlolo = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardrehihi[i] = new double[deg+1];
         forwardrelohi[i] = new double[deg+1];
         forwardrehilo[i] = new double[deg+1];
         forwardrelolo[i] = new double[deg+1];
         forwardimhihi[i] = new double[deg+1];
         forwardimlohi[i] = new double[deg+1];
         forwardimhilo[i] = new double[deg+1];
         forwardimlolo[i] = new double[deg+1];
         backwardrehihi[i] = new double[deg+1];
         backwardrelohi[i] = new double[deg+1];
         backwardrehilo[i] = new double[deg+1];
         backwardrelolo[i] = new double[deg+1];
         backwardimhihi[i] = new double[deg+1];
         backwardimlohi[i] = new double[deg+1];
         backwardimhilo[i] = new double[deg+1];
         backwardimlolo[i] = new double[deg+1];
         crossrehihi[i] = new double[deg+1];
         crossrelohi[i] = new double[deg+1];
         crossrehilo[i] = new double[deg+1];
         crossrelolo[i] = new double[deg+1];
         crossimhihi[i] = new double[deg+1];
         crossimlohi[i] = new double[deg+1];
         crossimhilo[i] = new double[deg+1];
         crossimlolo[i] = new double[deg+1];
      }
      forwardrehihi[nvr-2] = new double[deg+1];
      forwardrelohi[nvr-2] = new double[deg+1];
      forwardrehilo[nvr-2] = new double[deg+1];
      forwardrelolo[nvr-2] = new double[deg+1];
      forwardimhihi[nvr-2] = new double[deg+1];
      forwardimlohi[nvr-2] = new double[deg+1];
      forwardimhilo[nvr-2] = new double[deg+1];
      forwardimlolo[nvr-2] = new double[deg+1];
      forwardrehihi[nvr-1] = new double[deg+1];
      forwardrelohi[nvr-1] = new double[deg+1];
      forwardrehilo[nvr-1] = new double[deg+1];
      forwardrelolo[nvr-1] = new double[deg+1];
      forwardimhihi[nvr-1] = new double[deg+1];
      forwardimlohi[nvr-1] = new double[deg+1];
      forwardimhilo[nvr-1] = new double[deg+1];
      forwardimlolo[nvr-1] = new double[deg+1];

      CPU_cmplx4_speel(nvr,deg,idx,
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         inputrehihi,inputrelohi,inputrehilo,inputrelolo,
         inputimhihi,inputimlohi,inputimhilo,inputimlolo,
         forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
         forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
         backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
         backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,
         crossrehihi,crossrelohi,crossrehilo,crossrelolo,
         crossimhihi,crossimlohi,crossimhilo,crossimlolo);

      for(int i=0; i<deg+1; i++)          // assign value of the monomial
      {
         outputrehihi[dim][i] = forwardrehihi[nvr-1][i];
         outputrelohi[dim][i] = forwardrelohi[nvr-1][i];
         outputrehilo[dim][i] = forwardrehilo[nvr-1][i];
         outputrelolo[dim][i] = forwardrelolo[nvr-1][i];
         outputimhihi[dim][i] = forwardimhihi[nvr-1][i];
         outputimlohi[dim][i] = forwardimlohi[nvr-1][i];
         outputimhilo[dim][i] = forwardimhilo[nvr-1][i];
         outputimlolo[dim][i] = forwardimlolo[nvr-1][i];
      }

      if(nvr > 2)
      {
         int ix = idx[nvr-1];        // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++)
         {
            outputrehihi[ix][i] = forwardrehihi[nvr-2][i];
            outputrelohi[ix][i] = forwardrelohi[nvr-2][i];
            outputrehilo[ix][i] = forwardrehilo[nvr-2][i];
            outputrelolo[ix][i] = forwardrelolo[nvr-2][i];
            outputimhihi[ix][i] = forwardimhihi[nvr-2][i];
            outputimlohi[ix][i] = forwardimlohi[nvr-2][i];
            outputimhilo[ix][i] = forwardimhilo[nvr-2][i];
            outputimlolo[ix][i] = forwardimlolo[nvr-2][i];
         }

         ix = idx[0];                  // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputrehihi[ix][i] = backwardrehihi[nvr-3][i];
            outputrelohi[ix][i] = backwardrelohi[nvr-3][i];
            outputrehilo[ix][i] = backwardrehilo[nvr-3][i];
            outputrelolo[ix][i] = backwardrelolo[nvr-3][i];
            outputimhihi[ix][i] = backwardimhihi[nvr-3][i];
            outputimlohi[ix][i] = backwardimlohi[nvr-3][i];
            outputimhilo[ix][i] = backwardimhilo[nvr-3][i];
            outputimlolo[ix][i] = backwardimlolo[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];               // derivative with respect to x[k]
            for(int i=0; i<deg+1; i++)
            {
               outputrehihi[ix][i] = crossrehihi[k-1][i];
               outputrelohi[ix][i] = crossrelohi[k-1][i];
               outputrehilo[ix][i] = crossrehilo[k-1][i];
               outputrelolo[ix][i] = crossrelolo[k-1][i];
               outputimhihi[ix][i] = crossimhihi[k-1][i];
               outputimlohi[ix][i] = crossimlohi[k-1][i];
               outputimhilo[ix][i] = crossimhilo[k-1][i];
               outputimlolo[ix][i] = crossimlolo[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardimhihi[i]);  free(forwardimlohi[i]);
         free(forwardimhilo[i]);  free(forwardimlolo[i]);
         free(forwardrehihi[i]);  free(forwardrelohi[i]);
         free(forwardrehilo[i]);  free(forwardrelolo[i]);
         free(backwardrehihi[i]); free(backwardrelohi[i]);
         free(backwardrehilo[i]); free(backwardrelolo[i]);
         free(backwardimhihi[i]); free(backwardimlohi[i]);
         free(backwardimhilo[i]); free(backwardimlolo[i]);
         free(crossrehihi[i]);    free(crossrelohi[i]);
         free(crossrehilo[i]);    free(crossrelolo[i]);
         free(crossimhihi[i]);    free(crossimlohi[i]);
         free(crossimhilo[i]);    free(crossimlolo[i]);
      }
      free(forwardrehihi[nvr-2]); free(forwardrelohi[nvr-2]);
      free(forwardrehilo[nvr-2]); free(forwardrelolo[nvr-2]);
      free(forwardrehihi[nvr-1]); free(forwardrelohi[nvr-1]);
      free(forwardrehilo[nvr-1]); free(forwardrelolo[nvr-1]);
      free(forwardrehihi);        free(forwardrelohi);
      free(forwardrehilo);        free(forwardrelolo);
      free(forwardimhihi[nvr-2]); free(forwardimlohi[nvr-2]);
      free(forwardimhilo[nvr-2]); free(forwardimlolo[nvr-2]);
      free(forwardimhihi[nvr-1]); free(forwardimlohi[nvr-1]);
      free(forwardimhilo[nvr-1]); free(forwardimlolo[nvr-1]);
      free(forwardimhihi);        free(forwardimlohi);
      free(forwardimhilo);        free(forwardimlolo);
      free(backwardrehihi);       free(backwardrelohi);
      free(backwardrehilo);       free(backwardrelolo);
      free(backwardimhihi);       free(backwardimlohi);
      free(backwardimhilo);       free(backwardimlolo);
      free(crossrehihi);          free(crossrelohi);
      free(crossrehilo);          free(crossrelolo);
      free(crossimhihi);          free(crossimlohi);
      free(crossimhilo);          free(crossimlolo);
   }
}
