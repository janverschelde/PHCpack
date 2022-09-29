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
   int ix1 = idx[0];
   int ix2;
                                                           // f[0] = cff*x[0]
   CPU_cmplx2_product(deg,
      cffrehi,cffrelo,cffimhi,cffimlo,
      inputrehi[ix1],inputrelo[ix1],inputimhi[ix1],inputimlo[ix1],
      forwardrehi[0],forwardrelo[0],forwardimhi[0],forwardimlo[0]);
   for(int i=1; i<nvr; i++)
   {                                                    // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_cmplx2_product(deg,
         forwardrehi[i-1],forwardrelo[i-1],forwardimhi[i-1],forwardimlo[i-1],
         inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
         forwardrehi[i],forwardrelo[i],forwardimhi[i],forwardimlo[i]);
   }
   if(nvr > 2)
   {                                                  // b[0] = x[n-1]*x[n-2]
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];
      CPU_cmplx2_product(deg,
         inputrehi[ix1],inputrelo[ix1],inputimhi[ix1],inputimlo[ix1],
         inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
         backwardrehi[0],backwardrelo[0],backwardimhi[0],backwardimlo[0]);
      for(int i=1; i<nvr-2; i++)
      {                                             // b[i] = b[i-1]*x[x-2-i]
         ix2 = idx[nvr-2-i];
         CPU_cmplx2_product(deg,
            backwardrehi[i-1],backwardrelo[i-1],
            backwardimhi[i-1],backwardimlo[i-1],
            inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
            backwardrehi[i],backwardrelo[i],backwardimhi[i],backwardimlo[i]);
      }
                                                       // b[n-3] = b[n-3]*cff
      CPU_cmplx2_product(deg,
         backwardrehi[nvr-3],backwardrelo[nvr-3],
         backwardimhi[nvr-3],backwardimlo[nvr-3],
         cffrehi,cffrelo,cffimhi,cffimlo,
         crossrehi[0],crossrelo[0],crossimhi[0],crossimlo[0]); 
                                                       // cross is work space
      for(int i=0; i<=deg; i++)
      {
         backwardrehi[nvr-3][i] = crossrehi[0][i];
         backwardrelo[nvr-3][i] = crossrelo[0][i];
         backwardimhi[nvr-3][i] = crossimhi[0][i];
         backwardimlo[nvr-3][i] = crossimlo[0][i];
      }
      if(nvr == 3)
      {                                                   // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_cmplx2_product(deg,
            forwardrehi[0],forwardrelo[0],forwardimhi[0],forwardimlo[0],
            inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
            crossrehi[0],crossrelo[0],crossimhi[0],crossimlo[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                            // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_cmplx2_product(deg,
               forwardrehi[i],forwardrelo[i],forwardimhi[i],forwardimlo[i],
               backwardrehi[ix2],backwardrelo[ix2],
               backwardimhi[ix2],backwardimlo[ix2],
               crossrehi[i],crossrelo[i],crossimhi[i],crossimlo[i]);
         }
                                                    // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_cmplx2_product(deg,
            forwardrehi[nvr-3],forwardrelo[nvr-3],
            forwardimhi[nvr-3],forwardimlo[nvr-3],
            inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
            crossrehi[nvr-3],crossrelo[nvr-3],
            crossimhi[nvr-3],crossimlo[nvr-3]);
      }
   }
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
      // CPU_dbl2_product(deg,outputhi[dim],outputlo[dim],cffhi,cfflo,
      //                  outputhi[dim],outputlo[dim]); // wrong!
      double *acchi = new double[deg+1];
      double *acclo = new double[deg+1];
      CPU_dbl2_product(deg,outputhi[dim],outputlo[dim],
                       cffhi,cfflo,acchi,acclo); 
      for(int i=0; i<=deg; i++)
      {
         outputhi[dim][i] = acchi[i];
         outputlo[dim][i] = acclo[i];
      }
      free(acchi);
      free(acclo);

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
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_cmplx2_product(deg,
         inputrehi[ix],inputrelo[ix],inputimhi[ix],inputimlo[ix],
         cffrehi,cffrelo,cffimhi,cffimlo,
         outputrehi[dim],outputrelo[dim],outputimhi[dim],outputimlo[dim]);

      for(int i=0; i<=deg; i++) 
      {
         outputrehi[ix][i] = cffrehi[i]; outputimhi[ix][i] = cffimhi[i];
         outputrelo[ix][i] = cffrelo[i]; outputimlo[ix][i] = cffimlo[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0];
      int ix2 = idx[1];

      CPU_cmplx2_product(deg,
         inputrehi[ix1],inputrelo[ix1],inputimhi[ix1],inputimlo[ix1],
         inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
         outputrehi[dim],outputrelo[dim],outputimhi[dim],outputimlo[dim]);
      // CPU_cmplx2_product(deg, // wrong!
      //    outputrehi[dim],outputrelo[dim],outputimhi[dim],outputimlo[dim],
      //    cffrehi,cffrelo,cffimhi,cffimlo,
      //    outputrehi[dim],outputrelo[dim],outputimhi[dim],outputimlo[dim]);
      double *accrehi = new double[deg+1];
      double *accrelo = new double[deg+1];
      double *accimhi = new double[deg+1];
      double *accimlo = new double[deg+1];
      CPU_cmplx2_product(deg,
         outputrehi[dim],outputrelo[dim],outputimhi[dim],outputimlo[dim],
         cffrehi,cffrelo,cffimhi,cffimlo,accrehi,accrelo,accimhi,accimlo);
      for(int i=0; i<=deg; i++)
      {
         outputrehi[dim][i] = accrehi[i];
         outputrelo[dim][i] = accrelo[i];
         outputimhi[dim][i] = accimhi[i];
         outputimlo[dim][i] = accimlo[i];
      }
      free(accrehi); free(accrelo);
      free(accimhi); free(accimlo);

      CPU_cmplx2_product(deg,
         cffrehi,cffrelo,cffimhi,cffimlo,
         inputrehi[ix1],inputrelo[ix1],inputimhi[ix1],inputimlo[ix1],
         outputrehi[ix2],outputrelo[ix2],outputimhi[ix2],outputimlo[ix2]);
      CPU_cmplx2_product(deg,
         cffrehi,cffrelo,cffimhi,cffimlo,
         inputrehi[ix2],inputrelo[ix2],inputimhi[ix2],inputimlo[ix2],
         outputrehi[ix1],outputrelo[ix1],outputimhi[ix1],outputimlo[ix1]);
   }
   else
   {
      double **forwardrehi = new double*[nvr];
      double **forwardrelo = new double*[nvr];
      double **forwardimhi = new double*[nvr];
      double **forwardimlo = new double*[nvr];
      double **backwardrehi = new double*[nvr-2];
      double **backwardrelo = new double*[nvr-2];
      double **backwardimhi = new double*[nvr-2];
      double **backwardimlo = new double*[nvr-2];
      double **crossrehi = new double*[nvr-2];
      double **crossrelo = new double*[nvr-2];
      double **crossimhi = new double*[nvr-2];
      double **crossimlo = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardrehi[i] = new double[deg+1];
         forwardrelo[i] = new double[deg+1];
         forwardimhi[i] = new double[deg+1];
         forwardimlo[i] = new double[deg+1];
         backwardrehi[i] = new double[deg+1];
         backwardrelo[i] = new double[deg+1];
         backwardimhi[i] = new double[deg+1];
         backwardimlo[i] = new double[deg+1];
         crossrehi[i] = new double[deg+1];
         crossrelo[i] = new double[deg+1];
         crossimhi[i] = new double[deg+1];
         crossimlo[i] = new double[deg+1];
      }
      forwardrehi[nvr-2] = new double[deg+1];
      forwardrelo[nvr-2] = new double[deg+1];
      forwardimhi[nvr-2] = new double[deg+1];
      forwardimlo[nvr-2] = new double[deg+1];
      forwardrehi[nvr-1] = new double[deg+1];
      forwardrelo[nvr-1] = new double[deg+1];
      forwardimhi[nvr-1] = new double[deg+1];
      forwardimlo[nvr-1] = new double[deg+1];

      CPU_cmplx2_speel(nvr,deg,idx,
         cffrehi,cffrelo,cffimhi,cffimlo,
         inputrehi,inputrelo,inputimhi,inputimlo,
         forwardrehi,forwardrelo,forwardimhi,forwardimlo,
         backwardrehi,backwardrelo,backwardimhi,backwardimlo,
         crossrehi,crossrelo,crossimhi,crossimlo);

      for(int i=0; i<deg+1; i++)          // assign value of the monomial
      {
         outputrehi[dim][i] = forwardrehi[nvr-1][i];
         outputrelo[dim][i] = forwardrelo[nvr-1][i];
         outputimhi[dim][i] = forwardimhi[nvr-1][i];
         outputimlo[dim][i] = forwardimlo[nvr-1][i];
      }

      if(nvr > 2)
      {
         int ix = idx[nvr-1];        // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++)
         {
            outputrehi[ix][i] = forwardrehi[nvr-2][i];
            outputrelo[ix][i] = forwardrelo[nvr-2][i];
            outputimhi[ix][i] = forwardimhi[nvr-2][i];
            outputimlo[ix][i] = forwardimlo[nvr-2][i];
         }

         ix = idx[0];                  // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputrehi[ix][i] = backwardrehi[nvr-3][i];
            outputrelo[ix][i] = backwardrelo[nvr-3][i];
            outputimhi[ix][i] = backwardimhi[nvr-3][i];
            outputimlo[ix][i] = backwardimlo[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];               // derivative with respect to x[k]
            for(int i=0; i<deg+1; i++)
            {
               outputrehi[ix][i] = crossrehi[k-1][i];
               outputrelo[ix][i] = crossrelo[k-1][i];
               outputimhi[ix][i] = crossimhi[k-1][i];
               outputimlo[ix][i] = crossimlo[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardimhi[i]);  free(forwardimlo[i]);
         free(forwardrehi[i]);  free(forwardrelo[i]);
         free(backwardrehi[i]); free(backwardrelo[i]);
         free(backwardimhi[i]); free(backwardimlo[i]);
         free(crossrehi[i]);    free(crossrelo[i]);
         free(crossimhi[i]);    free(crossimlo[i]);
      }
      free(forwardrehi[nvr-2]); free(forwardrelo[nvr-2]);
      free(forwardrehi[nvr-1]); free(forwardrelo[nvr-1]);
      free(forwardrehi);        free(forwardrelo);
      free(forwardimhi[nvr-2]); free(forwardimlo[nvr-2]);
      free(forwardimhi[nvr-1]); free(forwardimlo[nvr-1]);
      free(forwardimhi);        free(forwardimlo);
      free(backwardrehi);       free(backwardrelo);
      free(backwardimhi);       free(backwardimlo);
      free(crossrehi);          free(crossrelo);
      free(crossimhi);          free(crossimlo);
   }
}
