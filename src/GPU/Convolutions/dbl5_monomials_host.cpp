/* The file dbl5_monomials_host.cpp defines functions specified
 * in dbl5_monomials_host.h. */

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
#include "dbl5_convolutions_host.h"
#include "dbl5_monomials_host.h"

void CPU_dbl5_speel
 ( int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi, double **inputrg,
   double **inputpk, double **forwardtb, double **forwardix,
   double **forwardmi, double **forwardrg, double **forwardpk,
   double **backwardtb, double **backwardix, double **backwardmi,
   double **backwardrg, double **backwardpk, double **crosstb,
   double **crossix, double **crossmi, double **crossrg, double **crosspk )
{
   int ix1 = idx[0];
   int ix2;
                                                   // f[0] = cff*x[0] 
   CPU_dbl5_product(deg,cfftb,cffix,cffmi,cffrg,cffpk,
      inputtb[ix1],inputix[ix1],inputmi[ix1],inputrg[ix1],inputpk[ix1],
      forwardtb[0],forwardix[0],forwardmi[0],forwardrg[0],forwardpk[0]);

   for(int i=1; i<nvr; i++)
   {                                               // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_dbl5_product(deg,
         forwardtb[i-1],forwardix[i-1],forwardmi[i-1],
         forwardrg[i-1],forwardpk[i-1],
         inputtb[ix2],inputix[ix2],inputmi[ix2],inputrg[ix2],inputpk[ix2],
         forwardtb[i],forwardix[i],forwardmi[i],forwardrg[i],forwardpk[i]);
   }

   if(nvr > 2)
   {
      ix1 = idx[nvr-1];                            // b[0] = x[n-1]*x[n-2]
      ix2 = idx[nvr-2];
      CPU_dbl5_product(deg,
         inputtb[ix1],inputix[ix1],inputmi[ix1],inputrg[ix1],inputpk[ix1],
         inputtb[ix2],inputix[ix2],inputmi[ix2],inputrg[ix2],inputpk[ix2],
         backwardtb[0],backwardix[0],backwardmi[0],
         backwardrg[0],backwardpk[0]);

      for(int i=1; i<nvr-2; i++)
      {                                            // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];
         CPU_dbl5_product
            (deg,backwardtb[i-1],backwardix[i-1],backwardmi[i-1],
                 backwardrg[i-1],backwardpk[i-1],
                 inputtb[ix2],inputix[ix2],inputmi[ix2],
                 inputrg[ix2],inputpk[ix2],
                 backwardtb[i],backwardix[i],backwardmi[i],
                 backwardrg[i],backwardpk[i]);
      }
                                                   // b[n-3] = cff*b[n-3]
      CPU_dbl5_product
         (deg,backwardtb[nvr-3],backwardix[nvr-3],backwardmi[nvr-3],
              backwardrg[nvr-3],backwardpk[nvr-3],
              cfftb,cffix,cffmi,cffrg,cffpk,
              crosstb[0],crossix[0],crossmi[0],crossrg[0],crosspk[0]);
      // cross[0] is work space, cannot write into backward[nvr-3]
      for(int i=0; i<=deg; i++)
      {
         backwardtb[nvr-3][i] = crosstb[0][i];
         backwardix[nvr-3][i] = crossix[0][i];
         backwardmi[nvr-3][i] = crossmi[0][i];
         backwardrg[nvr-3][i] = crossrg[0][i];
         backwardpk[nvr-3][i] = crosspk[0][i];
      }
      if(nvr == 3)
      {                                            // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_dbl5_product(deg,
            forwardtb[0],forwardix[0],forwardmi[0],forwardrg[0],forwardpk[0],
            inputtb[ix2],inputix[ix2],inputmi[ix2],inputrg[ix2],inputpk[ix2],
            crosstb[0],crossix[0],crossmi[0],crossrg[0],crosspk[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                         // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_dbl5_product
               (deg,forwardtb[i],forwardix[i],forwardmi[i],
                    forwardrg[i],forwardpk[i],
                    backwardtb[ix2],backwardix[ix2],backwardmi[ix2],
                    backwardrg[ix2],backwardpk[ix2],
                    crosstb[i],crossix[i],crossmi[i],crossrg[i],crosspk[i]);
         }
         ix2 = idx[nvr-1];                         // c[n-3] = f[n-3]*x[n-1]
         CPU_dbl5_product
            (deg,forwardtb[nvr-3],forwardix[nvr-3],forwardmi[nvr-3],
                 forwardrg[nvr-3],forwardpk[nvr-3],
                 inputtb[ix2],inputix[ix2],inputmi[ix2],
                 inputrg[ix2],inputpk[ix2],
                 crosstb[nvr-3],crossix[nvr-3],crossmi[nvr-3],
                 crossrg[nvr-3],crosspk[nvr-3]);
      }
   }
}

void CPU_cmplx5_speel
 ( int nvr, int deg, int *idx, double *cffretb, double *cffreix,
   double *cffremi, double *cffrerg, double *cffrepk, double *cffimtb,
   double *cffimix, double *cffimmi, double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **forwardretb, double **forwardreix, double **forwardremi,
   double **forwardrerg, double **forwardrepk,
   double **forwardimtb, double **forwardimix, double **forwardimmi,
   double **forwardimrg, double **forwardimpk,
   double **backwardretb, double **backwardreix, double **backwardremi,
   double **backwardrerg, double **backwardrepk,
   double **backwardimtb, double **backwardimix, double **backwardimmi,
   double **backwardimrg, double **backwardimpk,
   double **crossretb, double **crossreix, double **crossremi,
   double **crossrerg, double **crossrepk,
   double **crossimtb, double **crossimix, double **crossimmi,
   double **crossimrg, double **crossimpk )
{
   int ix1 = idx[0];
   int ix2;
                                                           // f[0] = cff*x[0]
   CPU_cmplx5_product(deg,
      cffretb,cffreix,cffremi,cffrerg,cffrepk,
      cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
      inputretb[ix1],inputreix[ix1],inputremi[ix1],
      inputrerg[ix1],inputrepk[ix1],
      inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
      inputimrg[ix1],inputimpk[ix1],
      forwardretb[0],forwardreix[0],forwardremi[0],
      forwardrerg[0],forwardrepk[0],
      forwardimtb[0],forwardimix[0],forwardimmi[0],
      forwardimrg[0],forwardimpk[0]);
   for(int i=1; i<nvr; i++)
   {                                                    // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_cmplx5_product(deg,
         forwardretb[i-1],forwardreix[i-1],forwardremi[i-1],
         forwardrerg[i-1],forwardrepk[i-1],
         forwardimtb[i-1],forwardimix[i-1],forwardimmi[i-1],
         forwardimrg[i-1],forwardimpk[i-1],
         inputretb[ix2],inputreix[ix2],inputremi[ix2],
         inputrerg[ix2],inputrepk[ix2],
         inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
         inputimrg[ix2],inputimpk[ix2],
         forwardretb[i],forwardreix[i],forwardremi[i],
         forwardrerg[i],forwardrepk[i],
         forwardimtb[i],forwardimix[i],forwardimmi[i],
         forwardimrg[i],forwardimpk[i]);
   }
   if(nvr > 2)
   {                                                  // b[0] = x[n-1]*x[n-2]
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];
      CPU_cmplx5_product(deg,
         inputretb[ix1],inputreix[ix1],inputremi[ix1],
         inputrerg[ix1],inputrepk[ix1],
         inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
         inputimrg[ix1],inputimpk[ix1],
         inputretb[ix2],inputreix[ix2],inputremi[ix2],
         inputrerg[ix2],inputrepk[ix2],
         inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
         inputimrg[ix2],inputimpk[ix2],
         backwardretb[0],backwardreix[0],backwardremi[0],
         backwardrerg[0],backwardrepk[0],
         backwardimtb[0],backwardimix[0],backwardimmi[0],
         backwardimrg[0],backwardimpk[0]);
      for(int i=1; i<nvr-2; i++)
      {                                             // b[i] = b[i-1]*x[x-2-i]
         ix2 = idx[nvr-2-i];
         CPU_cmplx5_product(deg,
            backwardretb[i-1],backwardreix[i-1],backwardremi[i-1],
            backwardrerg[i-1],backwardrepk[i-1],
            backwardimtb[i-1],backwardimix[i-1],backwardimmi[i-1],
            backwardimrg[i-1],backwardimpk[i-1],
            inputretb[ix2],inputreix[ix2],inputremi[ix2],
            inputrerg[ix2],inputrepk[ix2],
            inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
            inputimrg[ix2],inputimpk[ix2],
            backwardretb[i],backwardreix[i],backwardremi[i],
            backwardrerg[i],backwardrepk[i],
            backwardimtb[i],backwardimix[i],backwardimmi[i],
            backwardimrg[i],backwardimpk[i]);
      }
                                                       // b[n-3] = b[n-3]*cff
      CPU_cmplx5_product(deg,
         backwardretb[nvr-3],backwardreix[nvr-3],backwardremi[nvr-3],
         backwardrerg[nvr-3],backwardrepk[nvr-3],
         backwardimtb[nvr-3],backwardimix[nvr-3],backwardimmi[nvr-3],
         backwardimrg[nvr-3],backwardimpk[nvr-3],
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         crossretb[0],crossreix[0],crossremi[0],crossrerg[0],crossrepk[0],
         crossimtb[0],crossimix[0],crossimmi[0],crossimrg[0],crossimpk[0]);
                                                       // cross is work space
      for(int i=0; i<=deg; i++)
      {
         backwardretb[nvr-3][i] = crossretb[0][i];
         backwardreix[nvr-3][i] = crossreix[0][i];
         backwardremi[nvr-3][i] = crossremi[0][i];
         backwardrerg[nvr-3][i] = crossrerg[0][i];
         backwardrepk[nvr-3][i] = crossrepk[0][i];
         backwardimtb[nvr-3][i] = crossimtb[0][i];
         backwardimix[nvr-3][i] = crossimix[0][i];
         backwardimmi[nvr-3][i] = crossimmi[0][i];
         backwardimrg[nvr-3][i] = crossimrg[0][i];
         backwardimpk[nvr-3][i] = crossimpk[0][i];
      }
      if(nvr == 3)
      {                                                   // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_cmplx5_product(deg,
            forwardretb[0],forwardreix[0],forwardremi[0],
            forwardrerg[0],forwardrepk[0],
            forwardimtb[0],forwardimix[0],forwardimmi[0],
            forwardimrg[0],forwardimpk[0],
            inputretb[ix2],inputreix[ix2],inputremi[ix2],
            inputrerg[ix2],inputrepk[ix2],
            inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
            inputimrg[ix2],inputimpk[ix2],
            crossretb[0],crossreix[0],crossremi[0],crossrerg[0],crossrepk[0],
            crossimtb[0],crossimix[0],crossimmi[0],crossimrg[0],crossimpk[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                            // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_cmplx5_product(deg,
               forwardretb[i],forwardreix[i],forwardremi[i],
               forwardrerg[i],forwardrepk[i],
               forwardimtb[i],forwardimix[i],forwardimmi[i],
               forwardimrg[i],forwardimpk[i],
               backwardretb[ix2],backwardreix[ix2],backwardremi[ix2],
               backwardrerg[ix2],backwardrepk[ix2],
               backwardimtb[ix2],backwardimix[ix2],backwardimmi[ix2],
               backwardimrg[ix2],backwardimpk[ix2],
               crossretb[i],crossreix[i],crossremi[i],
               crossrerg[i],crossrepk[i],
               crossimtb[i],crossimix[i],crossimmi[i],
               crossimrg[i],crossimpk[i]);
         }
                                                    // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_cmplx5_product(deg,
            forwardretb[nvr-3],forwardreix[nvr-3],forwardremi[nvr-3],
            forwardrerg[nvr-3],forwardrepk[nvr-3],
            forwardimtb[nvr-3],forwardimix[nvr-3],forwardimmi[nvr-3],
            forwardimrg[nvr-3],forwardimpk[nvr-3],
            inputretb[ix2],inputreix[ix2],inputremi[ix2],
            inputrerg[ix2],inputrepk[ix2],
            inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
            inputimrg[ix2],inputimpk[ix2],
            crossretb[nvr-3],crossreix[nvr-3],crossremi[nvr-3],
            crossrerg[nvr-3],crossrepk[nvr-3],
            crossimtb[nvr-3],crossimix[nvr-3],crossimmi[nvr-3],
            crossimrg[nvr-3],crossimpk[nvr-3]);
      }
   }
}

void CPU_dbl5_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi, double **inputrg,
   double **inputpk, double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_dbl5_product(deg,
         inputtb[ix],inputix[ix],inputmi[ix],inputrg[ix],inputpk[ix],
         cfftb,cffix,cffmi,cffrg,cffpk,outputtb[dim],outputix[dim],
         outputmi[dim],outputrg[dim],outputpk[dim]);

      for(int i=0; i<=deg; i++)
      {
         outputtb[ix][i] = cfftb[i];
         outputix[ix][i] = cffix[i];
         outputmi[ix][i] = cffmi[i];
         outputrg[ix][i] = cffrg[i];
         outputpk[ix][i] = cffpk[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_dbl5_product(deg,
         inputtb[ix1],inputix[ix1],inputmi[ix1],inputrg[ix1],inputpk[ix1],
         inputtb[ix2],inputix[ix2],inputmi[ix2],inputrg[ix2],inputpk[ix2],
         outputtb[dim],outputix[dim],outputmi[dim],
         outputrg[dim],outputpk[dim]);
      CPU_dbl5_product(deg,
         outputtb[dim],outputix[dim],outputmi[dim],outputrg[dim],outputpk[dim],
         cfftb,cffix,cffmi,cffrg,cffpk,
         outputtb[dim],outputix[dim],outputmi[dim],
         outputrg[dim],outputpk[dim]);

      CPU_dbl5_product(deg,cfftb,cffix,cffmi,cffrg,cffpk,
         inputtb[ix1],inputix[ix1],inputmi[ix1],inputrg[ix1],inputpk[ix1],
         outputtb[ix2],outputix[ix2],outputmi[ix2],
         outputrg[ix2],outputpk[ix2]);
      CPU_dbl5_product(deg,cfftb,cffix,cffmi,cffrg,cffpk,
         inputtb[ix2],inputix[ix2],inputmi[ix2],inputrg[ix2],inputpk[ix2],
         outputtb[ix1],outputix[ix1],outputmi[ix1],
         outputrg[ix1],outputpk[ix1]);
   }
   else
   {
      double **forwardtb = new double*[nvr];
      double **forwardix = new double*[nvr];
      double **forwardmi = new double*[nvr];
      double **forwardrg = new double*[nvr];
      double **forwardpk = new double*[nvr];
      double **backwardtb = new double*[nvr-2];
      double **backwardix = new double*[nvr-2];
      double **backwardmi = new double*[nvr-2];
      double **backwardrg = new double*[nvr-2];
      double **backwardpk = new double*[nvr-2];
      double **crosstb = new double*[nvr-2];
      double **crossix = new double*[nvr-2];
      double **crossmi = new double*[nvr-2];
      double **crossrg = new double*[nvr-2];
      double **crosspk = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardtb[i] = new double[deg+1];
         forwardix[i] = new double[deg+1];
         forwardmi[i] = new double[deg+1];
         forwardrg[i] = new double[deg+1];
         forwardpk[i] = new double[deg+1];
         backwardtb[i] = new double[deg+1];
         backwardix[i] = new double[deg+1];
         backwardmi[i] = new double[deg+1];
         backwardrg[i] = new double[deg+1];
         backwardpk[i] = new double[deg+1];
         crosstb[i] = new double[deg+1];
         crossix[i] = new double[deg+1];
         crossmi[i] = new double[deg+1];
         crossrg[i] = new double[deg+1];
         crosspk[i] = new double[deg+1];
      }
      forwardtb[nvr-2] = new double[deg+1];
      forwardix[nvr-2] = new double[deg+1];
      forwardmi[nvr-2] = new double[deg+1];
      forwardrg[nvr-2] = new double[deg+1];
      forwardpk[nvr-2] = new double[deg+1];
      forwardtb[nvr-1] = new double[deg+1];
      forwardix[nvr-1] = new double[deg+1];
      forwardmi[nvr-1] = new double[deg+1];
      forwardrg[nvr-1] = new double[deg+1];
      forwardpk[nvr-1] = new double[deg+1];

      CPU_dbl5_speel(nvr,deg,idx,
         cfftb,cffix,cffmi,cffrg,cffpk,
         inputtb,inputix,inputmi,inputrg,inputpk,
         forwardtb,forwardix,forwardmi,forwardrg,forwardpk,
         backwardtb,backwardix,backwardmi,backwardrg,backwardpk,
         crosstb,crossix,crossmi,crossrg,crosspk);

      for(int i=0; i<deg+1; i++)     // assign the value of the monomial
      {
         outputtb[dim][i] = forwardtb[nvr-1][i];
         outputix[dim][i] = forwardix[nvr-1][i];
         outputmi[dim][i] = forwardmi[nvr-1][i];
         outputrg[dim][i] = forwardrg[nvr-1][i];
         outputpk[dim][i] = forwardpk[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];       // derivative with respect to x[n-1]
         for(int i=0; i<deg+1; i++)
         {
            outputtb[ix][i] = forwardtb[nvr-2][i];
            outputix[ix][i] = forwardix[nvr-2][i];
            outputmi[ix][i] = forwardmi[nvr-2][i];
            outputrg[ix][i] = forwardrg[nvr-2][i];
            outputpk[ix][i] = forwardpk[nvr-2][i];
         }
         ix = idx[0];               // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputtb[ix][i] = backwardtb[nvr-3][i];
            outputix[ix][i] = backwardix[nvr-3][i];
            outputmi[ix][i] = backwardmi[nvr-3][i];
            outputrg[ix][i] = backwardrg[nvr-3][i];
            outputpk[ix][i] = backwardpk[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {                          // derivative with respect to x[k]
            ix = idx[k];
            for(int i=0; i<deg+1; i++)
            {
               outputtb[ix][i] = crosstb[k-1][i];
               outputix[ix][i] = crossix[k-1][i];
               outputmi[ix][i] = crossmi[k-1][i];
               outputrg[ix][i] = crossrg[k-1][i];
               outputpk[ix][i] = crosspk[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardtb[i]); free(backwardtb[i]); free(crosstb[i]);
         free(forwardix[i]); free(backwardix[i]); free(crossix[i]);
         free(forwardmi[i]); free(backwardmi[i]); free(crossmi[i]);
         free(forwardrg[i]); free(backwardrg[i]); free(crossrg[i]);
         free(forwardpk[i]); free(backwardpk[i]); free(crosspk[i]);
      }
      free(forwardtb[nvr-2]); free(forwardtb[nvr-1]);
      free(forwardix[nvr-2]); free(forwardix[nvr-1]);
      free(forwardmi[nvr-2]); free(forwardmi[nvr-1]);
      free(forwardrg[nvr-2]); free(forwardrg[nvr-1]);
      free(forwardpk[nvr-2]); free(forwardpk[nvr-1]);
      free(forwardtb); free(backwardtb); free(crosstb);
      free(forwardix); free(backwardix); free(crossix);
      free(forwardmi); free(backwardmi); free(crossmi);
      free(forwardrg); free(backwardrg); free(crossrg);
      free(forwardpk); free(backwardpk); free(crosspk);
   }
}

void CPU_cmplx5_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cffretb, double *cffreix,
   double *cffremi, double *cffrerg, double *cffrepk, double *cffimtb,
   double *cffimix, double *cffimmi, double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_cmplx5_product(deg,
         inputretb[ix],inputreix[ix],inputremi[ix],
         inputrerg[ix],inputrepk[ix],
         inputimtb[ix],inputimix[ix],inputimmi[ix],
         inputimrg[ix],inputimpk[ix],
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         outputretb[dim],outputreix[dim],outputremi[dim],
         outputrerg[dim],outputrepk[dim],
         outputimtb[dim],outputimix[dim],outputimmi[dim],
         outputimrg[dim],outputimpk[dim]);

      for(int i=0; i<=deg; i++) 
      {
         outputretb[ix][i] = cffretb[i]; outputreix[ix][i] = cffreix[i];
         outputremi[ix][i] = cffremi[i]; outputrerg[ix][i] = cffrerg[i];
         outputrepk[ix][i] = cffrepk[i];
         outputimtb[ix][i] = cffimtb[i]; outputimix[ix][i] = cffimix[i];
         outputimmi[ix][i] = cffimmi[i]; outputimrg[ix][i] = cffimrg[i];
         outputimpk[ix][i] = cffimpk[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0];
      int ix2 = idx[1];

      CPU_cmplx5_product(deg,
         inputretb[ix1],inputreix[ix1],inputremi[ix1],
         inputrerg[ix1],inputrepk[ix1],
         inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
         inputimrg[ix1],inputimpk[ix1],
         inputretb[ix2],inputreix[ix2],inputremi[ix2],
         inputrerg[ix2],inputrepk[ix2],
         inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
         inputimrg[ix2],inputimpk[ix2],
         outputretb[dim],outputreix[dim],outputremi[dim],
         outputrerg[dim],outputrepk[dim],
         outputimtb[dim],outputimix[dim],outputimmi[dim],
         outputimrg[dim],outputimpk[dim]);
      CPU_cmplx5_product(deg,
         outputretb[dim],outputreix[dim],outputremi[dim],
         outputrerg[dim],outputrepk[dim],
         outputimtb[dim],outputimix[dim],outputimmi[dim],
         outputimrg[dim],outputimpk[dim],
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         outputretb[dim],outputreix[dim],outputremi[dim],
         outputrerg[dim],outputrepk[dim],
         outputimtb[dim],outputimix[dim],outputimmi[dim],
         outputimrg[dim],outputimpk[dim]);

      CPU_cmplx5_product(deg,
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         inputretb[ix1],inputreix[ix1],inputremi[ix1],
         inputrerg[ix1],inputrepk[ix1],
         inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
         inputimrg[ix1],inputimpk[ix1],
         outputretb[ix2],outputreix[ix2],outputremi[ix2],
         outputrerg[ix2],outputrepk[ix2],
         outputimtb[ix2],outputimix[ix2],outputimmi[ix2],
         outputimrg[ix2],outputimpk[ix2]);
      CPU_cmplx5_product(deg,
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         inputretb[ix2],inputreix[ix2],inputremi[ix2],
         inputrerg[ix2],inputrepk[ix2],
         inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
         inputimrg[ix1],inputimpk[ix1],
         outputretb[ix1],outputreix[ix1],outputremi[ix1],
         outputrerg[ix1],outputrepk[ix1],
         outputimtb[ix1],outputimix[ix1],outputimmi[ix1],
         outputimrg[ix1],outputimpk[ix1]);
   }
   else
   {
      double **forwardretb = new double*[nvr];
      double **forwardreix = new double*[nvr];
      double **forwardremi = new double*[nvr];
      double **forwardrerg = new double*[nvr];
      double **forwardrepk = new double*[nvr];
      double **forwardimtb = new double*[nvr];
      double **forwardimix = new double*[nvr];
      double **forwardimmi = new double*[nvr];
      double **forwardimrg = new double*[nvr];
      double **forwardimpk = new double*[nvr];
      double **backwardretb = new double*[nvr-2];
      double **backwardreix = new double*[nvr-2];
      double **backwardremi = new double*[nvr-2];
      double **backwardrerg = new double*[nvr-2];
      double **backwardrepk = new double*[nvr-2];
      double **backwardimtb = new double*[nvr-2];
      double **backwardimix = new double*[nvr-2];
      double **backwardimmi = new double*[nvr-2];
      double **backwardimrg = new double*[nvr-2];
      double **backwardimpk = new double*[nvr-2];
      double **crossretb = new double*[nvr-2];
      double **crossreix = new double*[nvr-2];
      double **crossremi = new double*[nvr-2];
      double **crossrerg = new double*[nvr-2];
      double **crossrepk = new double*[nvr-2];
      double **crossimtb = new double*[nvr-2];
      double **crossimix = new double*[nvr-2];
      double **crossimmi = new double*[nvr-2];
      double **crossimrg = new double*[nvr-2];
      double **crossimpk = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardretb[i] = new double[deg+1];
         forwardreix[i] = new double[deg+1];
         forwardremi[i] = new double[deg+1];
         forwardrerg[i] = new double[deg+1];
         forwardrepk[i] = new double[deg+1];
         forwardimtb[i] = new double[deg+1];
         forwardimix[i] = new double[deg+1];
         forwardimmi[i] = new double[deg+1];
         forwardimrg[i] = new double[deg+1];
         forwardimpk[i] = new double[deg+1];
         backwardretb[i] = new double[deg+1];
         backwardreix[i] = new double[deg+1];
         backwardremi[i] = new double[deg+1];
         backwardrerg[i] = new double[deg+1];
         backwardrepk[i] = new double[deg+1];
         backwardimtb[i] = new double[deg+1];
         backwardimix[i] = new double[deg+1];
         backwardimmi[i] = new double[deg+1];
         backwardimrg[i] = new double[deg+1];
         backwardimpk[i] = new double[deg+1];
         crossretb[i] = new double[deg+1];
         crossreix[i] = new double[deg+1];
         crossremi[i] = new double[deg+1];
         crossrerg[i] = new double[deg+1];
         crossrepk[i] = new double[deg+1];
         crossimtb[i] = new double[deg+1];
         crossimix[i] = new double[deg+1];
         crossimmi[i] = new double[deg+1];
         crossimrg[i] = new double[deg+1];
         crossimpk[i] = new double[deg+1];
      }
      forwardretb[nvr-2] = new double[deg+1];
      forwardreix[nvr-2] = new double[deg+1];
      forwardremi[nvr-2] = new double[deg+1];
      forwardrerg[nvr-2] = new double[deg+1];
      forwardrepk[nvr-2] = new double[deg+1];
      forwardimtb[nvr-2] = new double[deg+1];
      forwardimix[nvr-2] = new double[deg+1];
      forwardimmi[nvr-2] = new double[deg+1];
      forwardimrg[nvr-2] = new double[deg+1];
      forwardimpk[nvr-2] = new double[deg+1];
      forwardretb[nvr-1] = new double[deg+1];
      forwardreix[nvr-1] = new double[deg+1];
      forwardremi[nvr-1] = new double[deg+1];
      forwardrerg[nvr-1] = new double[deg+1];
      forwardrepk[nvr-1] = new double[deg+1];
      forwardimtb[nvr-1] = new double[deg+1];
      forwardimix[nvr-1] = new double[deg+1];
      forwardimmi[nvr-1] = new double[deg+1];
      forwardimrg[nvr-1] = new double[deg+1];
      forwardimpk[nvr-1] = new double[deg+1];

      CPU_cmplx5_speel(nvr,deg,idx,
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         inputretb,inputreix,inputremi,inputrerg,inputrepk,
         inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
         forwardretb,forwardreix,forwardremi,forwardrerg,forwardrepk,
         forwardimtb,forwardimix,forwardimmi,forwardimrg,forwardimpk,
         backwardretb,backwardreix,backwardremi,backwardrerg,backwardrepk,
         backwardimtb,backwardimix,backwardimmi,backwardimrg,backwardimpk,
         crossretb,crossreix,crossremi,crossrerg,crossrepk,
         crossimtb,crossimix,crossimmi,crossimrg,crossimpk);

      for(int i=0; i<deg+1; i++)          // assign value of the monomial
      {
         outputretb[dim][i] = forwardretb[nvr-1][i];
         outputreix[dim][i] = forwardreix[nvr-1][i];
         outputremi[dim][i] = forwardremi[nvr-1][i];
         outputrerg[dim][i] = forwardrerg[nvr-1][i];
         outputrepk[dim][i] = forwardrepk[nvr-1][i];
         outputimtb[dim][i] = forwardimtb[nvr-1][i];
         outputimix[dim][i] = forwardimix[nvr-1][i];
         outputimmi[dim][i] = forwardimmi[nvr-1][i];
         outputimrg[dim][i] = forwardimrg[nvr-1][i];
         outputimpk[dim][i] = forwardimpk[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];        // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++)
         {
            outputretb[ix][i] = forwardretb[nvr-2][i];
            outputreix[ix][i] = forwardreix[nvr-2][i];
            outputremi[ix][i] = forwardremi[nvr-2][i];
            outputrerg[ix][i] = forwardrerg[nvr-2][i];
            outputrepk[ix][i] = forwardrepk[nvr-2][i];
            outputimtb[ix][i] = forwardimtb[nvr-2][i];
            outputimix[ix][i] = forwardimix[nvr-2][i];
            outputimmi[ix][i] = forwardimmi[nvr-2][i];
            outputimrg[ix][i] = forwardimrg[nvr-2][i];
            outputimpk[ix][i] = forwardimpk[nvr-2][i];
         }

         ix = idx[0];                  // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputretb[ix][i] = backwardretb[nvr-3][i];
            outputreix[ix][i] = backwardreix[nvr-3][i];
            outputremi[ix][i] = backwardremi[nvr-3][i];
            outputrerg[ix][i] = backwardrerg[nvr-3][i];
            outputrepk[ix][i] = backwardrepk[nvr-3][i];
            outputimtb[ix][i] = backwardimtb[nvr-3][i];
            outputimix[ix][i] = backwardimix[nvr-3][i];
            outputimmi[ix][i] = backwardimmi[nvr-3][i];
            outputimrg[ix][i] = backwardimrg[nvr-3][i];
            outputimpk[ix][i] = backwardimpk[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];               // derivative with respect to x[k]
            for(int i=0; i<deg+1; i++)
            {
               outputretb[ix][i] = crossretb[k-1][i];
               outputreix[ix][i] = crossreix[k-1][i];
               outputremi[ix][i] = crossremi[k-1][i];
               outputrerg[ix][i] = crossrerg[k-1][i];
               outputrepk[ix][i] = crossrepk[k-1][i];
               outputimtb[ix][i] = crossimtb[k-1][i];
               outputimix[ix][i] = crossimix[k-1][i];
               outputimmi[ix][i] = crossimmi[k-1][i];
               outputimrg[ix][i] = crossimrg[k-1][i];
               outputimpk[ix][i] = crossimpk[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardretb[i]);  free(forwardreix[i]);  free(forwardremi[i]);
         free(forwardrerg[i]);  free(forwardrepk[i]);
         free(forwardimtb[i]);  free(forwardimix[i]);  free(forwardimmi[i]);
         free(forwardimrg[i]);  free(forwardimpk[i]);
         free(backwardretb[i]); free(backwardreix[i]); free(backwardremi[i]);
         free(backwardrerg[i]); free(backwardrepk[i]);
         free(backwardimtb[i]); free(backwardimix[i]); free(backwardimmi[i]);
         free(backwardimrg[i]); free(backwardimpk[i]);
         free(crossretb[i]);    free(crossreix[i]);    free(crossremi[i]);
         free(crossrerg[i]);    free(crossrepk[i]);
         free(crossimtb[i]);    free(crossimix[i]);    free(crossimmi[i]);
         free(crossimrg[i]);    free(crossimpk[i]);
      }
      free(forwardretb[nvr-2]); free(forwardreix[nvr-2]);
      free(forwardremi[nvr-2]); free(forwardrerg[nvr-2]);
      free(forwardrepk[nvr-2]);
      free(forwardimtb[nvr-2]); free(forwardimix[nvr-2]);
      free(forwardimmi[nvr-2]); free(forwardimrg[nvr-2]);
      free(forwardimpk[nvr-2]);
      free(forwardretb[nvr-1]); free(forwardreix[nvr-1]);
      free(forwardremi[nvr-1]); free(forwardrerg[nvr-1]);
      free(forwardrepk[nvr-1]);
      free(forwardimtb[nvr-1]); free(forwardimix[nvr-1]);
      free(forwardimmi[nvr-1]); free(forwardimrg[nvr-1]);
      free(forwardimpk[nvr-1]);
      free(forwardretb);  free(forwardreix);  free(forwardremi);
      free(forwardrerg);  free(forwardrepk);
      free(forwardimtb);  free(forwardimix);  free(forwardimmi);
      free(forwardimrg);  free(forwardimpk);
      free(backwardretb); free(backwardreix); free(backwardremi);
      free(backwardrerg); free(backwardrepk);
      free(backwardimtb); free(backwardimix); free(backwardimmi);
      free(backwardimrg); free(backwardimpk);
      free(crossretb);    free(crossreix);    free(crossremi);
      free(crossrerg);    free(crossrepk);
      free(crossimtb);    free(crossimix);    free(crossimmi);
      free(crossimrg);    free(crossimpk);
   }
}
