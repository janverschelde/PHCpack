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
}
