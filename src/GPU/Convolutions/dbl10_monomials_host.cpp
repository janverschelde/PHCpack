/* The file dbl10_monomials_host.cpp defines functions specified
 * in dbl10_monomials_host.h. */

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
#include "dbl10_convolutions_host.h"
#include "dbl10_monomials_host.h"

void CPU_dbl10_speel
 ( int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **forwardrtb, double **forwardrix, double **forwardrmi,
   double **forwardrrg, double **forwardrpk,
   double **forwardltb, double **forwardlix, double **forwardlmi,
   double **forwardlrg, double **forwardlpk,
   double **backwardrtb, double **backwardrix, double **backwardrmi,
   double **backwardrrg, double **backwardrpk,
   double **backwardltb, double **backwardlix, double **backwardlmi,
   double **backwardlrg, double **backwardlpk,
   double **crossrtb, double **crossrix, double **crossrmi,
   double **crossrrg, double **crossrpk,
   double **crossltb, double **crosslix, double **crosslmi,
   double **crosslrg, double **crosslpk )
{
   int ix1 = idx[0];
   int ix2;
                                                   // f[0] = cff*x[0] 
   CPU_dbl10_product(deg,
      cffrtb,cffrix,cffrmi,cffrrg,cffrpk,cffltb,cfflix,cfflmi,cfflrg,cfflpk,
      inputrtb[ix1],inputrix[ix1],inputrmi[ix1],inputrrg[ix1],inputrpk[ix1],
      inputltb[ix1],inputlix[ix1],inputlmi[ix1],inputlrg[ix1],inputlpk[ix1],
      forwardrtb[0],forwardrix[0],forwardrmi[0],forwardrrg[0],forwardrpk[0],
      forwardltb[0],forwardlix[0],forwardlmi[0],forwardlrg[0],forwardlpk[0]);

   for(int i=1; i<nvr; i++)
   {                                               // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_dbl10_product(deg,
         forwardrtb[i-1],forwardrix[i-1],forwardrmi[i-1],
         forwardrrg[i-1],forwardrpk[i-1],
         forwardltb[i-1],forwardlix[i-1],forwardlmi[i-1],
         forwardlrg[i-1],forwardlpk[i-1],
         inputrtb[ix2],inputrix[ix2],inputrmi[ix2],inputrrg[ix2],inputrpk[ix2],
         inputltb[ix2],inputlix[ix2],inputlmi[ix2],inputlrg[ix2],inputlpk[ix2],
         forwardrtb[i],forwardrix[i],forwardrmi[i],forwardrrg[i],forwardrpk[i],
         forwardltb[i],forwardlix[i],forwardlmi[i],
         forwardlrg[i],forwardlpk[i]);
   }

   if(nvr > 2)
   {
      ix1 = idx[nvr-1];                            // b[0] = x[n-1]*x[n-2]
      ix2 = idx[nvr-2];
      CPU_dbl10_product(deg,
         inputrtb[ix1],inputrix[ix1],inputrmi[ix1],inputrrg[ix1],inputrpk[ix1],
         inputltb[ix1],inputlix[ix1],inputlmi[ix1],inputlrg[ix1],inputlpk[ix1],
         inputrtb[ix2],inputrix[ix2],inputrmi[ix2],inputrrg[ix2],inputrpk[ix2],
         inputltb[ix2],inputlix[ix2],inputlmi[ix2],inputlrg[ix2],inputlpk[ix2],
         backwardrtb[0],backwardrix[0],backwardrmi[0],
         backwardrrg[0],backwardrpk[0],
         backwardltb[0],backwardrix[0],backwardlmi[0],
         backwardlrg[0],backwardrpk[0]);

      for(int i=1; i<nvr-2; i++)
      {                                            // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];
         CPU_dbl10_product
            (deg,backwardrtb[i-1],backwardrix[i-1],backwardrmi[i-1],
                 backwardrrg[i-1],backwardrpk[i-1],
                 backwardltb[i-1],backwardlix[i-1],backwardlmi[i-1],
                 backwardlrg[i-1],backwardlpk[i-1],
                 inputrtb[ix2],inputrix[ix2],inputrmi[ix2],
                 inputrrg[ix2],inputrpk[ix2],
                 inputltb[ix2],inputlix[ix2],inputlmi[ix2],
                 inputlrg[ix2],inputlpk[ix2],
                 backwardrtb[i],backwardrix[i],backwardrmi[i],
                 backwardrrg[i],backwardrpk[i],
                 backwardltb[i],backwardlix[i],backwardlmi[i],
                 backwardlrg[i],backwardlpk[i]);
      }
                                                   // b[n-3] = cff*b[n-3]
      CPU_dbl10_product
         (deg,backwardrtb[nvr-3],backwardrix[nvr-3],backwardrmi[nvr-3],
              backwardrrg[nvr-3],backwardrpk[nvr-3],
              backwardltb[nvr-3],backwardlix[nvr-3],backwardlmi[nvr-3],
              backwardlrg[nvr-3],backwardlpk[nvr-3],
              cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
              cffltb,cfflix,cfflmi,cfflrg,cfflpk,
              crossrtb[0],crossrix[0],crossrmi[0],crossrrg[0],crossrpk[0],
              crossltb[0],crosslix[0],crosslmi[0],crosslrg[0],crosslpk[0]);
      // cross[0] is work space, cannot write into backward[nvr-3]
      for(int i=0; i<=deg; i++)
      {
         backwardrtb[nvr-3][i] = crossrtb[0][i];
         backwardrix[nvr-3][i] = crossrix[0][i];
         backwardrmi[nvr-3][i] = crossrmi[0][i];
         backwardrrg[nvr-3][i] = crossrrg[0][i];
         backwardrpk[nvr-3][i] = crossrpk[0][i];
         backwardltb[nvr-3][i] = crossltb[0][i];
         backwardlix[nvr-3][i] = crosslix[0][i];
         backwardlmi[nvr-3][i] = crosslmi[0][i];
         backwardlrg[nvr-3][i] = crosslrg[0][i];
         backwardlpk[nvr-3][i] = crosslpk[0][i];
      }
      if(nvr == 3)
      {                                            // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_dbl10_product(deg,
            forwardrtb[0],forwardrix[0],forwardrmi[0],
            forwardrrg[0],forwardrpk[0],
            forwardltb[0],forwardlix[0],forwardlmi[0],
            forwardlrg[0],forwardlpk[0],
            inputrtb[ix2],inputrix[ix2],inputrmi[ix2],
            inputrrg[ix2],inputrpk[ix2],
            inputltb[ix2],inputlix[ix2],inputlmi[ix2],
            inputlrg[ix2],inputlpk[ix2],
            crossrtb[0],crossrix[0],crossrmi[0],crossrrg[0],crossrpk[0],
            crossltb[0],crosslix[0],crosslmi[0],crosslrg[0],crosslpk[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                         // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_dbl10_product
               (deg,forwardrtb[i],forwardrix[i],forwardrmi[i],
                    forwardrrg[i],forwardrpk[i],
                    forwardltb[i],forwardlix[i],forwardlmi[i],
                    forwardlrg[i],forwardlpk[i],
                    backwardrtb[ix2],backwardrix[ix2],backwardrmi[ix2],
                    backwardrrg[ix2],backwardrpk[ix2],
                    backwardltb[ix2],backwardlix[ix2],backwardlmi[ix2],
                    backwardlrg[ix2],backwardlpk[ix2],
                    crossrtb[i],crossrix[i],crossrmi[i],
                    crossrrg[i],crossrpk[i],
                    crossltb[i],crosslix[i],crosslmi[i],
                    crosslrg[i],crosslpk[i]);
         }
         ix2 = idx[nvr-1];                         // c[n-3] = f[n-3]*x[n-1]
         CPU_dbl10_product
            (deg,forwardrtb[nvr-3],forwardrix[nvr-3],forwardrmi[nvr-3],
                 forwardrrg[nvr-3],forwardrpk[nvr-3],
                 forwardltb[nvr-3],forwardlix[nvr-3],forwardlmi[nvr-3],
                 forwardlrg[nvr-3],forwardlpk[nvr-3],
                 inputrtb[ix2],inputrix[ix2],inputrmi[ix2],
                 inputrrg[ix2],inputrpk[ix2],
                 inputltb[ix2],inputlix[ix2],inputlmi[ix2],
                 inputlrg[ix2],inputlpk[ix2],
                 crossrtb[nvr-3],crossrix[nvr-3],crossrmi[nvr-3],
                 crossrrg[nvr-3],crossrpk[nvr-3],
                 crossltb[nvr-3],crosslix[nvr-3],crosslmi[nvr-3],
                 crosslrg[nvr-3],crosslpk[nvr-3]);
      }
   }
}

void CPU_cmplx10_speel
 ( int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double **inputrertb, double **inputrerix, double **inputrermi,
   double **inputrerrg, double **inputrerpk,
   double **inputreltb, double **inputrelix, double **inputrelmi,
   double **inputrelrg, double **inputrelpk,
   double **inputimrtb, double **inputimrix, double **inputimrmi,
   double **inputimrrg, double **inputimrpk,
   double **inputimltb, double **inputimlix, double **inputimlmi,
   double **inputimlrg, double **inputimlpk,
   double **forwardrertb, double **forwardrerix, double **forwardrermi,
   double **forwardrerrg, double **forwardrerpk,
   double **forwardreltb, double **forwardrelix, double **forwardrelmi,
   double **forwardrelrg, double **forwardrelpk,
   double **forwardimrtb, double **forwardimrix, double **forwardimrmi,
   double **forwardimrrg, double **forwardimrpk,
   double **forwardimltb, double **forwardimlix, double **forwardimlmi,
   double **forwardimlrg, double **forwardimlpk,
   double **backwardrertb, double **backwardrerix, double **backwardrermi,
   double **backwardrerrg, double **backwardrerpk,
   double **backwardreltb, double **backwardrelix, double **backwardrelmi,
   double **backwardrelrg, double **backwardrelpk,
   double **backwardimrtb, double **backwardimrix, double **backwardimrmi,
   double **backwardimrrg, double **backwardimrpk,
   double **backwardimltb, double **backwardimlix, double **backwardimlmi,
   double **backwardimlrg, double **backwardimlpk,
   double **crossrertb, double **crossrerix, double **crossrermi,
   double **crossrerrg, double **crossrerpk,
   double **crossreltb, double **crossrelix, double **crossrelmi,
   double **crossrelrg, double **crossrelpk,
   double **crossimrtb, double **crossimrix, double **crossimrmi,
   double **crossimrrg, double **crossimrpk,
   double **crossimltb, double **crossimlix, double **crossimlmi,
   double **crossimlrg, double **crossimlpk )
{
}

void CPU_dbl10_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi,
   double *cffrrg, double *cffrpk,
   double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_dbl10_product(deg,
         inputrtb[ix],inputrix[ix],inputrmi[ix],inputrrg[ix],inputrpk[ix],
         inputltb[ix],inputlix[ix],inputlmi[ix],inputlrg[ix],inputlpk[ix],
         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
         cffltb,cfflix,cfflmi,cfflrg,cfflpk,
         outputrtb[dim],outputrix[dim],outputrmi[dim],
         outputrrg[dim],outputrpk[dim],
         outputltb[dim],outputlix[dim],outputlmi[dim],
         outputlrg[dim],outputlpk[dim]);

      for(int i=0; i<=deg; i++)
      {
         outputrtb[ix][i] = cffrtb[i];
         outputrix[ix][i] = cffrix[i];
         outputrmi[ix][i] = cffrmi[i];
         outputrrg[ix][i] = cffrrg[i];
         outputrpk[ix][i] = cffrpk[i];
         outputltb[ix][i] = cffltb[i];
         outputlix[ix][i] = cfflix[i];
         outputlmi[ix][i] = cfflmi[i];
         outputlrg[ix][i] = cfflrg[i];
         outputlpk[ix][i] = cfflpk[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_dbl10_product(deg,
         inputrtb[ix1],inputrix[ix1],inputrmi[ix1],inputrrg[ix1],inputrpk[ix1],
         inputltb[ix1],inputlix[ix1],inputlmi[ix1],inputlrg[ix1],inputlpk[ix1],
         inputrtb[ix2],inputrix[ix2],inputrmi[ix2],inputrrg[ix2],inputrpk[ix2],
         inputltb[ix2],inputlix[ix2],inputlmi[ix2],inputlrg[ix2],inputlpk[ix2],
         outputrtb[dim],outputrix[dim],outputrmi[dim],
         outputrrg[dim],outputrpk[dim],
         outputltb[dim],outputlix[dim],outputlmi[dim],
         outputlrg[dim],outputlpk[dim]);
      CPU_dbl10_product(deg,
         outputrtb[dim],outputrix[dim],outputrmi[dim],
         outputrrg[dim],outputrpk[dim],
         outputltb[dim],outputlix[dim],outputlmi[dim],
         outputlrg[dim],outputlpk[dim],
         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,cffltb,cfflix,cfflmi,cfflrg,cfflpk,
         outputrtb[dim],outputrix[dim],outputrmi[dim],
         outputrrg[dim],outputrpk[dim],
         outputltb[dim],outputlix[dim],outputlmi[dim],
         outputlrg[dim],outputlpk[dim]);

      CPU_dbl10_product(deg,
         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,cffltb,cfflix,cfflmi,cfflrg,cfflpk,
         inputrtb[ix1],inputrix[ix1],inputrmi[ix1],inputrrg[ix1],inputrpk[ix1],
         inputltb[ix1],inputlix[ix1],inputlmi[ix1],inputlrg[ix1],inputlpk[ix1],
         outputrtb[ix2],outputrix[ix2],outputrmi[ix2],
         outputrrg[ix2],outputrpk[ix2],
         outputltb[ix2],outputlix[ix2],outputlmi[ix2],
         outputlrg[ix2],outputlpk[ix2]);
      CPU_dbl10_product(deg,
         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,cffltb,cfflix,cfflmi,cfflrg,cfflpk,
         inputrtb[ix2],inputrix[ix2],inputrmi[ix2],inputrrg[ix2],inputrpk[ix2],
         inputltb[ix2],inputlix[ix2],inputlmi[ix2],inputlrg[ix2],inputlpk[ix2],
         outputrtb[ix1],outputrix[ix1],outputrmi[ix1],
         outputrrg[ix1],outputrpk[ix1],
         outputltb[ix1],outputlix[ix1],outputlmi[ix1],
         outputlrg[ix1],outputlpk[ix1]);
   }
   else
   {
      double **forwardrtb = new double*[nvr];
      double **forwardrix = new double*[nvr];
      double **forwardrmi = new double*[nvr];
      double **forwardrrg = new double*[nvr];
      double **forwardrpk = new double*[nvr];
      double **forwardltb = new double*[nvr];
      double **forwardlix = new double*[nvr];
      double **forwardlmi = new double*[nvr];
      double **forwardlrg = new double*[nvr];
      double **forwardlpk = new double*[nvr];
      double **backwardrtb = new double*[nvr-2];
      double **backwardrix = new double*[nvr-2];
      double **backwardrmi = new double*[nvr-2];
      double **backwardrrg = new double*[nvr-2];
      double **backwardrpk = new double*[nvr-2];
      double **backwardltb = new double*[nvr-2];
      double **backwardlix = new double*[nvr-2];
      double **backwardlmi = new double*[nvr-2];
      double **backwardlrg = new double*[nvr-2];
      double **backwardlpk = new double*[nvr-2];
      double **crossrtb = new double*[nvr-2];
      double **crossrix = new double*[nvr-2];
      double **crossrmi = new double*[nvr-2];
      double **crossrrg = new double*[nvr-2];
      double **crossrpk = new double*[nvr-2];
      double **crossltb = new double*[nvr-2];
      double **crosslix = new double*[nvr-2];
      double **crosslmi = new double*[nvr-2];
      double **crosslrg = new double*[nvr-2];
      double **crosslpk = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardrtb[i] = new double[deg+1];
         forwardrix[i] = new double[deg+1];
         forwardrmi[i] = new double[deg+1];
         forwardrrg[i] = new double[deg+1];
         forwardrpk[i] = new double[deg+1];
         forwardltb[i] = new double[deg+1];
         forwardlix[i] = new double[deg+1];
         forwardlmi[i] = new double[deg+1];
         forwardlrg[i] = new double[deg+1];
         forwardlpk[i] = new double[deg+1];
         backwardrtb[i] = new double[deg+1];
         backwardrix[i] = new double[deg+1];
         backwardrmi[i] = new double[deg+1];
         backwardrrg[i] = new double[deg+1];
         backwardrpk[i] = new double[deg+1];
         backwardltb[i] = new double[deg+1];
         backwardlix[i] = new double[deg+1];
         backwardlmi[i] = new double[deg+1];
         backwardlrg[i] = new double[deg+1];
         backwardlpk[i] = new double[deg+1];
         crossrtb[i] = new double[deg+1];
         crossrix[i] = new double[deg+1];
         crossrmi[i] = new double[deg+1];
         crossrrg[i] = new double[deg+1];
         crossrpk[i] = new double[deg+1];
         crossltb[i] = new double[deg+1];
         crosslix[i] = new double[deg+1];
         crosslmi[i] = new double[deg+1];
         crosslrg[i] = new double[deg+1];
         crosslpk[i] = new double[deg+1];
      }
      forwardrtb[nvr-2] = new double[deg+1];
      forwardrix[nvr-2] = new double[deg+1];
      forwardrmi[nvr-2] = new double[deg+1];
      forwardrrg[nvr-2] = new double[deg+1];
      forwardrpk[nvr-2] = new double[deg+1];
      forwardltb[nvr-2] = new double[deg+1];
      forwardlix[nvr-2] = new double[deg+1];
      forwardlmi[nvr-2] = new double[deg+1];
      forwardlrg[nvr-2] = new double[deg+1];
      forwardlpk[nvr-2] = new double[deg+1];
      forwardrtb[nvr-1] = new double[deg+1];
      forwardrix[nvr-1] = new double[deg+1];
      forwardrmi[nvr-1] = new double[deg+1];
      forwardrrg[nvr-1] = new double[deg+1];
      forwardrpk[nvr-1] = new double[deg+1];
      forwardltb[nvr-1] = new double[deg+1];
      forwardlix[nvr-1] = new double[deg+1];
      forwardlmi[nvr-1] = new double[deg+1];
      forwardlrg[nvr-1] = new double[deg+1];
      forwardlpk[nvr-1] = new double[deg+1];

      CPU_dbl10_speel(nvr,deg,idx,
         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
         cffltb,cfflix,cfflmi,cfflrg,cfflpk,
         inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
         inputltb,inputlix,inputlmi,inputlrg,inputlpk,
         forwardrtb,forwardrix,forwardrmi,forwardrrg,forwardrpk,
         forwardltb,forwardlix,forwardlmi,forwardlrg,forwardlpk,
         backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
         backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,
         crossrtb,crossrix,crossrmi,crossrrg,crossrpk,
         crossltb,crosslix,crosslmi,crosslrg,crosslpk);

      for(int i=0; i<deg+1; i++)     // assign the value of the monomial
      {
         outputrtb[dim][i] = forwardrtb[nvr-1][i];
         outputrix[dim][i] = forwardrix[nvr-1][i];
         outputrmi[dim][i] = forwardrmi[nvr-1][i];
         outputrrg[dim][i] = forwardrrg[nvr-1][i];
         outputrpk[dim][i] = forwardrpk[nvr-1][i];
         outputltb[dim][i] = forwardltb[nvr-1][i];
         outputlix[dim][i] = forwardlix[nvr-1][i];
         outputlmi[dim][i] = forwardlmi[nvr-1][i];
         outputlrg[dim][i] = forwardlrg[nvr-1][i];
         outputlpk[dim][i] = forwardlpk[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];       // derivative with respect to x[n-1]
         for(int i=0; i<deg+1; i++)
         {
            outputrtb[ix][i] = forwardrtb[nvr-2][i];
            outputrix[ix][i] = forwardrix[nvr-2][i];
            outputrmi[ix][i] = forwardrmi[nvr-2][i];
            outputrrg[ix][i] = forwardrrg[nvr-2][i];
            outputrpk[ix][i] = forwardrpk[nvr-2][i];
            outputltb[ix][i] = forwardltb[nvr-2][i];
            outputlix[ix][i] = forwardlix[nvr-2][i];
            outputlmi[ix][i] = forwardlmi[nvr-2][i];
            outputlrg[ix][i] = forwardlrg[nvr-2][i];
            outputlpk[ix][i] = forwardlpk[nvr-2][i];
         }
         ix = idx[0];               // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputrtb[ix][i] = backwardrtb[nvr-3][i];
            outputrix[ix][i] = backwardrix[nvr-3][i];
            outputrmi[ix][i] = backwardrmi[nvr-3][i];
            outputrrg[ix][i] = backwardrrg[nvr-3][i];
            outputrpk[ix][i] = backwardrpk[nvr-3][i];
            outputltb[ix][i] = backwardltb[nvr-3][i];
            outputlix[ix][i] = backwardlix[nvr-3][i];
            outputlmi[ix][i] = backwardlmi[nvr-3][i];
            outputlrg[ix][i] = backwardlrg[nvr-3][i];
            outputlpk[ix][i] = backwardlpk[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {                          // derivative with respect to x[k]
            ix = idx[k];
            for(int i=0; i<deg+1; i++)
            {
               outputrtb[ix][i] = crossrtb[k-1][i];
               outputrix[ix][i] = crossrix[k-1][i];
               outputrmi[ix][i] = crossrmi[k-1][i];
               outputrrg[ix][i] = crossrrg[k-1][i];
               outputrpk[ix][i] = crossrpk[k-1][i];
               outputltb[ix][i] = crossltb[k-1][i];
               outputlix[ix][i] = crosslix[k-1][i];
               outputlmi[ix][i] = crosslmi[k-1][i];
               outputlrg[ix][i] = crosslrg[k-1][i];
               outputlpk[ix][i] = crosslpk[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardrtb[i]); free(backwardrtb[i]); free(crossrtb[i]);
         free(forwardrix[i]); free(backwardrix[i]); free(crossrix[i]);
         free(forwardrmi[i]); free(backwardrmi[i]); free(crossrmi[i]);
         free(forwardrrg[i]); free(backwardrrg[i]); free(crossrrg[i]);
         free(forwardrpk[i]); free(backwardrpk[i]); free(crossrpk[i]);
         free(forwardltb[i]); free(backwardltb[i]); free(crossltb[i]);
         free(forwardlix[i]); free(backwardlix[i]); free(crosslix[i]);
         free(forwardlmi[i]); free(backwardlmi[i]); free(crosslmi[i]);
         free(forwardlrg[i]); free(backwardlrg[i]); free(crosslrg[i]);
         free(forwardlpk[i]); free(backwardlpk[i]); free(crosslpk[i]);
      }
      free(forwardrtb[nvr-2]); free(forwardrtb[nvr-1]);
      free(forwardrix[nvr-2]); free(forwardrix[nvr-1]);
      free(forwardrmi[nvr-2]); free(forwardrmi[nvr-1]);
      free(forwardrrg[nvr-2]); free(forwardrrg[nvr-1]);
      free(forwardrpk[nvr-2]); free(forwardrpk[nvr-1]);
      free(forwardltb[nvr-2]); free(forwardltb[nvr-1]);
      free(forwardlix[nvr-2]); free(forwardlix[nvr-1]);
      free(forwardlmi[nvr-2]); free(forwardlmi[nvr-1]);
      free(forwardlrg[nvr-2]); free(forwardlrg[nvr-1]);
      free(forwardlpk[nvr-2]); free(forwardlpk[nvr-1]);
      free(forwardrtb); free(backwardrtb); free(crossrtb);
      free(forwardrix); free(backwardrix); free(crossrix);
      free(forwardrmi); free(backwardrmi); free(crossrmi);
      free(forwardrrg); free(backwardrrg); free(crossrrg);
      free(forwardrpk); free(backwardrpk); free(crossrpk);
      free(forwardltb); free(backwardltb); free(crossltb);
      free(forwardlix); free(backwardlix); free(crosslix);
      free(forwardlmi); free(backwardlmi); free(crosslmi);
      free(forwardlrg); free(backwardlrg); free(crosslrg);
      free(forwardlpk); free(backwardlpk); free(crosslpk);
   }
}

void CPU_cmplx10_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cffretb, double *cffreix,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double **inputrertb, double **inputrerix, double **inputrermi,
   double **inputrerrg, double **inputrerpk,
   double **inputreltb, double **inputrelix, double **inputrelmi,
   double **inputrelrg, double **inputrelpk,
   double **inputimrtb, double **inputimrix, double **inputimrmi,
   double **inputimrrg, double **inputimrpk,
   double **inputimltb, double **inputimlix, double **inputimlmi,
   double **inputimlrg, double **inputimlpk,
   double *cffremi, double *cffrerg, double *cffrepk, double *cffimtb,
   double *cffimix, double *cffimmi, double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputrertb, double **outputrerix, double **outputrermi,
   double **outputrerrg, double **outputrerpk,
   double **outputreltb, double **outputrelix, double **outputrelmi,
   double **outputrelrg, double **outputrelpk,
   double **outputimrtb, double **outputimrix, double **outputimrmi,
   double **outputimrrg, double **outputimrpk,
   double **outputimltb, double **outputimlix, double **outputimlmi,
   double **outputimlrg, double **outputimlpk )
{
}
