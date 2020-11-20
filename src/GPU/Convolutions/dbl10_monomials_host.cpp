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
         backwardltb[0],backwardlix[0],backwardlmi[0],
         backwardlrg[0],backwardlpk[0]);

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
   int ix1 = idx[0];
   int ix2;
                                                           // f[0] = cff*x[0]
   CPU_cmplx10_product(deg,
      cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
      cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
      cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
      cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
      inputrertb[ix1],inputrerix[ix1],inputrermi[ix1],
      inputrerrg[ix1],inputrerpk[ix1],
      inputreltb[ix1],inputrelix[ix1],inputrelmi[ix1],
      inputrelrg[ix1],inputrelpk[ix1],
      inputimrtb[ix1],inputimrix[ix1],inputimrmi[ix1],
      inputimrrg[ix1],inputimrpk[ix1],
      inputimltb[ix1],inputimlix[ix1],inputimlmi[ix1],
      inputimlrg[ix1],inputimlpk[ix1],
      forwardrertb[0],forwardrerix[0],forwardrermi[0],
      forwardrerrg[0],forwardrerpk[0],
      forwardreltb[0],forwardrelix[0],forwardrelmi[0],
      forwardrelrg[0],forwardrelpk[0],
      forwardimrtb[0],forwardimrix[0],forwardimrmi[0],
      forwardimrrg[0],forwardimrpk[0],
      forwardimltb[0],forwardimlix[0],forwardimlmi[0],
      forwardimlrg[0],forwardimlpk[0]);

   for(int i=1; i<nvr; i++)
   {                                                    // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_cmplx10_product(deg,
         forwardrertb[i-1],forwardrerix[i-1],forwardrermi[i-1],
         forwardrerrg[i-1],forwardrerpk[i-1],
         forwardreltb[i-1],forwardrelix[i-1],forwardrelmi[i-1],
         forwardrelrg[i-1],forwardrelpk[i-1],
         forwardimrtb[i-1],forwardimrix[i-1],forwardimrmi[i-1],
         forwardimrrg[i-1],forwardimrpk[i-1],
         forwardimltb[i-1],forwardimlix[i-1],forwardimlmi[i-1],
         forwardimlrg[i-1],forwardimlpk[i-1],
         inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
         inputrerrg[ix2],inputrerpk[ix2],
         inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
         inputrelrg[ix2],inputrelpk[ix2],
         inputimrtb[ix2],inputimrix[ix2],inputimrmi[ix2],
         inputimrrg[ix2],inputimrpk[ix2],
         inputimltb[ix2],inputimlix[ix2],inputimlmi[ix2],
         inputimlrg[ix2],inputimlpk[ix2],
         forwardrertb[i],forwardrerix[i],forwardrermi[i],
         forwardrerrg[i],forwardrerpk[i],
         forwardreltb[i],forwardrelix[i],forwardrelmi[i],
         forwardrelrg[i],forwardrelpk[i],
         forwardimrtb[i],forwardimrix[i],forwardimrmi[i],
         forwardimrrg[i],forwardimrpk[i],
         forwardimltb[i],forwardimlix[i],forwardimlmi[i],
         forwardimlrg[i],forwardimlpk[i]);
   }
   if(nvr > 2)
   {                                                  // b[0] = x[n-1]*x[n-2]
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];
      CPU_cmplx10_product(deg,
         inputrertb[ix1],inputrerix[ix1],inputrermi[ix1],
         inputrerrg[ix1],inputrerpk[ix1],
         inputreltb[ix1],inputrelix[ix1],inputrelmi[ix1],
         inputrelrg[ix1],inputrelpk[ix1],
         inputimrtb[ix1],inputimrix[ix1],inputimrmi[ix1],
         inputimrrg[ix1],inputimrpk[ix1],
         inputimltb[ix1],inputimlix[ix1],inputimlmi[ix1],
         inputimlrg[ix1],inputimlpk[ix1],
         inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
         inputrerrg[ix2],inputrerpk[ix2],
         inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
         inputrelrg[ix2],inputrelpk[ix2],
         inputimrtb[ix2],inputimrix[ix2],inputimrmi[ix2],
         inputimrrg[ix2],inputimrpk[ix2],
         inputimltb[ix2],inputimlix[ix2],inputimlmi[ix2],
         inputimlrg[ix2],inputimlpk[ix2],
         backwardrertb[0],backwardrerix[0],backwardrermi[0],
         backwardrerrg[0],backwardrerpk[0],
         backwardreltb[0],backwardrelix[0],backwardrelmi[0],
         backwardrelrg[0],backwardrelpk[0],
         backwardimrtb[0],backwardimrix[0],backwardimrmi[0],
         backwardimrrg[0],backwardimrpk[0],
         backwardimltb[0],backwardimlix[0],backwardimlmi[0],
         backwardimlrg[0],backwardimlpk[0]);

      for(int i=1; i<nvr-2; i++)
      {                                             // b[i] = b[i-1]*x[x-2-i]
         ix2 = idx[nvr-2-i];
         CPU_cmplx10_product(deg,
            backwardrertb[i-1],backwardrerix[i-1],backwardrermi[i-1],
            backwardrerrg[i-1],backwardrerpk[i-1],
            backwardreltb[i-1],backwardrelix[i-1],backwardrelmi[i-1],
            backwardrelrg[i-1],backwardrelpk[i-1],
            backwardimrtb[i-1],backwardimrix[i-1],backwardimrmi[i-1],
            backwardimrrg[i-1],backwardimrpk[i-1],
            backwardimltb[i-1],backwardimlix[i-1],backwardimlmi[i-1],
            backwardimlrg[i-1],backwardimlpk[i-1],
            inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
            inputrerrg[ix2],inputrerpk[ix2],
            inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
            inputrelrg[ix2],inputrelpk[ix2],
            inputimrtb[ix2],inputimrix[ix2],inputimrmi[ix2],
            inputimrrg[ix2],inputimrpk[ix2],
            inputimltb[ix2],inputimlix[ix2],inputimlmi[ix2],
            inputimlrg[ix2],inputimlpk[ix2],
            backwardrertb[i],backwardrerix[i],backwardrermi[i],
            backwardrerrg[i],backwardrerpk[i],
            backwardreltb[i],backwardrelix[i],backwardrelmi[i],
            backwardrelrg[i],backwardrelpk[i],
            backwardimrtb[i],backwardimrix[i],backwardimrmi[i],
            backwardimrrg[i],backwardimrpk[i],
            backwardimltb[i],backwardimlix[i],backwardimlmi[i],
            backwardimlrg[i],backwardimlpk[i]);
      }
                                                       // b[n-3] = b[n-3]*cff
      CPU_cmplx10_product(deg,
         backwardrertb[nvr-3],backwardrerix[nvr-3],backwardrermi[nvr-3],
         backwardrerrg[nvr-3],backwardrerpk[nvr-3],
         backwardreltb[nvr-3],backwardrelix[nvr-3],backwardrelmi[nvr-3],
         backwardrelrg[nvr-3],backwardrelpk[nvr-3],
         backwardimrtb[nvr-3],backwardimrix[nvr-3],backwardimrmi[nvr-3],
         backwardimrrg[nvr-3],backwardimrpk[nvr-3],
         backwardimltb[nvr-3],backwardimlix[nvr-3],backwardimlmi[nvr-3],
         backwardimlrg[nvr-3],backwardimlpk[nvr-3],
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         crossrertb[0],crossrerix[0],crossrermi[0],
         crossrerrg[0],crossrerpk[0],
         crossreltb[0],crossrelix[0],crossrelmi[0],
         crossrelrg[0],crossrelpk[0],
         crossimrtb[0],crossimrix[0],crossimrmi[0],
         crossimrrg[0],crossimrpk[0],
         crossimltb[0],crossimlix[0],crossimlmi[0],
         crossimlrg[0],crossimlpk[0]);                 // cross is work space

      for(int i=0; i<=deg; i++)
      {
         backwardrertb[nvr-3][i] = crossrertb[0][i];
         backwardrerix[nvr-3][i] = crossrerix[0][i];
         backwardrermi[nvr-3][i] = crossrermi[0][i];
         backwardrerrg[nvr-3][i] = crossrerrg[0][i];
         backwardrerpk[nvr-3][i] = crossrerpk[0][i];
         backwardreltb[nvr-3][i] = crossreltb[0][i];
         backwardrelix[nvr-3][i] = crossrelix[0][i];
         backwardrelmi[nvr-3][i] = crossrelmi[0][i];
         backwardrelrg[nvr-3][i] = crossrelrg[0][i];
         backwardrelpk[nvr-3][i] = crossrelpk[0][i];
         backwardimrtb[nvr-3][i] = crossimrtb[0][i];
         backwardimrix[nvr-3][i] = crossimrix[0][i];
         backwardimrmi[nvr-3][i] = crossimrmi[0][i];
         backwardimrrg[nvr-3][i] = crossimrrg[0][i];
         backwardimrpk[nvr-3][i] = crossimrpk[0][i];
         backwardimltb[nvr-3][i] = crossimltb[0][i];
         backwardimlix[nvr-3][i] = crossimlix[0][i];
         backwardimlmi[nvr-3][i] = crossimlmi[0][i];
         backwardimlrg[nvr-3][i] = crossimlrg[0][i];
         backwardimlpk[nvr-3][i] = crossimlpk[0][i];
      }
      if(nvr == 3)
      {                                                   // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_cmplx10_product(deg,
            forwardrertb[0],forwardrerix[0],forwardrermi[0],
            forwardrerrg[0],forwardrerpk[0],
            forwardreltb[0],forwardrelix[0],forwardrelmi[0],
            forwardrelrg[0],forwardrelpk[0],
            forwardimrtb[0],forwardimrix[0],forwardimrmi[0],
            forwardimrrg[0],forwardimrpk[0],
            forwardimltb[0],forwardimlix[0],forwardimlmi[0],
            forwardimlrg[0],forwardimlpk[0],
            inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
            inputrerrg[ix2],inputrerpk[ix2],
            inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
            inputrelrg[ix2],inputrelpk[ix2],
            inputimrtb[ix2],inputimrix[ix2],inputimrmi[ix2],
            inputimrrg[ix2],inputimrpk[ix2],
            inputimltb[ix2],inputimlix[ix2],inputimlmi[ix2],
            inputimlrg[ix2],inputimlpk[ix2],
            crossrertb[0],crossrerix[0],crossrermi[0],
            crossrerrg[0],crossrerpk[0],
            crossreltb[0],crossrelix[0],crossrelmi[0],
            crossrelrg[0],crossrelpk[0],
            crossimrtb[0],crossimrix[0],crossimrmi[0],
            crossimrrg[0],crossimrpk[0],
            crossimltb[0],crossimlix[0],crossimlmi[0],
            crossimlrg[0],crossimlpk[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                            // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_cmplx10_product(deg,
               forwardrertb[i],forwardrerix[i],forwardrermi[i],
               forwardrerrg[i],forwardrerpk[i],
               forwardreltb[i],forwardrelix[i],forwardrelmi[i],
               forwardrelrg[i],forwardrelpk[i],
               forwardimrtb[i],forwardimrix[i],forwardimrmi[i],
               forwardimrrg[i],forwardimrpk[i],
               forwardimltb[i],forwardimlix[i],forwardimlmi[i],
               forwardimlrg[i],forwardimlpk[i],
               backwardrertb[ix2],backwardrerix[ix2],backwardrermi[ix2],
               backwardrerrg[ix2],backwardrerpk[ix2],
               backwardreltb[ix2],backwardrelix[ix2],backwardrelmi[ix2],
               backwardrelrg[ix2],backwardrelpk[ix2],
               backwardimrtb[ix2],backwardimrix[ix2],backwardimrmi[ix2],
               backwardimrrg[ix2],backwardimrpk[ix2],
               backwardimltb[ix2],backwardimlix[ix2],backwardimlmi[ix2],
               backwardimlrg[ix2],backwardimlpk[ix2],
               crossrertb[i],crossrerix[i],crossrermi[i],
               crossrerrg[i],crossrerpk[i],
               crossreltb[i],crossrelix[i],crossrelmi[i],
               crossrelrg[i],crossrelpk[i],
               crossimrtb[i],crossimrix[i],crossimrmi[i],
               crossimrrg[i],crossimrpk[i],
               crossimltb[i],crossimlix[i],crossimlmi[i],
               crossimlrg[i],crossimlpk[i]);
         }
                                                    // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_cmplx10_product(deg,
            forwardrertb[nvr-3],forwardrerix[nvr-3],forwardrermi[nvr-3],
            forwardrerrg[nvr-3],forwardrerpk[nvr-3],
            forwardreltb[nvr-3],forwardrelix[nvr-3],forwardrelmi[nvr-3],
            forwardrelrg[nvr-3],forwardrelpk[nvr-3],
            forwardimrtb[nvr-3],forwardimrix[nvr-3],forwardimrmi[nvr-3],
            forwardimrrg[nvr-3],forwardimrpk[nvr-3],
            forwardimltb[nvr-3],forwardimlix[nvr-3],forwardimlmi[nvr-3],
            forwardimlrg[nvr-3],forwardimlpk[nvr-3],
            inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
            inputrerrg[ix2],inputrerpk[ix2],
            inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
            inputrelrg[ix2],inputrelpk[ix2],
            inputimrtb[ix2],inputimrix[ix2],inputimrmi[ix2],
            inputimrrg[ix2],inputimrpk[ix2],
            inputimltb[ix2],inputimlix[ix2],inputimlmi[ix2],
            inputimlrg[ix2],inputimlpk[ix2],
            crossrertb[nvr-3],crossrerix[nvr-3],crossrermi[nvr-3],
            crossrerrg[nvr-3],crossrerpk[nvr-3],
            crossreltb[nvr-3],crossrelix[nvr-3],crossrelmi[nvr-3],
            crossrelrg[nvr-3],crossrelpk[nvr-3],
            crossimrtb[nvr-3],crossimrix[nvr-3],crossimrmi[nvr-3],
            crossimrrg[nvr-3],crossimrpk[nvr-3],
            crossimltb[nvr-3],crossimlix[nvr-3],crossimlmi[nvr-3],
            crossimlrg[nvr-3],crossimlpk[nvr-3]);
      }
   }
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
 ( int dim, int nvr, int deg, int *idx,
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
   double **outputrertb, double **outputrerix, double **outputrermi,
   double **outputrerrg, double **outputrerpk,
   double **outputreltb, double **outputrelix, double **outputrelmi,
   double **outputrelrg, double **outputrelpk,
   double **outputimrtb, double **outputimrix, double **outputimrmi,
   double **outputimrrg, double **outputimrpk,
   double **outputimltb, double **outputimlix, double **outputimlmi,
   double **outputimlrg, double **outputimlpk )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_cmplx10_product(deg,
         inputrertb[ix],inputrerix[ix],inputrermi[ix],
         inputrerrg[ix],inputrerpk[ix],
         inputreltb[ix],inputrelix[ix],inputrelmi[ix],
         inputrelrg[ix],inputrelpk[ix],
         inputimrtb[ix],inputimrix[ix],inputimrmi[ix],
         inputimrrg[ix],inputimrpk[ix],
         inputimltb[ix],inputimlix[ix],inputimlmi[ix],
         inputimlrg[ix],inputimlpk[ix],
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         outputrertb[dim],outputrerix[dim],outputrermi[dim],
         outputrerrg[dim],outputrerpk[dim],
         outputreltb[dim],outputrelix[dim],outputrelmi[dim],
         outputrelrg[dim],outputrelpk[dim],
         outputimrtb[dim],outputimrix[dim],outputimrmi[dim],
         outputimrrg[dim],outputimrpk[dim],
         outputimltb[dim],outputimlix[dim],outputimlmi[dim],
         outputimlrg[dim],outputimlpk[dim]);

      for(int i=0; i<=deg; i++) 
      {
         outputrertb[ix][i] = cffrertb[i]; outputrerix[ix][i] = cffrerix[i];
         outputrermi[ix][i] = cffrermi[i]; outputrerrg[ix][i] = cffrerrg[i];
         outputrerpk[ix][i] = cffrerpk[i];
         outputreltb[ix][i] = cffreltb[i]; outputrelix[ix][i] = cffrelix[i];
         outputrelmi[ix][i] = cffrelmi[i]; outputrelrg[ix][i] = cffrelrg[i];
         outputrelpk[ix][i] = cffrelpk[i];
         outputimrtb[ix][i] = cffimrtb[i]; outputimrix[ix][i] = cffimrix[i];
         outputimrmi[ix][i] = cffimrmi[i]; outputimrrg[ix][i] = cffimrrg[i];
         outputimrpk[ix][i] = cffimrpk[i];
         outputimltb[ix][i] = cffimltb[i]; outputimlix[ix][i] = cffimlix[i];
         outputimlmi[ix][i] = cffimlmi[i]; outputimlrg[ix][i] = cffimlrg[i];
         outputimlpk[ix][i] = cffimlpk[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0];
      int ix2 = idx[1];

      CPU_cmplx10_product(deg,
         inputrertb[ix1],inputrerix[ix1],inputrermi[ix1],
         inputrerrg[ix1],inputrerpk[ix1],
         inputreltb[ix1],inputrelix[ix1],inputrelmi[ix1],
         inputrelrg[ix1],inputrelpk[ix1],
         inputimrtb[ix1],inputimrix[ix1],inputimrmi[ix1],
         inputimrrg[ix1],inputimrpk[ix1],
         inputimltb[ix1],inputimlix[ix1],inputimlmi[ix1],
         inputimlrg[ix1],inputimlpk[ix1],
         inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
         inputrerrg[ix2],inputrerpk[ix2],
         inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
         inputrelrg[ix2],inputrelpk[ix2],
         inputimrtb[ix2],inputimrix[ix2],inputimrmi[ix2],
         inputimrrg[ix2],inputimrpk[ix2],
         inputimltb[ix2],inputimlix[ix2],inputimlmi[ix2],
         inputimlrg[ix2],inputimlpk[ix2],
         outputrertb[dim],outputrerix[dim],outputrermi[dim],
         outputrerrg[dim],outputrerpk[dim],
         outputreltb[dim],outputrelix[dim],outputrelmi[dim],
         outputrelrg[dim],outputrelpk[dim],
         outputimrtb[dim],outputimrix[dim],outputimrmi[dim],
         outputimrrg[dim],outputimrpk[dim],
         outputimltb[dim],outputimlix[dim],outputimlmi[dim],
         outputimlrg[dim],outputimlpk[dim]);
      CPU_cmplx10_product(deg,
         outputrertb[dim],outputrerix[dim],outputrermi[dim],
         outputrerrg[dim],outputrerpk[dim],
         outputreltb[dim],outputrelix[dim],outputrelmi[dim],
         outputrelrg[dim],outputrelpk[dim],
         outputimrtb[dim],outputimrix[dim],outputimrmi[dim],
         outputimrrg[dim],outputimrpk[dim],
         outputimltb[dim],outputimlix[dim],outputimlmi[dim],
         outputimlrg[dim],outputimlpk[dim],
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         outputrertb[dim],outputrerix[dim],outputrermi[dim],
         outputrerrg[dim],outputrerpk[dim],
         outputreltb[dim],outputrelix[dim],outputrelmi[dim],
         outputrelrg[dim],outputrelpk[dim],
         outputimrtb[dim],outputimrix[dim],outputimrmi[dim],
         outputimrrg[dim],outputimrpk[dim],
         outputimltb[dim],outputimlix[dim],outputimlmi[dim],
         outputimlrg[dim],outputimlpk[dim]);

      CPU_cmplx10_product(deg,
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         inputrertb[ix1],inputrerix[ix1],inputrermi[ix1],
         inputrerrg[ix1],inputrerpk[ix1],
         inputreltb[ix1],inputrelix[ix1],inputrelmi[ix1],
         inputrelrg[ix1],inputrelpk[ix1],
         inputimrtb[ix1],inputimrix[ix1],inputimrmi[ix1],
         inputimrrg[ix1],inputimrpk[ix1],
         inputimltb[ix1],inputimlix[ix1],inputimlmi[ix1],
         inputimlrg[ix1],inputimlpk[ix1],
         outputrertb[ix2],outputrerix[ix2],outputrermi[ix2],
         outputrerrg[ix2],outputrerpk[ix2],
         outputreltb[ix2],outputrelix[ix2],outputrelmi[ix2],
         outputrelrg[ix2],outputrelpk[ix2],
         outputimrtb[ix2],outputimrix[ix2],outputimrmi[ix2],
         outputimrrg[ix2],outputimrpk[ix2],
         outputimltb[ix2],outputimlix[ix2],outputimlmi[ix2],
         outputimlrg[ix2],outputimlpk[ix2]);
      CPU_cmplx10_product(deg,
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         inputrertb[ix2],inputrerix[ix2],inputrermi[ix2],
         inputrerrg[ix2],inputrerpk[ix2],
         inputreltb[ix2],inputrelix[ix2],inputrelmi[ix2],
         inputrelrg[ix2],inputrelpk[ix2],
         inputimrtb[ix1],inputimrix[ix1],inputimrmi[ix1],
         inputimrrg[ix1],inputimrpk[ix1],
         inputimltb[ix1],inputimlix[ix1],inputimlmi[ix1],
         inputimlrg[ix1],inputimlpk[ix1],
         outputrertb[ix1],outputrerix[ix1],outputrermi[ix1],
         outputrerrg[ix1],outputrerpk[ix1],
         outputreltb[ix1],outputrelix[ix1],outputrelmi[ix1],
         outputrelrg[ix1],outputrelpk[ix1],
         outputimrtb[ix1],outputimrix[ix1],outputimrmi[ix1],
         outputimrrg[ix1],outputimrpk[ix1],
         outputimltb[ix1],outputimlix[ix1],outputimlmi[ix1],
         outputimlrg[ix1],outputimlpk[ix1]);
   }
   else
   {
      double **forwardrertb = new double*[nvr];
      double **forwardrerix = new double*[nvr];
      double **forwardrermi = new double*[nvr];
      double **forwardrerrg = new double*[nvr];
      double **forwardrerpk = new double*[nvr];
      double **forwardreltb = new double*[nvr];
      double **forwardrelix = new double*[nvr];
      double **forwardrelmi = new double*[nvr];
      double **forwardrelrg = new double*[nvr];
      double **forwardrelpk = new double*[nvr];
      double **forwardimrtb = new double*[nvr];
      double **forwardimrix = new double*[nvr];
      double **forwardimrmi = new double*[nvr];
      double **forwardimrrg = new double*[nvr];
      double **forwardimrpk = new double*[nvr];
      double **forwardimltb = new double*[nvr];
      double **forwardimlix = new double*[nvr];
      double **forwardimlmi = new double*[nvr];
      double **forwardimlrg = new double*[nvr];
      double **forwardimlpk = new double*[nvr];
      double **backwardrertb = new double*[nvr-2];
      double **backwardrerix = new double*[nvr-2];
      double **backwardrermi = new double*[nvr-2];
      double **backwardrerrg = new double*[nvr-2];
      double **backwardrerpk = new double*[nvr-2];
      double **backwardreltb = new double*[nvr-2];
      double **backwardrelix = new double*[nvr-2];
      double **backwardrelmi = new double*[nvr-2];
      double **backwardrelrg = new double*[nvr-2];
      double **backwardrelpk = new double*[nvr-2];
      double **backwardimrtb = new double*[nvr-2];
      double **backwardimrix = new double*[nvr-2];
      double **backwardimrmi = new double*[nvr-2];
      double **backwardimrrg = new double*[nvr-2];
      double **backwardimrpk = new double*[nvr-2];
      double **backwardimltb = new double*[nvr-2];
      double **backwardimlix = new double*[nvr-2];
      double **backwardimlmi = new double*[nvr-2];
      double **backwardimlrg = new double*[nvr-2];
      double **backwardimlpk = new double*[nvr-2];
      double **crossrertb = new double*[nvr-2];
      double **crossrerix = new double*[nvr-2];
      double **crossrermi = new double*[nvr-2];
      double **crossrerrg = new double*[nvr-2];
      double **crossrerpk = new double*[nvr-2];
      double **crossreltb = new double*[nvr-2];
      double **crossrelix = new double*[nvr-2];
      double **crossrelmi = new double*[nvr-2];
      double **crossrelrg = new double*[nvr-2];
      double **crossrelpk = new double*[nvr-2];
      double **crossimrtb = new double*[nvr-2];
      double **crossimrix = new double*[nvr-2];
      double **crossimrmi = new double*[nvr-2];
      double **crossimrrg = new double*[nvr-2];
      double **crossimrpk = new double*[nvr-2];
      double **crossimltb = new double*[nvr-2];
      double **crossimlix = new double*[nvr-2];
      double **crossimlmi = new double*[nvr-2];
      double **crossimlrg = new double*[nvr-2];
      double **crossimlpk = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardrertb[i] = new double[deg+1];
         forwardrerix[i] = new double[deg+1];
         forwardrermi[i] = new double[deg+1];
         forwardrerrg[i] = new double[deg+1];
         forwardrerpk[i] = new double[deg+1];
         forwardreltb[i] = new double[deg+1];
         forwardrelix[i] = new double[deg+1];
         forwardrelmi[i] = new double[deg+1];
         forwardrelrg[i] = new double[deg+1];
         forwardrelpk[i] = new double[deg+1];
         forwardimrtb[i] = new double[deg+1];
         forwardimrix[i] = new double[deg+1];
         forwardimrmi[i] = new double[deg+1];
         forwardimrrg[i] = new double[deg+1];
         forwardimrpk[i] = new double[deg+1];
         forwardimltb[i] = new double[deg+1];
         forwardimlix[i] = new double[deg+1];
         forwardimlmi[i] = new double[deg+1];
         forwardimlrg[i] = new double[deg+1];
         forwardimlpk[i] = new double[deg+1];
         backwardrertb[i] = new double[deg+1];
         backwardrerix[i] = new double[deg+1];
         backwardrermi[i] = new double[deg+1];
         backwardrerrg[i] = new double[deg+1];
         backwardrerpk[i] = new double[deg+1];
         backwardreltb[i] = new double[deg+1];
         backwardrelix[i] = new double[deg+1];
         backwardrelmi[i] = new double[deg+1];
         backwardrelrg[i] = new double[deg+1];
         backwardrelpk[i] = new double[deg+1];
         backwardimrtb[i] = new double[deg+1];
         backwardimrix[i] = new double[deg+1];
         backwardimrmi[i] = new double[deg+1];
         backwardimrrg[i] = new double[deg+1];
         backwardimrpk[i] = new double[deg+1];
         backwardimltb[i] = new double[deg+1];
         backwardimlix[i] = new double[deg+1];
         backwardimlmi[i] = new double[deg+1];
         backwardimlrg[i] = new double[deg+1];
         backwardimlpk[i] = new double[deg+1];
         crossrertb[i] = new double[deg+1];
         crossrerix[i] = new double[deg+1];
         crossrermi[i] = new double[deg+1];
         crossrerrg[i] = new double[deg+1];
         crossrerpk[i] = new double[deg+1];
         crossreltb[i] = new double[deg+1];
         crossrelix[i] = new double[deg+1];
         crossrelmi[i] = new double[deg+1];
         crossrelrg[i] = new double[deg+1];
         crossrelpk[i] = new double[deg+1];
         crossimrtb[i] = new double[deg+1];
         crossimrix[i] = new double[deg+1];
         crossimrmi[i] = new double[deg+1];
         crossimrrg[i] = new double[deg+1];
         crossimrpk[i] = new double[deg+1];
         crossimltb[i] = new double[deg+1];
         crossimlix[i] = new double[deg+1];
         crossimlmi[i] = new double[deg+1];
         crossimlrg[i] = new double[deg+1];
         crossimlpk[i] = new double[deg+1];
      }
      forwardrertb[nvr-2] = new double[deg+1];
      forwardrerix[nvr-2] = new double[deg+1];
      forwardrermi[nvr-2] = new double[deg+1];
      forwardrerrg[nvr-2] = new double[deg+1];
      forwardrerpk[nvr-2] = new double[deg+1];
      forwardreltb[nvr-2] = new double[deg+1];
      forwardrelix[nvr-2] = new double[deg+1];
      forwardrelmi[nvr-2] = new double[deg+1];
      forwardrelrg[nvr-2] = new double[deg+1];
      forwardrelpk[nvr-2] = new double[deg+1];
      forwardimrtb[nvr-2] = new double[deg+1];
      forwardimrix[nvr-2] = new double[deg+1];
      forwardimrmi[nvr-2] = new double[deg+1];
      forwardimrrg[nvr-2] = new double[deg+1];
      forwardimrpk[nvr-2] = new double[deg+1];
      forwardimltb[nvr-2] = new double[deg+1];
      forwardimlix[nvr-2] = new double[deg+1];
      forwardimlmi[nvr-2] = new double[deg+1];
      forwardimlrg[nvr-2] = new double[deg+1];
      forwardimlpk[nvr-2] = new double[deg+1];
      forwardrertb[nvr-1] = new double[deg+1];
      forwardrerix[nvr-1] = new double[deg+1];
      forwardrermi[nvr-1] = new double[deg+1];
      forwardrerrg[nvr-1] = new double[deg+1];
      forwardrerpk[nvr-1] = new double[deg+1];
      forwardreltb[nvr-1] = new double[deg+1];
      forwardrelix[nvr-1] = new double[deg+1];
      forwardrelmi[nvr-1] = new double[deg+1];
      forwardrelrg[nvr-1] = new double[deg+1];
      forwardrelpk[nvr-1] = new double[deg+1];
      forwardimrtb[nvr-1] = new double[deg+1];
      forwardimrix[nvr-1] = new double[deg+1];
      forwardimrmi[nvr-1] = new double[deg+1];
      forwardimrrg[nvr-1] = new double[deg+1];
      forwardimrpk[nvr-1] = new double[deg+1];
      forwardimltb[nvr-1] = new double[deg+1];
      forwardimlix[nvr-1] = new double[deg+1];
      forwardimlmi[nvr-1] = new double[deg+1];
      forwardimlrg[nvr-1] = new double[deg+1];
      forwardimlpk[nvr-1] = new double[deg+1];

      CPU_cmplx10_speel(nvr,deg,idx,
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         inputrertb,inputrerix,inputrermi,inputrerrg,inputrerpk,
         inputreltb,inputrelix,inputrelmi,inputrelrg,inputrelpk,
         inputimrtb,inputimrix,inputimrmi,inputimrrg,inputimrpk,
         inputimltb,inputimlix,inputimlmi,inputimlrg,inputimlpk,
         forwardrertb,forwardrerix,forwardrermi,forwardrerrg,forwardrerpk,
         forwardreltb,forwardrelix,forwardrelmi,forwardrelrg,forwardrelpk,
         forwardimrtb,forwardimrix,forwardimrmi,forwardimrrg,forwardimrpk,
         forwardimltb,forwardimlix,forwardimlmi,forwardimlrg,forwardimlpk,
         backwardrertb,backwardrerix,backwardrermi,backwardrerrg,backwardrerpk,
         backwardreltb,backwardrelix,backwardrelmi,backwardrelrg,backwardrelpk,
         backwardimrtb,backwardimrix,backwardimrmi,backwardimrrg,backwardimrpk,
         backwardimltb,backwardimlix,backwardimlmi,backwardimlrg,backwardimlpk,
         crossrertb,crossrerix,crossrermi,crossrerrg,crossrerpk,
         crossreltb,crossrelix,crossrelmi,crossrelrg,crossrelpk,
         crossimrtb,crossimrix,crossimrmi,crossimrrg,crossimrpk,
         crossimltb,crossimlix,crossimlmi,crossimlrg,crossimlpk);

      for(int i=0; i<deg+1; i++)          // assign value of the monomial
      {
         outputrertb[dim][i] = forwardrertb[nvr-1][i];
         outputrerix[dim][i] = forwardrerix[nvr-1][i];
         outputrermi[dim][i] = forwardrermi[nvr-1][i];
         outputrerrg[dim][i] = forwardrerrg[nvr-1][i];
         outputrerpk[dim][i] = forwardrerpk[nvr-1][i];
         outputreltb[dim][i] = forwardreltb[nvr-1][i];
         outputrelix[dim][i] = forwardrelix[nvr-1][i];
         outputrelmi[dim][i] = forwardrelmi[nvr-1][i];
         outputrelrg[dim][i] = forwardrelrg[nvr-1][i];
         outputrelpk[dim][i] = forwardrelpk[nvr-1][i];
         outputimrtb[dim][i] = forwardimrtb[nvr-1][i];
         outputimrix[dim][i] = forwardimrix[nvr-1][i];
         outputimrmi[dim][i] = forwardimrmi[nvr-1][i];
         outputimrrg[dim][i] = forwardimrrg[nvr-1][i];
         outputimrpk[dim][i] = forwardimrpk[nvr-1][i];
         outputimltb[dim][i] = forwardimltb[nvr-1][i];
         outputimlix[dim][i] = forwardimlix[nvr-1][i];
         outputimlmi[dim][i] = forwardimlmi[nvr-1][i];
         outputimlrg[dim][i] = forwardimlrg[nvr-1][i];
         outputimlpk[dim][i] = forwardimlpk[nvr-1][i];
      }
      if(nvr > 2)
      {
         int ix = idx[nvr-1];        // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++)
         {
            outputrertb[ix][i] = forwardrertb[nvr-2][i];
            outputrerix[ix][i] = forwardrerix[nvr-2][i];
            outputrermi[ix][i] = forwardrermi[nvr-2][i];
            outputrerrg[ix][i] = forwardrerrg[nvr-2][i];
            outputrerpk[ix][i] = forwardrerpk[nvr-2][i];
            outputreltb[ix][i] = forwardreltb[nvr-2][i];
            outputrelix[ix][i] = forwardrelix[nvr-2][i];
            outputrelmi[ix][i] = forwardrelmi[nvr-2][i];
            outputrelrg[ix][i] = forwardrelrg[nvr-2][i];
            outputrelpk[ix][i] = forwardrelpk[nvr-2][i];
            outputimrtb[ix][i] = forwardimrtb[nvr-2][i];
            outputimrix[ix][i] = forwardimrix[nvr-2][i];
            outputimrmi[ix][i] = forwardimrmi[nvr-2][i];
            outputimrrg[ix][i] = forwardimrrg[nvr-2][i];
            outputimrpk[ix][i] = forwardimrpk[nvr-2][i];
            outputimltb[ix][i] = forwardimltb[nvr-2][i];
            outputimlix[ix][i] = forwardimlix[nvr-2][i];
            outputimlmi[ix][i] = forwardimlmi[nvr-2][i];
            outputimlrg[ix][i] = forwardimlrg[nvr-2][i];
            outputimlpk[ix][i] = forwardimlpk[nvr-2][i];
         }

         ix = idx[0];                  // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputrertb[ix][i] = backwardrertb[nvr-3][i];
            outputrerix[ix][i] = backwardrerix[nvr-3][i];
            outputrermi[ix][i] = backwardrermi[nvr-3][i];
            outputrerrg[ix][i] = backwardrerrg[nvr-3][i];
            outputrerpk[ix][i] = backwardrerpk[nvr-3][i];
            outputreltb[ix][i] = backwardreltb[nvr-3][i];
            outputrelix[ix][i] = backwardrelix[nvr-3][i];
            outputrelmi[ix][i] = backwardrelmi[nvr-3][i];
            outputrelrg[ix][i] = backwardrelrg[nvr-3][i];
            outputrelpk[ix][i] = backwardrelpk[nvr-3][i];
            outputimrtb[ix][i] = backwardimrtb[nvr-3][i];
            outputimrix[ix][i] = backwardimrix[nvr-3][i];
            outputimrmi[ix][i] = backwardimrmi[nvr-3][i];
            outputimrrg[ix][i] = backwardimrrg[nvr-3][i];
            outputimrpk[ix][i] = backwardimrpk[nvr-3][i];
            outputimltb[ix][i] = backwardimltb[nvr-3][i];
            outputimlix[ix][i] = backwardimlix[nvr-3][i];
            outputimlmi[ix][i] = backwardimlmi[nvr-3][i];
            outputimlrg[ix][i] = backwardimlrg[nvr-3][i];
            outputimlpk[ix][i] = backwardimlpk[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];               // derivative with respect to x[k]
            for(int i=0; i<deg+1; i++)
            {
               outputrertb[ix][i] = crossrertb[k-1][i];
               outputrerix[ix][i] = crossrerix[k-1][i];
               outputrermi[ix][i] = crossrermi[k-1][i];
               outputrerrg[ix][i] = crossrerrg[k-1][i];
               outputrerpk[ix][i] = crossrerpk[k-1][i];
               outputreltb[ix][i] = crossreltb[k-1][i];
               outputrelix[ix][i] = crossrelix[k-1][i];
               outputrelmi[ix][i] = crossrelmi[k-1][i];
               outputrelrg[ix][i] = crossrelrg[k-1][i];
               outputrelpk[ix][i] = crossrelpk[k-1][i];
               outputimrtb[ix][i] = crossimrtb[k-1][i];
               outputimrix[ix][i] = crossimrix[k-1][i];
               outputimrmi[ix][i] = crossimrmi[k-1][i];
               outputimrrg[ix][i] = crossimrrg[k-1][i];
               outputimrpk[ix][i] = crossimrpk[k-1][i];
               outputimltb[ix][i] = crossimltb[k-1][i];
               outputimlix[ix][i] = crossimlix[k-1][i];
               outputimlmi[ix][i] = crossimlmi[k-1][i];
               outputimlrg[ix][i] = crossimlrg[k-1][i];
               outputimlpk[ix][i] = crossimlpk[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardrertb[i]);  free(forwardrerix[i]);
         free(forwardrermi[i]);  free(forwardrerrg[i]);
         free(forwardrerpk[i]);
         free(forwardreltb[i]);  free(forwardrelix[i]);
         free(forwardrelmi[i]);  free(forwardrelrg[i]);
         free(forwardrelpk[i]);
         free(forwardimrtb[i]);  free(forwardimrix[i]);
         free(forwardimrmi[i]);  free(forwardimrrg[i]);
         free(forwardimrpk[i]);
         free(forwardimltb[i]);  free(forwardimlix[i]);
         free(forwardimlmi[i]);  free(forwardimlrg[i]);
         free(forwardimlpk[i]);
         free(backwardrertb[i]); free(backwardrerix[i]);
         free(backwardrermi[i]); free(backwardrerrg[i]);
         free(backwardrerpk[i]);
         free(backwardreltb[i]); free(backwardrelix[i]);
         free(backwardrelmi[i]); free(backwardrelrg[i]);
         free(backwardrelpk[i]);
         free(backwardimrtb[i]); free(backwardimrix[i]);
         free(backwardimrmi[i]); free(backwardimrrg[i]);
         free(backwardimrpk[i]);
         free(backwardimltb[i]); free(backwardimlix[i]);
         free(backwardimlmi[i]); free(backwardimlrg[i]);
         free(backwardimlpk[i]);
         free(crossrertb[i]);    free(crossrerix[i]);    free(crossrermi[i]);
         free(crossrerrg[i]);    free(crossrerpk[i]);
         free(crossreltb[i]);    free(crossrelix[i]);    free(crossrelmi[i]);
         free(crossrelrg[i]);    free(crossrelpk[i]);
         free(crossimrtb[i]);    free(crossimrix[i]);    free(crossimrmi[i]);
         free(crossimrrg[i]);    free(crossimrpk[i]);
         free(crossimltb[i]);    free(crossimlix[i]);    free(crossimlmi[i]);
         free(crossimlrg[i]);    free(crossimlpk[i]);
      }
      free(forwardrertb[nvr-2]); free(forwardrerix[nvr-2]);
      free(forwardrermi[nvr-2]); free(forwardrerrg[nvr-2]);
      free(forwardrerpk[nvr-2]);
      free(forwardreltb[nvr-2]); free(forwardrelix[nvr-2]);
      free(forwardrelmi[nvr-2]); free(forwardrelrg[nvr-2]);
      free(forwardrelpk[nvr-2]);
      free(forwardimrtb[nvr-2]); free(forwardimrix[nvr-2]);
      free(forwardimrmi[nvr-2]); free(forwardimrrg[nvr-2]);
      free(forwardimrpk[nvr-2]);
      free(forwardimltb[nvr-2]); free(forwardimlix[nvr-2]);
      free(forwardimlmi[nvr-2]); free(forwardimlrg[nvr-2]);
      free(forwardimlpk[nvr-2]);
      free(forwardrertb[nvr-1]); free(forwardrerix[nvr-1]);
      free(forwardrermi[nvr-1]); free(forwardrerrg[nvr-1]);
      free(forwardrerpk[nvr-1]);
      free(forwardreltb[nvr-1]); free(forwardrelix[nvr-1]);
      free(forwardrelmi[nvr-1]); free(forwardrelrg[nvr-1]);
      free(forwardrelpk[nvr-1]);
      free(forwardimrtb[nvr-1]); free(forwardimrix[nvr-1]);
      free(forwardimrmi[nvr-1]); free(forwardimrrg[nvr-1]);
      free(forwardimrpk[nvr-1]);
      free(forwardimltb[nvr-1]); free(forwardimlix[nvr-1]);
      free(forwardimlmi[nvr-1]); free(forwardimlrg[nvr-1]);
      free(forwardimlpk[nvr-1]);
      free(forwardrertb);  free(forwardrerix);  free(forwardrermi);
      free(forwardrerrg);  free(forwardrerpk);
      free(forwardreltb);  free(forwardrelix);  free(forwardrelmi);
      free(forwardrelrg);  free(forwardrelpk);
      free(forwardimrtb);  free(forwardimrix);  free(forwardimrmi);
      free(forwardimrrg);  free(forwardimrpk);
      free(forwardimltb);  free(forwardimlix);  free(forwardimlmi);
      free(forwardimlrg);  free(forwardimlpk);
      free(backwardrertb); free(backwardrerix); free(backwardrermi);
      free(backwardrerrg); free(backwardrerpk);
      free(backwardreltb); free(backwardrelix); free(backwardrelmi);
      free(backwardrelrg); free(backwardrelpk);
      free(backwardimrtb); free(backwardimrix); free(backwardimrmi);
      free(backwardimrrg); free(backwardimrpk);
      free(backwardimltb); free(backwardimlix); free(backwardimlmi);
      free(backwardimlrg); free(backwardimlpk);
      free(crossrertb);    free(crossrerix);    free(crossrermi);
      free(crossrerrg);    free(crossrerpk);
      free(crossreltb);    free(crossrelix);    free(crossrelmi);
      free(crossrelrg);    free(crossrelpk);
      free(crossimrtb);    free(crossimrix);    free(crossimrmi);
      free(crossimrrg);    free(crossimrpk);
      free(crossimltb);    free(crossimlix);    free(crossimlmi);
      free(crossimlrg);    free(crossimlpk);
   }
}
