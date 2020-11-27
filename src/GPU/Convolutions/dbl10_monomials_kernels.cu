// The file dbl10_monomials_kernels.cu defines the kernels specified
// in dbl10_monomials_kernels.h.

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

#include "dbl10_convolutions_kernels.cu"
#include "dbl10_monomials_kernels.h"

void GPU_dbl10_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double *inputrtb, double *inputrix, double *inputrmi, double *inputrrg,
   double *inputrpk, double *inputltb, double *inputlix, double *inputlmi,
   double *inputlrg, double *inputlpk,
   double *forwardrtb, double *forwardrix, double *forwardrmi,
   double *forwardrrg, double *forwardrpk,
   double *forwardltb, double *forwardlix, double *forwardlmi,
   double *forwardlrg, double *forwardlpk,
   double *backwardrtb, double *backwardrix, double *backwardrmi,
   double *backwardrrg, double *backwardrpk,
   double *backwardltb, double *backwardlix, double *backwardlmi,
   double *backwardlrg, double *backwardlpk,
   double *crossrtb, double *crossrix, double *crossrmi, double *crossrrg,
   double *crossrpk, double *crossltb, double *crosslix, double *crosslmi,
   double *crosslrg, double *crosslpk )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl10_padded_convolute<<<1,BS>>>
      (cffrtb,cffrix,cffrmi,cffrrg,cffrpk,cffltb,cfflix,cfflmi,cfflrg,cfflpk,
       &inputrtb[ix1],&inputrix[ix1],&inputrmi[ix1],
       &inputrrg[ix1],&inputrpk[ix1],
       &inputltb[ix1],&inputlix[ix1],&inputlmi[ix1],
       &inputlrg[ix1],&inputlpk[ix1],
       forwardrtb,forwardrix,forwardrmi,forwardrrg,forwardrpk,
       forwardltb,forwardlix,forwardlmi,forwardlrg,forwardlpk,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl10_padded_convolute<<<1,BS>>>
         (&forwardrtb[ix1],&forwardrix[ix1],&forwardrmi[ix1],
          &forwardrrg[ix1],&forwardrpk[ix1],
          &forwardltb[ix1],&forwardlix[ix1],&forwardlmi[ix1],
          &forwardlrg[ix1],&forwardlpk[ix1],
          &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
          &inputrrg[ix2],&inputrpk[ix2],
          &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
          &inputlrg[ix2],&inputlpk[ix2],
          &forwardrtb[ix3],&forwardrix[ix3],&forwardrmi[ix3],
          &forwardrrg[ix3],&forwardrpk[ix3],
          &forwardltb[ix3],&forwardlix[ix3],&forwardlmi[ix3],
          &forwardlrg[ix3],&forwardlpk[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl10_padded_convolute<<<1,BS>>>
         (&inputrtb[ix1],&inputrix[ix1],&inputrmi[ix1],
          &inputrrg[ix1],&inputrpk[ix1],
          &inputltb[ix1],&inputlix[ix1],&inputlmi[ix1],
          &inputlrg[ix1],&inputlpk[ix1],
          &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
          &inputrrg[ix2],&inputrpk[ix2],
          &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
          &inputlrg[ix2],&inputlpk[ix2],
          backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
          backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl10_padded_convolute<<<1,BS>>>
            (&backwardrtb[ix1],&backwardrix[ix1],&backwardrmi[ix1],
             &backwardrrg[ix1],&backwardrpk[ix1],
             &backwardltb[ix1],&backwardlix[ix1],&backwardlmi[ix1],
             &backwardlrg[ix1],&backwardlpk[ix1],
             &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
             &inputrrg[ix2],&inputrpk[ix2],
             &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
             &inputlrg[ix2],&inputlpk[ix2],
             &backwardrtb[ix3],&backwardrix[ix3],&backwardrmi[ix3],
             &backwardrrg[ix3],&backwardrpk[ix3],
             &backwardltb[ix3],&backwardlix[ix3],&backwardlmi[ix3],
             &backwardlrg[ix3],&backwardlpk[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl10_padded_convolute<<<1,BS>>>
         (&backwardrtb[ix3],&backwardrix[ix3],&backwardrmi[ix3],
          &backwardrrg[ix3],&backwardrpk[ix3],
          &backwardltb[ix3],&backwardlix[ix3],&backwardlmi[ix3],
          &backwardlrg[ix3],&backwardlpk[ix3],
          cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
          cffltb,cfflix,cfflmi,cfflrg,cfflpk,
          &backwardrtb[ix2],&backwardrix[ix2],&backwardrmi[ix2],
          &backwardrrg[ix2],&backwardrpk[ix2],
          &backwardltb[ix2],&backwardlix[ix2],&backwardlmi[ix2],
          &backwardlrg[ix2],&backwardlpk[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl10_padded_convolute<<<1,BS>>>
            (forwardrtb,forwardrix,forwardrmi,forwardrrg,forwardrpk,
             forwardltb,forwardlix,forwardlmi,forwardlrg,forwardlpk,
             &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
             &inputrrg[ix2],&inputrpk[ix2],
             &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
             &inputlrg[ix2],&inputlpk[ix2],
             crossrtb,crossrix,crossrmi,crossrrg,crossrpk,
             crossltb,crosslix,crosslmi,crosslrg,crosslpk,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl10_padded_convolute<<<1,BS>>>
               (&forwardrtb[ix1],&forwardrix[ix1],&forwardrmi[ix1],
                &forwardrrg[ix1],&forwardrpk[ix1],
                &forwardltb[ix1],&forwardlix[ix1],&forwardlmi[ix1],
                &forwardlrg[ix1],&forwardlpk[ix1],
                &backwardrtb[ix2],&backwardrix[ix2],&backwardrmi[ix2],
                &backwardrrg[ix2],&backwardrpk[ix2],
                &backwardltb[ix2],&backwardlix[ix2],&backwardlmi[ix2],
                &backwardlrg[ix2],&backwardlpk[ix2],
                &crossrtb[ix1],&crossrix[ix1],&crossrmi[ix1],
                &crossrrg[ix1],&crossrpk[ix1],
                &crossltb[ix1],&crosslix[ix1],&crosslmi[ix1],
                &crosslrg[ix1],&crosslpk[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl10_padded_convolute<<<1,BS>>>
            (&forwardrtb[ix1],&forwardrix[ix1],&forwardrmi[ix1],
             &forwardrrg[ix1],&forwardrpk[ix1],
             &forwardltb[ix1],&forwardlix[ix1],&forwardlmi[ix1],
             &forwardlrg[ix1],&forwardlpk[ix1],
             &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
             &inputrrg[ix2],&inputrpk[ix2],
             &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
             &inputlrg[ix2],&inputlpk[ix2],
             &crossrtb[ix1],&crossrix[ix1],&crossrmi[ix1],
             &crossrrg[ix1],&crossrpk[ix1],
             &crossltb[ix1],&crosslix[ix1],&crosslmi[ix1],
             &crosslrg[ix1],&crosslpk[ix1],deg1);
      }
   }
}

void GPU_cmplx10_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double *inputrertb, double *inputrerix, double *inputrermi,
   double *inputrerrg, double *inputrerpk,
   double *inputreltb, double *inputrelix, double *inputrelmi,
   double *inputrelrg, double *inputrelpk,
   double *inputimrtb, double *inputimrix, double *inputimrmi,
   double *inputimrrg, double *inputimrpk,
   double *inputimltb, double *inputimlix, double *inputimlmi,
   double *inputimlrg, double *inputimlpk,
   double *forwardrertb, double *forwardrerix, double *forwardrermi,
   double *forwardrerrg, double *forwardrerpk,
   double *forwardreltb, double *forwardrelix, double *forwardrelmi,
   double *forwardrelrg, double *forwardrelpk,
   double *forwardimrtb, double *forwardimrix, double *forwardimrmi,
   double *forwardimrrg, double *forwardimrpk,
   double *forwardimltb, double *forwardimlix, double *forwardimlmi,
   double *forwardimlrg, double *forwardimlpk,
   double *backwardrertb, double *backwardrerix, double *backwardrermi,
   double *backwardrerrg, double *backwardrerpk,
   double *backwardreltb, double *backwardrelix, double *backwardrelmi,
   double *backwardrelrg, double *backwardrelpk,
   double *backwardimrtb, double *backwardimrix, double *backwardimrmi,
   double *backwardimrrg, double *backwardimrpk,
   double *backwardimltb, double *backwardimlix, double *backwardimlmi,
   double *backwardimlrg, double *backwardimlpk,
   double *crossrertb, double *crossrerix, double *crossrermi,
   double *crossrerrg, double *crossrerpk,
   double *crossreltb, double *crossrelix, double *crossrelmi,
   double *crossrelrg, double *crossrelpk,
   double *crossimrtb, double *crossimrix, double *crossimrmi,
   double *crossimrrg, double *crossimrpk,
   double *crossimltb, double *crossimlix, double *crossimlmi,
   double *crossimlrg, double *crossimlpk )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   cmplx10_padded_convolute<<<1,BS>>>
      (cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
       cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
       cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
       cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
       &inputrertb[ix1],&inputrerix[ix1],&inputrermi[ix1],
       &inputrerrg[ix1],&inputrerpk[ix1],
       &inputreltb[ix1],&inputrelix[ix1],&inputrelmi[ix1],
       &inputrelrg[ix1],&inputrelpk[ix1],
       &inputimrtb[ix1],&inputimrix[ix1],&inputimrmi[ix1],
       &inputimrrg[ix1],&inputimrpk[ix1],
       &inputimltb[ix1],&inputimlix[ix1],&inputimlmi[ix1],
       &inputimlrg[ix1],&inputimlpk[ix1],
       forwardrertb,forwardrerix,forwardrermi,forwardrerrg,forwardrerpk,
       forwardreltb,forwardrelix,forwardrelmi,forwardrelrg,forwardrelpk,
       forwardimrtb,forwardimrix,forwardimrmi,forwardimrrg,forwardimrpk,
       forwardimltb,forwardimlix,forwardimlmi,forwardimlrg,forwardimlpk,deg1); 

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-i]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      cmplx10_padded_convolute<<<1,BS>>>
         (&forwardrertb[ix1],&forwardrerix[ix1],&forwardrermi[ix1],
          &forwardrerrg[ix1],&forwardrerpk[ix1],
          &forwardreltb[ix1],&forwardrelix[ix1],&forwardrelmi[ix1],
          &forwardrelrg[ix1],&forwardrelpk[ix1],
          &forwardimrtb[ix1],&forwardimrix[ix1],&forwardimrmi[ix1],
          &forwardimrrg[ix1],&forwardimrpk[ix1],
          &forwardimltb[ix1],&forwardimlix[ix1],&forwardimlmi[ix1],
          &forwardimlrg[ix1],&forwardimlpk[ix1],
          &inputrertb[ix2],&inputrerix[ix2],&inputrermi[ix2],
          &inputrerrg[ix2],&inputrerpk[ix2],
          &inputreltb[ix2],&inputrelix[ix2],&inputrelmi[ix2],
          &inputrelrg[ix2],&inputrelpk[ix2],
          &inputimrtb[ix2],&inputimrix[ix2],&inputimrmi[ix2],
          &inputimrrg[ix2],&inputimrpk[ix2],
          &inputimltb[ix2],&inputimlix[ix2],&inputimlmi[ix2],
          &inputimlrg[ix2],&inputimlpk[ix2],
          &forwardrertb[ix3],&forwardrerix[ix3],&forwardrermi[ix3],
          &forwardrerrg[ix3],&forwardrerpk[ix3],
          &forwardreltb[ix3],&forwardrelix[ix3],&forwardrelmi[ix3],
          &forwardrelrg[ix3],&forwardrelpk[ix3],
          &forwardimrtb[ix3],&forwardimrix[ix3],&forwardimrmi[ix3],
          &forwardimrrg[ix3],&forwardimrpk[ix3],
          &forwardimltb[ix3],&forwardimlix[ix3],&forwardimlmi[ix3],
          &forwardimlrg[ix3],&forwardimlpk[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      cmplx10_padded_convolute<<<1,BS>>>
         (&inputrertb[ix1],&inputrerix[ix1],&inputrermi[ix1],
          &inputrerrg[ix1],&inputrerpk[ix1],
          &inputreltb[ix1],&inputrelix[ix1],&inputrelmi[ix1],
          &inputrelrg[ix1],&inputrelpk[ix1],
          &inputimrtb[ix1],&inputimrix[ix1],&inputimrmi[ix1],
          &inputimrrg[ix1],&inputimrpk[ix1],
          &inputimltb[ix1],&inputimlix[ix1],&inputimlmi[ix1],
          &inputimlrg[ix1],&inputimlpk[ix1],
          &inputrertb[ix2],&inputrerix[ix2],&inputrermi[ix2],
          &inputrerrg[ix2],&inputrerpk[ix2],
          &inputreltb[ix2],&inputrelix[ix2],&inputrelmi[ix2],
          &inputrelrg[ix2],&inputrelpk[ix2],
          &inputimrtb[ix2],&inputimrix[ix2],&inputimrmi[ix2],
          &inputimrrg[ix2],&inputimrpk[ix2],
          &inputimltb[ix2],&inputimlix[ix2],&inputimlmi[ix2],
          &inputimlrg[ix2],&inputimlpk[ix2],
          backwardrertb,backwardrerix,backwardrermi,
          backwardrerrg,backwardrerpk,
          backwardreltb,backwardrelix,backwardrelmi,
          backwardrelrg,backwardrelpk,
          backwardimrtb,backwardimrix,backwardimrmi,
          backwardimrrg,backwardimrpk,
          backwardimltb,backwardimlix,backwardimlmi,
          backwardimlrg,backwardimlpk,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         cmplx10_padded_convolute<<<1,BS>>>
            (&backwardrertb[ix1],&backwardrerix[ix1],&backwardrermi[ix1],
             &backwardrerrg[ix1],&backwardrerpk[ix1],
             &backwardreltb[ix1],&backwardrelix[ix1],&backwardrelmi[ix1],
             &backwardrelrg[ix1],&backwardrelpk[ix1],
             &backwardimrtb[ix1],&backwardimrix[ix1],&backwardimrmi[ix1],
             &backwardimrrg[ix1],&backwardimrpk[ix1],
             &backwardimltb[ix1],&backwardimlix[ix1],&backwardimlmi[ix1],
             &backwardimlrg[ix1],&backwardimlpk[ix1],
             &inputrertb[ix2],&inputrerix[ix2],&inputrermi[ix2],
             &inputrerrg[ix2],&inputrerpk[ix2],
             &inputreltb[ix2],&inputrelix[ix2],&inputrelmi[ix2],
             &inputrelrg[ix2],&inputrelpk[ix2],
             &inputimrtb[ix2],&inputimrix[ix2],&inputimrmi[ix2],
             &inputimrrg[ix2],&inputimrpk[ix2],
             &inputimltb[ix2],&inputimlix[ix2],&inputimlmi[ix2],
             &inputimlrg[ix2],&inputimlpk[ix2],
             &backwardrertb[ix3],&backwardrerix[ix3],&backwardrermi[ix3],
             &backwardrerrg[ix3],&backwardrerpk[ix3],
             &backwardreltb[ix3],&backwardrelix[ix3],&backwardrelmi[ix3],
             &backwardrelrg[ix3],&backwardrelpk[ix3],
             &backwardimrtb[ix3],&backwardimrix[ix3],&backwardimrmi[ix3],
             &backwardimrrg[ix3],&backwardimrpk[ix3],
             &backwardimltb[ix3],&backwardimlix[ix3],&backwardimlmi[ix3],
             &backwardimlrg[ix3],&backwardimlpk[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      cmplx10_padded_convolute<<<1,BS>>>
         (&backwardrertb[ix3],&backwardrerix[ix3],&backwardrermi[ix3],
          &backwardrerrg[ix3],&backwardrerpk[ix3],
          &backwardreltb[ix3],&backwardrelix[ix3],&backwardrelmi[ix3],
          &backwardrelrg[ix3],&backwardrelpk[ix3],
          &backwardimrtb[ix3],&backwardimrix[ix3],&backwardimrmi[ix3],
          &backwardimrrg[ix3],&backwardimrpk[ix3],
          &backwardimltb[ix3],&backwardimlix[ix3],&backwardimlmi[ix3],
          &backwardimlrg[ix3],&backwardimlpk[ix3],
          cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
          cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
          cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
          cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
          &backwardrertb[ix2],&backwardrerix[ix2],&backwardrermi[ix2],
          &backwardrerrg[ix2],&backwardrerpk[ix2],
          &backwardreltb[ix2],&backwardrelix[ix2],&backwardrelmi[ix2],
          &backwardrelrg[ix2],&backwardrelpk[ix2],
          &backwardimrtb[ix2],&backwardimrix[ix2],&backwardimrmi[ix2],
          &backwardimrrg[ix2],&backwardimrpk[ix2],
          &backwardimltb[ix2],&backwardimlix[ix2],&backwardimlmi[ix2],
          &backwardimlrg[ix2],&backwardimlpk[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         cmplx10_padded_convolute<<<1,BS>>>
            (forwardrertb,forwardrerix,forwardrermi,forwardrerrg,forwardrerpk,
             forwardreltb,forwardrelix,forwardrelmi,forwardrelrg,forwardrelpk,
             forwardimrtb,forwardimrix,forwardimrmi,forwardimrrg,forwardimrpk,
             forwardimltb,forwardimlix,forwardimlmi,forwardimlrg,forwardimlpk,
             &inputrertb[ix2],&inputrerix[ix2],&inputrermi[ix2],
             &inputrerrg[ix2],&inputrerpk[ix2],
             &inputreltb[ix2],&inputrelix[ix2],&inputrelmi[ix2],
             &inputrelrg[ix2],&inputrelpk[ix2],
             &inputimrtb[ix2],&inputimrix[ix2],&inputimrmi[ix2],
             &inputimrrg[ix2],&inputimrpk[ix2],
             &inputimltb[ix2],&inputimlix[ix2],&inputimlmi[ix2],
             &inputimlrg[ix2],&inputimlpk[ix2],
             crossrertb,crossrerix,crossrermi,crossrerrg,crossrerpk,
             crossreltb,crossrelix,crossrelmi,crossrelrg,crossrelpk,
             crossimrtb,crossimrix,crossimrmi,crossimrrg,crossimrpk,
             crossimltb,crossimlix,crossimlmi,crossimlrg,crossimlpk,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            cmplx10_padded_convolute<<<1,BS>>>
               (&forwardrertb[ix1],&forwardrerix[ix1],&forwardrermi[ix1],
                &forwardrerrg[ix1],&forwardrerpk[ix1],
                &forwardreltb[ix1],&forwardrelix[ix1],&forwardrelmi[ix1],
                &forwardrelrg[ix1],&forwardrelpk[ix1],
                &forwardimrtb[ix1],&forwardimrix[ix1],&forwardimrmi[ix1],
                &forwardimrrg[ix1],&forwardimrpk[ix1],
                &forwardimltb[ix1],&forwardimlix[ix1],&forwardimlmi[ix1],
                &forwardimlrg[ix1],&forwardimlpk[ix1],
                &backwardrertb[ix2],&backwardrerix[ix2],&backwardrermi[ix2],
                &backwardrerrg[ix2],&backwardrerpk[ix2],
                &backwardreltb[ix2],&backwardrelix[ix2],&backwardrelmi[ix2],
                &backwardrelrg[ix2],&backwardrelpk[ix2],
                &backwardimrtb[ix2],&backwardimrix[ix2],&backwardimrmi[ix2],
                &backwardimrrg[ix2],&backwardimrpk[ix2],
                &backwardimltb[ix2],&backwardimlix[ix2],&backwardimlmi[ix2],
                &backwardimlrg[ix2],&backwardimlpk[ix2],
                &crossrertb[ix1],&crossrerix[ix1],&crossrermi[ix1],
                &crossrerrg[ix1],&crossrerpk[ix1],
                &crossreltb[ix1],&crossrelix[ix1],&crossrelmi[ix1],
                &crossrelrg[ix1],&crossrelpk[ix1],
                &crossimrtb[ix1],&crossimrix[ix1],&crossimrmi[ix1],
                &crossimrrg[ix1],&crossimrpk[ix1],
                &crossimltb[ix1],&crossimlix[ix1],&crossimlmi[ix1],
                &crossimlrg[ix1],&crossimlpk[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         cmplx10_padded_convolute<<<1,BS>>>
            (&forwardrertb[ix1],&forwardrerix[ix1],&forwardrermi[ix1],
             &forwardrerrg[ix1],&forwardrerpk[ix1],
             &forwardreltb[ix1],&forwardrelix[ix1],&forwardrelmi[ix1],
             &forwardrelrg[ix1],&forwardrelpk[ix1],
             &forwardimrtb[ix1],&forwardimrix[ix1],&forwardimrmi[ix1],
             &forwardimrrg[ix1],&forwardimrpk[ix1],
             &forwardimltb[ix1],&forwardimlix[ix1],&forwardimlmi[ix1],
             &forwardimlrg[ix1],&forwardimlpk[ix1],
             &inputrertb[ix2],&inputrerix[ix2],&inputrermi[ix2],
             &inputrerrg[ix2],&inputrerpk[ix2],
             &inputreltb[ix2],&inputrelix[ix2],&inputrelmi[ix2],
             &inputrelrg[ix2],&inputrelpk[ix2],
             &inputimrtb[ix2],&inputimrix[ix2],&inputimrmi[ix2],
             &inputimrrg[ix2],&inputimrpk[ix2],
             &inputimltb[ix2],&inputimlix[ix2],&inputimlmi[ix2],
             &inputimlrg[ix2],&inputimlpk[ix2],
             &crossrertb[ix1],&crossrerix[ix1],&crossrermi[ix1],
             &crossrerrg[ix1],&crossrerpk[ix1],
             &crossreltb[ix1],&crossrelix[ix1],&crossrelmi[ix1],
             &crossrelrg[ix1],&crossrelpk[ix1],
             &crossimrtb[ix1],&crossimrix[ix1],&crossimrmi[ix1],
             &crossimrrg[ix1],&crossimrpk[ix1],
             &crossimltb[ix1],&crossimlix[ix1],&crossimlmi[ix1],
             &crossimlrg[ix1],&crossimlpk[ix1],deg1);
      }
   }
}

void GPU_dbl10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
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
   const int deg1 = deg+1;      // length of all vectors
   double *inputrtb_d;          // inputrtb_d is inputrtb on the device
   double *inputrix_d;          // inputrix_d is inputrix on the device
   double *inputrmi_d;          // inputrmi_d is inputrmi on the device
   double *inputrrg_d;          // inputrrg_d is inputrrg on the device
   double *inputrpk_d;          // inputrpk_d is inputrpk on the device
   double *inputltb_d;          // inputltb_d is inputltb on the device
   double *inputlix_d;          // inputlix_d is inputlix on the device
   double *inputlmi_d;          // inputlmi_d is inputlmi on the device
   double *inputlrg_d;          // inputlrg_d is inputlrg on the device
   double *inputlpk_d;          // inputlpk_d is inputlpk on the device
   double *forwardrtb_d;        // highest forward products on the device
   double *forwardrix_d;        // second highest forward products
   double *forwardrmi_d;        // third highest forward products
   double *forwardrrg_d;        // fourth highest forward products
   double *forwardrpk_d;        // fifth highest forward products
   double *forwardltb_d;        // fifth lowest forward products
   double *forwardlix_d;        // fourth lowest forward products
   double *forwardlmi_d;        // third lowest forward products
   double *forwardlrg_d;        // second lowest forward products
   double *forwardlpk_d;        // lowest forward products
   double *backwardrtb_d;       // highest backward products on the device
   double *backwardrix_d;       // second highest backward products
   double *backwardrmi_d;       // third highest backward products
   double *backwardrrg_d;       // fourth highest backward products
   double *backwardrpk_d;       // fifth highest backward products
   double *backwardltb_d;       // fifth lowest backward products
   double *backwardlix_d;       // fourth lowest backward products
   double *backwardlmi_d;       // third lowest backward products
   double *backwardlrg_d;       // second lowest backward products
   double *backwardlpk_d;       // lowest backward products
   double *crossrtb_d;          // highest cross products on the device
   double *crossrix_d;          // second highest cross products
   double *crossrmi_d;          // third highest cross products
   double *crossrrg_d;          // fourth highest cross products
   double *crossrpk_d;          // fifth highest cross products
   double *crossltb_d;          // fifth lowest cross products
   double *crosslix_d;          // fourth lowest cross products
   double *crosslmi_d;          // third lowest cross products
   double *crosslrg_d;          // second lowest cross products
   double *crosslpk_d;          // lowest cross products
   double *cffrtb_d;            // cffrtb_d is cffrtb on device
   double *cffrix_d;            // cffrix_d is cffrix on device
   double *cffrmi_d;            // cffrmi_d is cffrmi on device
   double *cffrrg_d;            // cffrrg_d is cffrrg on device
   double *cffrpk_d;            // cffrpk_d is cffrpk on device
   double *cffltb_d;            // cffltb_d is cffltb on device
   double *cfflix_d;            // cfflix_d is cfflix on device
   double *cfflmi_d;            // cfflmi_d is cfflmi on device
   double *cfflrg_d;            // cfflrg_d is cfflrg on device
   double *cfflpk_d;            // cfflpk_d is cfflpk on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);

   cudaMalloc((void**)&cffrtb_d,szcff);
   cudaMalloc((void**)&cffrix_d,szcff);
   cudaMalloc((void**)&cffrmi_d,szcff);
   cudaMalloc((void**)&cffrrg_d,szcff);
   cudaMalloc((void**)&cffrpk_d,szcff);
   cudaMalloc((void**)&cffltb_d,szcff);
   cudaMalloc((void**)&cfflix_d,szcff);
   cudaMalloc((void**)&cfflmi_d,szcff);
   cudaMalloc((void**)&cfflrg_d,szcff);
   cudaMalloc((void**)&cfflpk_d,szcff);
   cudaMalloc((void**)&inputrtb_d,szdim);
   cudaMalloc((void**)&inputrix_d,szdim);
   cudaMalloc((void**)&inputrmi_d,szdim);
   cudaMalloc((void**)&inputrrg_d,szdim);
   cudaMalloc((void**)&inputrpk_d,szdim);
   cudaMalloc((void**)&inputltb_d,szdim);
   cudaMalloc((void**)&inputlix_d,szdim);
   cudaMalloc((void**)&inputlmi_d,szdim);
   cudaMalloc((void**)&inputlrg_d,szdim);
   cudaMalloc((void**)&inputlpk_d,szdim);
   cudaMalloc((void**)&forwardrtb_d,sznvr);
   cudaMalloc((void**)&forwardrix_d,sznvr);
   cudaMalloc((void**)&forwardrmi_d,sznvr);
   cudaMalloc((void**)&forwardrrg_d,sznvr);
   cudaMalloc((void**)&forwardrpk_d,sznvr);
   cudaMalloc((void**)&forwardltb_d,sznvr);
   cudaMalloc((void**)&forwardlix_d,sznvr);
   cudaMalloc((void**)&forwardlmi_d,sznvr);
   cudaMalloc((void**)&forwardlrg_d,sznvr);
   cudaMalloc((void**)&forwardlpk_d,sznvr);
   cudaMalloc((void**)&backwardrtb_d,sznvr1);
   cudaMalloc((void**)&backwardrix_d,sznvr1);
   cudaMalloc((void**)&backwardrmi_d,sznvr1);
   cudaMalloc((void**)&backwardrrg_d,sznvr1);
   cudaMalloc((void**)&backwardrpk_d,sznvr1);
   cudaMalloc((void**)&backwardltb_d,sznvr1);
   cudaMalloc((void**)&backwardlix_d,sznvr1);
   cudaMalloc((void**)&backwardlmi_d,sznvr1);
   cudaMalloc((void**)&backwardlrg_d,sznvr1);
   cudaMalloc((void**)&backwardlpk_d,sznvr1);
   cudaMalloc((void**)&crossrtb_d,sznvr2);
   cudaMalloc((void**)&crossrix_d,sznvr2);
   cudaMalloc((void**)&crossrmi_d,sznvr2);
   cudaMalloc((void**)&crossrrg_d,sznvr2);
   cudaMalloc((void**)&crossrpk_d,sznvr2);
   cudaMalloc((void**)&crossltb_d,sznvr2);
   cudaMalloc((void**)&crosslix_d,sznvr2);
   cudaMalloc((void**)&crosslmi_d,sznvr2);
   cudaMalloc((void**)&crosslrg_d,sznvr2);
   cudaMalloc((void**)&crosslpk_d,sznvr2);

   double *inputrtb_h = new double[dim*(deg1)];
   double *inputrix_h = new double[dim*(deg1)];
   double *inputrmi_h = new double[dim*(deg1)];
   double *inputrrg_h = new double[dim*(deg1)];
   double *inputrpk_h = new double[dim*(deg1)];
   double *inputltb_h = new double[dim*(deg1)];
   double *inputlix_h = new double[dim*(deg1)];
   double *inputlmi_h = new double[dim*(deg1)];
   double *inputlrg_h = new double[dim*(deg1)];
   double *inputlpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrtb_h[ix] = inputrtb[i][j];
         inputrix_h[ix] = inputrix[i][j];
         inputrmi_h[ix] = inputrmi[i][j];
         inputrrg_h[ix] = inputrrg[i][j];
         inputrpk_h[ix] = inputrpk[i][j];
         inputltb_h[ix] = inputltb[i][j];
         inputlix_h[ix] = inputlix[i][j];
         inputlmi_h[ix] = inputlmi[i][j];
         inputlrg_h[ix] = inputlrg[i][j];
         inputlpk_h[ix++] = inputlpk[i][j];
      }

   cudaMemcpy(cffrtb_d,cffrtb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrix_d,cffrix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrmi_d,cffrmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrrg_d,cffrrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrpk_d,cffrpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffltb_d,cffltb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflix_d,cfflix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflmi_d,cfflmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflrg_d,cfflrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflpk_d,cfflpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrtb_d,inputrtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrix_d,inputrix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrmi_d,inputrmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrrg_d,inputrrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrpk_d,inputrpk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputltb_d,inputltb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlix_d,inputlix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlmi_d,inputlmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlrg_d,inputlrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlpk_d,inputlpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl10_speel
         (BS,nvr,deg,idx,
          cffrtb_d,cffrix_d,cffrmi_d,cffrrg_d,cffrpk_d,
          cffltb_d,cfflix_d,cfflmi_d,cfflrg_d,cfflpk_d,
          inputrtb_d,inputrix_d,inputrmi_d,inputrrg_d,inputrpk_d,
          inputltb_d,inputlix_d,inputlmi_d,inputlrg_d,inputlpk_d,
          forwardrtb_d,forwardrix_d,forwardrmi_d,forwardrrg_d,forwardrpk_d,
          forwardltb_d,forwardlix_d,forwardlmi_d,forwardlrg_d,forwardlpk_d,
          backwardrtb_d,backwardrix_d,backwardrmi_d,backwardrrg_d,
          backwardrpk_d,backwardltb_d,backwardlix_d,backwardlmi_d,
          backwardlrg_d,backwardlpk_d,
          crossrtb_d,crossrix_d,crossrmi_d,crossrrg_d,crossrpk_d,
          crossltb_d,crosslix_d,crosslmi_d,crosslrg_d,crosslpk_d);
   }
   double *forwardrtb_h = new double[(deg1)*nvr];
   double *forwardrix_h = new double[(deg1)*nvr];
   double *forwardrmi_h = new double[(deg1)*nvr];
   double *forwardrrg_h = new double[(deg1)*nvr];
   double *forwardrpk_h = new double[(deg1)*nvr];
   double *forwardltb_h = new double[(deg1)*nvr];
   double *forwardlix_h = new double[(deg1)*nvr];
   double *forwardlmi_h = new double[(deg1)*nvr];
   double *forwardlrg_h = new double[(deg1)*nvr];
   double *forwardlpk_h = new double[(deg1)*nvr];
   double *backwardrtb_h = new double[(deg1)*(nvr-1)];
   double *backwardrix_h = new double[(deg1)*(nvr-1)];
   double *backwardrmi_h = new double[(deg1)*(nvr-1)];
   double *backwardrrg_h = new double[(deg1)*(nvr-1)];
   double *backwardrpk_h = new double[(deg1)*(nvr-1)];
   double *backwardltb_h = new double[(deg1)*(nvr-1)];
   double *backwardlix_h = new double[(deg1)*(nvr-1)];
   double *backwardlmi_h = new double[(deg1)*(nvr-1)];
   double *backwardlrg_h = new double[(deg1)*(nvr-1)];
   double *backwardlpk_h = new double[(deg1)*(nvr-1)];
   double *crossrtb_h = new double[(deg1)*(nvr-2)];
   double *crossrix_h = new double[(deg1)*(nvr-2)];
   double *crossrmi_h = new double[(deg1)*(nvr-2)];
   double *crossrrg_h = new double[(deg1)*(nvr-2)];
   double *crossrpk_h = new double[(deg1)*(nvr-2)];
   double *crossltb_h = new double[(deg1)*(nvr-2)];
   double *crosslix_h = new double[(deg1)*(nvr-2)];
   double *crosslmi_h = new double[(deg1)*(nvr-2)];
   double *crosslrg_h = new double[(deg1)*(nvr-2)];
   double *crosslpk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrtb_h,forwardrtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrix_h,forwardrix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrmi_h,forwardrmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrrg_h,forwardrrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrpk_h,forwardrpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardltb_h,forwardltb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlix_h,forwardlix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlmi_h,forwardlmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlrg_h,forwardlrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlpk_h,forwardlpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrtb_h,backwardrtb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrix_h,backwardrix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrmi_h,backwardrmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrrg_h,backwardrrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrpk_h,backwardrpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardltb_h,backwardltb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlix_h,backwardlix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlmi_h,backwardlmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlrg_h,backwardlrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlpk_h,backwardlpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrtb_h,crossrtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrix_h,crossrix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrmi_h,crossrmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrrg_h,crossrrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrpk_h,crossrpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossltb_h,crossltb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslix_h,crosslix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslmi_h,crosslmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslrg_h,crosslrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslpk_h,crosslpk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputrtb[dim][i] = forwardrtb_h[offset+i];
      outputrix[dim][i] = forwardrix_h[offset+i];
      outputrmi[dim][i] = forwardrmi_h[offset+i];
      outputrrg[dim][i] = forwardrrg_h[offset+i];
      outputrpk[dim][i] = forwardrpk_h[offset+i];
      outputltb[dim][i] = forwardltb_h[offset+i];
      outputlix[dim][i] = forwardlix_h[offset+i];
      outputlmi[dim][i] = forwardlmi_h[offset+i];
      outputlrg[dim][i] = forwardlrg_h[offset+i];
      outputlpk[dim][i] = forwardlpk_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputrtb[ix][i] = forwardrtb_h[offset+i];
      outputrix[ix][i] = forwardrix_h[offset+i];
      outputrmi[ix][i] = forwardrmi_h[offset+i];
      outputrrg[ix][i] = forwardrrg_h[offset+i];
      outputrpk[ix][i] = forwardrpk_h[offset+i];
      outputltb[ix][i] = forwardltb_h[offset+i];
      outputlix[ix][i] = forwardlix_h[offset+i];
      outputlmi[ix][i] = forwardlmi_h[offset+i];
      outputlrg[ix][i] = forwardlrg_h[offset+i];
      outputlpk[ix][i] = forwardlpk_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputrtb[ix][i] = backwardrtb_h[offset+i];
      outputrix[ix][i] = backwardrix_h[offset+i];
      outputrmi[ix][i] = backwardrmi_h[offset+i];
      outputrrg[ix][i] = backwardrrg_h[offset+i];
      outputrpk[ix][i] = backwardrpk_h[offset+i];
      outputltb[ix][i] = backwardltb_h[offset+i];
      outputlix[ix][i] = backwardlix_h[offset+i];
      outputlmi[ix][i] = backwardlmi_h[offset+i];
      outputlrg[ix][i] = backwardlrg_h[offset+i];
      outputlpk[ix][i] = backwardlpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrtb[ix][i] = crossrtb_h[offset+i];
         outputrix[ix][i] = crossrix_h[offset+i];
         outputrmi[ix][i] = crossrmi_h[offset+i];
         outputrrg[ix][i] = crossrrg_h[offset+i];
         outputrpk[ix][i] = crossrpk_h[offset+i];
         outputltb[ix][i] = crossltb_h[offset+i];
         outputlix[ix][i] = crosslix_h[offset+i];
         outputlmi[ix][i] = crosslmi_h[offset+i];
         outputlrg[ix][i] = crosslrg_h[offset+i];
         outputlpk[ix][i] = crosslpk_h[offset+i];
      }
   }
}

void GPU_cmplx10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi,
   double *cffrerrg, double *cffrerpk,
   double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi,
   double *cffimrrg, double *cffimrpk,
   double *cffimltb, double *cffimlix, double *cffimlmi,
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
   const int deg1 = deg+1;    // length of all vectors
   double *inputrertb_d;      // inputrertb_d is inputrertb on the device
   double *inputrerix_d;      // inputrerix_d is inputrerix on the device
   double *inputrermi_d;      // inputrermi_d is inputrermi on the device
   double *inputrerrg_d;      // inputrerrg_d is inputrerrg on the device
   double *inputrerpk_d;      // inputrerpk_d is inputrerpk on the device
   double *inputreltb_d;      // inputreltb_d is inputreltb on the device
   double *inputrelix_d;      // inputrelix_d is inputrelix on the device
   double *inputrelmi_d;      // inputrelmi_d is inputrelmi on the device
   double *inputrelrg_d;      // inputrelrg_d is inputrelrg on the device
   double *inputrelpk_d;      // inputrelpk_d is inputrelpk on the device
   double *inputimrtb_d;      // inputimrtb_d is inputrertb on the device
   double *inputimrix_d;      // inputimrix_d is inputrerix on the device
   double *inputimrmi_d;      // inputimrmi_d is inputrermi on the device
   double *inputimrrg_d;      // inputimrrg_d is inputrerrg on the device
   double *inputimrpk_d;      // inputimrpk_d is inputrerpk on the device
   double *inputimltb_d;      // inputimltb_d is inputreltb on the device
   double *inputimlix_d;      // inputimlix_d is inputrelix on the device
   double *inputimlmi_d;      // inputimlmi_d is inputrelmi on the device
   double *inputimlrg_d;      // inputimlrg_d is inputrelrg on the device
   double *inputimlpk_d;      // inputimlpk_d is inputrelpk on the device
   double *forwardrertb_d;
   double *forwardrerix_d;
   double *forwardrermi_d;
   double *forwardrerrg_d;
   double *forwardrerpk_d;
   double *forwardreltb_d;
   double *forwardrelix_d;
   double *forwardrelmi_d;
   double *forwardrelrg_d;
   double *forwardrelpk_d;
   double *forwardimrtb_d;
   double *forwardimrix_d;
   double *forwardimrmi_d;
   double *forwardimrrg_d;
   double *forwardimrpk_d;
   double *forwardimltb_d;
   double *forwardimlix_d;
   double *forwardimlmi_d;
   double *forwardimlrg_d;
   double *forwardimlpk_d;
   double *backwardrertb_d;
   double *backwardrerix_d;
   double *backwardrermi_d;
   double *backwardrerrg_d;
   double *backwardrerpk_d;
   double *backwardreltb_d;
   double *backwardrelix_d;
   double *backwardrelmi_d;
   double *backwardrelrg_d;
   double *backwardrelpk_d;
   double *backwardimrtb_d;
   double *backwardimrix_d;
   double *backwardimrmi_d;
   double *backwardimrrg_d;
   double *backwardimrpk_d;
   double *backwardimltb_d;
   double *backwardimlix_d;
   double *backwardimlmi_d;
   double *backwardimlrg_d;
   double *backwardimlpk_d;
   double *crossrertb_d;
   double *crossrerix_d;
   double *crossrermi_d;
   double *crossrerrg_d;
   double *crossrerpk_d;
   double *crossreltb_d;
   double *crossrelix_d;
   double *crossrelmi_d;
   double *crossrelrg_d;
   double *crossrelpk_d;
   double *crossimrtb_d;
   double *crossimrix_d;
   double *crossimrmi_d;
   double *crossimrrg_d;
   double *crossimrpk_d;
   double *crossimltb_d;
   double *crossimlix_d;
   double *crossimlmi_d;
   double *crossimlrg_d;
   double *crossimlpk_d;
   double *cffrertb_d;         // cffrertb_d is cffrertb on the device
   double *cffrerix_d;         // cffrerix_d is cffrerix on the device
   double *cffrermi_d;         // cffrermi_d is cffrermi on the device
   double *cffrerrg_d;         // cffrermi_d is cffrerrg on the device
   double *cffrerpk_d;         // cffrerpk_d is cffrerpk on the device
   double *cffreltb_d;         // cffreltb_d is cffreltb on the device
   double *cffrelix_d;         // cffrelix_d is cffrelix on the device
   double *cffrelmi_d;         // cffrelmi_d is cffrelmi on the device
   double *cffrelrg_d;         // cffrelmi_d is cffrelrg on the device
   double *cffrelpk_d;         // cffrelpk_d is cffrelpk on the device
   double *cffimrtb_d;         // cffimrtb_d is cffimrtb on the device
   double *cffimrix_d;         // cffimrix_d is cffimrix on the device
   double *cffimrmi_d;         // cffimrmi_d is cffimrmi on the device
   double *cffimrrg_d;         // cffimrrg_d is cffimrrg on the device
   double *cffimrpk_d;         // cffimrpk_d is cffimrpk on the device
   double *cffimltb_d;         // cffimltb_d is cffimltb on the device
   double *cffimlix_d;         // cffimlix_d is cffimlix on the device
   double *cffimlmi_d;         // cffimlmi_d is cffimlmi on the device
   double *cffimlrg_d;         // cffimlrg_d is cffimlrg on the device
   double *cffimlpk_d;         // cffimlpk_d is cffimlpk on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&cffrertb_d,szcff);
   cudaMalloc((void**)&cffrerix_d,szcff);
   cudaMalloc((void**)&cffrermi_d,szcff);
   cudaMalloc((void**)&cffrerrg_d,szcff);
   cudaMalloc((void**)&cffrerpk_d,szcff);
   cudaMalloc((void**)&cffreltb_d,szcff);
   cudaMalloc((void**)&cffrelix_d,szcff);
   cudaMalloc((void**)&cffrelmi_d,szcff);
   cudaMalloc((void**)&cffrelrg_d,szcff);
   cudaMalloc((void**)&cffrelpk_d,szcff);
   cudaMalloc((void**)&cffimrtb_d,szcff);
   cudaMalloc((void**)&cffimrix_d,szcff);
   cudaMalloc((void**)&cffimrmi_d,szcff);
   cudaMalloc((void**)&cffimrrg_d,szcff);
   cudaMalloc((void**)&cffimrpk_d,szcff);
   cudaMalloc((void**)&cffimltb_d,szcff);
   cudaMalloc((void**)&cffimlix_d,szcff);
   cudaMalloc((void**)&cffimlmi_d,szcff);
   cudaMalloc((void**)&cffimlrg_d,szcff);
   cudaMalloc((void**)&cffimlpk_d,szcff);
   cudaMalloc((void**)&inputrertb_d,szdim);
   cudaMalloc((void**)&inputrerix_d,szdim);
   cudaMalloc((void**)&inputrermi_d,szdim);
   cudaMalloc((void**)&inputrerrg_d,szdim);
   cudaMalloc((void**)&inputrerpk_d,szdim);
   cudaMalloc((void**)&inputreltb_d,szdim);
   cudaMalloc((void**)&inputrelix_d,szdim);
   cudaMalloc((void**)&inputrelmi_d,szdim);
   cudaMalloc((void**)&inputrelrg_d,szdim);
   cudaMalloc((void**)&inputrelpk_d,szdim);
   cudaMalloc((void**)&inputimrtb_d,szdim);
   cudaMalloc((void**)&inputimrix_d,szdim);
   cudaMalloc((void**)&inputimrmi_d,szdim);
   cudaMalloc((void**)&inputimrrg_d,szdim);
   cudaMalloc((void**)&inputimrpk_d,szdim);
   cudaMalloc((void**)&inputimltb_d,szdim);
   cudaMalloc((void**)&inputimlix_d,szdim);
   cudaMalloc((void**)&inputimlmi_d,szdim);
   cudaMalloc((void**)&inputimlrg_d,szdim);
   cudaMalloc((void**)&inputimlpk_d,szdim);
   cudaMalloc((void**)&forwardrertb_d,sznvr);
   cudaMalloc((void**)&forwardrerix_d,sznvr);
   cudaMalloc((void**)&forwardrermi_d,sznvr);
   cudaMalloc((void**)&forwardrerrg_d,sznvr);
   cudaMalloc((void**)&forwardrerpk_d,sznvr);
   cudaMalloc((void**)&forwardreltb_d,sznvr);
   cudaMalloc((void**)&forwardrelix_d,sznvr);
   cudaMalloc((void**)&forwardrelmi_d,sznvr);
   cudaMalloc((void**)&forwardrelrg_d,sznvr);
   cudaMalloc((void**)&forwardrelpk_d,sznvr);
   cudaMalloc((void**)&forwardimrtb_d,sznvr);
   cudaMalloc((void**)&forwardimrix_d,sznvr);
   cudaMalloc((void**)&forwardimrmi_d,sznvr);
   cudaMalloc((void**)&forwardimrrg_d,sznvr);
   cudaMalloc((void**)&forwardimrpk_d,sznvr);
   cudaMalloc((void**)&forwardimltb_d,sznvr);
   cudaMalloc((void**)&forwardimlix_d,sznvr);
   cudaMalloc((void**)&forwardimlmi_d,sznvr);
   cudaMalloc((void**)&forwardimlrg_d,sznvr);
   cudaMalloc((void**)&forwardimlpk_d,sznvr);
   cudaMalloc((void**)&backwardrertb_d,sznvr1);
   cudaMalloc((void**)&backwardrerix_d,sznvr1);
   cudaMalloc((void**)&backwardrermi_d,sznvr1);
   cudaMalloc((void**)&backwardrerrg_d,sznvr1);
   cudaMalloc((void**)&backwardrerpk_d,sznvr1);
   cudaMalloc((void**)&backwardreltb_d,sznvr1);
   cudaMalloc((void**)&backwardrelix_d,sznvr1);
   cudaMalloc((void**)&backwardrelmi_d,sznvr1);
   cudaMalloc((void**)&backwardrelrg_d,sznvr1);
   cudaMalloc((void**)&backwardrelpk_d,sznvr1);
   cudaMalloc((void**)&backwardimrtb_d,sznvr1);
   cudaMalloc((void**)&backwardimrix_d,sznvr1);
   cudaMalloc((void**)&backwardimrmi_d,sznvr1);
   cudaMalloc((void**)&backwardimrrg_d,sznvr1);
   cudaMalloc((void**)&backwardimrpk_d,sznvr1);
   cudaMalloc((void**)&backwardimltb_d,sznvr1);
   cudaMalloc((void**)&backwardimlix_d,sznvr1);
   cudaMalloc((void**)&backwardimlmi_d,sznvr1);
   cudaMalloc((void**)&backwardimlrg_d,sznvr1);
   cudaMalloc((void**)&backwardimlpk_d,sznvr1);
   cudaMalloc((void**)&crossrertb_d,sznvr2);
   cudaMalloc((void**)&crossrerix_d,sznvr2);
   cudaMalloc((void**)&crossrermi_d,sznvr2);
   cudaMalloc((void**)&crossrerrg_d,sznvr2);
   cudaMalloc((void**)&crossrerpk_d,sznvr2);
   cudaMalloc((void**)&crossreltb_d,sznvr2);
   cudaMalloc((void**)&crossrelix_d,sznvr2);
   cudaMalloc((void**)&crossrelmi_d,sznvr2);
   cudaMalloc((void**)&crossrelrg_d,sznvr2);
   cudaMalloc((void**)&crossrelpk_d,sznvr2);
   cudaMalloc((void**)&crossimrtb_d,sznvr2);
   cudaMalloc((void**)&crossimrix_d,sznvr2);
   cudaMalloc((void**)&crossimrmi_d,sznvr2);
   cudaMalloc((void**)&crossimrrg_d,sznvr2);
   cudaMalloc((void**)&crossimrpk_d,sznvr2);
   cudaMalloc((void**)&crossimltb_d,sznvr2);
   cudaMalloc((void**)&crossimlix_d,sznvr2);
   cudaMalloc((void**)&crossimlmi_d,sznvr2);
   cudaMalloc((void**)&crossimlrg_d,sznvr2);
   cudaMalloc((void**)&crossimlpk_d,sznvr2);

   double *inputrertb_h = new double[dim*(deg1)];
   double *inputrerix_h = new double[dim*(deg1)];
   double *inputrermi_h = new double[dim*(deg1)];
   double *inputrerrg_h = new double[dim*(deg1)];
   double *inputrerpk_h = new double[dim*(deg1)];
   double *inputreltb_h = new double[dim*(deg1)];
   double *inputrelix_h = new double[dim*(deg1)];
   double *inputrelmi_h = new double[dim*(deg1)];
   double *inputrelrg_h = new double[dim*(deg1)];
   double *inputrelpk_h = new double[dim*(deg1)];
   double *inputimrtb_h = new double[dim*(deg1)];
   double *inputimrix_h = new double[dim*(deg1)];
   double *inputimrmi_h = new double[dim*(deg1)];
   double *inputimrrg_h = new double[dim*(deg1)];
   double *inputimrpk_h = new double[dim*(deg1)];
   double *inputimltb_h = new double[dim*(deg1)];
   double *inputimlix_h = new double[dim*(deg1)];
   double *inputimlmi_h = new double[dim*(deg1)];
   double *inputimlrg_h = new double[dim*(deg1)];
   double *inputimlpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrertb_h[ix] = inputrertb[i][j];
         inputrerix_h[ix] = inputrerix[i][j];
         inputrermi_h[ix] = inputrermi[i][j];
         inputrerrg_h[ix] = inputrerrg[i][j];
         inputrerpk_h[ix] = inputrerpk[i][j];
         inputreltb_h[ix] = inputreltb[i][j];
         inputrelix_h[ix] = inputrelix[i][j];
         inputrelmi_h[ix] = inputrelmi[i][j];
         inputrelrg_h[ix] = inputrelrg[i][j];
         inputrelpk_h[ix] = inputrelpk[i][j];
         inputimrtb_h[ix] = inputimrtb[i][j];
         inputimrix_h[ix] = inputimrix[i][j];
         inputimrmi_h[ix] = inputimrmi[i][j];
         inputimrrg_h[ix] = inputimrrg[i][j];
         inputimrpk_h[ix] = inputimrpk[i][j];
         inputimltb_h[ix] = inputimltb[i][j];
         inputimlix_h[ix] = inputimlix[i][j];
         inputimlmi_h[ix] = inputimlmi[i][j];
         inputimlrg_h[ix] = inputimlrg[i][j];
         inputimlpk_h[ix++] = inputimlpk[i][j];
      }

   cudaMemcpy(cffrertb_d,cffrertb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrerix_d,cffrerix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrermi_d,cffrermi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrerrg_d,cffrerrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrerpk_d,cffrerpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffreltb_d,cffreltb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelix_d,cffrelix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelmi_d,cffrelmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelrg_d,cffrelrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelpk_d,cffrelpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrtb_d,cffimrtb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrix_d,cffimrix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrmi_d,cffimrmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrrg_d,cffimrrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrpk_d,cffimrpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimltb_d,cffimltb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlix_d,cffimlix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlmi_d,cffimlmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlrg_d,cffimlrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlpk_d,cffimlpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrertb_d,inputrertb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrerix_d,inputrerix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrermi_d,inputrermi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrerrg_d,inputrerrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrerpk_d,inputrerpk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputreltb_d,inputreltb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelix_d,inputrelix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelmi_d,inputrelmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelrg_d,inputrelrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelpk_d,inputrelpk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrtb_d,inputimrtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrix_d,inputimrix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrmi_d,inputimrmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrrg_d,inputimrrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrpk_d,inputimrpk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimltb_d,inputimltb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlix_d,inputimlix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlmi_d,inputimlmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlrg_d,inputimlrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlpk_d,inputimlpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx10_speel(BS,nvr,deg,idx,
         cffrertb_d,cffrerix_d,cffrermi_d,cffrerrg_d,cffrerpk_d,
         cffreltb_d,cffrelix_d,cffrelmi_d,cffrelrg_d,cffrelpk_d,
         cffimrtb_d,cffimrix_d,cffimrmi_d,cffimrrg_d,cffimrpk_d,
         cffimltb_d,cffimlix_d,cffimlmi_d,cffimlrg_d,cffimlpk_d,
         inputrertb_d,inputrerix_d,inputrermi_d,inputrerrg_d,inputrerpk_d,
         inputreltb_d,inputrelix_d,inputrelmi_d,inputrelrg_d,inputrelpk_d,
         inputimrtb_d,inputimrix_d,inputimrmi_d,inputimrrg_d,inputimrpk_d,
         inputimltb_d,inputimlix_d,inputimlmi_d,inputimlrg_d,inputimlpk_d,
         forwardrertb_d,forwardrerix_d,forwardrermi_d,
         forwardrerrg_d,forwardrerpk_d,
         forwardreltb_d,forwardrelix_d,forwardrelmi_d,
         forwardrelrg_d,forwardrelpk_d,
         forwardimrtb_d,forwardimrix_d,forwardimrmi_d,
         forwardimrrg_d,forwardimrpk_d,
         forwardimltb_d,forwardimlix_d,forwardimlmi_d,
         forwardimlrg_d,forwardimlpk_d,
         backwardrertb_d,backwardrerix_d,backwardrermi_d,
         backwardrerrg_d,backwardrerpk_d,
         backwardreltb_d,backwardrelix_d,backwardrelmi_d,
         backwardrelrg_d,backwardrelpk_d,
         backwardimrtb_d,backwardimrix_d,backwardimrmi_d,
         backwardimrrg_d,backwardimrpk_d,
         backwardimltb_d,backwardimlix_d,backwardimlmi_d,
         backwardimlrg_d,backwardimlpk_d,
         crossrertb_d,crossrerix_d,crossrermi_d,crossrerrg_d,crossrerpk_d,
         crossreltb_d,crossrelix_d,crossrelmi_d,crossrelrg_d,crossrelpk_d,
         crossimrtb_d,crossimrix_d,crossimrmi_d,crossimrrg_d,crossimrpk_d,
         crossimltb_d,crossimlix_d,crossimlmi_d,crossimlrg_d,crossimlpk_d);
   }
   double *forwardrertb_h = new double[(deg1)*nvr];
   double *forwardrerix_h = new double[(deg1)*nvr];
   double *forwardrermi_h = new double[(deg1)*nvr];
   double *forwardrerrg_h = new double[(deg1)*nvr];
   double *forwardrerpk_h = new double[(deg1)*nvr];
   double *forwardreltb_h = new double[(deg1)*nvr];
   double *forwardrelix_h = new double[(deg1)*nvr];
   double *forwardrelmi_h = new double[(deg1)*nvr];
   double *forwardrelrg_h = new double[(deg1)*nvr];
   double *forwardrelpk_h = new double[(deg1)*nvr];
   double *forwardimrtb_h = new double[(deg1)*nvr];
   double *forwardimrix_h = new double[(deg1)*nvr];
   double *forwardimrmi_h = new double[(deg1)*nvr];
   double *forwardimrrg_h = new double[(deg1)*nvr];
   double *forwardimrpk_h = new double[(deg1)*nvr];
   double *forwardimltb_h = new double[(deg1)*nvr];
   double *forwardimlix_h = new double[(deg1)*nvr];
   double *forwardimlmi_h = new double[(deg1)*nvr];
   double *forwardimlrg_h = new double[(deg1)*nvr];
   double *forwardimlpk_h = new double[(deg1)*nvr];
   double *backwardrertb_h = new double[(deg1)*(nvr-1)];
   double *backwardrerix_h = new double[(deg1)*(nvr-1)];
   double *backwardrermi_h = new double[(deg1)*(nvr-1)];
   double *backwardrerrg_h = new double[(deg1)*(nvr-1)];
   double *backwardrerpk_h = new double[(deg1)*(nvr-1)];
   double *backwardreltb_h = new double[(deg1)*(nvr-1)];
   double *backwardrelix_h = new double[(deg1)*(nvr-1)];
   double *backwardrelmi_h = new double[(deg1)*(nvr-1)];
   double *backwardrelrg_h = new double[(deg1)*(nvr-1)];
   double *backwardrelpk_h = new double[(deg1)*(nvr-1)];
   double *backwardimrtb_h = new double[(deg1)*(nvr-1)];
   double *backwardimrix_h = new double[(deg1)*(nvr-1)];
   double *backwardimrmi_h = new double[(deg1)*(nvr-1)];
   double *backwardimrrg_h = new double[(deg1)*(nvr-1)];
   double *backwardimrpk_h = new double[(deg1)*(nvr-1)];
   double *backwardimltb_h = new double[(deg1)*(nvr-1)];
   double *backwardimlix_h = new double[(deg1)*(nvr-1)];
   double *backwardimlmi_h = new double[(deg1)*(nvr-1)];
   double *backwardimlrg_h = new double[(deg1)*(nvr-1)];
   double *backwardimlpk_h = new double[(deg1)*(nvr-1)];
   double *crossrertb_h = new double[(deg1)*(nvr-2)];
   double *crossrerix_h = new double[(deg1)*(nvr-2)];
   double *crossrermi_h = new double[(deg1)*(nvr-2)];
   double *crossrerrg_h = new double[(deg1)*(nvr-2)];
   double *crossrerpk_h = new double[(deg1)*(nvr-2)];
   double *crossreltb_h = new double[(deg1)*(nvr-2)];
   double *crossrelix_h = new double[(deg1)*(nvr-2)];
   double *crossrelmi_h = new double[(deg1)*(nvr-2)];
   double *crossrelrg_h = new double[(deg1)*(nvr-2)];
   double *crossrelpk_h = new double[(deg1)*(nvr-2)];
   double *crossimrtb_h = new double[(deg1)*(nvr-2)];
   double *crossimrix_h = new double[(deg1)*(nvr-2)];
   double *crossimrmi_h = new double[(deg1)*(nvr-2)];
   double *crossimrrg_h = new double[(deg1)*(nvr-2)];
   double *crossimrpk_h = new double[(deg1)*(nvr-2)];
   double *crossimltb_h = new double[(deg1)*(nvr-2)];
   double *crossimlix_h = new double[(deg1)*(nvr-2)];
   double *crossimlmi_h = new double[(deg1)*(nvr-2)];
   double *crossimlrg_h = new double[(deg1)*(nvr-2)];
   double *crossimlpk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrertb_h,forwardrertb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrerix_h,forwardrerix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrermi_h,forwardrermi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrerrg_h,forwardrerrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrerpk_h,forwardrerpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardreltb_h,forwardreltb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelix_h,forwardrelix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelmi_h,forwardrelmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelrg_h,forwardrelrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelpk_h,forwardrelpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrtb_h,forwardimrtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrix_h,forwardimrix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrmi_h,forwardimrmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrrg_h,forwardimrrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrpk_h,forwardimrpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimltb_h,forwardimltb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlix_h,forwardimlix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlmi_h,forwardimlmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlrg_h,forwardimlrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlpk_h,forwardimlpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrertb_h,backwardrertb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrerix_h,backwardrerix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrermi_h,backwardrermi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrerrg_h,backwardrerrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrerpk_h,backwardrerpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardreltb_h,backwardreltb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelix_h,backwardrelix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelmi_h,backwardrelmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelrg_h,backwardrelrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelpk_h,backwardrelpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrtb_h,backwardimrtb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrix_h,backwardimrix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrmi_h,backwardimrmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrrg_h,backwardimrrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrpk_h,backwardimrpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimltb_h,backwardimltb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlix_h,backwardimlix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlmi_h,backwardimlmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlrg_h,backwardimlrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlpk_h,backwardimlpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrertb_h,crossrertb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrerix_h,crossrerix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrermi_h,crossrermi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrerrg_h,crossrerrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrerpk_h,crossrerpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossreltb_h,crossreltb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelix_h,crossrelix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelmi_h,crossrelmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelrg_h,crossrelrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelpk_h,crossrelpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrtb_h,crossimrtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrix_h,crossimrix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrmi_h,crossimrmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrrg_h,crossimrrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrpk_h,crossimrpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimltb_h,crossimltb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlix_h,crossimlix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlmi_h,crossimlmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlrg_h,crossimlrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlpk_h,crossimlpk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputrertb[dim][i] = forwardrertb_h[offset+i];
      outputrerix[dim][i] = forwardrerix_h[offset+i];
      outputrermi[dim][i] = forwardrermi_h[offset+i];
      outputrerrg[dim][i] = forwardrerrg_h[offset+i];
      outputrerpk[dim][i] = forwardrerpk_h[offset+i];
      outputreltb[dim][i] = forwardreltb_h[offset+i];
      outputrelix[dim][i] = forwardrelix_h[offset+i];
      outputrelmi[dim][i] = forwardrelmi_h[offset+i];
      outputrelrg[dim][i] = forwardrelrg_h[offset+i];
      outputrelpk[dim][i] = forwardrelpk_h[offset+i];
      outputimrtb[dim][i] = forwardimrtb_h[offset+i];
      outputimrix[dim][i] = forwardimrix_h[offset+i];
      outputimrmi[dim][i] = forwardimrmi_h[offset+i];
      outputimrrg[dim][i] = forwardimrrg_h[offset+i];
      outputimrpk[dim][i] = forwardimrpk_h[offset+i];
      outputimltb[dim][i] = forwardimltb_h[offset+i];
      outputimlix[dim][i] = forwardimlix_h[offset+i];
      outputimlmi[dim][i] = forwardimlmi_h[offset+i];
      outputimlrg[dim][i] = forwardimlrg_h[offset+i];
      outputimlpk[dim][i] = forwardimlpk_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputrertb[ix][i] = forwardrertb_h[offset+i];
      outputrerix[ix][i] = forwardrerix_h[offset+i];
      outputrermi[ix][i] = forwardrermi_h[offset+i];
      outputrerrg[ix][i] = forwardrerrg_h[offset+i];
      outputrerpk[ix][i] = forwardrerpk_h[offset+i];
      outputreltb[ix][i] = forwardreltb_h[offset+i];
      outputrelix[ix][i] = forwardrelix_h[offset+i];
      outputrelmi[ix][i] = forwardrelmi_h[offset+i];
      outputrelrg[ix][i] = forwardrelrg_h[offset+i];
      outputrelpk[ix][i] = forwardrelpk_h[offset+i];
      outputimrtb[ix][i] = forwardimrtb_h[offset+i];
      outputimrix[ix][i] = forwardimrix_h[offset+i];
      outputimrmi[ix][i] = forwardimrmi_h[offset+i];
      outputimrrg[ix][i] = forwardimrrg_h[offset+i];
      outputimrpk[ix][i] = forwardimrpk_h[offset+i];
      outputimltb[ix][i] = forwardimltb_h[offset+i];
      outputimlix[ix][i] = forwardimlix_h[offset+i];
      outputimlmi[ix][i] = forwardimlmi_h[offset+i];
      outputimlrg[ix][i] = forwardimlrg_h[offset+i];
      outputimlpk[ix][i] = forwardimlpk_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputrertb[ix][i] = backwardrertb_h[offset+i];
      outputrerix[ix][i] = backwardrerix_h[offset+i];
      outputrermi[ix][i] = backwardrermi_h[offset+i];
      outputrerrg[ix][i] = backwardrerrg_h[offset+i];
      outputrerpk[ix][i] = backwardrerpk_h[offset+i];
      outputreltb[ix][i] = backwardreltb_h[offset+i];
      outputrelix[ix][i] = backwardrelix_h[offset+i];
      outputrelmi[ix][i] = backwardrelmi_h[offset+i];
      outputrelrg[ix][i] = backwardrelrg_h[offset+i];
      outputrelpk[ix][i] = backwardrelpk_h[offset+i];
      outputimrtb[ix][i] = backwardimrtb_h[offset+i];
      outputimrix[ix][i] = backwardimrix_h[offset+i];
      outputimrmi[ix][i] = backwardimrmi_h[offset+i];
      outputimrrg[ix][i] = backwardimrrg_h[offset+i];
      outputimrpk[ix][i] = backwardimrpk_h[offset+i];
      outputimltb[ix][i] = backwardimltb_h[offset+i];
      outputimlix[ix][i] = backwardimlix_h[offset+i];
      outputimlmi[ix][i] = backwardimlmi_h[offset+i];
      outputimlrg[ix][i] = backwardimlrg_h[offset+i];
      outputimlpk[ix][i] = backwardimlpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrertb[ix][i] = crossrertb_h[offset+i];
         outputrerix[ix][i] = crossrerix_h[offset+i];
         outputrermi[ix][i] = crossrermi_h[offset+i];
         outputrerrg[ix][i] = crossrerrg_h[offset+i];
         outputrerpk[ix][i] = crossrerpk_h[offset+i];
         outputreltb[ix][i] = crossreltb_h[offset+i];
         outputrelix[ix][i] = crossrelix_h[offset+i];
         outputrelmi[ix][i] = crossrelmi_h[offset+i];
         outputrelrg[ix][i] = crossrelrg_h[offset+i];
         outputrelpk[ix][i] = crossrelpk_h[offset+i];
         outputimrtb[ix][i] = crossimrtb_h[offset+i];
         outputimrix[ix][i] = crossimrix_h[offset+i];
         outputimrmi[ix][i] = crossimrmi_h[offset+i];
         outputimrrg[ix][i] = crossimrrg_h[offset+i];
         outputimrpk[ix][i] = crossimrpk_h[offset+i];
         outputimltb[ix][i] = crossimltb_h[offset+i];
         outputimlix[ix][i] = crossimlix_h[offset+i];
         outputimlmi[ix][i] = crossimlmi_h[offset+i];
         outputimlrg[ix][i] = crossimlrg_h[offset+i];
         outputimlpk[ix][i] = crossimlpk_h[offset+i];
      }
   }
}
