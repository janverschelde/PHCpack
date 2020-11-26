// The file dbl8_monomials_kernels.cu defines kernels specified
// in dbl8_monomials_kernels.h.

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

#include "dbl8_convolutions_kernels.cu"
#include "dbl8_monomials_kernels.h"

void GPU_dbl8_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo,
   double *inputhihihi, double *inputlohihi,
   double *inputhilohi, double *inputlolohi,
   double *inputhihilo, double *inputlohilo,
   double *inputhilolo, double *inputlololo,
   double *forwardhihihi, double *forwardlohihi,
   double *forwardhilohi, double *forwardlolohi,
   double *forwardhihilo, double *forwardlohilo,
   double *forwardhilolo, double *forwardlololo,
   double *backwardhihihi, double *backwardlohihi,
   double *backwardhilohi, double *backwardlolohi,
   double *backwardhihilo, double *backwardlohilo,
   double *backwardhilolo, double *backwardlololo,
   double *crosshihihi, double *crosslohihi,
   double *crosshilohi, double *crosslolohi,
   double *crosshihilo, double *crosslohilo,
   double *crosshilolo, double *crosslololo )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl8_padded_convolute<<<1,BS>>>
      (cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       &inputhihihi[ix1],&inputlohihi[ix1],
       &inputhilohi[ix1],&inputlolohi[ix1],
       &inputhihilo[ix1],&inputlohilo[ix1],
       &inputhilolo[ix1],&inputlololo[ix1],
       forwardhihihi,forwardlohihi,forwardhilohi,forwardlolohi,
       forwardhihilo,forwardlohilo,forwardhilolo,forwardlololo,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl8_padded_convolute<<<1,BS>>>
         (&forwardhihihi[ix1],&forwardlohihi[ix1],
          &forwardhilohi[ix1],&forwardlolohi[ix1],
          &forwardhihilo[ix1],&forwardlohilo[ix1],
          &forwardhilolo[ix1],&forwardlololo[ix1],
          &inputhihihi[ix2],&inputlohihi[ix2],
          &inputhilohi[ix2],&inputlolohi[ix2],
          &inputhihilo[ix2],&inputlohilo[ix2],
          &inputhilolo[ix2],&inputlololo[ix2],
          &forwardhihihi[ix3],&forwardlohihi[ix3],
          &forwardhilohi[ix3],&forwardlolohi[ix3],
          &forwardhihilo[ix3],&forwardlohilo[ix3],
          &forwardhilolo[ix3],&forwardlololo[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl8_padded_convolute<<<1,BS>>>
         (&inputhihihi[ix1],&inputlohihi[ix1],
          &inputhilohi[ix1],&inputlolohi[ix1],
          &inputhihilo[ix1],&inputlohilo[ix1],
          &inputhilolo[ix1],&inputlololo[ix1],
          &inputhihihi[ix2],&inputlohihi[ix2],
          &inputhilohi[ix2],&inputlolohi[ix2],
          &inputhihilo[ix2],&inputlohilo[ix2],
          &inputhilolo[ix2],&inputlololo[ix2],
          backwardhihihi,backwardlohihi,backwardhilohi,backwardlolohi,
          backwardhihilo,backwardlohilo,backwardhilolo,backwardlololo,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl8_padded_convolute<<<1,BS>>>
            (&backwardhihihi[ix1],&backwardlohihi[ix1],
             &backwardhilohi[ix1],&backwardlolohi[ix1],
             &backwardhihilo[ix1],&backwardlohilo[ix1],
             &backwardhilolo[ix1],&backwardlololo[ix1],
             &inputhihihi[ix2],&inputlohihi[ix2],
             &inputhilohi[ix2],&inputlolohi[ix2],
             &inputhihilo[ix2],&inputlohilo[ix2],
             &inputhilolo[ix2],&inputlololo[ix2],
             &backwardhihihi[ix3],&backwardlohihi[ix3],
             &backwardhilohi[ix3],&backwardlolohi[ix3],
             &backwardhihilo[ix3],&backwardlohilo[ix3],
             &backwardhilolo[ix3],&backwardlololo[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl8_padded_convolute<<<1,BS>>>
         (&backwardhihihi[ix3],&backwardlohihi[ix3],
          &backwardhilohi[ix3],&backwardlolohi[ix3],
          &backwardhihilo[ix3],&backwardlohilo[ix3],
          &backwardhilolo[ix3],&backwardlololo[ix3],
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          &backwardhihihi[ix2],&backwardlohihi[ix2],
          &backwardhilohi[ix2],&backwardlolohi[ix2],
          &backwardhihilo[ix2],&backwardlohilo[ix2],
          &backwardhilolo[ix2],&backwardlololo[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl8_padded_convolute<<<1,BS>>>
            (forwardhihihi,forwardlohihi,forwardhilohi,forwardlolohi,
             forwardhihilo,forwardlohilo,forwardhilolo,forwardlololo,
             &inputhihihi[ix2],&inputlohihi[ix2],
             &inputhilohi[ix2],&inputlolohi[ix2],
             &inputhihilo[ix2],&inputlohilo[ix2],
             &inputhilolo[ix2],&inputlololo[ix2],
             crosshihihi,crosslohihi,crosshilohi,crosslolohi,
             crosshihilo,crosslohilo,crosshilolo,crosslololo,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl8_padded_convolute<<<1,BS>>>
               (&forwardhihihi[ix1],&forwardlohihi[ix1],
                &forwardhilohi[ix1],&forwardlolohi[ix1],
                &forwardhihilo[ix1],&forwardlohilo[ix1],
                &forwardhilolo[ix1],&forwardlololo[ix1],
                &backwardhihihi[ix2],&backwardlohihi[ix2],
                &backwardhilohi[ix2],&backwardlolohi[ix2],
                &backwardhihilo[ix2],&backwardlohilo[ix2],
                &backwardhilolo[ix2],&backwardlololo[ix2],
                &crosshihihi[ix1],&crosslohihi[ix1],
                &crosshilohi[ix1],&crosslolohi[ix1],
                &crosshihilo[ix1],&crosslohilo[ix1],
                &crosshilolo[ix1],&crosslololo[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl8_padded_convolute<<<1,BS>>>
            (&forwardhihihi[ix1],&forwardlohihi[ix1],
             &forwardhilohi[ix1],&forwardlolohi[ix1],
             &forwardhihilo[ix1],&forwardlohilo[ix1],
             &forwardhilolo[ix1],&forwardlololo[ix1],
             &inputhihihi[ix2],&inputlohihi[ix2],
             &inputhilohi[ix2],&inputlolohi[ix2],
             &inputhihilo[ix2],&inputlohilo[ix2],
             &inputhilolo[ix2],&inputlololo[ix2],
             &crosshihihi[ix1],&crosslohihi[ix1],
             &crosshilohi[ix1],&crosslolohi[ix1],
             &crosshihilo[ix1],&crosslohilo[ix1],
             &crosshilolo[ix1],&crosslololo[ix1],deg1);
      }
   }
}

void GPU_cmplx8_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double *inputrehihihi, double *inputrelohihi,
   double *inputrehilohi, double *inputrelolohi,
   double *inputrehihilo, double *inputrelohilo,
   double *inputrehilolo, double *inputrelololo,
   double *inputimhihihi, double *inputimlohihi,
   double *inputimhilohi, double *inputimlolohi,
   double *inputimhihilo, double *inputimlohilo,
   double *inputimhilolo, double *inputimlololo,
   double *forwardrehihihi, double *forwardrelohihi,
   double *forwardrehilohi, double *forwardrelolohi,
   double *forwardrehihilo, double *forwardrelohilo,
   double *forwardrehilolo, double *forwardrelololo,
   double *forwardimhihihi, double *forwardimlohihi,
   double *forwardimhilohi, double *forwardimlolohi,
   double *forwardimhihilo, double *forwardimlohilo,
   double *forwardimhilolo, double *forwardimlololo,
   double *backwardrehihihi, double *backwardrelohihi,
   double *backwardrehilohi, double *backwardrelolohi,
   double *backwardrehihilo, double *backwardrelohilo,
   double *backwardrehilolo, double *backwardrelololo,
   double *backwardimhihihi, double *backwardimlohihi,
   double *backwardimhilohi, double *backwardimlolohi,
   double *backwardimhihilo, double *backwardimlohilo,
   double *backwardimhilolo, double *backwardimlololo,
   double *crossrehihihi, double *crossrelohihi,
   double *crossrehilohi, double *crossrelolohi,
   double *crossrehihilo, double *crossrelohilo,
   double *crossrehilolo, double *crossrelololo,
   double *crossimhihihi, double *crossimlohihi,
   double *crossimhilohi, double *crossimlolohi,
   double *crossimhihilo, double *crossimlohilo,
   double *crossimhilolo, double *crossimlololo )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   cmplx8_padded_convolute<<<1,BS>>>
      (cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       &inputrehihihi[ix1],&inputrelohihi[ix1],
       &inputrehilohi[ix1],&inputrelolohi[ix1],
       &inputrehihilo[ix1],&inputrelohilo[ix1],
       &inputrehilolo[ix1],&inputrelololo[ix1],
       &inputimhihihi[ix1],&inputimlohihi[ix1],
       &inputimhilohi[ix1],&inputimlolohi[ix1],
       &inputimhihilo[ix1],&inputimlohilo[ix1],
       &inputimhilolo[ix1],&inputimlololo[ix1],
       forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
       forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
       forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
       forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
       deg1); 

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-i]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      cmplx8_padded_convolute<<<1,BS>>>
         (&forwardrehihihi[ix1],&forwardrelohihi[ix1],
          &forwardrehilohi[ix1],&forwardrelolohi[ix1],
          &forwardrehihilo[ix1],&forwardrelohilo[ix1],
          &forwardrehilolo[ix1],&forwardrelololo[ix1],
          &forwardimhihihi[ix1],&forwardimlohihi[ix1],
          &forwardimhilohi[ix1],&forwardimlolohi[ix1],
          &forwardimhihilo[ix1],&forwardimlohilo[ix1],
          &forwardimhilolo[ix1],&forwardimlololo[ix1],
          &inputrehihihi[ix2],&inputrelohihi[ix2],
          &inputrehilohi[ix2],&inputrelolohi[ix2],
          &inputrehihilo[ix2],&inputrelohilo[ix2],
          &inputrehilolo[ix2],&inputrelololo[ix2],
          &inputimhihihi[ix2],&inputimlohihi[ix2],
          &inputimhilohi[ix2],&inputimlolohi[ix2],
          &inputimhihilo[ix2],&inputimlohilo[ix2],
          &inputimhilolo[ix2],&inputimlololo[ix2],
          &forwardrehihihi[ix3],&forwardrelohihi[ix3],
          &forwardrehilohi[ix3],&forwardrelolohi[ix3],
          &forwardrehihilo[ix3],&forwardrelohilo[ix3],
          &forwardrehilolo[ix3],&forwardrelololo[ix3],
          &forwardimhihihi[ix3],&forwardimlohihi[ix3],
          &forwardimhilohi[ix3],&forwardimlolohi[ix3],
          &forwardimhihilo[ix3],&forwardimlohilo[ix3],
          &forwardimhilolo[ix3],&forwardimlololo[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      cmplx8_padded_convolute<<<1,BS>>>
         (&inputrehihihi[ix1],&inputrelohihi[ix1],
          &inputrehilohi[ix1],&inputrelolohi[ix1],
          &inputrehihilo[ix1],&inputrelohilo[ix1],
          &inputrehilolo[ix1],&inputrelololo[ix1],
          &inputimhihihi[ix1],&inputimlohihi[ix1],
          &inputimhilohi[ix1],&inputimlolohi[ix1],
          &inputimhihilo[ix1],&inputimlohilo[ix1],
          &inputimhilolo[ix1],&inputimlololo[ix1],
          &inputrehihihi[ix2],&inputrelohihi[ix2],
          &inputrehilohi[ix2],&inputrelolohi[ix2],
          &inputrehihilo[ix2],&inputrelohilo[ix2],
          &inputrehilolo[ix2],&inputrelololo[ix2],
          &inputimhihihi[ix2],&inputimlohihi[ix2],
          &inputimhilohi[ix2],&inputimlolohi[ix2],
          &inputimhihilo[ix2],&inputimlohilo[ix2],
          &inputimhilolo[ix2],&inputimlololo[ix2],
          backwardrehihihi,backwardrelohihi,backwardrehilohi,backwardrelolohi,
          backwardrehihilo,backwardrelohilo,backwardrehilolo,backwardrelololo,
          backwardimhihihi,backwardimlohihi,backwardimhilohi,backwardimlolohi,
          backwardimhihilo,backwardimlohilo,backwardimhilolo,backwardimlololo,
          deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         cmplx8_padded_convolute<<<1,BS>>>
            (&backwardrehihihi[ix1],&backwardrelohihi[ix1],
             &backwardrehilohi[ix1],&backwardrelolohi[ix1],
             &backwardrehihilo[ix1],&backwardrelohilo[ix1],
             &backwardrehilolo[ix1],&backwardrelololo[ix1],
             &backwardimhihihi[ix1],&backwardimlohihi[ix1],
             &backwardimhilohi[ix1],&backwardimlolohi[ix1],
             &backwardimhihilo[ix1],&backwardimlohilo[ix1],
             &backwardimhilolo[ix1],&backwardimlololo[ix1],
             &inputrehihihi[ix2],&inputrelohihi[ix2],
             &inputrehilohi[ix2],&inputrelolohi[ix2],
             &inputrehihilo[ix2],&inputrelohilo[ix2],
             &inputrehilolo[ix2],&inputrelololo[ix2],
             &inputimhihihi[ix2],&inputimlohihi[ix2],
             &inputimhilohi[ix2],&inputimlolohi[ix2],
             &inputimhihilo[ix2],&inputimlohilo[ix2],
             &inputimhilolo[ix2],&inputimlololo[ix2],
             &backwardrehihihi[ix3],&backwardrelohihi[ix3],
             &backwardrehilohi[ix3],&backwardrelolohi[ix3],
             &backwardrehihilo[ix3],&backwardrelohilo[ix3],
             &backwardrehilolo[ix3],&backwardrelololo[ix3],
             &backwardimhihihi[ix3],&backwardimlohihi[ix3],
             &backwardimhilohi[ix3],&backwardimlolohi[ix3],
             &backwardimhihilo[ix3],&backwardimlohilo[ix3],
             &backwardimhilolo[ix3],&backwardimlololo[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      cmplx8_padded_convolute<<<1,BS>>>
         (&backwardrehihihi[ix3],&backwardrelohihi[ix3],
          &backwardrehilohi[ix3],&backwardrelolohi[ix3],
          &backwardrehihilo[ix3],&backwardrelohilo[ix3],
          &backwardrehilolo[ix3],&backwardrelololo[ix3],
          &backwardimhihihi[ix3],&backwardimlohihi[ix3],
          &backwardimhilohi[ix3],&backwardimlolohi[ix3],
          &backwardimhihilo[ix3],&backwardimlohilo[ix3],
          &backwardimhilolo[ix3],&backwardimlololo[ix3],
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
          &backwardrehihihi[ix2],&backwardrelohihi[ix2],
          &backwardrehilohi[ix2],&backwardrelolohi[ix2],
          &backwardrehihilo[ix2],&backwardrelohilo[ix2],
          &backwardrehilolo[ix2],&backwardrelololo[ix2],
          &backwardimhihihi[ix2],&backwardimlohihi[ix2],
          &backwardimhilohi[ix2],&backwardimlolohi[ix2],
          &backwardimhihilo[ix2],&backwardimlohilo[ix2],
          &backwardimhilolo[ix2],&backwardimlololo[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         cmplx8_padded_convolute<<<1,BS>>>
            (forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
             forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
             forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
             forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
             &inputrehihihi[ix2],&inputrelohihi[ix2],
             &inputrehilohi[ix2],&inputrelolohi[ix2],
             &inputrehihilo[ix2],&inputrelohilo[ix2],
             &inputrehilolo[ix2],&inputrelololo[ix2],
             &inputimhihihi[ix2],&inputimlohihi[ix2],
             &inputimhilohi[ix2],&inputimlolohi[ix2],
             &inputimhihilo[ix2],&inputimlohilo[ix2],
             &inputimhilolo[ix2],&inputimlololo[ix2],
             crossrehihihi,crossrelohihi,crossrehilohi,crossrelolohi,
             crossrehihilo,crossrelohilo,crossrehilolo,crossrelololo,
             crossimhihihi,crossimlohihi,crossimhilohi,crossimlolohi,
             crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            cmplx8_padded_convolute<<<1,BS>>>
               (&forwardrehihihi[ix1],&forwardrelohihi[ix1],
                &forwardrehilohi[ix1],&forwardrelolohi[ix1],
                &forwardrehihilo[ix1],&forwardrelohilo[ix1],
                &forwardrehilolo[ix1],&forwardrelololo[ix1],
                &forwardimhihihi[ix1],&forwardimlohihi[ix1],
                &forwardimhilohi[ix1],&forwardimlolohi[ix1],
                &forwardimhihilo[ix1],&forwardimlohilo[ix1],
                &forwardimhilolo[ix1],&forwardimlololo[ix1],
                &backwardrehihihi[ix2],&backwardrelohihi[ix2],
                &backwardrehilohi[ix2],&backwardrelolohi[ix2],
                &backwardrehihilo[ix2],&backwardrelohilo[ix2],
                &backwardrehilolo[ix2],&backwardrelololo[ix2],
                &backwardimhihihi[ix2],&backwardimlohihi[ix2],
                &backwardimhilohi[ix2],&backwardimlolohi[ix2],
                &backwardimhihilo[ix2],&backwardimlohilo[ix2],
                &backwardimhilolo[ix2],&backwardimlololo[ix2],
                &crossrehihihi[ix1],&crossrelohihi[ix1],
                &crossrehilohi[ix1],&crossrelolohi[ix1],
                &crossrehihilo[ix1],&crossrelohilo[ix1],
                &crossrehilolo[ix1],&crossrelololo[ix1],
                &crossimhihihi[ix1],&crossimlohihi[ix1],
                &crossimhilohi[ix1],&crossimlolohi[ix1],
                &crossimhihilo[ix1],&crossimlohilo[ix1],
                &crossimhilolo[ix1],&crossimlololo[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         cmplx8_padded_convolute<<<1,BS>>>
            (&forwardrehihihi[ix1],&forwardrelohihi[ix1],
             &forwardrehilohi[ix1],&forwardrelolohi[ix1],
             &forwardrehihilo[ix1],&forwardrelohilo[ix1],
             &forwardrehilolo[ix1],&forwardrelololo[ix1],
             &forwardimhihihi[ix1],&forwardimlohihi[ix1],
             &forwardimhilohi[ix1],&forwardimlolohi[ix1],
             &forwardimhihilo[ix1],&forwardimlohilo[ix1],
             &forwardimhilolo[ix1],&forwardimlololo[ix1],
             &inputrehihihi[ix2],&inputrelohihi[ix2],
             &inputrehilohi[ix2],&inputrelolohi[ix2],
             &inputrehihilo[ix2],&inputrelohilo[ix2],
             &inputrehilolo[ix2],&inputrelololo[ix2],
             &inputimhihihi[ix2],&inputimlohihi[ix2],
             &inputimhilohi[ix2],&inputimlolohi[ix2],
             &inputimhihilo[ix2],&inputimlohilo[ix2],
             &inputimhilolo[ix2],&inputimlololo[ix2],
             &crossrehihihi[ix1],&crossrelohihi[ix1],
             &crossrehilohi[ix1],&crossrelolohi[ix1],
             &crossrehihilo[ix1],&crossrelohilo[ix1],
             &crossrehilolo[ix1],&crossrelololo[ix1],
             &crossimhihihi[ix1],&crossimlohihi[ix1],
             &crossimhilohi[ix1],&crossimlolohi[ix1],
             &crossimhihilo[ix1],&crossimlohilo[ix1],
             &crossimhilolo[ix1],&crossimlololo[ix1],deg1);
      }
   }
}

void GPU_dbl8_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
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
   const int deg1 = deg+1;   // length of all vectors
   double *inputhihihi_d;    // inputhihihi_d is inputhihihi on the device
   double *inputlohihi_d;    // inputlohihi_d is inputlohihi on the device
   double *inputhilohi_d;    // inputhilohi_d is inputhilohi on the device
   double *inputlolohi_d;    // inputlolohi_d is inputlolohi on the device
   double *inputhihilo_d;    // inputhihilo_d is inputhihilo on the device
   double *inputlohilo_d;    // inputlohilo_d is inputlohilo on the device
   double *inputhilolo_d;    // inputhilolo_d is inputhilolo on the device
   double *inputlololo_d;    // inputlololo_d is inputlololo on the device
   double *forwardhihihi_d;  // highest forward products on the device
   double *forwardlohihi_d;  // second highest forward products on the device
   double *forwardhilohi_d;  // third highest forward products on the device
   double *forwardlolohi_d;  // fourth highest forward products on the device
   double *forwardhihilo_d;  // fourth lowest forward products on the device
   double *forwardlohilo_d;  // third lowest forward products on the device
   double *forwardhilolo_d;  // second lowest forward products on the device
   double *forwardlololo_d;  // lowest forward products on the device
   double *backwardhihihi_d; // highest backward products on the device
   double *backwardlohihi_d; // second highest backward products on the device
   double *backwardhilohi_d; // third highest backward products on the device
   double *backwardlolohi_d; // fourth highest backward products on the device
   double *backwardhihilo_d; // fourth lowest backward products on the device
   double *backwardlohilo_d; // third lowest backward products on the device
   double *backwardhilolo_d; // second lowest backward products on the device
   double *backwardlololo_d; // lowest backward products on the device
   double *crosshihihi_d;    // highest cross products on the device
   double *crosslohihi_d;    // second highest cross products on the device
   double *crosshilohi_d;    // third highest cross products on the device
   double *crosslolohi_d;    // fourth highest cross products on the device
   double *crosshihilo_d;    // fourth lowest cross products on the device
   double *crosslohilo_d;    // third lowest cross products on the device
   double *crosshilolo_d;    // second lowest cross products on the device
   double *crosslololo_d;    // lowest cross products on the device
   double *cffhihihi_d;      // cffhihihi_d is cffhihihi on device
   double *cfflohihi_d;      // cfflohihi_d is cfflohihi on device
   double *cffhilohi_d;      // cffhilohi_d is cffhilohi on device
   double *cfflolohi_d;      // cfflolohi_d is cfflolohi on device
   double *cffhihilo_d;      // cffhihilo_d is cffhihilo on device
   double *cfflohilo_d;      // cfflohilo_d is cfflohilo on device
   double *cffhilolo_d;      // cffhilolo_d is cffhilolo on device
   double *cfflololo_d;      // cfflololo_d is cfflololo on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);

   cudaMalloc((void**)&cffhihihi_d,szcff);
   cudaMalloc((void**)&cfflohihi_d,szcff);
   cudaMalloc((void**)&cffhilohi_d,szcff);
   cudaMalloc((void**)&cfflolohi_d,szcff);
   cudaMalloc((void**)&cffhihilo_d,szcff);
   cudaMalloc((void**)&cfflohilo_d,szcff);
   cudaMalloc((void**)&cffhilolo_d,szcff);
   cudaMalloc((void**)&cfflololo_d,szcff);
   cudaMalloc((void**)&inputhihihi_d,szdim);
   cudaMalloc((void**)&inputlohihi_d,szdim);
   cudaMalloc((void**)&inputhilohi_d,szdim);
   cudaMalloc((void**)&inputlolohi_d,szdim);
   cudaMalloc((void**)&inputhihilo_d,szdim);
   cudaMalloc((void**)&inputlohilo_d,szdim);
   cudaMalloc((void**)&inputhilolo_d,szdim);
   cudaMalloc((void**)&inputlololo_d,szdim);
   cudaMalloc((void**)&forwardhihihi_d,sznvr);
   cudaMalloc((void**)&forwardlohihi_d,sznvr);
   cudaMalloc((void**)&forwardhilohi_d,sznvr);
   cudaMalloc((void**)&forwardlolohi_d,sznvr);
   cudaMalloc((void**)&forwardhihilo_d,sznvr);
   cudaMalloc((void**)&forwardlohilo_d,sznvr);
   cudaMalloc((void**)&forwardhilolo_d,sznvr);
   cudaMalloc((void**)&forwardlololo_d,sznvr);
   cudaMalloc((void**)&backwardhihihi_d,sznvr1);
   cudaMalloc((void**)&backwardlohihi_d,sznvr1);
   cudaMalloc((void**)&backwardhilohi_d,sznvr1);
   cudaMalloc((void**)&backwardlolohi_d,sznvr1);
   cudaMalloc((void**)&backwardhihilo_d,sznvr1);
   cudaMalloc((void**)&backwardlohilo_d,sznvr1);
   cudaMalloc((void**)&backwardhilolo_d,sznvr1);
   cudaMalloc((void**)&backwardlololo_d,sznvr1);
   cudaMalloc((void**)&crosshihihi_d,sznvr1);
   cudaMalloc((void**)&crosslohihi_d,sznvr1);
   cudaMalloc((void**)&crosshilohi_d,sznvr1);
   cudaMalloc((void**)&crosslolohi_d,sznvr1);
   cudaMalloc((void**)&crosshihilo_d,sznvr1);
   cudaMalloc((void**)&crosslohilo_d,sznvr1);
   cudaMalloc((void**)&crosshilolo_d,sznvr1);
   cudaMalloc((void**)&crosslololo_d,sznvr1);

   double *inputhihihi_h = new double[dim*(deg1)];
   double *inputlohihi_h = new double[dim*(deg1)];
   double *inputhilohi_h = new double[dim*(deg1)];
   double *inputlolohi_h = new double[dim*(deg1)];
   double *inputhihilo_h = new double[dim*(deg1)];
   double *inputlohilo_h = new double[dim*(deg1)];
   double *inputhilolo_h = new double[dim*(deg1)];
   double *inputlololo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputhihihi_h[ix] = inputhihihi[i][j];
         inputlohihi_h[ix] = inputlohihi[i][j];
         inputhilohi_h[ix] = inputhilohi[i][j];
         inputlolohi_h[ix] = inputlolohi[i][j];
         inputhihilo_h[ix] = inputhihilo[i][j];
         inputlohilo_h[ix] = inputlohilo[i][j];
         inputhilolo_h[ix] = inputhilolo[i][j];
         inputlololo_h[ix++] = inputlololo[i][j];
      }

   cudaMemcpy(cffhihihi_d,cffhihihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflohihi_d,cfflohihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhilohi_d,cffhilohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflolohi_d,cfflolohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhihilo_d,cffhihilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflohilo_d,cfflohilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhilolo_d,cffhilolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflololo_d,cfflololo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhihihi_d,inputhihihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlohihi_d,inputlohihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhilohi_d,inputhilohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlolohi_d,inputlolohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhihilo_d,inputhihilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlohilo_d,inputlohilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhilolo_d,inputhilolo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlololo_d,inputlololo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl8_speel
         (BS,nvr,deg,idx,
          cffhihihi_d,cfflohihi_d,cffhilohi_d,cfflolohi_d,
          cffhihilo_d,cfflohilo_d,cffhilolo_d,cfflololo_d,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          forwardhihihi_d,forwardlohihi_d,forwardhilohi_d,forwardlolohi_d,
          forwardhihilo_d,forwardlohilo_d,forwardhilolo_d,forwardlololo_d,
          backwardhihihi_d,backwardlohihi_d,backwardhilohi_d,backwardlolohi_d,
          backwardhihilo_d,backwardlohilo_d,backwardhilolo_d,backwardlololo_d,
          crosshihihi_d,crosslohihi_d,crosshilohi_d,crosslolohi_d,
          crosshihilo_d,crosslohilo_d,crosshilolo_d,crosslololo_d);
   }
   double *forwardhihihi_h = new double[(deg1)*nvr];
   double *forwardlohihi_h = new double[(deg1)*nvr];
   double *forwardhilohi_h = new double[(deg1)*nvr];
   double *forwardlolohi_h = new double[(deg1)*nvr];
   double *forwardhihilo_h = new double[(deg1)*nvr];
   double *forwardlohilo_h = new double[(deg1)*nvr];
   double *forwardhilolo_h = new double[(deg1)*nvr];
   double *forwardlololo_h = new double[(deg1)*nvr];
   double *backwardhihihi_h = new double[(deg1)*(nvr-1)];
   double *backwardlohihi_h = new double[(deg1)*(nvr-1)];
   double *backwardhilohi_h = new double[(deg1)*(nvr-1)];
   double *backwardlolohi_h = new double[(deg1)*(nvr-1)];
   double *backwardhihilo_h = new double[(deg1)*(nvr-1)];
   double *backwardlohilo_h = new double[(deg1)*(nvr-1)];
   double *backwardhilolo_h = new double[(deg1)*(nvr-1)];
   double *backwardlololo_h = new double[(deg1)*(nvr-1)];
   double *crosshihihi_h = new double[(deg1)*(nvr-2)];
   double *crosslohihi_h = new double[(deg1)*(nvr-2)];
   double *crosshilohi_h = new double[(deg1)*(nvr-2)];
   double *crosslolohi_h = new double[(deg1)*(nvr-2)];
   double *crosshihilo_h = new double[(deg1)*(nvr-2)];
   double *crosslohilo_h = new double[(deg1)*(nvr-2)];
   double *crosshilolo_h = new double[(deg1)*(nvr-2)];
   double *crosslololo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardhihihi_h,forwardhihihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlohihi_h,forwardlohihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhilohi_h,forwardhilohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlolohi_h,forwardlolohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhihilo_h,forwardhihilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlohilo_h,forwardlohilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhilolo_h,forwardhilolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlololo_h,forwardlololo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhihihi_h,backwardhihihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlohihi_h,backwardlohihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhilohi_h,backwardhilohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlolohi_h,backwardlolohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhihilo_h,backwardhihilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlohilo_h,backwardlohilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhilolo_h,backwardhilolo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlololo_h,backwardlololo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshihihi_h,crosshihihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslohihi_h,crosslohihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshilohi_h,crosshilohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslolohi_h,crosslolohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshihilo_h,crosshihilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslohilo_h,crosslohilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshilolo_h,crosshilolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslololo_h,crosslololo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputhihihi[dim][i] = forwardhihihi_h[offset+i];
      outputlohihi[dim][i] = forwardlohihi_h[offset+i];
      outputhilohi[dim][i] = forwardhilohi_h[offset+i];
      outputlolohi[dim][i] = forwardlolohi_h[offset+i];
      outputhihilo[dim][i] = forwardhihilo_h[offset+i];
      outputlohilo[dim][i] = forwardlohilo_h[offset+i];
      outputhilolo[dim][i] = forwardhilolo_h[offset+i];
      outputlololo[dim][i] = forwardlololo_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihihi[ix][i] = forwardhihihi_h[offset+i];
      outputlohihi[ix][i] = forwardlohihi_h[offset+i];
      outputhilohi[ix][i] = forwardhilohi_h[offset+i];
      outputlolohi[ix][i] = forwardlolohi_h[offset+i];
      outputhihilo[ix][i] = forwardhihilo_h[offset+i];
      outputlohilo[ix][i] = forwardlohilo_h[offset+i];
      outputhilolo[ix][i] = forwardhilolo_h[offset+i];
      outputlololo[ix][i] = forwardlololo_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihihi[ix][i] = backwardhihihi_h[offset+i];
      outputlohihi[ix][i] = backwardlohihi_h[offset+i];
      outputhilohi[ix][i] = backwardhilohi_h[offset+i];
      outputlolohi[ix][i] = backwardlolohi_h[offset+i];
      outputhihilo[ix][i] = backwardhihilo_h[offset+i];
      outputlohilo[ix][i] = backwardlohilo_h[offset+i];
      outputhilolo[ix][i] = backwardhilolo_h[offset+i];
      outputlololo[ix][i] = backwardlololo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputhihihi[ix][i] = crosshihihi_h[offset+i];
         outputlohihi[ix][i] = crosslohihi_h[offset+i];
         outputhilohi[ix][i] = crosshilohi_h[offset+i];
         outputlolohi[ix][i] = crosslolohi_h[offset+i];
         outputhihilo[ix][i] = crosshihilo_h[offset+i];
         outputlohilo[ix][i] = crosslohilo_h[offset+i];
         outputhilolo[ix][i] = crosshilolo_h[offset+i];
         outputlololo[ix][i] = crosslololo_h[offset+i];
      }
   }
}

void GPU_cmplx8_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
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
   const int deg1 = deg+1;      // length of all vectors
   double *inputrehihihi_d;  // inputrehihihi_d is inputrehihihi on the device
   double *inputrelohihi_d;  // inputrelohihi_d is inputrelohihi on the device
   double *inputrehilohi_d;  // inputrehilohi_d is inputrehilohi on the device
   double *inputrelolohi_d;  // inputrelolohi_d is inputrelolohi on the device
   double *inputrehihilo_d;  // inputrehihilo_d is inputrehihilo on the device
   double *inputrelohilo_d;  // inputrelohilo_d is inputrelohilo on the device
   double *inputrehilolo_d;  // inputrehilolo_d is inputrehilolo on the device
   double *inputrelololo_d;  // inputrelololo_d is inputrelololo on the device
   double *inputimhihihi_d;  // inputimhihihi_d is inputrehihihi on the device
   double *inputimlohihi_d;  // inputimlohihi_d is inputrelohihi on the device
   double *inputimhilohi_d;  // inputimhilohi_d is inputrehilohi on the device
   double *inputimlolohi_d;  // inputimlolohi_d is inputrelolohi on the device
   double *inputimhihilo_d;  // inputimhihilo_d is inputrehihilo on the device
   double *inputimlohilo_d;  // inputimlohilo_d is inputrelohilo on the device
   double *inputimhilolo_d;  // inputimhilolo_d is inputrehilolo on the device
   double *inputimlololo_d;  // inputimlololo_d is inputrelololo on the device
   double *forwardrehihihi_d;
   double *forwardrelohihi_d;
   double *forwardrehilohi_d;
   double *forwardrelolohi_d;
   double *forwardrehihilo_d;
   double *forwardrelohilo_d;
   double *forwardrehilolo_d;
   double *forwardrelololo_d;
   double *forwardimhihihi_d;
   double *forwardimlohihi_d;
   double *forwardimhilohi_d;
   double *forwardimlolohi_d;
   double *forwardimhihilo_d;
   double *forwardimlohilo_d;
   double *forwardimhilolo_d;
   double *forwardimlololo_d;
   double *backwardrehihihi_d;
   double *backwardrelohihi_d;
   double *backwardrehilohi_d;
   double *backwardrelolohi_d;
   double *backwardrehihilo_d;
   double *backwardrelohilo_d;
   double *backwardrehilolo_d;
   double *backwardrelololo_d;
   double *backwardimhihihi_d;
   double *backwardimlohihi_d;
   double *backwardimhilohi_d;
   double *backwardimlolohi_d;
   double *backwardimhihilo_d;
   double *backwardimlohilo_d;
   double *backwardimhilolo_d;
   double *backwardimlololo_d;
   double *crossrehihihi_d;
   double *crossrelohihi_d;
   double *crossrehilohi_d;
   double *crossrelolohi_d;
   double *crossrehihilo_d;
   double *crossrelohilo_d;
   double *crossrehilolo_d;
   double *crossrelololo_d;
   double *crossimhihihi_d;
   double *crossimlohihi_d;
   double *crossimhilohi_d;
   double *crossimlolohi_d;
   double *crossimhihilo_d;
   double *crossimlohilo_d;
   double *crossimhilolo_d;
   double *crossimlololo_d;
   double *cffrehihihi_d;     // cffrehihihi_d is cffrehihihi on the device
   double *cffrelohihi_d;     // cffrelohihi_d is cffrelohihi on the device
   double *cffrehilohi_d;     // cffrehilohi_d is cffrehilohi on the device
   double *cffrelolohi_d;     // cffrelolohi_d is cffrelolohi on the device
   double *cffrehihilo_d;     // cffrehihilo_d is cffrehihilo on the device
   double *cffrelohilo_d;     // cffrelohilo_d is cffrelohilo on the device
   double *cffrehilolo_d;     // cffrehilolo_d is cffrehilolo on the device
   double *cffrelololo_d;     // cffrelololo_d is cffrelololo on the device
   double *cffimhihihi_d;     // cffimhihihi_d is cffimhihihi on the device
   double *cffimlohihi_d;     // cffimlohihi_d is cffimlohihi on the device
   double *cffimhilohi_d;     // cffimhilohi_d is cffimhilohi on the device
   double *cffimlolohi_d;     // cffimlolohi_d is cffimlolohi on the device
   double *cffimhihilo_d;     // cffimhihilo_d is cffimhihilo on the device
   double *cffimlohilo_d;     // cffimlohilo_d is cffimlohilo on the device
   double *cffimhilolo_d;     // cffimhilolo_d is cffimhilolo on the device
   double *cffimlololo_d;     // cffimlololo_d is cffimlololo on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&cffrehihihi_d,szcff);
   cudaMalloc((void**)&cffrelohihi_d,szcff);
   cudaMalloc((void**)&cffrehilohi_d,szcff);
   cudaMalloc((void**)&cffrelolohi_d,szcff);
   cudaMalloc((void**)&cffrehihilo_d,szcff);
   cudaMalloc((void**)&cffrelohilo_d,szcff);
   cudaMalloc((void**)&cffrehilolo_d,szcff);
   cudaMalloc((void**)&cffrelololo_d,szcff);
   cudaMalloc((void**)&cffimhihihi_d,szcff);
   cudaMalloc((void**)&cffimlohihi_d,szcff);
   cudaMalloc((void**)&cffimhilohi_d,szcff);
   cudaMalloc((void**)&cffimlolohi_d,szcff);
   cudaMalloc((void**)&cffimhihilo_d,szcff);
   cudaMalloc((void**)&cffimlohilo_d,szcff);
   cudaMalloc((void**)&cffimhilolo_d,szcff);
   cudaMalloc((void**)&cffimlololo_d,szcff);
   cudaMalloc((void**)&inputrehihihi_d,szdim);
   cudaMalloc((void**)&inputrelohihi_d,szdim);
   cudaMalloc((void**)&inputrehilohi_d,szdim);
   cudaMalloc((void**)&inputrelolohi_d,szdim);
   cudaMalloc((void**)&inputrehihilo_d,szdim);
   cudaMalloc((void**)&inputrelohilo_d,szdim);
   cudaMalloc((void**)&inputrehilolo_d,szdim);
   cudaMalloc((void**)&inputrelololo_d,szdim);
   cudaMalloc((void**)&inputimhihihi_d,szdim);
   cudaMalloc((void**)&inputimlohihi_d,szdim);
   cudaMalloc((void**)&inputimhilohi_d,szdim);
   cudaMalloc((void**)&inputimlolohi_d,szdim);
   cudaMalloc((void**)&inputimhihilo_d,szdim);
   cudaMalloc((void**)&inputimlohilo_d,szdim);
   cudaMalloc((void**)&inputimhilolo_d,szdim);
   cudaMalloc((void**)&inputimlololo_d,szdim);
   cudaMalloc((void**)&forwardrehihihi_d,sznvr);
   cudaMalloc((void**)&forwardrelohihi_d,sznvr);
   cudaMalloc((void**)&forwardrehilohi_d,sznvr);
   cudaMalloc((void**)&forwardrelolohi_d,sznvr);
   cudaMalloc((void**)&forwardrehihilo_d,sznvr);
   cudaMalloc((void**)&forwardrelohilo_d,sznvr);
   cudaMalloc((void**)&forwardrehilolo_d,sznvr);
   cudaMalloc((void**)&forwardrelololo_d,sznvr);
   cudaMalloc((void**)&forwardimhihihi_d,sznvr);
   cudaMalloc((void**)&forwardimlohihi_d,sznvr);
   cudaMalloc((void**)&forwardimhilohi_d,sznvr);
   cudaMalloc((void**)&forwardimlolohi_d,sznvr);
   cudaMalloc((void**)&forwardimhihilo_d,sznvr);
   cudaMalloc((void**)&forwardimlohilo_d,sznvr);
   cudaMalloc((void**)&forwardimhilolo_d,sznvr);
   cudaMalloc((void**)&forwardimlololo_d,sznvr);
   cudaMalloc((void**)&backwardrehihihi_d,sznvr1);
   cudaMalloc((void**)&backwardrelohihi_d,sznvr1);
   cudaMalloc((void**)&backwardrehilohi_d,sznvr1);
   cudaMalloc((void**)&backwardrelolohi_d,sznvr1);
   cudaMalloc((void**)&backwardrehihilo_d,sznvr1);
   cudaMalloc((void**)&backwardrelohilo_d,sznvr1);
   cudaMalloc((void**)&backwardrehilolo_d,sznvr1);
   cudaMalloc((void**)&backwardrelololo_d,sznvr1);
   cudaMalloc((void**)&backwardimhihihi_d,sznvr1);
   cudaMalloc((void**)&backwardimlohihi_d,sznvr1);
   cudaMalloc((void**)&backwardimhilohi_d,sznvr1);
   cudaMalloc((void**)&backwardimlolohi_d,sznvr1);
   cudaMalloc((void**)&backwardimhihilo_d,sznvr1);
   cudaMalloc((void**)&backwardimlohilo_d,sznvr1);
   cudaMalloc((void**)&backwardimhilolo_d,sznvr1);
   cudaMalloc((void**)&backwardimlololo_d,sznvr1);
   cudaMalloc((void**)&crossrehihihi_d,sznvr2);
   cudaMalloc((void**)&crossrelohihi_d,sznvr2);
   cudaMalloc((void**)&crossrehilohi_d,sznvr2);
   cudaMalloc((void**)&crossrelolohi_d,sznvr2);
   cudaMalloc((void**)&crossrehihilo_d,sznvr2);
   cudaMalloc((void**)&crossrelohilo_d,sznvr2);
   cudaMalloc((void**)&crossrehilolo_d,sznvr2);
   cudaMalloc((void**)&crossrelololo_d,sznvr2);
   cudaMalloc((void**)&crossimhihihi_d,sznvr2);
   cudaMalloc((void**)&crossimlohihi_d,sznvr2);
   cudaMalloc((void**)&crossimhilohi_d,sznvr2);
   cudaMalloc((void**)&crossimlolohi_d,sznvr2);
   cudaMalloc((void**)&crossimhihilo_d,sznvr2);
   cudaMalloc((void**)&crossimlohilo_d,sznvr2);
   cudaMalloc((void**)&crossimhilolo_d,sznvr2);
   cudaMalloc((void**)&crossimlololo_d,sznvr2);

   double *inputrehihihi_h = new double[dim*(deg1)];
   double *inputrelohihi_h = new double[dim*(deg1)];
   double *inputrehilohi_h = new double[dim*(deg1)];
   double *inputrelolohi_h = new double[dim*(deg1)];
   double *inputrehihilo_h = new double[dim*(deg1)];
   double *inputrelohilo_h = new double[dim*(deg1)];
   double *inputrehilolo_h = new double[dim*(deg1)];
   double *inputrelololo_h = new double[dim*(deg1)];
   double *inputimhihihi_h = new double[dim*(deg1)];
   double *inputimlohihi_h = new double[dim*(deg1)];
   double *inputimhilohi_h = new double[dim*(deg1)];
   double *inputimlolohi_h = new double[dim*(deg1)];
   double *inputimhihilo_h = new double[dim*(deg1)];
   double *inputimlohilo_h = new double[dim*(deg1)];
   double *inputimhilolo_h = new double[dim*(deg1)];
   double *inputimlololo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrehihihi_h[ix] = inputrehihihi[i][j];
         inputrelohihi_h[ix] = inputrelohihi[i][j];
         inputrehilohi_h[ix] = inputrehilohi[i][j];
         inputrelolohi_h[ix] = inputrelolohi[i][j];
         inputrehihilo_h[ix] = inputrehihilo[i][j];
         inputrelohilo_h[ix] = inputrelohilo[i][j];
         inputrehilolo_h[ix] = inputrehilolo[i][j];
         inputrelololo_h[ix] = inputrelololo[i][j];
         inputimhihihi_h[ix] = inputimhihihi[i][j];
         inputimlohihi_h[ix] = inputimlohihi[i][j];
         inputimhilohi_h[ix] = inputimhilohi[i][j];
         inputimlolohi_h[ix] = inputimlolohi[i][j];
         inputimhihilo_h[ix] = inputimhihilo[i][j];
         inputimlohilo_h[ix] = inputimlohilo[i][j];
         inputimhilolo_h[ix] = inputimhilolo[i][j];
         inputimlololo_h[ix++] = inputimlololo[i][j];
      }

   cudaMemcpy(cffrehihihi_d,cffrehihihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelohihi_d,cffrelohihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehilohi_d,cffrehilohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelolohi_d,cffrelolohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehihilo_d,cffrehihilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelohilo_d,cffrelohilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehilolo_d,cffrehilolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelololo_d,cffrelololo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhihihi_d,cffimhihihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlohihi_d,cffimlohihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhilohi_d,cffimhilohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlolohi_d,cffimlolohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhihilo_d,cffimhihilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlohilo_d,cffimlohilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhilolo_d,cffimhilolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlololo_d,cffimlololo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehihihi_d,inputrehihihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelohihi_d,inputrelohihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehilohi_d,inputrehilohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelolohi_d,inputrelolohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehihilo_d,inputrehihilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelohilo_d,inputrelohilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehilolo_d,inputrehilolo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelololo_d,inputrelololo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhihihi_d,inputimhihihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlohihi_d,inputimlohihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhilohi_d,inputimhilohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlolohi_d,inputimlolohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhihilo_d,inputimhihilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlohilo_d,inputimlohilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhilolo_d,inputimhilolo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlololo_d,inputimlololo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx8_speel
         (BS,nvr,deg,idx,
          cffrehihihi_d,cffrelohihi_d,cffrehilohi_d,cffrelolohi_d,
          cffrehihilo_d,cffrelohilo_d,cffrehilolo_d,cffrelololo_d,
          cffimhihihi_d,cffimlohihi_d,cffimhilohi_d,cffimlolohi_d,
          cffimhihilo_d,cffimlohilo_d,cffimhilolo_d,cffimlololo_d,
          inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
          inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
          inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
          inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
          forwardrehihihi_d,forwardrelohihi_d,
          forwardrehilohi_d,forwardrelolohi_d,
          forwardrehihilo_d,forwardrelohilo_d,
          forwardrehilolo_d,forwardrelololo_d,
          forwardimhihihi_d,forwardimlohihi_d,
          forwardimhilohi_d,forwardimlolohi_d,
          forwardimhihilo_d,forwardimlohilo_d,
          forwardimhilolo_d,forwardimlololo_d,
          backwardrehihihi_d,backwardrelohihi_d,
          backwardrehilohi_d,backwardrelolohi_d,
          backwardrehihilo_d,backwardrelohilo_d,
          backwardrehilolo_d,backwardrelololo_d,
          backwardimhihihi_d,backwardimlohihi_d,
          backwardimhilohi_d,backwardimlolohi_d,
          backwardimhihilo_d,backwardimlohilo_d,
          backwardimhilolo_d,backwardimlololo_d,
          crossrehihihi_d,crossrelohihi_d,crossrehilohi_d,crossrelolohi_d,
          crossrehihilo_d,crossrelohilo_d,crossrehilolo_d,crossrelololo_d,
          crossimhihihi_d,crossimlohihi_d,crossimhilohi_d,crossimlolohi_d,
          crossimhihilo_d,crossimlohilo_d,crossimhilolo_d,crossimlololo_d);
   }
   double *forwardrehihihi_h = new double[(deg1)*nvr];
   double *forwardrelohihi_h = new double[(deg1)*nvr];
   double *forwardrehilohi_h = new double[(deg1)*nvr];
   double *forwardrelolohi_h = new double[(deg1)*nvr];
   double *forwardrehihilo_h = new double[(deg1)*nvr];
   double *forwardrelohilo_h = new double[(deg1)*nvr];
   double *forwardrehilolo_h = new double[(deg1)*nvr];
   double *forwardrelololo_h = new double[(deg1)*nvr];
   double *forwardimhihihi_h = new double[(deg1)*nvr];
   double *forwardimlohihi_h = new double[(deg1)*nvr];
   double *forwardimhilohi_h = new double[(deg1)*nvr];
   double *forwardimlolohi_h = new double[(deg1)*nvr];
   double *forwardimhihilo_h = new double[(deg1)*nvr];
   double *forwardimlohilo_h = new double[(deg1)*nvr];
   double *forwardimhilolo_h = new double[(deg1)*nvr];
   double *forwardimlololo_h = new double[(deg1)*nvr];
   double *backwardrehihihi_h = new double[(deg1)*(nvr-1)];
   double *backwardrelohihi_h = new double[(deg1)*(nvr-1)];
   double *backwardrehilohi_h = new double[(deg1)*(nvr-1)];
   double *backwardrelolohi_h = new double[(deg1)*(nvr-1)];
   double *backwardrehihilo_h = new double[(deg1)*(nvr-1)];
   double *backwardrelohilo_h = new double[(deg1)*(nvr-1)];
   double *backwardrehilolo_h = new double[(deg1)*(nvr-1)];
   double *backwardrelololo_h = new double[(deg1)*(nvr-1)];
   double *backwardimhihihi_h = new double[(deg1)*(nvr-1)];
   double *backwardimlohihi_h = new double[(deg1)*(nvr-1)];
   double *backwardimhilohi_h = new double[(deg1)*(nvr-1)];
   double *backwardimlolohi_h = new double[(deg1)*(nvr-1)];
   double *backwardimhihilo_h = new double[(deg1)*(nvr-1)];
   double *backwardimlohilo_h = new double[(deg1)*(nvr-1)];
   double *backwardimhilolo_h = new double[(deg1)*(nvr-1)];
   double *backwardimlololo_h = new double[(deg1)*(nvr-1)];
   double *crossrehihihi_h = new double[(deg1)*(nvr-2)];
   double *crossrelohihi_h = new double[(deg1)*(nvr-2)];
   double *crossrehilohi_h = new double[(deg1)*(nvr-2)];
   double *crossrelolohi_h = new double[(deg1)*(nvr-2)];
   double *crossrehihilo_h = new double[(deg1)*(nvr-2)];
   double *crossrelohilo_h = new double[(deg1)*(nvr-2)];
   double *crossrehilolo_h = new double[(deg1)*(nvr-2)];
   double *crossrelololo_h = new double[(deg1)*(nvr-2)];
   double *crossimhihihi_h = new double[(deg1)*(nvr-2)];
   double *crossimlohihi_h = new double[(deg1)*(nvr-2)];
   double *crossimhilohi_h = new double[(deg1)*(nvr-2)];
   double *crossimlolohi_h = new double[(deg1)*(nvr-2)];
   double *crossimhihilo_h = new double[(deg1)*(nvr-2)];
   double *crossimlohilo_h = new double[(deg1)*(nvr-2)];
   double *crossimhilolo_h = new double[(deg1)*(nvr-2)];
   double *crossimlololo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrehihihi_h,
              forwardrehihihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelohihi_h,
              forwardrelohihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrehilohi_h,
              forwardrehilohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelolohi_h,
              forwardrelolohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrehihilo_h,
              forwardrehihilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelohilo_h,
              forwardrelohilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrehilolo_h,
              forwardrehilolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelololo_h,
              forwardrelololo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhihihi_h,
              forwardimhihihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlohihi_h,
              forwardimlohihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhilohi_h,
              forwardimhilohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlolohi_h,
              forwardimlolohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhihilo_h,
              forwardimhihilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlohilo_h,
              forwardimlohilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhilolo_h,
              forwardimhilolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlololo_h,
              forwardimlololo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehihihi_h,
              backwardrehihihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelohihi_h,
              backwardrelohihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehilohi_h,
              backwardrehilohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelolohi_h,
              backwardrelolohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehihilo_h,
              backwardrehihilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelohilo_h,
              backwardrelohilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehilolo_h,
              backwardrehilolo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelololo_h,
              backwardrelololo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhihihi_h,
              backwardimhihihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlohihi_h,
              backwardimlohihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhilohi_h,
              backwardimhilohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlolohi_h,
              backwardimlolohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhihilo_h,
              backwardimhihilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlohilo_h,
              backwardimlohilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhilolo_h,
              backwardimhilolo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlololo_h,
              backwardimlololo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehihihi_h,crossrehihihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelohihi_h,crossrelohihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehilohi_h,crossrehilohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelolohi_h,crossrelolohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehihilo_h,crossrehihilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelohilo_h,crossrelohilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehilolo_h,crossrehilolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelololo_h,crossrelololo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhihihi_h,crossimhihihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlohihi_h,crossimlohihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhilohi_h,crossimhilohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlolohi_h,crossimlolohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhihilo_h,crossimhihilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlohilo_h,crossimlohilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhilolo_h,crossimhilolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlololo_h,crossimlololo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputrehihihi[dim][i] = forwardrehihihi_h[offset+i];
      outputrelohihi[dim][i] = forwardrelohihi_h[offset+i];
      outputrehilohi[dim][i] = forwardrehilohi_h[offset+i];
      outputrelolohi[dim][i] = forwardrelolohi_h[offset+i];
      outputrehihilo[dim][i] = forwardrehihilo_h[offset+i];
      outputrelohilo[dim][i] = forwardrelohilo_h[offset+i];
      outputrehilolo[dim][i] = forwardrehilolo_h[offset+i];
      outputrelololo[dim][i] = forwardrelololo_h[offset+i];
      outputimhihihi[dim][i] = forwardimhihihi_h[offset+i];
      outputimlohihi[dim][i] = forwardimlohihi_h[offset+i];
      outputimhilohi[dim][i] = forwardimhilohi_h[offset+i];
      outputimlolohi[dim][i] = forwardimlolohi_h[offset+i];
      outputimhihilo[dim][i] = forwardimhihilo_h[offset+i];
      outputimlohilo[dim][i] = forwardimlohilo_h[offset+i];
      outputimhilolo[dim][i] = forwardimhilolo_h[offset+i];
      outputimlololo[dim][i] = forwardimlololo_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputrehihihi[ix][i] = forwardrehihihi_h[offset+i];
      outputrelohihi[ix][i] = forwardrelohihi_h[offset+i];
      outputrehilohi[ix][i] = forwardrehilohi_h[offset+i];
      outputrelolohi[ix][i] = forwardrelolohi_h[offset+i];
      outputrehihilo[ix][i] = forwardrehihilo_h[offset+i];
      outputrelohilo[ix][i] = forwardrelohilo_h[offset+i];
      outputrehilolo[ix][i] = forwardrehilolo_h[offset+i];
      outputrelololo[ix][i] = forwardrelololo_h[offset+i];
      outputimhihihi[ix][i] = forwardimhihihi_h[offset+i];
      outputimlohihi[ix][i] = forwardimlohihi_h[offset+i];
      outputimhilohi[ix][i] = forwardimhilohi_h[offset+i];
      outputimlolohi[ix][i] = forwardimlolohi_h[offset+i];
      outputimhihilo[ix][i] = forwardimhihilo_h[offset+i];
      outputimlohilo[ix][i] = forwardimlohilo_h[offset+i];
      outputimhilolo[ix][i] = forwardimhilolo_h[offset+i];
      outputimlololo[ix][i] = forwardimlololo_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputrehihihi[ix][i] = backwardrehihihi_h[offset+i];
      outputrelohihi[ix][i] = backwardrelohihi_h[offset+i];
      outputrehilohi[ix][i] = backwardrehilohi_h[offset+i];
      outputrelolohi[ix][i] = backwardrelolohi_h[offset+i];
      outputrehihilo[ix][i] = backwardrehihilo_h[offset+i];
      outputrelohilo[ix][i] = backwardrelohilo_h[offset+i];
      outputrehilolo[ix][i] = backwardrehilolo_h[offset+i];
      outputrelololo[ix][i] = backwardrelololo_h[offset+i];
      outputimhihihi[ix][i] = backwardimhihihi_h[offset+i];
      outputimlohihi[ix][i] = backwardimlohihi_h[offset+i];
      outputimhilohi[ix][i] = backwardimhilohi_h[offset+i];
      outputimlolohi[ix][i] = backwardimlolohi_h[offset+i];
      outputimhihilo[ix][i] = backwardimhihilo_h[offset+i];
      outputimlohilo[ix][i] = backwardimlohilo_h[offset+i];
      outputimhilolo[ix][i] = backwardimhilolo_h[offset+i];
      outputimlololo[ix][i] = backwardimlololo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrehihihi[ix][i] = crossrehihihi_h[offset+i];
         outputrelohihi[ix][i] = crossrelohihi_h[offset+i];
         outputrehilohi[ix][i] = crossrehilohi_h[offset+i];
         outputrelolohi[ix][i] = crossrelolohi_h[offset+i];
         outputrehihilo[ix][i] = crossrehihilo_h[offset+i];
         outputrelohilo[ix][i] = crossrelohilo_h[offset+i];
         outputrehilolo[ix][i] = crossrehilolo_h[offset+i];
         outputrelololo[ix][i] = crossrelololo_h[offset+i];
         outputimhihihi[ix][i] = crossimhihihi_h[offset+i];
         outputimlohihi[ix][i] = crossimlohihi_h[offset+i];
         outputimhilohi[ix][i] = crossimhilohi_h[offset+i];
         outputimlolohi[ix][i] = crossimlolohi_h[offset+i];
         outputimhihilo[ix][i] = crossimhihilo_h[offset+i];
         outputimlohilo[ix][i] = crossimlohilo_h[offset+i];
         outputimhilolo[ix][i] = crossimhilolo_h[offset+i];
         outputimlololo[ix][i] = crossimlololo_h[offset+i];
      }
   }
}
