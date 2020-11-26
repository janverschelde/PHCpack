// The file dbl4_monomials_kernels.cu defines kernels specified
// in dbl4_monomials_kernels.h.

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

#include "dbl4_convolutions_kernels.cu"
#include "dbl4_monomials_kernels.h"

void GPU_dbl4_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double *inputhihi, double *inputlohi, double *inputhilo, double *inputlolo,
   double *forwardhihi, double *forwardlohi,
   double *forwardhilo, double *forwardlolo,
   double *backwardhihi, double *backwardlohi,
   double *backwardhilo, double *backwardlolo,
   double *crosshihi, double *crosslohi,
   double *crosshilo, double *crosslolo )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl4_padded_convolute<<<1,BS>>>
      (cffhihi,cfflohi,cffhilo,cfflolo,
       &inputhihi[ix1],&inputlohi[ix1],&inputhilo[ix1],&inputlolo[ix1],
       forwardhihi,forwardlohi,forwardhilo,forwardlolo,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl4_padded_convolute<<<1,BS>>>
         (&forwardhihi[ix1],&forwardlohi[ix1],
          &forwardhilo[ix1],&forwardlolo[ix1],
          &inputhihi[ix2],&inputlohi[ix2],
          &inputhilo[ix2],&inputlolo[ix2],
          &forwardhihi[ix3],&forwardlohi[ix3],
          &forwardhilo[ix3],&forwardlolo[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl4_padded_convolute<<<1,BS>>>
         (&inputhihi[ix1],&inputlohi[ix1],&inputhilo[ix1],&inputlolo[ix1],
          &inputhihi[ix2],&inputlohi[ix2],&inputhilo[ix2],&inputlolo[ix2],
          backwardhihi,backwardlohi,backwardhilo,backwardlolo,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl4_padded_convolute<<<1,BS>>>
            (&backwardhihi[ix1],&backwardlohi[ix1],
             &backwardhilo[ix1],&backwardlolo[ix1],
             &inputhihi[ix2],&inputlohi[ix2],&inputhilo[ix2],&inputlolo[ix2],
             &backwardhihi[ix3],&backwardlohi[ix3],
             &backwardhilo[ix3],&backwardlolo[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl4_padded_convolute<<<1,BS>>>
         (&backwardhihi[ix3],&backwardlohi[ix3],
          &backwardhilo[ix3],&backwardlolo[ix3],
          cffhihi,cfflohi,cffhilo,cfflolo,
          &backwardhihi[ix2],&backwardlohi[ix2],
          &backwardhilo[ix2],&backwardlolo[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl4_padded_convolute<<<1,BS>>>
            (forwardhihi,forwardlohi,forwardhilo,forwardlolo,
             &inputhihi[ix2],&inputlohi[ix2],&inputhilo[ix2],&inputlolo[ix2],
             crosshihi,crosslohi,crosshilo,crosslolo,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl4_padded_convolute<<<1,BS>>>
               (&forwardhihi[ix1],&forwardlohi[ix1],
                &forwardhilo[ix1],&forwardlolo[ix1],
                &backwardhihi[ix2],&backwardlohi[ix2],
                &backwardhilo[ix2],&backwardlolo[ix2],
                &crosshihi[ix1],&crosslohi[ix1],
                &crosshilo[ix1],&crosslolo[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl4_padded_convolute<<<1,BS>>>
            (&forwardhihi[ix1],&forwardlohi[ix1],
             &forwardhilo[ix1],&forwardlolo[ix1],
             &inputhihi[ix2],&inputlohi[ix2],&inputhilo[ix2],&inputlolo[ix2],
             &crosshihi[ix1],&crosslohi[ix1],&crosshilo[ix1],&crosslolo[ix1],
             deg1);
      }
   }
}

void GPU_cmplx4_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double *inputrehihi, double *inputrelohi,
   double *inputrehilo, double *inputrelolo,
   double *inputimhihi, double *inputimlohi,
   double *inputimhilo, double *inputimlolo,
   double *forwardrehihi, double *forwardrelohi,
   double *forwardrehilo, double *forwardrelolo,
   double *forwardimhihi, double *forwardimlohi,
   double *forwardimhilo, double *forwardimlolo,
   double *backwardrehihi, double *backwardrelohi,
   double *backwardrehilo, double *backwardrelolo,
   double *backwardimhihi, double *backwardimlohi,
   double *backwardimhilo, double *backwardimlolo,
   double *crossrehihi, double *crossrelohi,
   double *crossrehilo, double *crossrelolo,
   double *crossimhihi, double *crossimlohi,
   double *crossimhilo, double *crossimlolo )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   cmplx4_padded_convolute<<<1,BS>>>
      (cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       &inputrehihi[ix1],&inputrelohi[ix1],&inputrehilo[ix1],&inputrelolo[ix1],
       &inputimhihi[ix1],&inputimlohi[ix1],&inputimhilo[ix1],&inputimlolo[ix1],
       forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
       forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,deg1); 

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-i]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      cmplx4_padded_convolute<<<1,BS>>>
         (&forwardrehihi[ix1],&forwardrelohi[ix1],
          &forwardrehilo[ix1],&forwardrelolo[ix1],
          &forwardimhihi[ix1],&forwardimlohi[ix1],
          &forwardimhilo[ix1],&forwardimlolo[ix1],
          &inputrehihi[ix2],&inputrelohi[ix2],
          &inputrehilo[ix2],&inputrelolo[ix2],
          &inputimhihi[ix2],&inputimlohi[ix2],
          &inputimhilo[ix2],&inputimlolo[ix2],
          &forwardrehihi[ix3],&forwardrelohi[ix3],
          &forwardrehilo[ix3],&forwardrelolo[ix3],
          &forwardimhihi[ix3],&forwardimlohi[ix3],
          &forwardimhilo[ix3],&forwardimlolo[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      cmplx4_padded_convolute<<<1,BS>>>
         (&inputrehihi[ix1],&inputrelohi[ix1],
          &inputrehilo[ix1],&inputrelolo[ix1],
          &inputimhihi[ix1],&inputimlohi[ix1],
          &inputimhilo[ix1],&inputimlolo[ix1],
          &inputrehihi[ix2],&inputrelohi[ix2],
          &inputrehilo[ix2],&inputrelolo[ix2],
          &inputimhihi[ix2],&inputimlohi[ix2],
          &inputimhilo[ix2],&inputimlolo[ix2],
          backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
          backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         cmplx4_padded_convolute<<<1,BS>>>
            (&backwardrehihi[ix1],&backwardrelohi[ix1],
             &backwardrehilo[ix1],&backwardrelolo[ix1],
             &backwardimhihi[ix1],&backwardimlohi[ix1],
             &backwardimhilo[ix1],&backwardimlolo[ix1],
             &inputrehihi[ix2],&inputrelohi[ix2],
             &inputrehilo[ix2],&inputrelolo[ix2],
             &inputimhihi[ix2],&inputimlohi[ix2],
             &inputimhilo[ix2],&inputimlolo[ix2],
             &backwardrehihi[ix3],&backwardrelohi[ix3],
             &backwardrehilo[ix3],&backwardrelolo[ix3],
             &backwardimhihi[ix3],&backwardimlohi[ix3],
             &backwardimhilo[ix3],&backwardimlolo[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      cmplx4_padded_convolute<<<1,BS>>>
         (&backwardrehihi[ix3],&backwardrelohi[ix3],
          &backwardrehilo[ix3],&backwardrelolo[ix3],
          &backwardimhihi[ix3],&backwardimlohi[ix3],
          &backwardimhilo[ix3],&backwardimlolo[ix3],
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,
          &backwardrehihi[ix2],&backwardrelohi[ix2],
          &backwardrehilo[ix2],&backwardrelolo[ix2],
          &backwardimhihi[ix2],&backwardimlohi[ix2],
          &backwardimhilo[ix2],&backwardimlolo[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         cmplx4_padded_convolute<<<1,BS>>>
            (forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
             forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
             &inputrehihi[ix2],&inputrelohi[ix2],
             &inputrehilo[ix2],&inputrelolo[ix2],
             &inputimhihi[ix2],&inputimlohi[ix2],
             &inputimhilo[ix2],&inputimlolo[ix2],
             crossrehihi,crossrelohi,crossrehilo,crossrelolo,
             crossimhihi,crossimlohi,crossimhilo,crossimlolo,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            cmplx4_padded_convolute<<<1,BS>>>
               (&forwardrehihi[ix1],&forwardrelohi[ix1],
                &forwardrehilo[ix1],&forwardrelolo[ix1],
                &forwardimhihi[ix1],&forwardimlohi[ix1],
                &forwardimhilo[ix1],&forwardimlolo[ix1],
                &backwardrehihi[ix2],&backwardrelohi[ix2],
                &backwardrehilo[ix2],&backwardrelolo[ix2],
                &backwardimhihi[ix2],&backwardimlohi[ix2],
                &backwardimhilo[ix2],&backwardimlolo[ix2],
                &crossrehihi[ix1],&crossrelohi[ix1],
                &crossrehilo[ix1],&crossrelolo[ix1],
                &crossimhihi[ix1],&crossimlohi[ix1],
                &crossimhilo[ix1],&crossimlolo[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         cmplx4_padded_convolute<<<1,BS>>>
            (&forwardrehihi[ix1],&forwardrelohi[ix1],
             &forwardrehilo[ix1],&forwardrelolo[ix1],
             &forwardimhihi[ix1],&forwardimlohi[ix1],
             &forwardimhilo[ix1],&forwardimlolo[ix1],
             &inputrehihi[ix2],&inputrelohi[ix2],
             &inputrehilo[ix2],&inputrelolo[ix2],
             &inputimhihi[ix2],&inputimlohi[ix2],
             &inputimhilo[ix2],&inputimlolo[ix2],
             &crossrehihi[ix1],&crossrelohi[ix1],
             &crossrehilo[ix1],&crossrelolo[ix1],
             &crossimhihi[ix1],&crossimlohi[ix1],
             &crossimhilo[ix1],&crossimlolo[ix1],deg1);
      }
   }
}

void GPU_dbl4_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo )
{
   const int deg1 = deg+1;   // length of all vectors
   double *inputhihi_d;      // inputhihi_d is input on the device
   double *inputlohi_d;      // inputlohi_d is input on the device
   double *inputhilo_d;      // inputhilo_d is input on the device
   double *inputlolo_d;      // inputlolo_d is input on the device
   double *forwardhihi_d;    // highest forward products on the device
   double *forwardlohi_d;    // second highest forward products on the device
   double *forwardhilo_d;    // second lowest forward products on the device
   double *forwardlolo_d;    // lowest forward products on the device
   double *backwardhihi_d;   // highest backward products on the device
   double *backwardlohi_d;   // second highest backward products on the device
   double *backwardhilo_d;   // second lowest backward products on the device
   double *backwardlolo_d;   // lowest backward products on the device
   double *crosshihi_d;      // highest cross products on the device
   double *crosslohi_d;      // second highest cross products on the device
   double *crosshilo_d;      // second lowest cross products on the device
   double *crosslolo_d;      // lowest cross products on the device
   double *cffhihi_d;        // cffhihi_d is cffhihi on device
   double *cfflohi_d;        // cfflohi_d is cfflohi on device
   double *cffhilo_d;        // cffhilo_d is cffhilo on device
   double *cfflolo_d;        // cfflolo_d is cfflolo on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);

   cudaMalloc((void**)&cffhihi_d,szcff);
   cudaMalloc((void**)&cfflohi_d,szcff);
   cudaMalloc((void**)&cffhilo_d,szcff);
   cudaMalloc((void**)&cfflolo_d,szcff);
   cudaMalloc((void**)&inputhihi_d,szdim);
   cudaMalloc((void**)&inputlohi_d,szdim);
   cudaMalloc((void**)&inputhilo_d,szdim);
   cudaMalloc((void**)&inputlolo_d,szdim);
   cudaMalloc((void**)&forwardhihi_d,sznvr);
   cudaMalloc((void**)&forwardlohi_d,sznvr);
   cudaMalloc((void**)&forwardhilo_d,sznvr);
   cudaMalloc((void**)&forwardlolo_d,sznvr);
   cudaMalloc((void**)&backwardhihi_d,sznvr1);
   cudaMalloc((void**)&backwardlohi_d,sznvr1);
   cudaMalloc((void**)&backwardhilo_d,sznvr1);
   cudaMalloc((void**)&backwardlolo_d,sznvr1);
   cudaMalloc((void**)&crosshihi_d,sznvr2);
   cudaMalloc((void**)&crosslohi_d,sznvr2);
   cudaMalloc((void**)&crosshilo_d,sznvr2);
   cudaMalloc((void**)&crosslolo_d,sznvr2);

   double *inputhihi_h = new double[dim*(deg1)];
   double *inputlohi_h = new double[dim*(deg1)];
   double *inputhilo_h = new double[dim*(deg1)];
   double *inputlolo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputhihi_h[ix] = inputhihi[i][j];
         inputlohi_h[ix] = inputlohi[i][j];
         inputhilo_h[ix] = inputhilo[i][j];
         inputlolo_h[ix++] = inputlolo[i][j];
      }

   cudaMemcpy(cffhihi_d,cffhihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflohi_d,cfflohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhilo_d,cffhilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflolo_d,cfflolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhihi_d,inputhihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlohi_d,inputlohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhilo_d,inputhilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlolo_d,inputlolo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl4_speel
         (BS,nvr,deg,idx,cffhihi_d,cfflohi_d,cffhilo_d,cfflolo_d,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          forwardhihi_d,forwardlohi_d,forwardhilo_d,forwardlolo_d,
          backwardhihi_d,backwardlohi_d,backwardhilo_d,backwardlolo_d,
          crosshihi_d,crosslohi_d,crosshilo_d,crosslolo_d);
   }
   double *forwardhihi_h = new double[(deg1)*nvr];
   double *forwardlohi_h = new double[(deg1)*nvr];
   double *forwardhilo_h = new double[(deg1)*nvr];
   double *forwardlolo_h = new double[(deg1)*nvr];
   double *backwardhihi_h = new double[(deg1)*(nvr-1)];
   double *backwardlohi_h = new double[(deg1)*(nvr-1)];
   double *backwardhilo_h = new double[(deg1)*(nvr-1)];
   double *backwardlolo_h = new double[(deg1)*(nvr-1)];
   double *crosshihi_h = new double[(deg1)*(nvr-2)];
   double *crosslohi_h = new double[(deg1)*(nvr-2)];
   double *crosshilo_h = new double[(deg1)*(nvr-2)];
   double *crosslolo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardhihi_h,forwardhihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlohi_h,forwardlohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhilo_h,forwardhilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlolo_h,forwardlolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhihi_h,backwardhihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlohi_h,backwardlohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhilo_h,backwardhilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlolo_h,backwardlolo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshihi_h,crosshihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslohi_h,crosslohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshilo_h,crosshilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslolo_h,crosslolo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputhihi[dim][i] = forwardhihi_h[offset+i];
      outputlohi[dim][i] = forwardlohi_h[offset+i];
      outputhilo[dim][i] = forwardhilo_h[offset+i];
      outputlolo[dim][i] = forwardlolo_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihi[ix][i] = forwardhihi_h[offset+i];
      outputlohi[ix][i] = forwardlohi_h[offset+i];
      outputhilo[ix][i] = forwardhilo_h[offset+i];
      outputlolo[ix][i] = forwardlolo_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihi[ix][i] = backwardhihi_h[offset+i];
      outputlohi[ix][i] = backwardlohi_h[offset+i];
      outputhilo[ix][i] = backwardhilo_h[offset+i];
      outputlolo[ix][i] = backwardlolo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputhihi[ix][i] = crosshihi_h[offset+i];
         outputlohi[ix][i] = crosslohi_h[offset+i];
         outputhilo[ix][i] = crosshilo_h[offset+i];
         outputlolo[ix][i] = crosslolo_h[offset+i];
      }
   }
}

void GPU_cmplx4_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
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
   const int deg1 = deg+1;      // length of all vectors
   double *inputrehihi_d;       // inputrehihi_d is inputrehihi on the device
   double *inputrelohi_d;       // inputrelohi_d is inputrelohi on the device
   double *inputrehilo_d;       // inputrehilo_d is inputrehilo on the device
   double *inputrelolo_d;       // inputrelolo_d is inputrelolo on the device
   double *inputimhihi_d;       // inputimhihi_d is inputrehihi on the device
   double *inputimlohi_d;       // inputimlohi_d is inputrelohi on the device
   double *inputimhilo_d;       // inputimhilo_d is inputrehilo on the device
   double *inputimlolo_d;       // inputimlolo_d is inputrelolo on the device
   double *forwardrehihi_d;
   double *forwardrelohi_d;
   double *forwardrehilo_d;
   double *forwardrelolo_d;
   double *forwardimhihi_d;
   double *forwardimlohi_d;
   double *forwardimhilo_d;
   double *forwardimlolo_d;
   double *backwardrehihi_d;
   double *backwardrelohi_d;
   double *backwardrehilo_d;
   double *backwardrelolo_d;
   double *backwardimhihi_d;
   double *backwardimlohi_d;
   double *backwardimhilo_d;
   double *backwardimlolo_d;
   double *crossrehihi_d;
   double *crossrelohi_d;
   double *crossrehilo_d;
   double *crossrelolo_d;
   double *crossimhihi_d;
   double *crossimlohi_d;
   double *crossimhilo_d;
   double *crossimlolo_d;
   double *cffrehihi_d;             // cffrehihi_d is cffrehihi on the device
   double *cffrelohi_d;             // cffrelohi_d is cffrelohi on the device
   double *cffrehilo_d;             // cffrehilo_d is cffrehilo on the device
   double *cffrelolo_d;             // cffrelolo_d is cffrelolo on the device
   double *cffimhihi_d;             // cffimhihi_d is cffimhihi on the device
   double *cffimlohi_d;             // cffimlohi_d is cffimlohi on the device
   double *cffimhilo_d;             // cffimhilo_d is cffimhilo on the device
   double *cffimlolo_d;             // cffimlolo_d is cffimlolo on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&cffrehihi_d,szcff);
   cudaMalloc((void**)&cffrelohi_d,szcff);
   cudaMalloc((void**)&cffrehilo_d,szcff);
   cudaMalloc((void**)&cffrelolo_d,szcff);
   cudaMalloc((void**)&cffimhihi_d,szcff);
   cudaMalloc((void**)&cffimlohi_d,szcff);
   cudaMalloc((void**)&cffimhilo_d,szcff);
   cudaMalloc((void**)&cffimlolo_d,szcff);
   cudaMalloc((void**)&inputrehihi_d,szdim);
   cudaMalloc((void**)&inputrelohi_d,szdim);
   cudaMalloc((void**)&inputrehilo_d,szdim);
   cudaMalloc((void**)&inputrelolo_d,szdim);
   cudaMalloc((void**)&inputimhihi_d,szdim);
   cudaMalloc((void**)&inputimlohi_d,szdim);
   cudaMalloc((void**)&inputimhilo_d,szdim);
   cudaMalloc((void**)&inputimlolo_d,szdim);
   cudaMalloc((void**)&forwardrehihi_d,sznvr);
   cudaMalloc((void**)&forwardrelohi_d,sznvr);
   cudaMalloc((void**)&forwardrehilo_d,sznvr);
   cudaMalloc((void**)&forwardrelolo_d,sznvr);
   cudaMalloc((void**)&forwardimhihi_d,sznvr);
   cudaMalloc((void**)&forwardimlohi_d,sznvr);
   cudaMalloc((void**)&forwardimhilo_d,sznvr);
   cudaMalloc((void**)&forwardimlolo_d,sznvr);
   cudaMalloc((void**)&backwardrehihi_d,sznvr1);
   cudaMalloc((void**)&backwardrelohi_d,sznvr1);
   cudaMalloc((void**)&backwardrehilo_d,sznvr1);
   cudaMalloc((void**)&backwardrelolo_d,sznvr1);
   cudaMalloc((void**)&backwardimhihi_d,sznvr1);
   cudaMalloc((void**)&backwardimlohi_d,sznvr1);
   cudaMalloc((void**)&backwardimhilo_d,sznvr1);
   cudaMalloc((void**)&backwardimlolo_d,sznvr1);
   cudaMalloc((void**)&crossrehihi_d,sznvr2);
   cudaMalloc((void**)&crossrelohi_d,sznvr2);
   cudaMalloc((void**)&crossrehilo_d,sznvr2);
   cudaMalloc((void**)&crossrelolo_d,sznvr2);
   cudaMalloc((void**)&crossimhihi_d,sznvr2);
   cudaMalloc((void**)&crossimlohi_d,sznvr2);
   cudaMalloc((void**)&crossimhilo_d,sznvr2);
   cudaMalloc((void**)&crossimlolo_d,sznvr2);

   double *inputrehihi_h = new double[dim*(deg1)];
   double *inputrelohi_h = new double[dim*(deg1)];
   double *inputrehilo_h = new double[dim*(deg1)];
   double *inputrelolo_h = new double[dim*(deg1)];
   double *inputimhihi_h = new double[dim*(deg1)];
   double *inputimlohi_h = new double[dim*(deg1)];
   double *inputimhilo_h = new double[dim*(deg1)];
   double *inputimlolo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrehihi_h[ix] = inputrehihi[i][j];
         inputrelohi_h[ix] = inputrelohi[i][j];
         inputrehilo_h[ix] = inputrehilo[i][j];
         inputrelolo_h[ix] = inputrelolo[i][j];
         inputimhihi_h[ix] = inputimhihi[i][j];
         inputimlohi_h[ix] = inputimlohi[i][j];
         inputimhilo_h[ix] = inputimhilo[i][j];
         inputimlolo_h[ix++] = inputimlolo[i][j];
      }

   cudaMemcpy(cffrehihi_d,cffrehihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelohi_d,cffrelohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehilo_d,cffrehilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelolo_d,cffrelolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhihi_d,cffimhihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlohi_d,cffimlohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhilo_d,cffimhilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlolo_d,cffimlolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehihi_d,inputrehihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelohi_d,inputrelohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehilo_d,inputrehilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelolo_d,inputrelolo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhihi_d,inputimhihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlohi_d,inputimlohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhilo_d,inputimhilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlolo_d,inputimlolo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx4_speel
         (BS,nvr,deg,idx,
          cffrehihi_d,cffrelohi_d,cffrehilo_d,cffrelolo_d,
          cffimhihi_d,cffimlohi_d,cffimhilo_d,cffimlolo_d,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
          forwardrehihi_d,forwardrelohi_d,forwardrehilo_d,forwardrelolo_d,
          forwardimhihi_d,forwardimlohi_d,forwardimhilo_d,forwardimlolo_d,
          backwardrehihi_d,backwardrelohi_d,backwardrehilo_d,backwardrelolo_d,
          backwardimhihi_d,backwardimlohi_d,backwardimhilo_d,backwardimlolo_d,
          crossrehihi_d,crossrelohi_d,crossrehilo_d,crossrelolo_d,
          crossimhihi_d,crossimlohi_d,crossimhilo_d,crossimlolo_d);
   }
   double *forwardrehihi_h = new double[(deg1)*nvr];
   double *forwardrelohi_h = new double[(deg1)*nvr];
   double *forwardrehilo_h = new double[(deg1)*nvr];
   double *forwardrelolo_h = new double[(deg1)*nvr];
   double *forwardimhihi_h = new double[(deg1)*nvr];
   double *forwardimlohi_h = new double[(deg1)*nvr];
   double *forwardimhilo_h = new double[(deg1)*nvr];
   double *forwardimlolo_h = new double[(deg1)*nvr];
   double *backwardrehihi_h = new double[(deg1)*(nvr-1)];
   double *backwardrelohi_h = new double[(deg1)*(nvr-1)];
   double *backwardrehilo_h = new double[(deg1)*(nvr-1)];
   double *backwardrelolo_h = new double[(deg1)*(nvr-1)];
   double *backwardimhihi_h = new double[(deg1)*(nvr-1)];
   double *backwardimlohi_h = new double[(deg1)*(nvr-1)];
   double *backwardimhilo_h = new double[(deg1)*(nvr-1)];
   double *backwardimlolo_h = new double[(deg1)*(nvr-1)];
   double *crossrehihi_h = new double[(deg1)*(nvr-2)];
   double *crossrelohi_h = new double[(deg1)*(nvr-2)];
   double *crossrehilo_h = new double[(deg1)*(nvr-2)];
   double *crossrelolo_h = new double[(deg1)*(nvr-2)];
   double *crossimhihi_h = new double[(deg1)*(nvr-2)];
   double *crossimlohi_h = new double[(deg1)*(nvr-2)];
   double *crossimhilo_h = new double[(deg1)*(nvr-2)];
   double *crossimlolo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrehihi_h,forwardrehihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelohi_h,forwardrelohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrehilo_h,forwardrehilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelolo_h,forwardrelolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhihi_h,forwardimhihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlohi_h,forwardimlohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhilo_h,forwardimhilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlolo_h,forwardimlolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehihi_h,backwardrehihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelohi_h,backwardrelohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehilo_h,backwardrehilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelolo_h,backwardrelolo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhihi_h,backwardimhihi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlohi_h,backwardimlohi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhilo_h,backwardimhilo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlolo_h,backwardimlolo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehihi_h,crossrehihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelohi_h,crossrelohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehilo_h,crossrehilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelolo_h,crossrelolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhihi_h,crossimhihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlohi_h,crossimlohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhilo_h,crossimhilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlolo_h,crossimlolo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputrehihi[dim][i] = forwardrehihi_h[offset+i];
      outputrelohi[dim][i] = forwardrelohi_h[offset+i];
      outputrehilo[dim][i] = forwardrehilo_h[offset+i];
      outputrelolo[dim][i] = forwardrelolo_h[offset+i];
      outputimhihi[dim][i] = forwardimhihi_h[offset+i];
      outputimlohi[dim][i] = forwardimlohi_h[offset+i];
      outputimhilo[dim][i] = forwardimhilo_h[offset+i];
      outputimlolo[dim][i] = forwardimlolo_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputrehihi[ix][i] = forwardrehihi_h[offset+i];
      outputrelohi[ix][i] = forwardrelohi_h[offset+i];
      outputrehilo[ix][i] = forwardrehilo_h[offset+i];
      outputrelolo[ix][i] = forwardrelolo_h[offset+i];
      outputimhihi[ix][i] = forwardimhihi_h[offset+i];
      outputimlohi[ix][i] = forwardimlohi_h[offset+i];
      outputimhilo[ix][i] = forwardimhilo_h[offset+i];
      outputimlolo[ix][i] = forwardimlolo_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputrehihi[ix][i] = backwardrehihi_h[offset+i];
      outputrelohi[ix][i] = backwardrelohi_h[offset+i];
      outputrehilo[ix][i] = backwardrehilo_h[offset+i];
      outputrelolo[ix][i] = backwardrelolo_h[offset+i];
      outputimhihi[ix][i] = backwardimhihi_h[offset+i];
      outputimlohi[ix][i] = backwardimlohi_h[offset+i];
      outputimhilo[ix][i] = backwardimhilo_h[offset+i];
      outputimlolo[ix][i] = backwardimlolo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrehihi[ix][i] = crossrehihi_h[offset+i];
         outputrelohi[ix][i] = crossrelohi_h[offset+i];
         outputrehilo[ix][i] = crossrehilo_h[offset+i];
         outputrelolo[ix][i] = crossrelolo_h[offset+i];
         outputimhihi[ix][i] = crossimhihi_h[offset+i];
         outputimlohi[ix][i] = crossimlohi_h[offset+i];
         outputimhilo[ix][i] = crossimhilo_h[offset+i];
         outputimlolo[ix][i] = crossimlolo_h[offset+i];
      }
   }
}
