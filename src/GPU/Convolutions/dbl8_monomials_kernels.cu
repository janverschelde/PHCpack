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
   dbl8_convolute<<<1,BS>>>
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
      dbl8_convolute<<<1,BS>>>
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
      dbl8_convolute<<<1,BS>>>
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
         dbl8_convolute<<<1,BS>>>
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
      dbl8_convolute<<<1,BS>>>
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
         dbl8_convolute<<<1,BS>>>
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
            dbl8_convolute<<<1,BS>>>
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
         dbl8_convolute<<<1,BS>>>
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
}
