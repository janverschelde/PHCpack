// Defines code of the functions in cmplx4_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector,
// in quad double precision,
// for vectors of small, medium, and large size.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "cmplx4_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim, int dimLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   int j = threadIdx.x;

   __shared__ double shvrehihi[qd_shmemsize];
   __shared__ double shvrelohi[qd_shmemsize];
   __shared__ double shvrehilo[qd_shmemsize];
   __shared__ double shvrelolo[qd_shmemsize];
   __shared__ double shvimhihi[qd_shmemsize];
   __shared__ double shvimlohi[qd_shmemsize];
   __shared__ double shvimhilo[qd_shmemsize];
   __shared__ double shvimlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];
   __shared__ double sumhihi[qd_shmemsize];
   __shared__ double sumlohi[qd_shmemsize];
   __shared__ double sumhilo[qd_shmemsize];
   __shared__ double sumlolo[qd_shmemsize];

   shvrehihi[j] = vrehihi[j]; // reading real parts into shared memory
   shvrelohi[j] = vrelohi[j];
   shvrehilo[j] = vrehilo[j];
   shvrelolo[j] = vrelolo[j];
   shvimhihi[j] = vimhihi[j]; // reading imaginary parts into shared memory
   shvimlohi[j] = vimlohi[j];
   shvimhilo[j] = vimhilo[j];
   shvimlolo[j] = vimlolo[j];

   qdg_sqr(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
            &sumhihi[j], &sumlohi[j], &sumhilo[j], &sumlolo[j]);
   qdg_sqr(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
            &prdhihi[j], &prdlohi[j], &prdhilo[j], &prdlolo[j]);
   qdg_inc(&sumhihi[j],&sumlohi[j],&sumhilo[j],&sumlolo[j],
            prdhihi[j], prdlohi[j], prdhilo[j], prdlolo[j]);

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            qdg_inc(&sumhihi[j],&sumlohi[j],&sumhilo[j],&sumlolo[j],
                     sumhihi[j+powTwo],sumlohi[j+powTwo],
                     sumhilo[j+powTwo],sumlolo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
      qdg_sqrt( sumhihi[0], sumlohi[0], sumhilo[0], sumlolo[0],
               &sumhihi[0],&sumlohi[0],&sumhilo[0],&sumlolo[0]); 
   if(j == 0)
   {
      *normhihi = sumhihi[0];
      *normlohi = sumlohi[0];
      *normhilo = sumhilo[0];
      *normlolo = sumlolo[0];
   }
   __syncthreads();
   qdg_div(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
             sumhihi[0],  sumlohi[0],  sumhilo[0],  sumlolo[j],
            &vrehihi[j], &vrelohi[j], &vrehilo[j], &vrelolo[j]);
   qdg_div(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
             sumhihi[0],  sumlohi[0],  sumhilo[0],  sumlolo[j],
            &vimhihi[j], &vimlohi[j], &vimhilo[j], &vimlolo[j]);
}

__global__ void medium_normalize_vector
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvrehihi[qd_shmemsize];
   __shared__ double shvrelohi[qd_shmemsize];
   __shared__ double shvrehilo[qd_shmemsize];
   __shared__ double shvrelolo[qd_shmemsize];
   __shared__ double shvimhihi[qd_shmemsize];
   __shared__ double shvimlohi[qd_shmemsize];
   __shared__ double shvimhilo[qd_shmemsize];
   __shared__ double shvimlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];
   __shared__ double acchihi[qd_shmemsize];
   __shared__ double acclohi[qd_shmemsize];
   __shared__ double acchilo[qd_shmemsize];
   __shared__ double acclolo[qd_shmemsize];
   __shared__ double sumshihi[maxrounds];
   __shared__ double sumslohi[maxrounds];
   __shared__ double sumshilo[maxrounds];
   __shared__ double sumslolo[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         acchihi[j] = 0.0;
         acclohi[j] = 0.0;
         acchilo[j] = 0.0;
         acclolo[j] = 0.0;
      }
      else
      {
         shvrehihi[j] = vrehihi[vBSind+j];  // reading into shared memory
         shvrelohi[j] = vrelohi[vBSind+j];
         shvrehilo[j] = vrehilo[vBSind+j];
         shvrelolo[j] = vrelolo[vBSind+j];
         shvimhihi[j] = vimhihi[vBSind+j]; 
         shvimlohi[j] = vimlohi[vBSind+j]; 
         shvimhilo[j] = vimhilo[vBSind+j]; 
         shvimlolo[j] = vimlolo[vBSind+j]; 

         qdg_sqr(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
                  &acchihi[j], &acclohi[j], &acchilo[j], &acclolo[j]);
         qdg_sqr(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
                  &prdhihi[j], &prdlohi[j], &prdhilo[j], &prdlolo[j]);
         qdg_inc(&acchihi[j],&acclohi[j],&acchilo[j],&acclolo[j],
                  prdhihi[j], prdlohi[j], prdhilo[j], prdlolo[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               qdg_inc(&acchihi[j],&acclohi[j],&acchilo[j],&acclolo[j],
                       acchihi[j+powTwo],acclohi[j+powTwo],
                       acchilo[j+powTwo],acclolo[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumshihi[i] = acchihi[0]; 
         sumslohi[i] = acclohi[0]; 
         sumshilo[i] = acchilo[0]; 
         sumslolo[i] = acclolo[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            qdg_inc(&sumshihi[j],&sumslohi[j],&sumshilo[j],&sumslolo[j],
                     sumshihi[j+powTwo],sumslohi[j+powTwo],
                     sumshilo[j+powTwo],sumslolo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
   {
      qdg_sqrt( sumshihi[0], sumslohi[0], sumshilo[0], sumslolo[0],
               &sumshihi[0],&sumslohi[0],&sumshilo[0],&sumslolo[0]);
      *normhihi = sumshihi[0];
      *normlohi = sumslohi[0];
      *normhilo = sumshilo[0];
      *normlolo = sumslolo[0];
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvrehihi[j] = vrehihi[vBSind+j];     // read into shared memory
         shvrelohi[j] = vrelohi[vBSind+j];
         shvrehilo[j] = vrehilo[vBSind+j];
         shvrelolo[j] = vrelolo[vBSind+j];
         shvimhihi[j] = vimhihi[vBSind+j];
         shvimlohi[j] = vimlohi[vBSind+j];
         shvimhilo[j] = vimhilo[vBSind+j];
         shvimlolo[j] = vimlolo[vBSind+j];
         // normalize vector
         qdg_div(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
                  sumshihi[0], sumslohi[0], sumshilo[0], sumslolo[0],
                  &vrehihi[vBSind+j],&vrelohi[vBSind+j],
                  &vrehilo[vBSind+j],&vrelolo[vBSind+j]);
         qdg_div(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
                  sumshihi[0], sumslohi[0], sumshilo[0], sumslolo[0],
                  &vimhihi[vBSind+j],&vimlohi[vBSind+j],
                  &vimhilo[vBSind+j],&vimlolo[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}
__global__ void large_sum_the_squares
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehihi[qd_shmemsize];
   __shared__ double shvrelohi[qd_shmemsize];
   __shared__ double shvrehilo[qd_shmemsize];
   __shared__ double shvrelolo[qd_shmemsize];
   __shared__ double shvimhihi[qd_shmemsize];
   __shared__ double shvimlohi[qd_shmemsize];
   __shared__ double shvimhilo[qd_shmemsize];
   __shared__ double shvimlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];
   __shared__ double acchihi[qd_shmemsize];
   __shared__ double acclohi[qd_shmemsize];
   __shared__ double acchilo[qd_shmemsize];
   __shared__ double acclolo[qd_shmemsize];

   shvrehihi[j] = vrehihi[k];
   shvrelohi[j] = vrelohi[k];
   shvrehilo[j] = vrehilo[k];
   shvrelolo[j] = vrelolo[k];
   shvimhihi[j] = vimhihi[k];
   shvimlohi[j] = vimlohi[k];
   shvimhilo[j] = vimhilo[k];
   shvimlolo[j] = vimlolo[k];

   qdg_sqr(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
            &acchihi[j], &acclohi[j], &acchilo[j], &acclolo[j]);
   qdg_sqr(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
            &prdhihi[j], &prdlohi[j], &prdhilo[j], &prdlolo[j]);
   qdg_inc(&acchihi[j],&acclohi[j],&acchilo[j],&acclolo[j],
            prdhihi[j], prdlohi[j], prdhilo[j], prdlolo[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)
            qdg_inc(&acchihi[j],&acclohi[j],&acchilo[j],&acclolo[j],
                     acchihi[j+powTwo],acclohi[j+powTwo],
                     acchilo[j+powTwo],acclolo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                               // thread 0 writes the sum
   {
      sumshihi[i] = acchihi[0];
      sumslohi[i] = acclohi[0];
      sumshilo[i] = acchilo[0];
      sumslolo[i] = acclolo[0];
   }
}
__global__ void large_normalize_vector
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int nbsums, int nbsumsLog2, int BS,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehihi[qd_shmemsize];
   __shared__ double shvrelohi[qd_shmemsize];
   __shared__ double shvrehilo[qd_shmemsize];
   __shared__ double shvrelolo[qd_shmemsize];
   __shared__ double shvimhihi[qd_shmemsize];
   __shared__ double shvimlohi[qd_shmemsize];
   __shared__ double shvimhilo[qd_shmemsize];
   __shared__ double shvimlolo[qd_shmemsize];

   if(j < nbsums)
   {
      shvrehihi[j] = sumshihi[j];
      shvrelohi[j] = sumslohi[j];
      shvrehilo[j] = sumshilo[j];
      shvrelolo[j] = sumslolo[j];
   }

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            qdg_inc(&shvrehihi[j],&shvrelohi[j],&shvrehilo[j],&shvrelolo[j],
                     shvrehihi[j+powTwo],shvrelohi[j+powTwo],
                     shvrehilo[j+powTwo],shvrelolo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                      // every thread 0 of all blocks
   if(j == 0)                            // compute the 2-norm and assigns
   {                                     // to the output parameter
      qdg_sqrt(shvrehihi[0],shvrelohi[0],shvrehilo[0],shvrelolo[0],
                normhihi,    normlohi,    normhilo,    normlolo);
   }
   __syncthreads();                    

   if(k < dim)
   {
      shvrehihi[j] = vrehihi[k];
      shvrelohi[j] = vrelohi[k];
      shvrehilo[j] = vrehilo[k];
      shvrelolo[j] = vrelolo[k];
      shvimhihi[j] = vimhihi[k];
      shvimlohi[j] = vimlohi[k];
      shvimhilo[j] = vimhilo[k];
      shvimlolo[j] = vimlolo[k];

      qdg_div( shvrehihi[j], shvrelohi[j], shvrehilo[j], shvrelolo[j],
              *normhihi,*normlohi,*normhilo,*normlolo,
              &shvrehihi[j],&shvrelohi[j],&shvrehilo[j],&shvrelolo[j]);
      qdg_div( shvimhihi[j], shvimlohi[j], shvimhilo[j], shvimlolo[j],
              *normhihi,*normlohi,*normhilo,*normlolo,
              &shvimhihi[j],&shvimlohi[j],&shvimhilo[j],&shvimlolo[j]);

      vrehihi[k] = shvrehihi[j];
      vrelohi[k] = shvrelohi[j];
      vrehilo[k] = shvrehilo[j];
      vrelolo[k] = shvrelolo[j];
      vimhihi[k] = shvimhihi[j];
      vimlohi[k] = shvimlohi[j];
      vimhilo[k] = shvimhilo[j];
      vimlolo[k] = shvimlolo[j];
   }
}

void GPU_norm
 ( double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   int dim, int freq, int BS,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo,
   int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vrehihi_d;                 // highest real parts on device
   double* vrelohi_d;                 // 2nd highest real parts on device
   double* vrehilo_d;                 // 2nd lowest real parts on device
   double* vrelolo_d;                 // lowest real parts on device
   double* vimhihi_d;                 // highest imaginary parts on device
   double* vimlohi_d;                 // 2nd highest imaginary parts on device
   double* vimhilo_d;                 // 2nd lowest imaginary parts on device
   double* vimlolo_d;                 // lowest imaginary parts on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vrehihi_d,size);
   cudaMalloc((void**)&vrelohi_d,size);
   cudaMalloc((void**)&vrehilo_d,size);
   cudaMalloc((void**)&vrelolo_d,size);
   cudaMalloc((void**)&vimhihi_d,size);
   cudaMalloc((void**)&vimlohi_d,size);
   cudaMalloc((void**)&vimhilo_d,size);
   cudaMalloc((void**)&vimlolo_d,size);
   cudaMemcpy(vrehihi_d,vrehihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelohi_d,vrelohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrehilo_d,vrehilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelolo_d,vrelolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhihi_d,vimhihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlohi_d,vimlohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhilo_d,vimhilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlolo_d,vimlolo_h,size,cudaMemcpyHostToDevice);
   double* normhihi_d;
   double* normlohi_d;
   double* normhilo_d;
   double* normlolo_d;
   cudaMalloc((void**)&normhihi_d,sizeof(double));
   cudaMalloc((void**)&normlohi_d,sizeof(double));
   cudaMalloc((void**)&normhilo_d,sizeof(double));
   cudaMalloc((void**)&normlolo_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vrehihi_d,vrelohi_d,vrehilo_d,vrelolo_d,
             vimhihi_d,vimlohi_d,vimhilo_d,vimlolo_d,
             dim,BSLog2,normhihi_d,normlohi_d,normhilo_d,normlolo_d);
   }
   else if(blocked == 0)
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vrehihi_d,vrelohi_d,vrehilo_d,vrelolo_d,
             vimhihi_d,vimlohi_d,vimhilo_d,vimlolo_d,dim,
             rf,rfLog2,BS,BSLog2,normhihi_d,normlohi_d,normhilo_d,normlolo_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumshihi_d; // highest parts of sums of squares for each block
      double* sumslohi_d; // 2nd highest parts of sums of squares
      double* sumshilo_d; // 2nd lowest parts of sums of squares for each block
      double* sumslolo_d; // lowest parts of sums of squares for each block
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumshihi_d,sums_size);
      cudaMalloc((void**)&sumslohi_d,sums_size);
      cudaMalloc((void**)&sumshilo_d,sums_size);
      cudaMalloc((void**)&sumslolo_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vrehihi_d,vrelohi_d,vrehilo_d,vrelolo_d,
             vimhihi_d,vimlohi_d,vimhilo_d,vimlolo_d,dim,
             sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vrehihi_d,vrelohi_d,vrehilo_d,vrelolo_d,
             vimhihi_d,vimlohi_d,vimhilo_d,vimlolo_d,dim,
             sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,
             nblocks,nblocksLog2,BS,
             normhihi_d,normlohi_d,normhilo_d,normlolo_d);
      }
   }
   cudaMemcpy(vrehihi_h,vrehihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelohi_h,vrelohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrehilo_h,vrehilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelolo_h,vrelolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhihi_h,vimhihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlohi_h,vimlohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhilo_h,vimhilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlolo_h,vimlolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhihi,normhihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlohi,normlohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhilo,normhilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlolo,normlolo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
