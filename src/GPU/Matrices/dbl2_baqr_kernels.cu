/* The file dbl2_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl2_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_baqr_kernels.h"
#include "double_double_functions.h"
#include "dbl_baqr_flopcounts.h"

using namespace std;

__global__ void dbl2_small_house
 ( double *x0hi, double *x0lo, double *x1hi, double *x1lo,
   int dim, int dimLog2,
   double *vhi, double *vlo, double *betahi, double *betalo )
{
   const int j = threadIdx.x;

   __shared__ double shvhi[dd_shmemsize];
   __shared__ double shvlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];

   bool stopflag = false;
   double acchi,acclo,muhi,mulo,v0hi,v0lo,v0p2hi,v0p2lo;

   shvhi[j] = x1hi[j];          // reading of vector into shared memory
   shvlo[j] = x1lo[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   ddg_sqr(shvhi[j],shvlo[j],&prdhi[j],&prdlo[j]);

   vhi[j+1] = shvhi[j];         // copies x to v, in case beta is zero
   vlo[j+1] = shvlo[j];
   if(j == 0) vhi[0] = 1.0;
   if(j == 0) vlo[0] = 0.0;

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) // prd[j] = prd[j] + prd[j+powTwo];
            ddg_inc(&prdhi[j],&prdlo[j],prdhi[j+powTwo],prdlo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {
      if((prdhi[0] == 0.0) && (prdlo[0] == 0.0))  // prd[0] is sigma of house
      {
         *betahi = 0.0; *betalo = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   if(j == 0)                              // thread zero sets beta
   {
      // mu = sqrt((*x0)*(*x0) + prd[0]);
      ddg_sqr(*x0hi,*x0lo,&acchi,&acclo);
      ddg_inc(&acchi,&acclo,prdhi[0],prdlo[0]);
      ddg_sqrt(acchi,acclo,&muhi,&mulo);
      if(*x0hi <= 0.0)
      {
         // v0 = *x0 - mu;
         ddg_sub(*x0hi,*x0lo,muhi,mulo,&v0hi,&v0lo);
      }
      else
      {
         // v0 = -prd[0]/(*x0 + mu);
         ddg_add(*x0hi,*x0lo,muhi,mulo,&acchi,&acclo);
         ddg_div(prdhi[0],prdlo[0],acchi,acclo,&v0hi,&v0lo);
         ddg_minus(&v0hi,&v0lo);
      }
      // v0p2 = v0*v0;
      ddg_sqr(v0hi,v0lo,&v0p2hi,&v0p2lo);
      // *beta = 2.0*v0p2/(prd[0] + v0p2);
      ddg_add(prdhi[0],prdlo[0],v0p2hi,v0p2lo,&acchi,&acclo);
      ddg_div(v0p2hi,v0p2lo,acchi,acclo,betahi,betalo);
      ddg_mlt_d(betahi,betalo,2.0);
      prdhi[0] = v0hi;
      prdlo[0] = v0lo;                     // v0 needed for normalization
   }
   __syncthreads();
   // shv[j] = shv[j]/prd[0];
   vhi[j+1] = 0.0;
   vlo[j+1] = 0.0;
   if(1.0 + *betahi + *betalo != 1.0)
   {
      ddg_div(shvhi[j],shvlo[j],prdhi[0],prdlo[0],&acchi,&acclo);
      vhi[j+1] = acchi;
      vlo[j+1] = acclo;
   }
   if(j == 0) vhi[0] = 1.0;
   if(j == 0) vlo[0] = 0.0;
}

__global__ void cmplx2_small_house
 ( double *x0rehi, double *x0relo, double *x0imhi, double *x0imlo,
   double *x1rehi, double *x1relo, double *x1imhi, double *x1imlo,
   int dim, int dimLog2,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double *betahi, double *betalo )
{
   const int j = threadIdx.x;

   __shared__ double shvrehi[cdd_shmemsize];
   __shared__ double shvrelo[cdd_shmemsize];
   __shared__ double shvimhi[cdd_shmemsize];
   __shared__ double shvimlo[cdd_shmemsize];
   __shared__ double prdhi[cdd_shmemsize];
   __shared__ double prdlo[cdd_shmemsize];
   __shared__ double v0parts[4];

   bool stopflag = false;
   double muhi,mulo,v0rehi,v0relo,v0imhi,v0imlo;
   double x0radhi,x0radlo,sqrx0hi,sqrx0lo,sqrv0hi,sqrv0lo;
   double inv0rehi,inv0relo,inv0imhi,inv0imlo;
   double zrehi,zrelo,zimhi,zimlo,acchi,acclo;

   shvrehi[j] = x1rehi[j];      // reading of vector into shared memory
   shvrelo[j] = x1relo[j];
   shvimhi[j] = x1imhi[j];
   shvimlo[j] = x1imlo[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   // prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];
   ddg_sqr(shvrehi[j],shvrelo[j],&prdhi[j],&prdlo[j]);
   ddg_sqr(shvimhi[j],shvimlo[j],&acchi,&acclo);
   ddg_inc(&prdhi[j],&prdlo[j],acchi,acclo);

   vrehi[j+1] = shvrehi[j];     // copies x to v, in case beta is zero
   vrelo[j+1] = shvrelo[j];
   vimhi[j+1] = shvimhi[j];
   vimlo[j+1] = shvimlo[j];
   if(j == 0)
   {
      vrehi[0] = 1.0;
      vrelo[0] = 0.0;
      vimhi[0] = 0.0;
      vimlo[0] = 0.0;
   }
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)     // prd[j] = prd[j] + prd[j+powTwo];
            ddg_inc(&prdhi[j],&prdlo[j],prdhi[j+powTwo],prdlo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {
      if((prdhi[0] == 0.0) && (prdlo[0] == 0.0))  // prd[0] is sigma of house
      {
         *betahi = 0.0; *betalo = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   if(j == 0)                              // thread zero sets beta
   {
      // sqrx0 = (*x0re)*(*x0re) + (*x0im)*(*x0im);
      ddg_sqr(*x0rehi,*x0relo,&sqrx0hi,&sqrx0lo);
      ddg_sqr(*x0imhi,*x0imlo,&acchi,&acclo);
      ddg_inc(&sqrx0hi,&sqrx0lo,acchi,acclo);
      // x0rad = sqrt(sqrx0);
      ddg_sqrt(sqrx0hi,sqrx0lo,&x0radhi,&x0radlo);
      // mu = sqrt(sqrx0 + prd[0]);
      ddg_inc(&sqrx0hi,&sqrx0lo,prdhi[0],prdlo[0]);
      ddg_sqrt(sqrx0hi,sqrx0lo,&muhi,&mulo);

      if((x0radhi == 0.0) && (x0radlo == 0.0))
      {
         v0rehi = muhi;
         v0relo = mulo;
         ddg_minus(&v0rehi,&v0relo);
         v0imhi = 0.0;
         v0imlo = 0.0;
      }
      else
      {
         // mu = mu/x0rad;
         ddg_div(muhi,mulo,x0radhi,x0radlo,&acchi,&acclo);
         muhi = acchi;
         mulo = acclo;
         // v0re = (*x0re) - mu*(*x0re);
         ddg_mul(muhi,mulo,x0rehi[0],x0relo[0],&acchi,&acclo);
         ddg_sub(x0rehi[0],x0relo[0],acchi,acclo,&v0rehi,&v0relo);
         // v0im = (*x0im) - mu*(*x0im);
         ddg_mul(muhi,mulo,x0imhi[0],x0imlo[0],&acchi,&acclo);
         ddg_sub(x0imhi[0],x0imlo[0],acchi,acclo,&v0imhi,&v0imlo);
      }
      // sqrv0 = v0re*v0re + v0im*v0im;
      ddg_sqr(v0rehi,v0relo,&sqrv0hi,&sqrv0lo);
      ddg_sqr(v0imhi,v0imlo,&acchi,&acclo);
      ddg_inc(&sqrv0hi,&sqrv0lo,acchi,acclo);
      // *beta = 2.0*sqrv0/(prd[0] + sqrv0);
      ddg_add(prdhi[0],prdlo[0],sqrv0hi,sqrv0lo,&acchi,&acclo);
      ddg_div(sqrv0hi,sqrv0lo,acchi,acclo,betahi,betalo);
      ddg_mlt_d(betahi,betalo,2.0);

      prdhi[0] = sqrv0hi;                 // sqrv0 needed for normalization
      prdlo[0] = sqrv0lo;
      v0parts[0] = v0rehi;                // share v0rehi with all threads
      v0parts[1] = v0relo;                // share v0relo with all threads
      v0parts[2] = v0imhi;                // share v0imhi with all threads
      v0parts[3] = v0imlo;                // share v0imlo with all threads
   }
   __syncthreads(); // important synchronization!
   // inv0re = v0parts[0]/prd[0];               // real part of 1/v[0]
   vrehi[j+1] = 0.0;
   vrelo[j+1] = 0.0; // assign just in case of a zero beta ...
   vimhi[j+1] = 0.0;
   vimlo[j+1] = 0.0;
   if(1.0 + prdhi[0] + prdlo[0] != 1.0)
   {
      ddg_div(v0parts[0],v0parts[1],prdhi[0],prdlo[0],&inv0rehi,&inv0relo);
      // inv0im = -v0parts[1]/prd[0];              // imag part of 1/v[0]
      ddg_div(v0parts[2],v0parts[3],prdhi[0],prdlo[0],&inv0imhi,&inv0imlo);
      ddg_minus(&inv0imhi,&inv0imlo);
      // zre = shvre[j]*inv0re - shvim[j]*inv0im;  // real part of v[j]/v[0]
      ddg_mul(shvrehi[j],shvrelo[j],inv0rehi,inv0relo,&zrehi,&zrelo);
      ddg_mul(shvimhi[j],shvimlo[j],inv0imhi,inv0imlo,&acchi,&acclo);
      ddg_dec(&zrehi,&zrelo,acchi,acclo);
      // zim = shvim[j]*inv0re + shvre[j]*inv0im;  // imag part of v[j]/v[0]
      ddg_mul(shvimhi[j],shvimlo[j],inv0rehi,inv0relo,&zimhi,&zimlo);
      ddg_mul(shvrehi[j],shvrelo[j],inv0imhi,inv0imlo,&acchi,&acclo);
      ddg_inc(&zimhi,&zimlo,acchi,acclo);
      vrehi[j+1] = zrehi;
      vrelo[j+1] = zrelo;
      vimhi[j+1] = zimhi;
      vimlo[j+1] = zimlo;
   }
   if(j == 0)
   {
      vrehi[0] = 1.0;
      vrelo[0] = 0.0;
      vimhi[0] = 0.0;
      vimlo[0] = 0.0;
   }
}

__global__ void dbl2_large_sum_of_squares
 ( double *vhi, double *vlo, double *sumshi, double *sumslo,
   int dim, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhi[inner_dd_shmemsize];
   __shared__ double shvlo[inner_dd_shmemsize];
   __shared__ double prdhi[inner_dd_shmemsize];
   __shared__ double prdlo[inner_dd_shmemsize];

   shvhi[j] = vhi[k];
   shvlo[j] = vlo[k];
   if(k >= dim)
   {
      shvhi[j] = 0.0;
      shvlo[j] = 0.0;
   }
   ddg_sqr(shvhi[j],shvlo[j],&prdhi[j],&prdlo[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS) 
            ddg_inc(&prdhi[j],&prdlo[j],prdhi[j+powTwo],prdlo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                              // thread 0 writes the sum
   {
      sumshi[i] = prdhi[0];
      sumslo[i] = prdlo[0];
   }
}

__global__ void cmplx2_large_sum_of_squares
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double *sumshi, double *sumslo, int dim, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehi[inner_dd_shmemsize];
   __shared__ double shvrelo[inner_dd_shmemsize];
   __shared__ double shvimhi[inner_dd_shmemsize];
   __shared__ double shvimlo[inner_dd_shmemsize];
   __shared__ double prdhi[inner_dd_shmemsize];
   __shared__ double prdlo[inner_dd_shmemsize];

   shvrehi[j] = vrehi[k];
   shvrelo[j] = vrelo[k];
   shvimhi[j] = vimhi[k];
   shvimlo[j] = vimlo[k];

   if(k >= dim)
   {
      shvrehi[j] = 0.0;
      shvrelo[j] = 0.0;
      shvimhi[j] = 0.0;
      shvimlo[j] = 0.0;
   }
   double acchi,acclo;

   ddg_sqr(shvrehi[j],shvrelo[j],&prdhi[j],&prdlo[j]);
   ddg_sqr(shvimhi[j],shvimlo[j],&acchi,&acclo);
   ddg_inc(&prdhi[j],&prdlo[j],acchi,acclo);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int k=0; k < BSLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)     // prd[j] = prd[j] + prd[j+powTwo];
            ddg_inc(&prdhi[j],&prdlo[j],prdhi[j+powTwo],prdlo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)                              // thread 0 writes the sum
   {
      sumshi[i] = prdhi[0];
      sumslo[i] = prdlo[0];
   }
}

__global__ void dbl2_sum_accumulator
 ( double *sumshi, double *sumslo, int nbsums, int nbsumsLog2,
   double *acchi, double *acclo )
{
   const int j = threadIdx.x;

   __shared__ double shvhi[outer_dd_shmemsize];
   __shared__ double shvlo[outer_dd_shmemsize];

   shvhi[j] = sumshi[j];
   shvlo[j] = sumslo[j];

   if(j >= nbsums)
   {
      shvhi[j] = 0.0;
      shvlo[j] = 0.0;
   }
   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            ddg_inc(&shvhi[j],&shvlo[j],shvhi[j+powTwo],shvlo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();
   if(j == 0)
   {
      *acchi = shvhi[0];
      *acclo = shvlo[0];
   }
}

__global__ void dbl2_normalize
 ( int dim, int szt, double *xhi, double *xlo, double *v0hi, double *v0lo,
   double *vhi, double *vlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvhi[inner_dd_shmemsize];
   __shared__ double shvlo[inner_dd_shmemsize];

   shvhi[tdx] = xhi[idx];
   shvlo[tdx] = xlo[idx];
   __syncthreads();

   double resulthi,resultlo;

   // shv[j] = shv[j]/v0;
   ddg_div(shvhi[tdx],shvlo[tdx],v0hi[0],v0lo[0],&resulthi,&resultlo);

   __syncthreads();
   if(idx < dim)
   {
      vhi[idx] = resulthi;    
      vlo[idx] = resultlo;    
   }
}

__global__ void cmplx2_normalize
 ( int dim, int szt,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *inv0rehi, double *inv0relo, double *inv0imhi, double *inv0imlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvrehi[inner_dd_shmemsize];
   __shared__ double shvrelo[inner_dd_shmemsize];
   __shared__ double shvimhi[inner_dd_shmemsize];
   __shared__ double shvimlo[inner_dd_shmemsize];

   shvrehi[tdx] = xrehi[idx];
   shvrelo[tdx] = xrelo[idx];
   shvimhi[tdx] = ximhi[idx];
   shvimlo[tdx] = ximlo[idx];
   __syncthreads();

   double resultrehi,resultrelo,resultimhi,resultimlo;
   double acchi,acclo;

   // shv[j] = shv[j]/v0;

   // resultre = vre[i]*inv0re - vim[i]*inv0im;
   ddg_mul(shvrehi[tdx],shvrelo[tdx],*inv0rehi,*inv0relo,
           &resultrehi,&resultrelo);
   ddg_mul(shvimhi[tdx],shvimlo[tdx],*inv0imhi,*inv0imlo,&acchi,&acclo);
   ddg_dec(&resultrehi,&resultrelo,acchi,acclo);
   // zim = vim[i]*inv0re + vre[i]*inv0im;
   ddg_mul(shvimhi[tdx],shvimlo[tdx],*inv0rehi,*inv0relo,
           &resultimhi,&resultimlo);
   ddg_mul(shvrehi[tdx],shvrelo[tdx],*inv0imhi,*inv0imlo,&acchi,&acclo);
   ddg_inc(&resultimhi,&resultimlo,acchi,acclo);

   __syncthreads();
   if(idx < dim)
   {
      vrehi[idx] = resultrehi;    
      vrelo[idx] = resultrelo;    
      vimhi[idx] = resultimhi;    
      vimlo[idx] = resultimlo;    
   }
}

__global__ void dbl2_small_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *Rhi, double *Rlo,
   double *vhi, double *vlo, double *betahi, double *betalo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   int Rcolidx;
   double whi,wlo,Rtdxhi,Rtdxlo,acchi,acclo;

   __shared__ double shvhi[dd_shmemsize]; // slice of v
   __shared__ double shvlo[dd_shmemsize]; 

   shvhi[tdx] = vhi[tdx];
   shvlo[tdx] = vlo[tdx];
   __syncthreads();
   whi = 0.0;
   wlo = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdxhi = Rhi[Rcolidx];
      Rtdxlo = Rlo[Rcolidx];
      // w = w + Rtdx*shv[i];
      ddg_mul(Rtdxhi,Rtdxlo,shvhi[i],shvlo[i],&acchi,&acclo);
      ddg_inc(&whi,&wlo,acchi,acclo);
   }
   // w = (*beta)*w;
   // ddg_mlt(&whi,&wlo,*betahi,*betalo); <-- this does not work!
   ddg_mul(*betahi,*betalo,whi,wlo,&acchi,&acclo);
   whi = acchi;
   wlo = acclo;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdxhi = Rhi[Rcolidx];
      Rtdxlo = Rlo[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w;
      ddg_mul(shvhi[i],shvlo[i],whi,wlo,&acchi,&acclo);
      ddg_dec(&Rtdxhi,&Rtdxlo,acchi,acclo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rhi[Rcolidx] = Rtdxhi;
         Rlo[Rcolidx] = Rtdxlo;
      }
      __syncthreads();
   }
}

__global__ void cmplx2_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehi, double *Rrelo, double *Rimhi, double *Rimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double *betahi, double *betalo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   int Rcolidx;
   double w_rehi,w_relo,w_imhi,w_imlo;
   double Rtdx_rehi,Rtdx_relo,Rtdx_imhi,Rtdx_imlo;
   double acchi,acclo,bthi,btlo;

   __shared__ double shvrehi[cdd_shmemsize]; // slice of vrehi
   __shared__ double shvrelo[cdd_shmemsize]; // slice of vrelo
   __shared__ double shvimhi[cdd_shmemsize]; // slice of vimhi
   __shared__ double shvimlo[cdd_shmemsize]; // slice of vimlo

   shvrehi[tdx] = vrehi[tdx];
   shvrelo[tdx] = vrelo[tdx];
   shvimhi[tdx] = vimhi[tdx];
   shvimlo[tdx] = vimlo[tdx];
   __syncthreads();
   w_rehi = 0.0;
   w_relo = 0.0;
   w_imhi = 0.0;
   w_imlo = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx_rehi = Rrehi[Rcolidx];
      Rtdx_relo = Rrelo[Rcolidx];
      Rtdx_imhi = Rimhi[Rcolidx];
      Rtdx_imlo = Rimlo[Rcolidx];
      // w = w + Rtdx*shv[i]; beware of the Hermitian transpose!
      // w_re = w_re + Rtdx_re*shvre[i] + Rtdx_im*shvim[i];
      ddg_mul(Rtdx_rehi,Rtdx_relo,shvrehi[i],shvrelo[i],&acchi,&acclo);
      ddg_inc(&w_rehi,&w_relo,acchi,acclo);
      ddg_mul(Rtdx_imhi,Rtdx_imlo,shvimhi[i],shvimlo[i],&acchi,&acclo);
      ddg_inc(&w_rehi,&w_relo,acchi,acclo);
      // w_im = w_im - Rtdx_im*shvre[i] + Rtdx_re*shvim[i];
      ddg_mul(Rtdx_imhi,Rtdx_imlo,shvrehi[i],shvrelo[i],&acchi,&acclo);
      ddg_dec(&w_imhi,&w_imlo,acchi,acclo);
      ddg_mul(Rtdx_rehi,Rtdx_relo,shvimhi[i],shvimlo[i],&acchi,&acclo);
      ddg_inc(&w_imhi,&w_imlo,acchi,acclo);
   }
   bthi = *betahi;
   btlo = *betalo;
   // w_re = acc*w_re;
   ddg_mul(w_rehi,w_relo,bthi,btlo,&acchi,&acclo);
   w_rehi = acchi;
   w_relo = acclo;
   // w_im = acc*w_im;
   ddg_mul(w_imhi,w_imlo,bthi,btlo,&acchi,&acclo);
   w_imhi = acchi;
   w_imlo = acclo;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx_rehi = Rrehi[Rcolidx];
      Rtdx_relo = Rrelo[Rcolidx];
      Rtdx_imhi = Rimhi[Rcolidx];
      Rtdx_imlo = Rimlo[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w; beware of the Hermitian transpose!
      // Rtdx_re = Rtdx_re - (shvre[i]*w_re + shvim[i]*w_im);
      ddg_mul(shvrehi[i],shvrelo[i],w_rehi,w_relo,&acchi,&acclo);
      ddg_dec(&Rtdx_rehi,&Rtdx_relo,acchi,acclo);
      ddg_mul(shvimhi[i],shvimlo[i],w_imhi,w_imlo,&acchi,&acclo);
      ddg_dec(&Rtdx_rehi,&Rtdx_relo,acchi,acclo);
      // Rtdx_im = Rtdx_im - (shvim[i]*w_re - shvre[i]*w_im);
      ddg_mul(shvimhi[i],shvimlo[i],w_rehi,w_relo,&acchi,&acclo);
      ddg_dec(&Rtdx_imhi,&Rtdx_imlo,acchi,acclo);
      ddg_mul(shvrehi[i],shvrelo[i],w_imhi,w_imlo,&acchi,&acclo);
      ddg_inc(&Rtdx_imhi,&Rtdx_imlo,acchi,acclo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rrehi[Rcolidx] = Rtdx_rehi;
         Rrelo[Rcolidx] = Rtdx_relo;
         Rimhi[Rcolidx] = Rtdx_imhi;
         Rimlo[Rcolidx] = Rtdx_imlo;
      }
   }
}

__global__ void dbl2_small_betaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhi, double *Rlo, double *vhi, double *vlo,
   double *betahi, double *betalo, double *whi, double *wlo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   double resulthi = 0.0;
   double resultlo = 0.0;
   const double mybetahi = *betahi;
   const double mybetalo = *betalo;
   double Rtdxhi,Rtdxlo,acchi,acclo;
   int Rcolidx;

   __shared__ double shvhi[dd_shmemsize]; // slice of v
   __shared__ double shvlo[dd_shmemsize]; // low doubles of v

   shvhi[tdx] = vhi[tdx];
   shvlo[tdx] = vlo[tdx];
   __syncthreads();

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdxhi = Rhi[Rcolidx];
      Rtdxlo = Rlo[Rcolidx];
      // result = result + Rtdx*shv[i];
      ddg_mul(Rtdxhi,Rtdxlo,shvhi[i],shvlo[i],&acchi,&acclo);
      ddg_inc(&resulthi,&resultlo,acchi,acclo);
   }
   // result = mybeta*result;
   ddg_mul(mybetahi,mybetalo,resulthi,resultlo,&acchi,&acclo);
   whi[tdx] = acchi;
   wlo[tdx] = acclo;
}

__global__ void dbl2_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhi, double *Rlo, double *vhi, double *vlo,
   double *RTdotvhi, double *RTdotvlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   const double Vvalhi = vhi[vdx];
   const double Vvallo = vlo[vdx];
   const double Rvalhi = Rhi[Rdx];
   const double Rvallo = Rlo[Rdx];
   // double result = Rval*Vval;
   double resulthi,resultlo;

   ddg_mul(Rvalhi,Rvallo,Vvalhi,Vvallo,&resulthi,&resultlo);

   RTdotvhi[idx] = resulthi;
   RTdotvlo[idx] = resultlo;
}

__global__ void cmplx2_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehi, double *Rrelo, double *Rimhi, double *Rimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double *RHdotvrehi, double *RHdotvrelo,
   double *RHdotvimhi, double *RHdotvimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   const double Vvalrehi = vrehi[vdx];
   const double Vvalrelo = vrelo[vdx];
   const double Vvalimhi = vimhi[vdx];
   const double Vvalimlo = vimlo[vdx];
   const double Rvalrehi = Rrehi[Rdx];
   const double Rvalrelo = Rrelo[Rdx];
   const double Rvalimhi = Rimhi[Rdx];
   const double Rvalimlo = Rimlo[Rdx];
   // double result = Rval*Vval;
   double resultrehi,resultrelo,resultimhi,resultimlo,acchi,acclo;

   ddg_mul(Rvalrehi,Rvalrelo,Vvalrehi,Vvalrelo,&resultrehi,&resultrelo);
   ddg_mul(Rvalimhi,Rvalimlo,Vvalimhi,Vvalimlo,&acchi,&acclo);
   ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
   ddg_mul(Rvalrehi,Rvalrelo,Vvalimhi,Vvalimlo,&resultimhi,&resultimlo);
   ddg_mul(Rvalimhi,Rvalimlo,Vvalrehi,Vvalrelo,&acchi,&acclo);
   ddg_dec(&resultimhi,&resultimlo,acchi,acclo);

   RHdotvrehi[idx] = resultrehi;
   RHdotvrelo[idx] = resultrelo;
   RHdotvimhi[idx] = resultimhi;
   RHdotvimlo[idx] = resultimlo;
}

__global__ void dbl2_sum_betaRTdotv
 ( int nrows, double *betahi, double *betalo,
   double *RTdotvhi, double *RTdotvlo, double *whi, double *wlo )
{
   const int tdx = threadIdx.x;  // tdx sums elements on row tdx
   const int offset = tdx*nrows; // number of rows before current row
   int idx;

   double resulthi = 0.0;
   double resultlo = 0.0;
   double Rvalhi,Rvallo;

   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      Rvalhi = RTdotvhi[idx];
      Rvallo = RTdotvlo[idx];
      // result = result + Rval;
      ddg_inc(&resulthi,&resultlo,Rvalhi,Rvallo);
   }
   Rvalhi = *betahi;
   Rvallo = *betalo;
   // w[tdx] = Rval*result;
   ddg_mul(Rvalhi,Rvallo,resulthi,resultlo,&whi[tdx],&wlo[tdx]);
}

__global__ void cmplx2_sum_betaRHdotv
 ( int nrows, double *betahi, double *betalo,
   double *RTdotvrehi, double *RTdotvrelo,
   double *RTdotvimhi, double *RTdotvimlo,
   double *wrehi, double *wrelo, double *wimhi, double *wimlo )
{
   const int tdx = threadIdx.x;  // tdx sums elements on row tdx
   const int offset = tdx*nrows; // number of rows before current row
   int idx;

   double resultrehi = 0.0;
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;
   double Rvalrehi,Rvalrelo,Rvalimhi,Rvalimlo;

   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      Rvalrehi = RTdotvrehi[idx];
      Rvalrelo = RTdotvrelo[idx];
      Rvalimhi = RTdotvimhi[idx];
      Rvalimlo = RTdotvimlo[idx];
      // result = result + Rval;
      ddg_inc(&resultrehi,&resultrelo,Rvalrehi,Rvalrelo);
      ddg_inc(&resultimhi,&resultimlo,Rvalimhi,Rvalimlo);
   }
   Rvalrehi = *betahi;
   Rvalrelo = *betalo;
   // w[tdx] = Rval*result;
   ddg_mul(Rvalrehi,Rvalrelo,resultrehi,resultrelo,
           &wrehi[tdx],&wrelo[tdx]);
   ddg_mul(Rvalrehi,Rvalrelo,resultimhi,resultimlo,
           &wimhi[tdx],&wimlo[tdx]);
}

__global__ void dbl2_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhi, double *Rlo, double *vhi, double *vlo,
   double *betahi, double *betalo, double *whi, double *wlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwhi[dd_shmemsize];  // values in beta*R^T*v
   __shared__ double shwlo[dd_shmemsize];  // are less in number than szt
   shwhi[tdx] = whi[tdx];
   shwlo[tdx] = wlo[tdx];
   __syncthreads();

   double Rwidxhi = Rhi[Ridx];         // number that tdx updates
   double Rwidxlo = Rlo[Ridx];
   double vValhi = vhi[rowidx];        // value in Householder vector
   double vVallo = vlo[rowidx];
   double wValhi = shwhi[colidx];      // value in beta*R^T*v
   double wVallo = shwlo[colidx];
   double acchi,acclo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   ddg_mul(vValhi,vVallo,wValhi,wVallo,&acchi,&acclo);
   ddg_dec(&Rwidxhi,&Rwidxlo,acchi,acclo);

   if(widx < bound)                    // if() takes care of padding
   {
      Rhi[Ridx] = Rwidxhi;
      Rlo[Ridx] = Rwidxlo;
   }
}

__global__ void cmplx2_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehi, double *Rrelo, double *Rimhi, double *Rimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double *betahi, double *betalo,
   double *wrehi, double *wrelo, double *wimhi, double *wimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwrehi[cdd_shmemsize];  // values in beta*R^H*v
   __shared__ double shwrelo[cdd_shmemsize];  // are less in number than szt
   __shared__ double shwimhi[cdd_shmemsize]; 
   __shared__ double shwimlo[cdd_shmemsize];
   shwrehi[tdx] = wrehi[tdx];
   shwrelo[tdx] = wrelo[tdx];
   shwimhi[tdx] = wimhi[tdx];
   shwimlo[tdx] = wimlo[tdx];
   __syncthreads();

   double Rwidxrehi = Rrehi[Ridx];         // number that tdx updates
   double Rwidxrelo = Rrelo[Ridx];
   double Rwidximhi = Rimhi[Ridx];
   double Rwidximlo = Rimlo[Ridx];
   double vValrehi = vrehi[rowidx];        // value in Householder vector
   double vValrelo = vrelo[rowidx];
   double vValimhi = vimhi[rowidx];
   double vValimlo = vimlo[rowidx];
   double wValrehi = shwrehi[colidx];      // value in beta*R^H*v
   double wValrelo = shwrelo[colidx];
   double wValimhi = shwimhi[colidx];
   double wValimlo = shwimlo[colidx];
   double acchi,acclo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w
   ddg_mul(vValrehi,vValrelo,wValrehi,wValrelo,&acchi,&acclo);
   ddg_dec(&Rwidxrehi,&Rwidxrelo,acchi,acclo);
   ddg_mul(vValimhi,vValimlo,wValimhi,wValimlo,&acchi,&acclo);
   ddg_dec(&Rwidxrehi,&Rwidxrelo,acchi,acclo);
   ddg_mul(vValimhi,vValimlo,wValrehi,wValrelo,&acchi,&acclo);
   ddg_dec(&Rwidximhi,&Rwidximlo,acchi,acclo);
   ddg_mul(vValrehi,vValrelo,wValimhi,wValimlo,&acchi,&acclo);
   ddg_inc(&Rwidximhi,&Rwidximlo,acchi,acclo);

   if(widx < bound)                    // if() takes care of padding
   {
      Rrehi[Ridx] = Rwidxrehi;
      Rrelo[Ridx] = Rwidxrelo;
      Rimhi[Ridx] = Rwidximhi;
      Rimlo[Ridx] = Rwidximlo;
   }
}

__global__ void dbl2_beta_times_V
 ( int nrows, int szt, double *Bhi, double *Blo,
   double *Vhi, double *Vlo, double *Whi, double *Wlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // thread tdx computes W[idx]
   double resulthi,resultlo;

   __shared__ double shvhi[dd_shmemsize]; // to store a slice of V
   __shared__ double shvlo[dd_shmemsize];

   shvhi[tdx] = Vhi[idx]; // thread tdx loads the data at the global index
   shvlo[tdx] = Vlo[idx];

   // result = -B[0]*shv[tdx];
   ddg_mul(-Bhi[0],-Blo[0],shvhi[tdx],shvlo[tdx],&resulthi,&resultlo);

   if(idx < nrows)
   {
      Whi[idx] = resulthi;
      Wlo[idx] = resultlo;
   }
}

__global__ void cmplx2_beta_times_V
 ( int nrows, int szt, double *Bhi, double *Blo,
   double *Vrehi, double *Vrelo, double *Vimhi, double *Vimlo,
   double *Wrehi, double *Wrelo, double *Wimhi, double *Wimlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // thread tdx computes W[idx]
   double resultrehi,resultrelo;
   double resultimhi,resultimlo;

   __shared__ double shvrehi[cdd_shmemsize]; // to store a slice of V
   __shared__ double shvrelo[cdd_shmemsize];
   __shared__ double shvimhi[cdd_shmemsize];
   __shared__ double shvimlo[cdd_shmemsize];

   shvrehi[tdx] = Vrehi[idx]; // thread tdx loads data at the global index
   shvrelo[tdx] = Vrelo[idx];
   shvimhi[tdx] = Vimhi[idx];
   shvimlo[tdx] = Vimlo[idx];

   // resultre = -B[0]*shvre[tdx];
   // resultim = -B[0]*shvim[tdx];
   ddg_mul(-Bhi[0],-Blo[0],shvrehi[tdx],shvrelo[tdx],&resultrehi,&resultrelo);
   ddg_mul(-Bhi[0],-Blo[0],shvimhi[tdx],shvimlo[tdx],&resultimhi,&resultimlo);

   if(idx < nrows)
   {
      Wrehi[idx] = resultrehi;
      Wrelo[idx] = resultrelo;
      Wimhi[idx] = resultimhi;
      Wimlo[idx] = resultimlo;
   }
}

__global__ void dbl2_initialize_WYT
 ( int dim, int szt, double *Vhi, double *Vlo,
   double *Whi, double *Wlo, double *WYThi, double *WYTlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYT
   const int col = idx % dim;         // column index in WYT

   const double Vvalhi = Vhi[col];
   const double Vvallo = Vlo[col];
   const double Wvalhi = Whi[row];
   const double Wvallo = Wlo[row];
   // const double result = Vval*Wval;
   double resulthi,resultlo;

   ddg_mul(Vvalhi,Vvallo,Wvalhi,Wvallo,&resulthi,&resultlo);

   if(idx < dim*dim)
   {
      WYThi[idx] = resulthi;
      WYTlo[idx] = resultlo;
   }
}

__global__ void cmplx2_initialize_WYH
 ( int dim, int szt,
   double *Vrehi, double *Vrelo, double *Vimhi, double *Vimlo,
   double *Wrehi, double *Wrelo, double *Wimhi, double *Wimlo,
   double *WYHrehi, double *WYHrelo, double *WYHimhi, double *WYHimlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYH
   const int col = idx % dim;         // column index in WYH

   const double Vvalrehi = Vrehi[col];
   const double Vvalrelo = Vrelo[col];
   const double Vvalimhi = Vimhi[col];
   const double Vvalimlo = Vimlo[col];
   const double Wvalrehi = Wrehi[row];
   const double Wvalrelo = Wrelo[row];
   const double Wvalimhi = Wimhi[row];
   const double Wvalimlo = Wimlo[row];
   // const double result = Vval*Wval;
   double resultrehi,resultrelo,resultimhi,resultimlo,acchi,acclo;

   // take the Hermitian transpose of V
   ddg_mul(Vvalrehi,Vvalrelo,Wvalrehi,Wvalrelo,&resultrehi,&resultrelo);
   ddg_mul(Vvalimhi,Vvalimlo,Wvalimhi,Wvalimlo,&acchi,&acclo);
   ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
   ddg_mul(Vvalrehi,Vvalrelo,Wvalimhi,Wvalimlo,&resultimhi,&resultimlo);
   ddg_mul(Vvalimhi,Vvalimlo,Wvalrehi,Wvalrelo,&acchi,&acclo);
   ddg_dec(&resultimhi,&resultimlo,acchi,acclo);

   if(idx < dim*dim)
   {
      WYHrehi[idx] = resultrehi;
      WYHrelo[idx] = resultrelo;
      WYHimhi[idx] = resultimhi;
      WYHimlo[idx] = resultimlo;
   }
}

__global__ void dbl2_update_WYT
 ( int dim, int szt, double *Vhi, double *Vlo, double *Whi, double *Wlo,
   double *WYThi, double *WYTlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYT
   const int col = idx % dim;         // column index in WYT

   const double Vvalhi = Vhi[col];
   const double Vvallo = Vlo[col];
   const double Wvalhi = Whi[row];
   const double Wvallo = Wlo[row];
   double resulthi = WYThi[idx];
   double resultlo = WYTlo[idx];
   double acchi,acclo;

   // result = result + Vval*Wval;

   ddg_mul(Vvalhi,Vvallo,Wvalhi,Wvallo,&acchi,&acclo);
   ddg_inc(&resulthi,&resultlo,acchi,acclo);
   
   if(idx < dim*dim)
   {
      WYThi[idx] = resulthi;
      WYTlo[idx] = resultlo;
   }
}

__global__ void cmplx2_update_WYH
 ( int dim, int szt,
   double *Vrehi, double *Vrelo, double *Vimhi, double *Vimlo,
   double *Wrehi, double *Wrelo, double *Wimhi, double *Wimlo,
   double *WYHrehi, double *WYHrelo, double *WYHimhi, double *WYHimlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYH
   const int col = idx % dim;         // column index in WYH

   const double Vvalrehi = Vrehi[col];
   const double Vvalrelo = Vrelo[col];
   const double Vvalimhi = Vimhi[col];
   const double Vvalimlo = Vimlo[col];
   const double Wvalrehi = Wrehi[row];
   const double Wvalrelo = Wrelo[row];
   const double Wvalimhi = Wimhi[row];
   const double Wvalimlo = Wimlo[row];
   // const double result = Vval*Wval;
   double resultrehi = WYHrehi[idx];
   double resultrelo = WYHrelo[idx];
   double resultimhi = WYHimhi[idx];
   double resultimlo = WYHimlo[idx];
   double acchi,acclo;

   // take the Hermitian transpose of V
   ddg_mul(Vvalrehi,Vvalrelo,Wvalrehi,Wvalrelo,&acchi,&acclo);
   ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
   ddg_mul(Vvalimhi,Vvalimlo,Wvalimhi,Wvalimlo,&acchi,&acclo);
   ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
   ddg_mul(Vvalrehi,Vvalrelo,Wvalimhi,Wvalimlo,&acchi,&acclo);
   ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
   ddg_mul(Vvalimhi,Vvalimlo,Wvalrehi,Wvalrelo,&acchi,&acclo);
   ddg_dec(&resultimhi,&resultimlo,acchi,acclo);

   if(idx < dim*dim)
   {
      WYHrehi[idx] = resultrehi;
      WYHrelo[idx] = resultrelo;
      WYHimhi[idx] = resultimhi;
      WYHimlo[idx] = resultimlo;
   }
}

__global__ void dbl2_beta_next_W
 ( int nrows, int szt, double *Bhi, double *Blo, double *Vhi, double *Vlo,
   double *Whi, double *Wlo, double *WYThi, double *WYTlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int WYToff = idx*nrows;      // start of idx row in YWT
   const double mybetahi = Bhi[0];
   const double mybetalo = Blo[0];
   int vdx,ydx;
   double resulthi,WYTvalhi,Vvalhi,acchi;
   double resultlo,WYTvallo,Vvallo,acclo;

   __shared__ double shVhi[dd_shmemsize];   // to store a slice of V
   __shared__ double shVlo[dd_shmemsize];

   shVhi[tdx] = Vhi[idx]; // thread tdx loads the data at the global index
   shVlo[tdx] = Vlo[idx];

   __syncthreads();
   resulthi = shVhi[tdx]; // thread tdx computes the value at index idx
   resultlo = shVlo[tdx];

   for(int i=0; i<nrows/szt; i++)
   {
      vdx = i*szt + tdx;                 // index in V and in YWT
      shVhi[tdx] = Vhi[vdx];             // threads load next szt values
      shVlo[tdx] = Vlo[vdx];

      __syncthreads();
      for(int j=0; j<szt; j++)           // multiply szt values with YWT
      {
         ydx = WYToff + i*szt + j;       // WYT is stored row by row
         WYTvalhi = WYThi[ydx];
         WYTvallo = WYTlo[ydx];
         Vvalhi = shVhi[j];
         Vvallo = shVlo[j];
         // result = result + WYTval*Vvalue;
         ddg_mul(Vvalhi,Vvallo,WYTvalhi,WYTvallo,&acchi,&acclo);
         ddg_inc(&resulthi,&resultlo,acchi,acclo);
      }
      __syncthreads();
   }
   int quot = nrows/szt;
   int rest = nrows - quot*szt;          // remainder to compute

   vdx = quot*szt + tdx;                 // next index to compute
   shVhi[tdx] = Vhi[vdx];
   shVlo[tdx] = Vlo[vdx];

   for(int j=0; j<rest; j++)            // rest < szt prevents overflow
   {
      __syncthreads();
      ydx = WYToff + quot*szt + j;
      WYTvalhi = WYThi[ydx];
      WYTvallo = WYTlo[ydx];
      Vvalhi = shVhi[j];
      Vvallo = shVlo[j];
      // result = result + WYTval*Vvalue;
      ddg_mul(Vvalhi,Vvallo,WYTvalhi,WYTvallo,&acchi,&acclo);
      ddg_inc(&resulthi,&resultlo,acchi,acclo);
   }
   // result = -mybeta*result;
   ddg_mul(-mybetahi,-mybetalo,resulthi,resultlo,&acchi,&acclo);

   if(idx < nrows) 
   {
      Whi[idx] = acchi;
      Wlo[idx] = acclo;
   }
}

__global__ void cmplx2_beta_next_W
 ( int nrows, int szt, double *Bhi, double *Blo,
   double *Vrehi, double *Vrelo, double *Vimhi, double *Vimlo,
   double *Wrehi, double *Wrelo, double *Wimhi, double *Wimlo,
   double *WYHrehi, double *WYHrelo, double *WYHimhi, double *WYHimlo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int WYHoff = idx*nrows;      // start of idx row in YWH
   const double mybetahi = Bhi[0];
   const double mybetalo = Blo[0];
   int vdx,ydx;
   double resultrehi,WYHvalrehi,Vvalrehi,acchi;
   double resultrelo,WYHvalrelo,Vvalrelo,acclo;
   double resultimhi,WYHvalimhi,Vvalimhi;
   double resultimlo,WYHvalimlo,Vvalimlo;

   __shared__ double shVrehi[cdd_shmemsize];   // to store a slice of V
   __shared__ double shVrelo[cdd_shmemsize];
   __shared__ double shVimhi[cdd_shmemsize];
   __shared__ double shVimlo[cdd_shmemsize];

   shVrehi[tdx] = Vrehi[idx]; // thread tdx loads data at the global index
   shVrelo[tdx] = Vrelo[idx];
   shVimhi[tdx] = Vimhi[idx];
   shVimlo[tdx] = Vimlo[idx];

   __syncthreads();
   resultrehi = shVrehi[tdx]; // thread tdx computes the value at index idx
   resultrelo = shVrelo[tdx];
   resultimhi = shVimhi[tdx];
   resultimlo = shVimlo[tdx];

   for(int i=0; i<nrows/szt; i++)
   {
      vdx = i*szt + tdx;                 // index in V and in YWT
      shVrehi[tdx] = Vrehi[vdx];         // threads load next szt values
      shVrelo[tdx] = Vrelo[vdx];
      shVimhi[tdx] = Vimhi[vdx];
      shVimlo[tdx] = Vimlo[vdx];

      __syncthreads();
      for(int j=0; j<szt; j++)           // multiply szt values with YWT
      {
         ydx = WYHoff + i*szt + j;       // WYH is stored row by row
         WYHvalrehi = WYHrehi[ydx];
         WYHvalrelo = WYHrelo[ydx];
         WYHvalimhi = WYHimhi[ydx];
         WYHvalimlo = WYHimlo[ydx];
         Vvalrehi = shVrehi[j];
         Vvalrelo = shVrelo[j];
         Vvalimhi = shVimhi[j];
         Vvalimlo = shVimlo[j];
         // result = result + WYHval*Vvalue;
         // take the Hermitian transpose of V
         ddg_mul(Vvalrehi,Vvalrelo,WYHvalrehi,WYHvalrelo,&acchi,&acclo);
         ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
         ddg_mul(Vvalimhi,Vvalimlo,WYHvalimhi,WYHvalimlo,&acchi,&acclo);
         ddg_dec(&resultrehi,&resultrelo,acchi,acclo);
         ddg_mul(Vvalimhi,Vvalimlo,WYHvalrehi,WYHvalrelo,&acchi,&acclo);
         ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
         ddg_mul(Vvalrehi,Vvalrelo,WYHvalimhi,WYHvalimlo,&acchi,&acclo);
         ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
      }
      __syncthreads();
   }
   int quot = nrows/szt;
   int rest = nrows - quot*szt;          // remainder to compute

   vdx = quot*szt + tdx;                 // next index to compute
   shVrehi[tdx] = Vrehi[vdx];
   shVrelo[tdx] = Vrelo[vdx];
   shVimhi[tdx] = Vimhi[vdx];
   shVimlo[tdx] = Vimlo[vdx];

   for(int j=0; j<rest; j++)            // rest < szt prevents overflow
   {
      __syncthreads();
      ydx = WYHoff + quot*szt + j;
      WYHvalrehi = WYHrehi[ydx];
      WYHvalrelo = WYHrelo[ydx];
      WYHvalimhi = WYHimhi[ydx];
      WYHvalimlo = WYHimlo[ydx];
      Vvalrehi = shVrehi[j];
      Vvalrelo = shVrelo[j];
      Vvalimhi = shVimhi[j];
      Vvalimlo = shVimlo[j];
      // result = result + WYTval*Vvalue;
      // take the Hermitian transpose of V
      ddg_mul(Vvalrehi,Vvalrelo,WYHvalrehi,WYHvalrelo,&acchi,&acclo);
      ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
      ddg_mul(Vvalimhi,Vvalimlo,WYHvalimhi,WYHvalimlo,&acchi,&acclo);
      ddg_dec(&resultrehi,&resultrelo,acchi,acclo);
      ddg_mul(Vvalimhi,Vvalimlo,WYHvalrehi,WYHvalrelo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
      ddg_mul(Vvalrehi,Vvalrelo,WYHvalimhi,WYHvalimlo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
   }
   // result = -mybeta*result;
   ddg_mul(-mybetahi,-mybetalo,resultrehi,resultrelo,&acchi,&acclo);

   if(idx < nrows) 
   {
      Wrehi[idx] = acchi;
      Wrelo[idx] = acclo;
   }
   ddg_mul(-mybetahi,-mybetalo,resultimhi,resultimlo,&acchi,&acclo);

   if(idx < nrows) 
   {
      Wimhi[idx] = acchi;
      Wimlo[idx] = acclo;
   }
}

__global__ void dbl2_small_WYT
 ( int nrows, int szt, double *Whi, double *Wlo, double *Yhi, double *Ylo,
   double *WYThi, double *WYTlo )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resulthi = 0.0;
   double resultlo = 0.0;
   double ahi,alo,bhi,blo,chi,clo;

   for(int k=0; k<szt; k++)
   {
      ahi = Whi[k*nrows + row];   // if(nrows == szt) then row = bdx
      alo = Wlo[k*nrows + row];
      bhi = Yhi[k*nrows + col];   // if(nrows == szt) then col = tdx
      blo = Ylo[k*nrows + col]; 
      // result = result + a*b;
      ddg_mul(ahi,alo,bhi,blo,&chi,&clo);
      ddg_inc(&resulthi,&resultlo,chi,clo);
   }
   __syncthreads();
   WYThi[offset] = resulthi;
   WYTlo[offset] = resultlo;
}

__global__ void cmplx2_small_WYH
 ( int nrows, int szt,
   double *Wrehi, double *Wrelo, double *Wimhi, double *Wimlo,
   double *Yrehi, double *Yrelo, double *Yimhi, double *Yimlo,
   double *WYTrehi, double *WYTrelo, double *WYTimhi, double *WYTimlo )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resultrehi = 0.0;
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;
   double a_rehi,a_relo,a_imhi,a_imlo;
   double b_rehi,b_relo,b_imhi,b_imlo,acchi,acclo;
   int Widx,Yidx;

   for(int k=0; k<szt; k++)
   {
      Widx = k*nrows + row;
      a_rehi = Wrehi[Widx];            // if(nrows == szt) then row = bdx
      a_relo = Wrelo[Widx];
      a_imhi = Wimhi[Widx]; 
      a_imlo = Wimlo[Widx]; 
      Yidx = k*nrows + col;
      b_rehi = Yrehi[Yidx];            // if(nrows == szt) then col = tdx
      b_relo = Yrelo[Yidx];
      b_imhi = Yimhi[Yidx];
      b_imlo = Yimlo[Yidx];
      // result = result + a*b; with Hermitian transpose of Y
      // resultre = resultre + a_re*b_re + a_im*b_im;
      ddg_mul(a_rehi,a_relo,b_rehi,b_relo,&acchi,&acclo);
      ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
      ddg_mul(a_imhi,a_imlo,b_imhi,b_imlo,&acchi,&acclo);
      ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
      // resultim = resultim + a_im*b_re - a_re*b_im;
      ddg_mul(a_imhi,a_imlo,b_rehi,b_relo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
      ddg_mul(a_rehi,a_relo,b_imhi,b_imlo,&acchi,&acclo);
      ddg_dec(&resultimhi,&resultimlo,acchi,acclo);
   }
   __syncthreads();
   WYTrehi[offset] = resultrehi;
   WYTrelo[offset] = resultrelo;
   WYTimhi[offset] = resultimhi;
   WYTimlo[offset] = resultimlo;
}

__global__ void dbl2_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhi, double *Qlo, double *WYThi, double *WYTlo,
   double *QWYThi, double *QWYTlo )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resulthi = 0.0;
   double resultlo = 0.0;
   double ahi,alo,bhi,blo,chi,clo;
   int idx;

   for(int k=0; k<rowdim; k++)       // run over rowdim, not just szt
   {                                 // coloff shifts by col*row elements
      idx = row*dim + coloff + k;
      ahi = Qhi[idx];                // row = bdx,
      alo = Qlo[idx];                // if dim == szt, coloff == 0
      idx = k*rowdim + col;
      bhi = WYThi[idx];              // if(dim == szt) then col = tdx
      blo = WYTlo[idx];              // if(dim == szt) then col = tdx
      // result = result + a*b;
      ddg_mul(ahi,alo,bhi,blo,&chi,&clo);
      ddg_inc(&resulthi,&resultlo,chi,clo);
   }
   __syncthreads();
   QWYThi[offset] = resulthi;        // no column offset in saving QWYT
   QWYTlo[offset] = resultlo;  
}

__global__ void cmplx2_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehi, double *Qrelo, double *Qimhi, double *Qimlo,
   double *WYTrehi, double *WYTrelo, double *WYTimhi, double *WYTimlo,
   double *QWYTrehi, double *QWYTrelo, double *QWYTimhi, double *QWYTimlo )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resultrehi = 0.0;
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;
   double a_rehi,a_relo,a_imhi,a_imlo;
   double b_rehi,b_relo,b_imhi,b_imlo;
   double acchi,acclo;
   int Qidx,WYTidx;

   for(int k=0; k<rowdim; k++)          // run over rowdim, not just szt
   {                                    // coloff shifts by col*row elements
      Qidx = row*dim + coloff + k;
      a_rehi = Qrehi[Qidx];             // row = bdx,
      a_relo = Qrelo[Qidx];
      a_imhi = Qimhi[Qidx];             // if dim == szt, coloff == 0
      a_imlo = Qimlo[Qidx];
      WYTidx = k*rowdim + col;
      b_rehi = WYTrehi[WYTidx];         // if(dim == szt) then col = tdx
      b_relo = WYTrelo[WYTidx];
      b_imhi = WYTimhi[WYTidx];
      b_imlo = WYTimlo[WYTidx];
      // result = result + a*b;
      // resultre = resultre + a_re*b_re - a_im*b_im;
      ddg_mul(a_rehi,a_relo,b_rehi,b_relo,&acchi,&acclo);
      ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
      ddg_mul(a_imhi,a_imlo,b_imhi,b_imlo,&acchi,&acclo);
      ddg_dec(&resultrehi,&resultrelo,acchi,acclo);
      // resultim = resultim + a_im*b_re + a_re*b_im;
      ddg_mul(a_imhi,a_imlo,b_rehi,b_relo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
      ddg_mul(a_rehi,a_relo,b_imhi,b_imlo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
   }
   __syncthreads();
   QWYTrehi[offset] = resultrehi;       // no column offset in saving QWYT
   QWYTrelo[offset] = resultrelo;
   QWYTimhi[offset] = resultimhi;
   QWYTimlo[offset] = resultimlo;
}

__global__ void dbl2_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWThi, double *YWTlo,
   double *Chi, double *Clo, double *YWTChi, double *YWTClo )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C

   double resulthi = 0.0;
   double resultlo = 0.0;
   double ahi,alo,bhi,blo,chi,clo;
   int idx;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      idx = row*rowdim + k;
      ahi = YWThi[idx];                // YWT is stored row by row
      alo = YWTlo[idx];
      idx = colCoff0 + k;
      bhi = Chi[idx];                  // but C is stored column by column
      blo = Clo[idx];
      // result = result + a*b;
      ddg_mul(ahi,alo,bhi,blo,&chi,&clo);
      ddg_inc(&resulthi,&resultlo,chi,clo);
   }
   __syncthreads();
   idx = (coloff + col)*nrows + (rowoff + row);
   YWTChi[idx] = resulthi;
   YWTClo[idx] = resultlo;
}

__global__ void cmplx2_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehi, double *YWTrelo, double *YWTimhi, double *YWTimlo,
   double *Crehi, double *Crelo, double *Cimhi, double *Cimlo,
   double *YWTCrehi, double *YWTCrelo, double *YWTCimhi, double *YWTCimlo )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C
   int idx;

   double resultrehi = 0.0;
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;
   double a_rehi,a_relo,a_imhi,a_imlo;
   double b_rehi,b_relo,b_imhi,b_imlo;
   double acchi,acclo;
   int YWTidx,Cidx;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      YWTidx = row*rowdim + k;
      a_rehi = YWTrehi[YWTidx];        // YWT is stored row by row
      a_relo = YWTrelo[YWTidx];
      a_imhi = YWTimhi[YWTidx];
      a_imlo = YWTimlo[YWTidx];
      Cidx = colCoff0 + k;
      b_rehi = Crehi[Cidx];            // but C is stored column by column
      b_relo = Crelo[Cidx];
      b_imhi = Cimhi[Cidx];
      b_imlo = Cimlo[Cidx];
      // result = result + a*b;
      // resultre = resultre + a_re*b_re - a_im*b_im;
      ddg_mul(a_rehi,a_relo,b_rehi,b_relo,&acchi,&acclo);
      ddg_inc(&resultrehi,&resultrelo,acchi,acclo);
      ddg_mul(a_imhi,a_imlo,b_imhi,b_imlo,&acchi,&acclo);
      ddg_dec(&resultrehi,&resultrelo,acchi,acclo);
      // resultim = resultim + a_im*b_re + a_re*b_im;
      ddg_mul(a_imhi,a_imlo,b_rehi,b_relo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
      ddg_mul(a_rehi,a_relo,b_imhi,b_imlo,&acchi,&acclo);
      ddg_inc(&resultimhi,&resultimlo,acchi,acclo);
   }
   __syncthreads();
   idx = (coloff + col)*nrows + (rowoff + row);
   YWTCrehi[idx] = resultrehi;
   YWTCrelo[idx] = resultrelo;
   YWTCimhi[idx] = resultimhi;
   YWTCimlo[idx] = resultimlo;
}

__global__ void dbl2_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhi, double *Qlo, double *QWYThi, double *QWYTlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double ahi,alo,bhi,blo;

   ahi = Qhi[idx1];       // row = bdx, if dim == szt, coloff == 0
   alo = Qlo[idx1];
   bhi = QWYThi[offset];  // if(dim == szt) then col = tdx
   blo = QWYTlo[offset];
   // a = a + b;
   ddg_inc(&ahi,&alo,bhi,blo);
   __syncthreads();
   Qhi[idx1] = ahi;
   Qlo[idx1] = alo;
}

__global__ void cmplx2_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehi, double *Qrelo, double *Qimhi, double *Qimlo,
   double *QWYTrehi, double *QWYTrelo, double *QWYTimhi, double *QWYTimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double a_rehi,a_relo,a_imhi,a_imlo;
   double b_rehi,b_relo,b_imhi,b_imlo;

   a_rehi = Qrehi[idx1];       // row = bdx, if dim == szt, coloff == 0
   a_relo = Qrelo[idx1];
   a_imhi = Qimhi[idx1];
   a_imlo = Qimlo[idx1];
   b_rehi = QWYTrehi[offset];  // if(dim == szt) then col = tdx
   b_relo = QWYTrelo[offset];
   b_imhi = QWYTimhi[offset];
   b_imlo = QWYTimlo[offset];
   // a_re = a_re + b_re;
   ddg_inc(&a_rehi,&a_relo,b_rehi,b_relo);
   // a_im = a_im + b_im;
   ddg_inc(&a_imhi,&a_imlo,b_imhi,b_imlo);

   __syncthreads();
   Qrehi[idx1] = a_rehi;
   Qrelo[idx1] = a_relo;
   Qimhi[idx1] = a_imhi;
   Qimlo[idx1] = a_imlo;
}

__global__ void dbl2_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhi, double *Rlo, double *YWTChi, double *YWTClo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double ahi,alo,bhi,blo;
   
   ahi = Rhi[idx];
   alo = Rlo[idx];
   bhi = YWTChi[idx];
   blo = YWTClo[idx];
   // a = a + b;
   ddg_inc(&ahi,&alo,bhi,blo);
  
   __syncthreads();
   Rhi[idx] = ahi;
   Rlo[idx] = alo;
}

__global__ void cmplx2_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehi, double *Rrelo, double *Rimhi, double *Rimlo, 
   double *YWTCrehi, double *YWTCrelo, double *YWTCimhi, double *YWTCimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double a_rehi,a_relo,a_imhi,a_imlo;
   double b_rehi,b_relo,b_imhi,b_imlo;
   
   a_rehi = Rrehi[idx];
   a_relo = Rrelo[idx];
   a_imhi = Rimhi[idx];
   a_imlo = Rimlo[idx];
   b_rehi = YWTCrehi[idx];
   b_relo = YWTCrelo[idx];
   b_imhi = YWTCimhi[idx];
   b_imlo = YWTCimlo[idx];
   // a_re = a_re + b_re;
   ddg_inc(&a_rehi,&a_relo,b_rehi,b_relo);
   // a_im = a_im + b_im;
   ddg_inc(&a_imhi,&a_imlo,b_imhi,b_imlo);
  
   __syncthreads();
   Rrehi[idx] = a_rehi;
   Rrelo[idx] = a_relo;
   Rimhi[idx] = a_imhi;
   Rimlo[idx] = a_imlo;
}

void GPU_dbl2_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *vhi_h, double *vlo_h, double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   const int nrLog2 = ceil(log2((double) nrows1));
   const int rowidx = colidx*(nrows+1);       // start of number in A_h
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   if(verbose)
   {
      cout << "nrows : " << nrows
           << "  nVrows : " << nVrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  nbt : " << nbt << endl;
      cout << "k : " << k 
           << "  L : " << L
           << "  nrows1 : " << nrows1
           << "  colidx : " << colidx
           << "  rowidx : " << rowidx << endl;
   }
   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vhi_h[i] = 0.0;
         vlo_h[i] = 0.0;
      }
      cudaMemcpy(&Vhi_d[L*nVrows],vhi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlo_d[L*nVrows],vlo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahi_h[L] = 0.0; vhi_h[0] = 1.0;
      betalo_h[L] = 0.0; vlo_h[0] = 0.0;
      cudaMemcpy(&betahi_d[L],&betahi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalo_d[L],&betalo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhi_d[L*nVrows+L],vhi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlo_d[L*nVrows+L],vlo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl2_small_house<<<1,nrows1>>>
         (&Ahi_d[rowidx],&Alo_d[rowidx],&Ahi_d[rowidx+1],&Alo_d[rowidx+1],
          nrows1,nrLog2,&Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
          &betahi_d[L],&betalo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   cudaMemcpy(&betahi_h[L],&betahi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalo_h[L],&betalo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(vhi_h,&Vhi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vlo_h,&Vlo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahi_h[L] << "  " << betalo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : " << vhi_h[i] << "  " << vlo_h[i] << endl;
   }
}

void GPU_cmplx2_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehi_h, double *Arelo_h, double *Aimhi_h, double *Aimlo_h,
   double *Arehi_d, double *Arelo_d, double *Aimhi_d, double *Aimlo_d,
   double *vrehi_h, double *vrelo_h, double *vimhi_h, double *vimlo_h,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   const int nrLog2 = ceil(log2((double) nrows1));
   const int rowidx = colidx*(nrows+1);       // start of number in A_h
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   if(verbose)
   {
      cout << "nrows : " << nrows
           << "  nVrows : " << nVrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  nbt : " << nbt << endl;
      cout << "k : " << k 
           << "  L : " << L
           << "  nrows1 : " << nrows1
           << "  colidx : " << colidx
           << "  rowidx : " << rowidx << endl;
   }
   if(L > 0)
   {
      for(int i=0; i<L; i++)   // insert zeros
      {
         vrehi_h[i] = 0.0; vrelo_h[i] = 0.0;
         vimhi_h[i] = 0.0; vimlo_h[i] = 0.0;
      }
      cudaMemcpy(&Vrehi_d[L*nVrows],vrehi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelo_d[L*nVrows],vrelo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhi_d[L*nVrows],vimhi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlo_d[L*nVrows],vimlo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahi_h[L] = 0.0; betalo_h[L] = 0.0;
      vrehi_h[0] = 1.0; vrelo_h[0] = 0.0;
      vimhi_h[0] = 0.0; vimlo_h[0] = 0.0;
      cudaMemcpy(&betahi_d[L],&betahi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalo_d[L],&betalo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehi_d[L*nVrows+L],vrehi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelo_d[L*nVrows+L],vrelo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhi_d[L*nVrows+L],vimhi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlo_d[L*nVrows+L],vimlo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      if(verbose)  // to verify the input column is correct ...
      {
         cout << "The column of A : " << endl;
         cudaMemcpy(&Arehi_h[rowidx],&Arehi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelo_h[rowidx],&Arelo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhi_h[rowidx],&Aimhi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlo_h[rowidx],&Aimlo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         for(int i=0; i<=nrows1; i++)
         {
            cout << "A[" << i << "]re : " 
                 << Arehi_h[i] << "  " << Arelo_h[i] << endl;
            cout << "A[" << i << "]im : " 
                 << Aimhi_h[i] << "  " << Aimlo_h[i] << endl;
         }
      }
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      cmplx2_small_house<<<1,nrows1>>>
         (&Arehi_d[rowidx],&Arelo_d[rowidx],
          &Aimhi_d[rowidx],&Aimlo_d[rowidx],
          &Arehi_d[rowidx+1],&Arelo_d[rowidx+1],
          &Aimhi_d[rowidx+1],&Aimlo_d[rowidx+1],nrows1,nrLog2,
          &Vrehi_d[L*nVrows+L],&Vrelo_d[L*nVrows+L],
          &Vimhi_d[L*nVrows+L],&Vimlo_d[L*nVrows+L],
          &betahi_d[L],&betalo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   cudaMemcpy(&betahi_h[L],&betahi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalo_h[L],&betalo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(vrehi_h,&Vrehi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelo_h,&Vrelo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhi_h,&Vimhi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlo_h,&Vimlo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahi_h[L] << "  " << betalo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
      {
         cout << "v[" << i << "]re : "
              << vrehi_h[i] << "  " << vrelo_h[i] << endl;
         cout << "v[" << i << "]im : "
              << vimhi_h[i] << "  " << vimlo_h[i] << endl;
      }
   }
}

void GPU_dbl2_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *vhi_h, double *vlo_h, double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *sumshi_h, double *sumslo_h, double *sumshi_d, double *sumslo_d,
   double *sigmahi_h, double *sigmalo_h, double *sigmahi_d, double *sigmalo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   // nrows1 = nrows - colidx - 1 = size of Householder vector
   const int nblocks = ceil(((double) nrows1)/szt); // sufficient threads
   const int nblLog2 = ceil(log2((double) nblocks));
   const int sztLog2 = ceil(log2((double) szt));
   const int rowidx = colidx*(nrows+1);         // start of number in A_h
   const int nVrows = nrows - k*szt;             // dimension of V matrix

   cudaEvent_t start,stop;            // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vhi_h[i] = 0.0;
         vlo_h[i] = 0.0;
      }
   }
   vhi_h[L] = 1.0;                    // set one on the diagonal
   vlo_h[L] = 0.0;
   cudaMemcpy(&Vhi_d[L*nVrows],vhi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlo_d[L*nVrows],vlo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);

   if(verbose)
   {
      cout << "-> launching " << nblocks << " blocks of "
           << szt << " threads to compute the sum of squares ..." << endl;
      cout << "nrows1 : " << nrows1 << "  rowidx : " << rowidx;
      cout << "  ceil(log2(#blocks)) : " << nblLog2;
      cout << "  ceil(log2(szt)) : " << sztLog2 << endl;
   }
   for(int i=0; i<nblocks; i++)
   {
      sumshi_h[i] = 0.0;
      sumslo_h[i] = 0.0;
   }
   cudaMemcpy(sumshi_d,sumshi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslo_d,sumslo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl2_large_sum_of_squares<<<nblocks,szt>>>
      (&Ahi_d[rowidx+1],&Alo_d[rowidx+1],sumshi_d,sumslo_d,
       nrows1,szt,sztLog2);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

   if(verbose)
   {
      cout << "-> launching 1 block of " << nblocks
           << " threads to accumulate the sums ..." << endl;
   }
   cudaEventRecord(start);
   dbl2_sum_accumulator<<<1,nblocks>>>
      (sumshi_d,sumslo_d,nblocks,nblLog2,sigmahi_d,sigmalo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahi_h,sigmahi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalo_h,sigmalo_d,sizeof(double),cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahi_h[0] == 0.0) && (sigmalo_h[0] == 0.0))
   {
      betahi_h[L] = 0.0; betalo_h[L] = 0.0; done = true;
      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      // const double x0hi = Ahi_h[rowidx];
      // const double x0lo = Alo_h[rowidx];
      double acchi,acclo,muhi,mulo,v0hi,v0lo,v0p2hi,v0p2lo;
      double x0hi,x0lo;

      cudaMemcpy(&x0hi,&Ahi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lo,&Alo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // mu = sqrt((*x0)*(*x0) + sigma[0]);
      ddf_sqr(x0hi,x0lo,&acchi,&acclo);
      ddf_inc(&acchi,&acclo,sigmahi_h[0],sigmalo_h[0]);
      ddf_sqrt(acchi,acclo,&muhi,&mulo);
      if(x0hi <= 0.0)
      {
         // v0 = *x0 - mu;
         ddf_sub(x0hi,x0lo,muhi,mulo,&v0hi,&v0lo);
      }
      else
      {
         // v0 = -sigma[0]/(*x0 + mu);
         ddf_add(x0hi,x0lo,muhi,mulo,&acchi,&acclo);
         ddf_div(sigmahi_h[0],sigmalo_h[0],acchi,acclo,&v0hi,&v0lo);
         ddf_minus(&v0hi,&v0lo);
      }
      // v0p2 = v0*v0;
      ddf_sqr(v0hi,v0lo,&v0p2hi,&v0p2lo);
      // *beta = 2.0*v0p2/(sigma[0] + v0p2);
      ddf_add(sigmahi_h[0],sigmalo_h[0],v0p2hi,v0p2lo,&acchi,&acclo);
      ddf_div(v0p2hi,v0p2lo,acchi,acclo,&betahi_h[L],&betalo_h[L]);
      ddf_mlt_d(&betahi_h[L],&betalo_h[L],2.0);
      sigmahi_h[0] = v0hi;
      sigmalo_h[0] = v0lo;                     // v0 needed for normalization
      // update the flop counts
      *add += 3;
      *mul += 3;
      *div += 2;
      *sqrtfun += 1;
   }
   if(verbose)
   {
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahi_h[L] << "  " << betalo_h[L] << endl;
   }
   cudaMemcpy(&betahi_d[L],&betahi_h[L],sizeof(double),cudaMemcpyHostToDevice);
   cudaMemcpy(&betalo_d[L],&betalo_h[L],sizeof(double),cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0hi, v0lo).
      cudaMemcpy(sigmahi_d,sigmahi_h,sizeof(double),cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalo_d,sigmalo_h,sizeof(double),cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
      cudaEventRecord(start);
      dbl2_normalize<<<nblocks,szt>>>
         (nrows1,szt,&Ahi_d[rowidx+1],&Alo_d[rowidx+1],sigmahi_d,sigmalo_d,
          &Vhi_d[L*nVrows+L+1],&Vlo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_normalize(nblocks,szt,div);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahi_h[L],&betahi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalo_h[L],&betalo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhi_h,&Vhi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vlo_h,&Vlo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahi_h[L] << "  " << betalo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : " << vhi_h[i] << "  " << vlo_h[i] << endl;
   }
}

void GPU_cmplx2_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehi_h, double *Arelo_h, double *Aimhi_h, double *Aimlo_h,
   double *Arehi_d, double *Arelo_d, double *Aimhi_d, double *Aimlo_d,
   double *vrehi_h, double *vrelo_h, double *vimhi_h, double *vimlo_h,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *sumshi_h, double *sumslo_h, double *sumshi_d, double *sumslo_d,
   double *sigmahi_h, double *sigmalo_h, double *sigmahi_d, double *sigmalo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   // nrows1 = nrows - colidx - 1 = size of Householder vector
   const int nblocks = ceil(((double) nrows1)/szt); // sufficient threads
   const int nblLog2 = ceil(log2((double) nblocks));
   const int sztLog2 = ceil(log2((double) szt));
   const int rowidx = colidx*(nrows+1);         // start of number in A_h
   const int nVrows = nrows - k*szt;             // dimension of V matrix

   cudaEvent_t start,stop;            // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vrehi_h[i] = 0.0;
         vrelo_h[i] = 0.0;
         vimhi_h[i] = 0.0;
         vimlo_h[i] = 0.0;
      }
   }
   vrehi_h[L] = 1.0;                    // set one on the diagonal
   vrelo_h[L] = 0.0;
   vimhi_h[L] = 0.0;
   vimlo_h[L] = 0.0;
   cudaMemcpy(&Vrehi_d[L*nVrows],vrehi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelo_d[L*nVrows],vrelo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhi_d[L*nVrows],vimhi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlo_d[L*nVrows],vimlo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);

   if(verbose)
   {
      cout << "-> launching " << nblocks << " blocks of "
           << szt << " threads to compute the sum of squares ..." << endl;
      cout << "nrows1 : " << nrows1 << "  rowidx : " << rowidx;
      cout << "  ceil(log2(#blocks)) : " << nblLog2;
      cout << "  ceil(log2(szt)) : " << sztLog2 << endl;
   }
   for(int i=0; i<nblocks; i++)
   {
      sumshi_h[i] = 0.0;
      sumslo_h[i] = 0.0;
   }
   cudaMemcpy(sumshi_d,sumshi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslo_d,sumslo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   cmplx2_large_sum_of_squares<<<nblocks,szt>>>
      (&Arehi_d[rowidx+1],&Arelo_d[rowidx+1],
       &Aimhi_d[rowidx+1],&Aimlo_d[rowidx+1],
       sumshi_d,sumslo_d,nrows1,szt,sztLog2);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

   if(verbose)
   {
      cout << "-> launching 1 block of " << nblocks
           << " threads to accumulate the sums ..." << endl;
   }
   cudaEventRecord(start);
   dbl2_sum_accumulator<<<1,nblocks>>>
      (sumshi_d,sumslo_d,nblocks,nblLog2,sigmahi_d,sigmalo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahi_h,sigmahi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalo_h,sigmalo_d,sizeof(double),cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahi_h[0] == 0.0) && (sigmalo_h[0] == 0.0))
   {
      betahi_h[L] = 0.0; betalo_h[L] = 0.0; done = true;
      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      double acchi,acclo,muhi,mulo,sqrv0hi,sqrv0lo;
      double x0rehi,x0relo,x0imhi,x0imlo,sqrx0hi,sqrx0lo;
      double v0rehi,v0relo,v0imhi,v0imlo,x0radhi,x0radlo;
      double inv0rehi,inv0relo,inv0imhi,inv0imlo;

      cudaMemcpy(&x0rehi,&Arehi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relo,&Arelo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhi,&Aimhi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlo,&Aimlo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      ddf_sqr(x0rehi,x0relo,&sqrx0hi,&sqrx0lo);
      ddf_sqr(x0imhi,x0imlo,&acchi,&acclo);
      ddf_inc(&sqrx0hi,&sqrx0lo,acchi,acclo);
      // x0rad = sqrt(sqrx0);
      ddf_sqrt(sqrx0hi,sqrx0lo,&x0radhi,&x0radlo);
      // mu = sqrt(sqrx0 + sigma); // norm of the vector x
      ddf_inc(&sqrx0hi,&sqrx0lo,*sigmahi_h,*sigmalo_h);
      ddf_sqrt(sqrx0hi,sqrx0lo,&muhi,&mulo);

      if((x0radhi == 0.0) && (x0radlo == 0.0))
      {
         v0rehi = -muhi; v0imhi = 0.0;
         v0relo = -mulo; v0imlo = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         // mu = mu/x0rad;
         ddf_div(muhi,mulo,x0radhi,x0radlo,&acchi,&acclo);
         muhi = acchi; mulo = acclo;
         // vre[0] = xre[0] - mu*xre[0];
         ddf_mul(muhi,mulo,x0rehi,x0relo,&acchi,&acclo);
         ddf_sub(x0rehi,x0relo,acchi,acclo,&v0rehi,&v0relo);
         // vim[0] = xim[0] - mu*xim[0];
         ddf_mul(muhi,mulo,x0imhi,x0imlo,&acchi,&acclo);
         ddf_sub(x0imhi,x0imlo,acchi,acclo,&v0imhi,&v0imlo);
      }
      // sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      ddf_sqr(v0rehi,v0relo,&sqrv0hi,&sqrv0lo);
      ddf_sqr(v0imhi,v0imlo,&acchi,&acclo);
      ddf_inc(&sqrv0hi,&sqrv0lo,acchi,acclo);
      // *beta = 2.0*sqrv0/(sigma + sqrv0);
      ddf_inc(sigmahi_h,sigmalo_h,sqrv0hi,sqrv0lo);
      ddf_div(sqrv0hi,sqrv0lo,*sigmahi_h,*sigmalo_h,
              &betahi_h[L],&betalo_h[L]);
      ddf_mlt_d(&betahi_h[L],&betalo_h[L],2.0);
      // inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      ddf_div(v0rehi,v0relo,sqrv0hi,sqrv0lo,&inv0rehi,&inv0relo);
      // inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]
      ddf_div(v0imhi,v0imlo,sqrv0hi,sqrv0lo,&inv0imhi,&inv0imlo);
      ddf_minus(&inv0imhi,&inv0imlo);
      *sigmahi_h = inv0rehi;
      *sigmalo_h = inv0relo;
      betahi_h[szt] = inv0imhi;
      betalo_h[szt] = inv0imlo;
      // update the flop counts
      *add += 6;
      *mul += 7;
      *div += 4;
      *sqrtfun += 2;
   }
   if(verbose)
   {
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahi_h[L] << "  " << betalo_h[L] << endl;
   }
   cudaMemcpy(&betahi_d[L],&betahi_h[L],sizeof(double),cudaMemcpyHostToDevice);
   cudaMemcpy(&betalo_d[L],&betalo_h[L],sizeof(double),cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0rehi, v0relo)
      cudaMemcpy(sigmahi_d,sigmahi_h,sizeof(double),cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalo_d,sigmalo_h,sizeof(double),cudaMemcpyHostToDevice);
      // (betahi_h[szt], betalo_h[szt]) has the values for (v0imhi, v0rimlo).
      cudaMemcpy(&betahi_d[szt],&betahi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalo_d[szt],&betalo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
      cudaEventRecord(start);
      cmplx2_normalize<<<nblocks,szt>>>
         (nrows1,szt,&Arehi_d[rowidx+1],&Arelo_d[rowidx+1],
                     &Aimhi_d[rowidx+1],&Aimlo_d[rowidx+1],
          sigmahi_d,sigmalo_d,&betahi_d[szt],&betalo_d[szt],
          &Vrehi_d[L*nVrows+L+1],&Vrelo_d[L*nVrows+L+1],
          &Vimhi_d[L*nVrows+L+1],&Vimlo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_normalize(nblocks,szt,add,mul);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahi_h[L],&betahi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalo_h[L],&betalo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehi_h,&Vrehi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelo_h,&Vrelo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhi_h,&Vimhi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlo_h,&Vimlo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahi_h[L] << "  " << betalo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
      {
         cout << "v[" << i << "]re : "
              << vrehi_h[i] << "  " << vrelo_h[i] << endl;
         cout << "v[" << i << "]im : "
              << vimhi_h[i] << "  " << vimlo_h[i] << endl;
      }
   }
}

void GPU_dbl2_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);           // 2nd argument: ncols -> szt
   // changed second argument ncols into szt
   // to avoid updating the next tile
   dbl2_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,Ahi_d,Alo_d,
       &Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],&betahi_d[L],&betalo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Ahi_h,Ahi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alo_h,Alo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahi_h[j*nrows+i] << "  "
                 << Alo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx2_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehi_h, double *Arelo_h, double *Aimhi_h, double *Aimlo_h,
   double *Arehi_d, double *Arelo_d, double *Aimhi_d, double *Aimlo_d,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);
   cmplx2_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,
       Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
       &Vrehi_d[L*nVrows+L],&Vrelo_d[L*nVrows+L],
       &Vimhi_d[L*nVrows+L],&Vimlo_d[L*nVrows+L],
       &betahi_d[L],&betalo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Arehi_h,Arehi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelo_h,Arelo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhi_h,Aimhi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlo_h,Aimlo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A_d[" << i << "][" << j << "]re : "
                 << Arehi_h[j*nrows+i] << "  "
                 << Arelo_h[j*nrows+i] << endl;
            cout << "A_d[" << i << "][" << j << "]im : "
                 << Aimhi_h[j*nrows+i] << "  "
                 << Aimlo_h[j*nrows+i] << endl;
         }
   }
}

void GPU_dbl2_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *RTdotvhi_h, double *RTdotvlo_h,
   double *RTdotvhi_d, double *RTdotvlo_d,
   double *whi_h, double *wlo_h, double *whi_d, double *wlo_d,
   double *RTvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix
   const int nhouse = nrows - colidx;  // length of Householder vector
   // total number of entries in R that will be modified
   const int RToffset = colidx*nrows;
   const int dimRTdotv = endcol - colidx;
   // total number of entries in R that will be modified
   const int sizenum = (nrows - colidx)*dimRTdotv;
   const int nbrblocks = (int) ceil(sizenum/((double) szt));

   // changed second argument ncols into endcol
   // to avoid updating the next tile
   // dbl_medium_betaRTv<<<nbrblocks,szt>>>
   //   (nrows,endcol,szt,colidx,A_d,&V_d[L*nVrows+L],&beta_d[L],w_d);
   // number of threads must be ncols - colidx, not endcol - colidx
   // dbl2_small_betaRTv<<<1,nrows-colidx>>> // nrows ...
   //   (nrows,endcol,szt,colidx,Ahi_d,Alo_d,
   //    &Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
   //    &betahi_d[L],&betalo_d[L],whi_d,wlo_d);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to compute RTdotv ..." << endl;
      cout << "   nhouse : " << nhouse << "  RToffset : " << RToffset
           << "  dimRTdotv : " << dimRTdotv << endl;
   }

   cudaEventRecord(start);
   dbl2_RTdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RToffset,dimRTdotv,Ahi_d,Alo_d,
       &Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],RTdotvhi_d,RTdotvlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   cudaEventRecord(start);
   dbl2_sum_betaRTdotv<<<1,dimRTdotv>>>
      (nhouse,&betahi_d[L],&betalo_d[L],RTdotvhi_d,RTdotvlo_d,whi_d,wlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   flopcount_dbl_RTdotv(nhouse,szt,mul);
   flopcount_dbl_sum_betaRTdotv(nhouse,dimRTdotv,add,mul);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to update " << sizenum << " numbers ..." << endl;
      cout << "   nrows : " << nrows << "  endcol : " << endcol
           << "  szt : " << szt << "  colidx : " << colidx << endl;
   }
   cudaEventRecord(start);
   dbl2_medium_subvbetaRTv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Ahi_d,Alo_d,&Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
       &betahi_d[L],&betalo_d[L],whi_d,wlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
   flopcount_dbl_medium_subvbetaRTv(nrows,endcol,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRTv = dimRTdotv*sizeof(double);
      const size_t szRTdotv = nVrows*szbRTv;

      cudaMemcpy(RTdotvhi_h,RTdotvhi_d,szRTdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlo_h,RTdotvlo_d,szRTdotv,cudaMemcpyDeviceToHost);
      cout << "the matrix R^T dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RTdotv[" << i << "][" << j << "] : "
                 << RTdotvhi_h[ix] << "  "
                 << RTdotvlo_h[ix] << endl;
            ix = ix + 1;
         }
      }

      cudaMemcpy(whi_h,whi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlo_h,wlo_d,szbRTv,cudaMemcpyDeviceToHost);
      cout << "the vector w = beta*R^T*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
         cout << "w[" << i << "] : "
              << whi_h[i] << "  " << wlo_h[i] << endl;

      cudaMemcpy(Ahi_h,Ahi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alo_h,Alo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahi_h[j*nrows+i] << "  "
                 << Alo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx2_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehi_h, double *Arelo_h, double *Aimhi_h, double *Aimlo_h,
   double *Arehi_d, double *Arelo_d, double *Aimhi_d, double *Aimlo_d,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *RHdotvrehi_h, double *RHdotvrelo_h,
   double *RHdotvimhi_h, double *RHdotvimlo_h,
   double *RHdotvrehi_d, double *RHdotvrelo_d,
   double *RHdotvimhi_d, double *RHdotvimlo_d,
   double *wrehi_h, double *wrelo_h, double *wimhi_h, double *wimlo_h,
   double *wrehi_d, double *wrelo_d, double *wimhi_d, double *wimlo_d,
   double *RHvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix
   const int nhouse = nrows - colidx;  // length of Householder vector
   // total number of entries in R that will be modified
   const int RHoffset = colidx*nrows;
   const int dimRHdotv = endcol - colidx;
   // total number of entries in R that will be modified
   const int sizenum = (nrows - colidx)*dimRHdotv;
   const int nbrblocks = (int) ceil(sizenum/((double) szt));

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to compute RHdotv ..." << endl;
      cout << "   nhouse : " << nhouse << "  RHoffset : " << RHoffset
           << "  dimRHdotv : " << dimRHdotv << endl;
   }

   cudaEventRecord(start);
   cmplx2_RHdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
       &Vrehi_d[L*nVrows+L],&Vrelo_d[L*nVrows+L],
       &Vimhi_d[L*nVrows+L],&Vimlo_d[L*nVrows+L],
       RHdotvrehi_d,RHdotvrelo_d,RHdotvimhi_d,RHdotvimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
   cudaEventRecord(start);
   cmplx2_sum_betaRHdotv<<<1,dimRHdotv>>>
      (nhouse,&betahi_d[L],&betalo_d[L],
       RHdotvrehi_d,RHdotvrelo_d,RHdotvimhi_d,RHdotvimlo_d,
       wrehi_d,wrelo_d,wimhi_d,wimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
   flopcount_cmplx_RHdotv(nhouse,szt,add,mul);
   flopcount_cmplx_sum_betaRHdotv(nhouse,dimRHdotv,add,mul);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to update " << sizenum << " numbers ..." << endl;
      cout << "   nrows : " << nrows << "  endcol : " << endcol
           << "  szt : " << szt << "  colidx : " << colidx << endl;
   }

   cudaEventRecord(start);
   cmplx2_medium_subvbetaRHv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
       &Vrehi_d[L*nVrows+L],&Vrelo_d[L*nVrows+L],
       &Vimhi_d[L*nVrows+L],&Vimlo_d[L*nVrows+L],
       &betahi_d[L],&betalo_d[L],wrehi_d,wrelo_d,wimhi_d,wimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
   flopcount_cmplx_medium_subvbetaRHv(nrows,endcol,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRHv = dimRHdotv*sizeof(double);
      const size_t szRHdotv = nVrows*szbRHv;

      cudaMemcpy(RHdotvrehi_h,RHdotvrehi_d,szRHdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelo_h,RHdotvrelo_d,szRHdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhi_h,RHdotvimhi_d,szRHdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlo_h,RHdotvimlo_d,szRHdotv,cudaMemcpyDeviceToHost);
      cout << "the matrix R^H dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RHdotv[" << i << "][" << j << "]re : "
                 << RHdotvrehi_h[ix] << "  "
                 << RHdotvrelo_h[ix] << endl;
            cout << "RHdotv[" << i << "][" << j << "]im : "
                 << RHdotvimhi_h[ix] << "  "
                 << RHdotvimlo_h[ix] << endl;
            ix = ix + 1;
         }
      }

      cudaMemcpy(wrehi_h,wrehi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelo_h,wrelo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhi_h,wimhi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlo_h,wimlo_d,szbRHv,cudaMemcpyDeviceToHost);
      cout << "the vector w = beta*R^H*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
      {
         cout << "w[" << i << "]re : "
              << wrehi_h[i] << "  " << wimlo_h[i] << endl;
         cout << "w[" << i << "]im : "
              << wimhi_h[i] << "  " << wimlo_h[i] << endl;
      }
      cudaMemcpy(Arehi_h,Arehi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelo_h,Arelo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhi_h,Aimhi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlo_h,Aimlo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A_d[" << i << "][" << j << "]re : "
                 << Arehi_h[j*nrows+i] << "  "
                 << Arelo_h[j*nrows+i] << endl;
            cout << "A_d[" << i << "][" << j << "]im : "
                 << Aimhi_h[j*nrows+i] << "  "
                 << Aimlo_h[j*nrows+i] << endl;
         }
   }
}

void GPU_dbl2_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhi_h, double *Vlo_h, double *Vhi_d, double *Vlo_d,
   double *Whi_h, double *Wlo_h, double *Whi_d, double *Wlo_d,
   double *WYThi_h, double *WYTlo_h, double *WYThi_d, double *WYTlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   dbl2_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahi_d,betalo_d,Vhi_d,Vlo_d,Whi_d,Wlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl2_initialize_WYT<<<nbrblocks2,szt>>>
      (rowdim,szt,Vhi_d,Vlo_d,Whi_d,Wlo_d,WYThi_d,WYTlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_initialize_WYT(rowdim,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      dbl2_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,&betahi_d[j],&betalo_d[j],&Vhi_d[j*rowdim],
          &Vlo_d[j*rowdim],&Whi_d[j*rowdim],&Wlo_d[j*rowdim],WYThi_d,WYTlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      dbl2_update_WYT<<<nbrblocks2,szt>>>
         (rowdim,szt,&Vhi_d[j*rowdim],&Vlo_d[j*rowdim],&Whi_d[j*rowdim],
          &Wlo_d[j*rowdim],WYThi_d,WYTlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_update_WYT(rowdim,add,mul);
   }
   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = rowdim*sizeof(double);
      const size_t szVandW = szt*szhouse;
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(betahi_h,betahi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalo_h,betalo_d,szbeta,cudaMemcpyDeviceToHost);
      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahi_h[j] << "  " << betalo_h[j] << endl;

      cudaMemcpy(Vhi_h,Vhi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlo_h,Vlo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "] : "
                 << Vhi_h[ix] << "  " << Vlo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Whi_h,Whi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlo_h,Wlo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "] : "
                 << Whi_h[ix] << "  " << Wlo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYThi_h,WYThi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlo_h,WYTlo_d,szmat,cudaMemcpyDeviceToHost);
      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThi_h[ix] << "  " << WYTlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx2_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehi_h, double *Vrelo_h, double *Vimhi_h, double *Vimlo_h,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *Wrehi_h, double *Wrelo_h, double *Wimhi_h, double *Wimlo_h,
   double *Wrehi_d, double *Wrelo_d, double *Wimhi_d, double *Wimlo_d,
   double *WYHrehi_h, double *WYHrelo_h, double *WYHimhi_h, double *WYHimlo_h,
   double *WYHrehi_d, double *WYHrelo_d, double *WYHimhi_d, double *WYHimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx2_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahi_d,betalo_d,
       Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,Wrehi_d,Wrelo_d,Wimhi_d,Wimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx2_initialize_WYH<<<nbrblocks2,szt>>>
      (rowdim,szt,Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,Wrehi_d,Wrelo_d,
       Wimhi_d,Wimlo_d,WYHrehi_d,WYHrelo_d,WYHimhi_d,WYHimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_initialize_WYH(rowdim,add,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      cmplx2_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,&betahi_d[j],&betalo_d[j],
          &Vrehi_d[j*rowdim],&Vrelo_d[j*rowdim],
          &Vimhi_d[j*rowdim],&Vimlo_d[j*rowdim],
          &Wrehi_d[j*rowdim],&Wrelo_d[j*rowdim],
          &Wimhi_d[j*rowdim],&Wimlo_d[j*rowdim],
          WYHrehi_d,WYHrelo_d,WYHimhi_d,WYHimlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      cmplx2_update_WYH<<<nbrblocks2,szt>>>
         (rowdim,szt,&Vrehi_d[j*rowdim],&Vrelo_d[j*rowdim],
                     &Vimhi_d[j*rowdim],&Vimlo_d[j*rowdim],
                     &Wrehi_d[j*rowdim],&Wrelo_d[j*rowdim],
                     &Wimhi_d[j*rowdim],&Wimlo_d[j*rowdim],
          WYHrehi_d,WYHrelo_d,WYHimhi_d,WYHimlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_update_WYH(rowdim,add,mul);
   }
   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = rowdim*sizeof(double);
      const size_t szVandW = szt*szhouse;
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(betahi_h,betahi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalo_h,betalo_d,szbeta,cudaMemcpyDeviceToHost);
      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahi_h[j] << "  " << betalo_h[j] << endl;

      cudaMemcpy(Vrehi_h,Vrehi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelo_h,Vrelo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhi_h,Vimhi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlo_h,Vimlo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "]re : "
                 << Vrehi_h[ix] << "  " << Vrelo_h[ix] << endl;
            cout << "V[" << i << "][" << j << "]im : "
                 << Vimhi_h[ix] << "  " << Vimlo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Wrehi_h,Wrehi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelo_h,Wrelo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhi_h,Wimhi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlo_h,Wimlo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehi_h[ix] << "  " << Wrelo_h[ix] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhi_h[ix] << "  " << Wimlo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYHrehi_h,WYHrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelo_h,WYHrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhi_h,WYHimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlo_h,WYHimlo_d,szmat,cudaMemcpyDeviceToHost);
      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYH[" << i << "][" << j << "]re : "
                 << WYHrehi_h[ix] << "  " << WYHrelo_h[ix] << endl;
            cout << "WYH[" << i << "][" << j << "]im : "
                 << WYHimhi_h[ix] << "  " << WYHimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl2_small_WYT
 ( int nrows, int szt,
   double *Whi_d, double *Wlo_d, double *Yhi_d, double *Ylo_d,
   double *WYThi_d, double *WYTlo_d, double *WYThi_h, double *WYTlo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   dbl2_small_WYT<<<nbrblocks,szt>>>
      (nrows,szt,Whi_d,Wlo_d,Yhi_d,Ylo_d,WYThi_d,WYTlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_dbl_small_WYT(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYThi_h,WYThi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlo_h,WYTlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThi_h[ix] << "  " << WYTlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx2_small_WYH
 ( int nrows, int szt,
   double *Wrehi_d, double *Wrelo_d, double *Wimhi_d, double *Wimlo_d,
   double *Yrehi_d, double *Yrelo_d, double *Yimhi_d, double *Yimlo_d,
   double *WYTrehi_d, double *WYTrelo_d, double *WYTimhi_d, double *WYTimlo_d,
   double *WYTrehi_h, double *WYTrelo_h, double *WYTimhi_h, double *WYTimlo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   cmplx2_small_WYH<<<nbrblocks,szt>>>
      (nrows,szt,Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
                 Yrehi_d,  Yrelo_d,  Yimhi_d,  Yimlo_d,
               WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_cmplx_small_WYH(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYTrehi_h,WYTrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelo_h,WYTrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhi_h,WYTimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlo_h,WYTimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "]re : "
                 << WYTrehi_h[ix] << "  " << WYTrelo_h[ix] << endl;
            cout << "WYT[" << i << "][" << j << "]im : "
                 << WYTimhi_h[ix] << "  " << WYTimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl2_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhi_d, double *Ylo_d, double *Whi_d, double *Wlo_d,
   double *YWThi_d, double *YWTlo_d, double *YWThi_h, double *YWTlo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl2_small_WYT<<<nbrblocks,szt>>>
      (rowdim,szt,Yhi_d,Ylo_d,Whi_d,Wlo_d,YWThi_d,YWTlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_WYT(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWThi_h,YWThi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlo_h,YWTlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThi_h[ix] << "  " << YWTlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx2_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehi_d, double *Yrelo_d, double *Yimhi_d, double *Yimlo_d,
   double *Wrehi_d, double *Wrelo_d, double *Wimhi_d, double *Wimlo_d,
   double *YWTrehi_d, double *YWTrelo_d, double *YWTimhi_d, double *YWTimlo_d,
   double *YWTrehi_h, double *YWTrelo_h, double *YWTimhi_h, double *YWTimlo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx2_small_WYH<<<nbrblocks,szt>>>
      (rowdim,szt,Yrehi_d,  Yrelo_d,  Yimhi_d,  Yimlo_d,
                  Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
                YWTrehi_d,YWTrelo_d,YWTimhi_d,YWTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_WYH(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWTrehi_h,YWTrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelo_h,YWTrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhi_h,YWTimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlo_h,YWTimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "]re : "
                 << YWTrehi_h[ix] << "  " << YWTrelo_h[ix] << endl;
            cout << "YWT[" << i << "][" << j << "]im : "
                 << YWTimhi_h[ix] << "  " << YWTimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl2_small_QWYT
 ( int dim, int szt, int idx, double *Qhi_d, double *Qlo_d,
   double *WYThi_d, double *WYTlo_d, double *QWYThi_d, double *QWYTlo_d,
   double *QWYThi_h, double *QWYTlo_h, double *Qhi_h, double *Qlo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhi_h,Qhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlo_h,Qlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhi_h[ix] << "  " << Qlo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   dbl2_small_QWYT<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qhi_d,Qlo_d,WYThi_d,WYTlo_d,QWYThi_d,QWYTlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_QWYT(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYThi_h,QWYThi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlo_h,QWYTlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYThi_h[ix] << "  " << QWYTlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx2_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehi_d, double *Qrelo_d, double *Qimhi_d, double *Qimlo_d,
   double *WYTrehi_d, double *WYTrelo_d, double *WYTimhi_d, double *WYTimlo_d,
   double *QWYTrehi_d, double *QWYTrelo_d,
   double *QWYTimhi_d, double *QWYTimlo_d,
   double *QWYTrehi_h, double *QWYTrelo_h,
   double *QWYTimhi_h, double *QWYTimlo_h,
   double *Qrehi_h, double *Qrelo_h, double *Qimhi_h, double *Qimlo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qrehi_h,Qrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelo_h,Qrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhi_h,Qimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlo_h,Qimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehi_h[ix] << "  " << Qrelo_h[ix] << endl;
            cout << "Q[" << i << "][" << j << "] : "
                 << Qimhi_h[ix] << "  " << Qimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   cmplx2_small_QWYH<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
        WYTrehi_d, WYTrelo_d, WYTimhi_d, WYTimlo_d,
       QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_QWYH(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYTrehi_h,QWYTrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelo_h,QWYTrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhi_h,QWYTimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlo_h,QWYTimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "]re : "
                 << QWYTrehi_h[ix] << "  " << QWYTrelo_h[ix] << endl;
            cout << "QWYT[" << i << "][" << j << "]im : "
                 << QWYTimhi_h[ix] << "  " << QWYTimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl2_small_YWTC
 ( int nrows, int ncols, int szt, int idx, double *YWThi_d, double *YWTlo_d,
   double *Chi_d, double *Clo_d, double *YWTChi_d, double *YWTClo_d,
   double *YWTChi_h, double *YWTClo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   if(verbose)
   {
      cout << "in GPU_dbl2_small_YWTC ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff
           << "  nbrblocks : " << nbrblocks << endl;

      double *Chi_h = new double[nrows*ncols];
      double *Clo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Chi_h,Chi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clo_h,Clo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
            cout << "C_h[" << i << "][" << j << "] : "
                 << Chi_h[j*nrows+i] << "  "
                 << Clo_h[j*nrows+i] << endl;

      free(Chi_h); free(Clo_h);
   }

   cudaEventRecord(start);
   dbl2_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
       YWThi_d,YWTlo_d,Chi_d,Clo_d,YWTChi_d,YWTClo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_YWTC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTChi_h,YWTChi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClo_h,YWTClo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "YWTC[" << i << "][" << j << "] : "
                 << YWTChi_h[j*nrows + i] << "  "
                 << YWTClo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx2_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehi_d, double *YWTrelo_d, double *YWTimhi_d, double *YWTimlo_d,
   double *Crehi_d, double *Crelo_d, double *Cimhi_d, double *Cimlo_d,
   double *YWTCrehi_d, double *YWTCrelo_d,
   double *YWTCimhi_d, double *YWTCimlo_d,
   double *YWTCrehi_h, double *YWTCrelo_h,
   double *YWTCimhi_h, double *YWTCimlo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   if(verbose)
   {
      cout << "in GPU_dbl_small_YWTC ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff
           << "  nbrblocks : " << nbrblocks << endl;

      double *Crehi_h = new double[nrows*ncols];
      double *Crelo_h = new double[nrows*ncols];
      double *Cimhi_h = new double[nrows*ncols];
      double *Cimlo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Crehi_h,Crehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelo_h,Crelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhi_h,Cimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlo_h,Cimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
         {
            cout << "C_h[" << i << "][" << j << "]re : "
                 << Crehi_h[j*nrows+i] << "  "
                 << Crelo_h[j*nrows+i] << endl;
            cout << "C_h[" << i << "][" << j << "]im : "
                 << Cimhi_h[j*nrows+i] << "  "
                 << Cimlo_h[j*nrows+i] << endl;
         }

      free(Crehi_h); free(Crelo_h); free(Cimhi_h); free(Cimlo_h);
   }
   cudaEventRecord(start);
   cmplx2_small_YWHC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWTrehi_d, YWTrelo_d, YWTimhi_d, YWTimlo_d,
          Crehi_d,   Crelo_d,   Cimhi_d,   Cimlo_d,
       YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_YWHC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTCrehi_h,YWTCrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelo_h,YWTCrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhi_h,YWTCimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlo_h,YWTCimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
         {
            cout << "YWTC[" << i << "][" << j << "]re : "
                 << YWTCrehi_h[j*nrows + i] << "  "
                 << YWTCrelo_h[j*nrows + i] << endl;
            cout << "YWTC[" << i << "][" << j << "]im : "
                 << YWTCimhi_h[j*nrows + i] << "  "
                 << YWTCimlo_h[j*nrows + i] << endl;
         }
   }
}

void GPU_dbl2_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhi_d, double *Qlo_d, double *QWYThi_d, double *QWYTlo_d,
   double *Qhi_h, double *Qlo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl2_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qhi_d,Qlo_d,QWYThi_d,QWYTlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhi_h,Qhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlo_h,Qlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhi_h[ix] << "  " << Qlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx2_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehi_d, double *Qrelo_d, double *Qimhi_d, double *Qimlo_d,
   double *QWYTrehi_d, double *QWYTrelo_d,
   double *QWYTimhi_d, double *QWYTimlo_d,
   double *Qrehi_h, double *Qrelo_h, double *Qimhi_h, double *Qimlo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx2_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
       QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qrehi_h,Qrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelo_h,Qrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhi_h,Qimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlo_h,Qimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehi_h[ix] << "  " << Qrelo_h[ix] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhi_h[ix] << "  " << Qimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl2_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx, double *Rhi_d, double *Rlo_d,
   double *YWTChi_d, double *YWTClo_d, double *Rhi_h, double *Rlo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   cudaEventRecord(start);
   dbl2_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,Rhi_d,Rlo_d,YWTChi_d,YWTClo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_R_add_YWTC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rhi_h,Rhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlo_h,Rlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhi_h[j*nrows + i] << "  "
                 << Rlo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx2_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehi_d, double *Rrelo_d, double *Rimhi_d, double *Rimlo_d,
   double *YWTCrehi_d, double *YWTCrelo_d,
   double *YWTCimhi_d, double *YWTCimlo_d,
   double *Rrehi_h, double *Rrelo_h, double *Rimhi_h, double *Rimlo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   cudaEventRecord(start);
   cmplx2_small_R_add_YWHC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
          Rrehi_d,   Rrelo_d,   Rimhi_d,   Rimlo_d,
       YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_R_add_YWHC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rrehi_h,Rrehi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelo_h,Rrelo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhi_h,Rimhi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlo_h,Rimlo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
         {
            cout << "R[" << i << "][" << j << "]re : "
                 << Rrehi_h[j*nrows + i] << "  "
                 << Rrelo_h[j*nrows + i] << endl;
            cout << "R[" << i << "][" << j << "]im : "
                 << Rimhi_h[j*nrows + i] << "  "
                 << Rimlo_h[j*nrows + i] << endl;
         }
   }
}

void GPU_dbl2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahi, double **Alo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;         // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Ahi_h = new double[dim];     // high doubles of A on the host
   double *Alo_h = new double[dim];     // low doubles of A on the host
   double *Ahi_d;                       // Ahi on the device
   double *Alo_d;                       // Alo on the device
   double *Qhi_h = new double[nrows2];  // high doubles of Q on the host
   double *Qlo_h = new double[nrows2];  // low doubles of Q on the host
   double *Qhi_d;                       // Qhi on the device
   double *Qlo_d;                       // Qlo on the device
   double *vhi_h = new double[nrows];   // high doubles of Householder vector
   double *vlo_h = new double[nrows];   // low doubles of Householder vector
   double *betahi_h = new double[szt];  // high doubles of beta on the host
   double *betalo_h = new double[szt];  // low doubles of beta on the host
   double *betahi_d;                      // betahi on the device
   double *betalo_d;                      // betalo on the device
   double *Vhi_h = new double[nrows*szt]; // high doubles of V matrix
   double *Vlo_h = new double[nrows*szt]; // low doubles of V matrix
   double *Vhi_d;                         // Vhi on the device
   double *Vlo_d;                         // Vlo on the device
   double *Whi_h = new double[nrows*szt]; // high doubes of W on the host
   double *Wlo_h = new double[nrows*szt]; // low doules of W on the host
   double *Whi_d;                         // Whi on the device
   double *Wlo_d;                         // Wlo on the device
   double *WYThi_h = new double[nrows2];  // high doubles of W*Y^T 
   double *WYTlo_h = new double[nrows2];  // low doubles of W*Y^T
   double *WYThi_d;                       // WYThi on the device
   double *WYTlo_d;                       // WYTlo on the device
   double *YWThi_h = new double[nrows2];  // high doubles of Y*W^T
   double *YWTlo_h = new double[nrows2];  // low doubles of Y*W^T 
   double *YWThi_d;                       // YWThi on the device
   double *YWTlo_d;                       // YWTlo on the device
   double *QWYThi_h = new double[nrows2]; // high doubles of Q*WY^T
   double *QWYTlo_h = new double[nrows2]; // low doubles of Q*WY^T
   double *QWYThi_d;                      // QWYThi on the device
   double *QWYTlo_d;                      // QWYTlo on the device
   double *YWTChi_h = new double[dim];    // YWT*C on the host
   double *YWTClo_h = new double[dim];    // YWT*C on the host
   double *YWTChi_d;                      // YWTChi on the device
   double *YWTClo_d;                      // YWTClo on the device
   double *RTdotvhi_h = new double[nrows2]; // high R^T dotted with v
   double *RTdotvlo_h = new double[nrows2]; // low R^T dotted with v
   double *RTdotvhi_d;                      // RTdotvhi on the device
   double *RTdotvlo_d;                      // RTdotvlo on the device
   double *bRTvhi_h = new double[nrows];  // high doubles of beta*R^T*v
   double *bRTvlo_h = new double[nrows];  // low doubles of beta*R^T*v
   double *bRTvhi_d;                      // bRTvhi_h on the device
   double *bRTvlo_d;                      // bRTvlo_h on the device
   double *sumshi_h = new double[nrows];  // subsums for large house
   double *sumslo_h = new double[nrows];  // low doubles of subsums
   double *sumshi_d;                      // sumshi on the device
   double *sumslo_d;                      // sumslo on the device
   double sigmahi_h,sigmalo_h;
   double *sigmahi_d;                     // sigmahi on the device
   double *sigmalo_d;                     // sigmalo on the device

   int ix = 0;                          // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Ahi_h[ix]   = Ahi[i][j];
         Alo_h[ix++] = Alo[i][j];
      }

   ix = 0;                              // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qhi_h[ix]   = 1.0;
            Qlo_h[ix++] = 0.0;
         }
         else
         {
            Qhi_h[ix]   = 0.0;
            Qlo_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Ahi_d,sznum);
   cudaMalloc((void**)&Alo_d,sznum);
   cudaMemcpy(Ahi_d,Ahi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alo_d,Alo_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&betahi_d,szbeta);
   cudaMalloc((void**)&betalo_d,szbeta);
   for(int i=0; i<szt; i++)
   {
      betahi_h[i] = 0.0;
      betalo_h[i] = 0.0;
   }
   cudaMemcpy(betahi_d,betahi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalo_d,betalo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);  // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vhi_d,szVandW + szpad); // padding only in allocation
   cudaMalloc((void**)&Vlo_d,szVandW + szpad);
   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vhi_h[ix] = 0.0; 
      Vlo_h[ix++] = 0.0; 
   }
   Vhi_h[--ix] = 1.0; // initialize last vector for square tiles
   cudaMemcpy(Vhi_d,Vhi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlo_d,Vlo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Whi_d,szVandW + szpad); // padding only in allocation
   cudaMalloc((void**)&Wlo_d,szVandW + szpad); 

   cudaMalloc((void**)&RTdotvhi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlo_d,szVandW + szpad);
   cudaMalloc((void**)&bRTvhi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshi_d,szhouse);
   cudaMalloc((void**)&sumslo_d,szhouse);
   cudaMalloc((void**)&sigmahi_d,sizeof(double));
   cudaMalloc((void**)&sigmalo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYThi_d,szWYT + szpad); // padding for W*Y^T product
   cudaMalloc((void**)&WYTlo_d,szWYT + szpad); 
   cudaMalloc((void**)&Qhi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlo_d,szWYT + szpad);
   cudaMemcpy(Qhi_d,Qhi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlo_d,Qlo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYThi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWThi_d,szYWT + szpad); // padding for Y*W^T product
   cudaMalloc((void**)&YWTlo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTChi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClo_d,sznum + szpad);

   *houselapms = 0.0; *RTvlapms = 0.0; *tileRlapms = 0.0; *vb2Wlapms = 0.0;
   *WYTlapms = 0.0; *QWYTlapms = 0.0; *Qaddlapms = 0.0;
   *YWTlapms = 0.0; *YWTClapms = 0.0; *Raddlapms = 0.0;
   *addcnt = 0; *mulcnt = 0; *divcnt = 0; *sqrtcnt = 0;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      int colidx,nrows1;

      for(int L=0; L<szt; L++)  // L runs over the columns in one block
      {
         colidx = k*szt + L;              // index of the current column
         nrows1 = nrows - colidx - 1;     // #rows in Householder vector - 1

         if(verbose)
            cout << "-> current column : " << colidx << endl
                 << "-> #nrows in Householder vector - 1 : "
                 << nrows1 << endl;

         if(nrows1 <= szt)
         {
            GPU_dbl2_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                Ahi_h,Alo_h,Ahi_d,Alo_d,vhi_h,vlo_h,Vhi_d,Vlo_d,
                betahi_h,betalo_h,betahi_d,betalo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahi_h[L] == 0.0) && (betalo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_dbl2_small_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,Ahi_h,Alo_h,Ahi_d,Alo_d,
                   Vhi_d,Vlo_d,betahi_h,betalo_h,betahi_d,betalo_d,
                   tileRlapms,addcnt,mulcnt,verbose);
            }
         }
         else
         {
            GPU_dbl2_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                Ahi_h,Alo_h,Ahi_d,Alo_d,vhi_h,vlo_h,Vhi_d,Vlo_d,
                betahi_h,betalo_h,betahi_d,betalo_d,
                sumshi_h,sumslo_h,sumshi_d,sumslo_d,
                &sigmahi_h,&sigmalo_h,sigmahi_d,sigmalo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahi_h[L] == 0.0) && (betalo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_dbl2_medium_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,Ahi_h,Alo_h,Ahi_d,Alo_d,
                   Vhi_d,Vlo_d,betahi_h,betalo_h,betahi_d,betalo_d,
                   RTdotvhi_h,RTdotvlo_h,RTdotvhi_d,RTdotvlo_d,
                   bRTvhi_h,bRTvlo_h,bRTvhi_d,bRTvlo_d,
                   RTvlapms,tileRlapms,addcnt,mulcnt,verbose);
            }
         }
      }
      GPU_dbl2_medium_VB_to_W
         (nrows,szt,szt,k,Vhi_h,Vlo_h,Vhi_d,Vlo_d,Whi_h,Wlo_h,Whi_d,Wlo_d,
          WYThi_h,WYTlo_h,WYThi_d,WYTlo_d,betahi_h,betalo_h,betahi_d,betalo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);
/*
      GPU_dbl2_small_WYT
         (nrows-k*szt,szt,Whi_d,Wlo_d,Vhi_d,Vlo_d,WYThi_d,WYTlo_d,
          WYThi_h,WYTlo_h,WYTlapms,verbose);
 */
      GPU_dbl2_small_QWYT
         (nrows,szt,k,Qhi_d,Qlo_d,WYThi_d,WYTlo_d,QWYThi_d,QWYTlo_d,
          QWYThi_h,QWYTlo_h,Qhi_h,Qlo_h,QWYTlapms,addcnt,mulcnt,verbose);
      GPU_dbl2_small_Qupdate
         (nrows,szt,k,Qhi_d,Qlo_d,QWYThi_d,QWYTlo_d,Qhi_h,Qlo_h,
          Qaddlapms,addcnt,verbose);
      if(k < nbt-1)                                           // update R
      {
         GPU_dbl2_small_YWT
            (nrows,szt,k,Vhi_d,Vlo_d,Whi_d,Wlo_d,YWThi_d,YWTlo_d,
             YWThi_h,YWTlo_h,YWTlapms,addcnt,mulcnt,verbose);
         GPU_dbl2_small_YWTC
            (nrows,ncols,szt,k,YWThi_d,YWTlo_d,Ahi_d,Alo_d,
             YWTChi_d,YWTClo_d,YWTChi_h,YWTClo_h,
             YWTClapms,addcnt,mulcnt,verbose);
         GPU_dbl2_small_R_add_YWTC
            (nrows,ncols,szt,k,Ahi_d,Alo_d,YWTChi_d,YWTClo_d,
             Ahi_h,Alo_h,Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qhi_h,Qhi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlo_h,Qlo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qhi[i][j] = Qhi_h[ix];
         Qlo[i][j] = Qlo_h[ix++];
      }

   cudaMemcpy(Ahi_h,Ahi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alo_h,Alo_d,sznum,cudaMemcpyDeviceToHost);
   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rhi[i][j] = Ahi_h[j*nrows+i];
         Rlo[i][j] = Alo_h[j*nrows+i];
      }

   free(Ahi_h); free(Alo_h); free(Qhi_h); free(Qlo_h); 
   free(vhi_h); free(vlo_h); free(Vhi_h); free(Vlo_h);
   free(Whi_h); free(Wlo_h); free(sumshi_h); free(sumslo_h);
   free(RTdotvhi_h); free(RTdotvlo_h); free(bRTvhi_h); free(bRTvlo_h);
   free(WYThi_h); free(QWYThi_h); free(YWThi_h); free(YWTChi_h);
   free(WYTlo_h); free(QWYTlo_h); free(YWTlo_h); free(YWTClo_h);
}

void GPU_cmplx2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;        // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Arehi_h = new double[dim];  // high doubles of the real parts of A
   double *Arelo_h = new double[dim];  // low doubles of the real parts of A 
   double *Aimhi_h = new double[dim];  // high doubles of the imag parts of A
   double *Aimlo_h = new double[dim];  // low doubles of the imag parts of A
   double *Arehi_d;                    // Arehi on the device
   double *Arelo_d;                    // Arelo on the device
   double *Aimhi_d;                    // Aimhi on the device
   double *Aimlo_d;                    // Aimlo on the device
   double *Qrehi_h = new double[nrows2]; // high doubles of real parts of Q 
   double *Qrelo_h = new double[nrows2]; // low doubles of real parts of Q
   double *Qimhi_h = new double[nrows2]; // high doubles of imag parts of Q
   double *Qimlo_h = new double[nrows2]; // low doubles of imag parts of Q
   double *Qrehi_d;                      // Qrehi on device
   double *Qrelo_d;                      // Qrelo on device
   double *Qimhi_d;                      // Qimhi on device
   double *Qimlo_d;                      // Qimlo on device
   double *vrehi_h = new double[nrows];  // high real parts of Householder v
   double *vrelo_h = new double[nrows];  // low real parts of Householder v
   double *vimhi_h = new double[nrows];  // high imag parts of Householder v
   double *vimlo_h = new double[nrows];  // low imag parts of Householder v
   double *betahi_h = new double[szt+1]; // high doubles of beta
   double *betalo_h = new double[szt+1]; // low doubles of beta
   double *betahi_d;                     // betahi on the device
   double *betalo_d;                     // betalo on the device
   double *Vrehi_h = new double[nrows*szt]; // high real parts of Householder
   double *Vrelo_h = new double[nrows*szt]; // low real parts of Householder
   double *Vimhi_h = new double[nrows*szt]; // high imag parts of Householder
   double *Vimlo_h = new double[nrows*szt]; // low imag parts of Householder
   double *Vrehi_d;                         // Vrehi on device
   double *Vrelo_d;                         // Vrelo on device
   double *Vimhi_d;                         // Vimhi on device
   double *Vimlo_d;                         // Vimlo on device
   double *Wrehi_h = new double[nrows*szt]; // high doubles of real parts of W
   double *Wrelo_h = new double[nrows*szt]; // low doubles of real parts of W
   double *Wimhi_h = new double[nrows*szt]; // high doubles of imag parts of W
   double *Wimlo_h = new double[nrows*szt]; // low doubles of imag parts of W
   double *Wrehi_d;                         // Wrehi on the device
   double *Wrelo_d;                         // Wrelo on the device
   double *Wimhi_d;                         // Wimhi on the device
   double *Wimlo_d;                         // Wimlo on the device
   double *WYTrehi_h = new double[nrows2];  // high real parts of W*Y^T
   double *WYTrelo_h = new double[nrows2];  // low real parts of W*Y^T
   double *WYTimhi_h = new double[nrows2];  // high imag parts of W*Y^T
   double *WYTimlo_h = new double[nrows2];  // low imag parts of W*Y^T
   double *WYTrehi_d;                       // WYTrehi on the device 
   double *WYTrelo_d;                       // WYTrelo on the device 
   double *WYTimhi_d;                       // WYTimhi on the device
   double *WYTimlo_d;                       // WYTimlo on the device
   double *YWTrehi_h = new double[nrows2];  // high real parts of Y*W^T
   double *YWTrelo_h = new double[nrows2];  // low real parts of Y*W^T
   double *YWTimhi_h = new double[nrows2];  // high imag parts of Y*W^T
   double *YWTimlo_h = new double[nrows2];  // low imag parts of Y*W^T
   double *YWTrehi_d;                       // YWTrehi on the device
   double *YWTrelo_d;                       // YWTrelo on the device
   double *YWTimhi_d;                       // YWTimhi on the device
   double *YWTimlo_d;                       // YWTimlo on the device
   double *QWYTrehi_h = new double[nrows2]; // high real parts of Q*WY^T
   double *QWYTrelo_h = new double[nrows2]; // low real parts of Q*WY^T
   double *QWYTimhi_h = new double[nrows2]; // high imag parts of Q*WY^T
   double *QWYTimlo_h = new double[nrows2]; // low imag parts of Q*WY^T
   double *QWYTrehi_d;                      // QWYTrehi on the device
   double *QWYTrelo_d;                      // QWYTrelo on the device
   double *QWYTimhi_d;                      // QWYTimhi on the device
   double *QWYTimlo_d;                      // QWYTimlo on the device
   double *YWTCrehi_h = new double[dim];    // high real parts of YWT*C
   double *YWTCrelo_h = new double[dim];    // low real parts of YWT*C
   double *YWTCimhi_h = new double[dim];    // high imag parts of YWT*C
   double *YWTCimlo_h = new double[dim];    // low imag parts of YWT*C
   double *YWTCrehi_d;                      // YWTCrehi on the device
   double *YWTCrelo_d;                      // YWTCrelo on the device
   double *YWTCimhi_d;                      // YWTCimhi on the device
   double *YWTCimlo_d;                      // YWTCimlo on the device
   double *RHdotvrehi_h = new double[nrows2]; // high real R^H dotted with v
   double *RHdotvrelo_h = new double[nrows2]; // low real R^H dotted with v
   double *RHdotvimhi_h = new double[nrows2]; // high imag R^H dotted with v
   double *RHdotvimlo_h = new double[nrows2]; // low imag R^H dotted with v
   double *RHdotvrehi_d;                      // RHdotvrehi on the device
   double *RHdotvrelo_d;                      // RHdotvrelo on the device
   double *RHdotvimhi_d;                      // RHdotvimhi on the device
   double *RHdotvimlo_d;                      // RHdotvimlo on the device
   double *bRHvrehi_h = new double[nrows];  // high real parts of beta*R^H*v
   double *bRHvrelo_h = new double[nrows];  // low real parts of beta*R^H*v
   double *bRHvimhi_h = new double[nrows];  // high imag parts of beta*R^H*v
   double *bRHvimlo_h = new double[nrows];  // low imag parts of beta*R^H*v
   double *bRHvrehi_d;                      // bRHvrehi on the device
   double *bRHvrelo_d;                      // bRHvrelo on the device
   double *bRHvimhi_d;                      // bRHvimhi on the device
   double *bRHvimlo_d;                      // bRHvimlo on the device
   double *sumshi_h = new double[nrows];  // subsums for large house
   double *sumslo_h = new double[nrows];  // low doubles of subsums
   double *sumshi_d;                      // sumshi on the device
   double *sumslo_d;                      // sumslo on the device
   double sigmahi_h,sigmalo_h;
   double *sigmahi_d;                     // sigmahi on the device
   double *sigmalo_d;                     // sigmalo on the device

   int ix = 0;                            // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Arehi_h[ix]   = Arehi[i][j];
         Arelo_h[ix]   = Arelo[i][j];
         Aimhi_h[ix]   = Aimhi[i][j];
         Aimlo_h[ix++] = Aimlo[i][j];
      }

   ix = 0;                                // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qrehi_h[ix]   = 1.0;
            Qrelo_h[ix]   = 0.0;
            Qimhi_h[ix]   = 0.0;
            Qimlo_h[ix++] = 0.0;
         }
         else
         {
            Qrehi_h[ix]   = 0.0;
            Qrelo_h[ix]   = 0.0;
            Qimhi_h[ix]   = 0.0;
            Qimlo_h[ix++] = 0.0;
         }
         // cout << "Q[" << ix-1 << "] : "
         //      << Qre_h[ix-1] << "  " << Qim_h[ix-1] << endl;
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Arehi_d,sznum);
   cudaMalloc((void**)&Arelo_d,sznum);
   cudaMalloc((void**)&Aimhi_d,sznum);
   cudaMalloc((void**)&Aimlo_d,sznum);
   cudaMemcpy(Arehi_d,Arehi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelo_d,Arelo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhi_d,Aimhi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlo_d,Aimlo_h,sznum,cudaMemcpyHostToDevice);
   // allocate one extra beta for use in the cmplx2_normalize
   const size_t szbeta = (szt+1)*sizeof(double);
   cudaMalloc((void**)&betahi_d,szbeta);
   cudaMalloc((void**)&betalo_d,szbeta);
   for(int i=0; i<szt; i++)
   {
      betahi_h[i] = 0.0;
      betalo_h[i] = 0.0;
   }
   cudaMemcpy(betahi_d,betahi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalo_d,betalo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);    // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vrehi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Vrelo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlo_d,szVandW + szpad);
   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vrehi_h[ix]   = 0.0; 
      Vrelo_h[ix]   = 0.0; 
      Vimhi_h[ix]   = 0.0; 
      Vimlo_h[ix++] = 0.0; 
   }
   Vrehi_h[--ix] = 1.0; // initialize last vector for square tiles
   cudaMemcpy(Vrehi_d,Vrehi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelo_d,Vrelo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhi_d,Vimhi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlo_d,Vimlo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Wrehi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Wrelo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlo_d,szVandW + szpad);

   cudaMalloc((void**)&RHdotvrehi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlo_d,szVandW + szpad);
   cudaMalloc((void**)&bRHvrehi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshi_d,szhouse);
   cudaMalloc((void**)&sumslo_d,szhouse);
   cudaMalloc((void**)&sigmahi_d,sizeof(double));
   cudaMalloc((void**)&sigmalo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYTrehi_d,szWYT + szpad); // padding for W*Y^T 
   cudaMalloc((void**)&WYTrelo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlo_d,szWYT + szpad);
   cudaMemcpy(Qrehi_d,Qrehi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelo_d,Qrelo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhi_d,Qimhi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlo_d,Qimlo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYTrehi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWTrehi_d,szYWT + szpad); // padding for Y*W^T
   cudaMalloc((void**)&YWTrelo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTCrehi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlo_d,sznum + szpad);

   *houselapms = 0.0; *RHvlapms = 0.0; *tileRlapms = 0.0; *vb2Wlapms = 0.0;
   *WYHlapms = 0.0; *QWYHlapms = 0.0; *Qaddlapms = 0.0;
   *YWHlapms = 0.0; *YWHClapms = 0.0; *Raddlapms = 0.0;
   *addcnt = 0; *mulcnt = 0; *divcnt = 0; *sqrtcnt = 0;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      int colidx,nrows1;

      for(int L=0; L<szt; L++)  // L runs over the columns in one block
      {
         colidx = k*szt + L;              // index of the current column
         nrows1 = nrows - colidx - 1;     // #rows in Householder vector - 1
         if(nrows1 <= szt)
         {
            GPU_cmplx2_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                Arehi_h,Arelo_h,Aimhi_h,Aimlo_h,
                Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
                vrehi_h,vrelo_h,vimhi_h,vimlo_h,
                Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
                betahi_h,betalo_h,betahi_d,betalo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahi_h[L] == 0.0) && (betalo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_cmplx2_small_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                   Arehi_h,Arelo_h,Aimhi_h,Aimlo_h,
                   Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
                   Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
                   betahi_h,betalo_h,betahi_d,betalo_d,
                   tileRlapms,addcnt,mulcnt,verbose);
            }
         }
         else
         {
            GPU_cmplx2_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                Arehi_h,Arelo_h,Aimhi_h,Aimlo_h,
                Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
                vrehi_h,vrelo_h,vimhi_h,vimlo_h,
                Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
                betahi_h,betalo_h,betahi_d,betalo_d,
                sumshi_h,sumslo_h,sumshi_d,sumslo_d,
                &sigmahi_h,&sigmalo_h,sigmahi_d,sigmalo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahi_h[L] == 0.0) && (betalo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_cmplx2_medium_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                   Arehi_h,Arelo_h,Aimhi_h,Aimlo_h,
                   Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
                   Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
                   betahi_h,betalo_h,betahi_d,betalo_d,
                   RHdotvrehi_h,RHdotvrelo_h,RHdotvimhi_h,RHdotvimlo_h,
                   RHdotvrehi_d,RHdotvrelo_d,RHdotvimhi_d,RHdotvimlo_d,
                   bRHvrehi_h,bRHvrelo_h,bRHvimhi_h,bRHvimlo_h,
                   bRHvrehi_d,bRHvrelo_d,bRHvimhi_d,bRHvimlo_d,
                   RHvlapms,tileRlapms,addcnt,mulcnt,verbose);
            }
         }
      }
      GPU_cmplx2_medium_VB_to_W
         (nrows,szt,szt,k,
          Vrehi_h,Vrelo_h,Vimhi_h,Vimlo_h,Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
          Wrehi_h,Wrelo_h,Wimhi_h,Wimlo_h,Wrehi_d,Wrelo_d,Wimhi_d,Wimlo_d,
          WYTrehi_h,WYTrelo_h,WYTimhi_h,WYTimlo_h,
          WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d,
          betahi_h,betalo_h,betahi_d,betalo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);
/*
      GPU_cmplx2_small_WYH(nrows-k*szt,szt,
            Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
            Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
          WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d,
          WYTrehi_h,WYTrelo_h,WYTimhi_h,WYTimlo_h,WYHlapms,verbose);
 */
      GPU_cmplx2_small_QWYH(nrows,szt,k,
             Qrehi_d,   Qrelo_d,   Qimhi_d,   Qimlo_d,
           WYTrehi_d, WYTrelo_d, WYTimhi_d, WYTimlo_d,
          QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d,
          QWYTrehi_h,QWYTrelo_h,QWYTimhi_h,QWYTimlo_h,
             Qrehi_h,   Qrelo_h,   Qimhi_h,   Qimlo_h,
          QWYHlapms,addcnt,mulcnt,verbose);
      GPU_cmplx2_small_Qupdate
         (nrows,szt,k,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
          QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Qaddlapms,addcnt,verbose);
      if(k < nbt-1)                              // update R
      {
         GPU_cmplx2_small_YWH
            (nrows,szt,k,
               Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
               Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
             YWTrehi_d,YWTrelo_d,YWTimhi_d,YWTimlo_d,
             YWTrehi_h,YWTrelo_h,YWTimhi_h,YWTimlo_h,
             YWHlapms,addcnt,mulcnt,verbose);
         GPU_cmplx2_small_YWHC
            (nrows,ncols,szt,k,
              YWTrehi_d, YWTrelo_d, YWTimhi_d, YWTimlo_d,
                Arehi_d,   Arelo_d,   Aimhi_d,   Aimlo_d,
             YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d,
             YWTCrehi_h,YWTCrelo_h,YWTCimhi_h,YWTCimlo_h,
             YWHClapms,addcnt,mulcnt,verbose);
         GPU_cmplx2_small_R_add_YWHC
            (nrows,ncols,szt,k,
                Arehi_d,   Arelo_d,   Aimhi_d,   Aimlo_d,
             YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d,
                Arehi_h,   Arelo_h,   Aimhi_h,   Aimlo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qrehi_h,Qrehi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelo_h,Qrelo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhi_h,Qimhi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlo_h,Qimlo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qrehi[i][j] = Qrehi_h[ix];
         Qrelo[i][j] = Qrelo_h[ix];
         Qimhi[i][j] = Qimhi_h[ix];
         Qimlo[i][j] = Qimlo_h[ix++];
      }

   cudaMemcpy(Arehi_h,Arehi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelo_h,Arelo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhi_h,Aimhi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlo_h,Aimlo_d,sznum,cudaMemcpyDeviceToHost);
   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rrehi[i][j] = Arehi_h[j*nrows+i];
         Rrelo[i][j] = Arelo_h[j*nrows+i];
         Rimhi[i][j] = Aimhi_h[j*nrows+i];
         Rimlo[i][j] = Aimlo_h[j*nrows+i];
      }

   free(Arehi_h); free(Aimhi_h); free(Qrehi_h); free(Qimhi_h);
   free(Arelo_h); free(Aimlo_h); free(Qrelo_h); free(Qimlo_h);
   free(vrehi_h); free(vimhi_h); free(Vrehi_h); free(Vimhi_h);
   free(vrelo_h); free(vimlo_h); free(Vrelo_h); free(Vimlo_h);
   free(RHdotvrehi_h); free(RHdotvrelo_h);
   free(RHdotvimhi_h); free(RHdotvimlo_h);
   free(bRHvrehi_h); free(bRHvrelo_h); free(bRHvimhi_h); free(bRHvimlo_h);
   free(sumshi_h); free(sumslo_h);
   free(Wrehi_h); free(Wimhi_h); free(Wrelo_h); free(Wimlo_h);
   free(WYTrehi_h); free(QWYTrehi_h); free(YWTrehi_h); free(YWTCrehi_h);
   free(WYTrelo_h); free(QWYTrelo_h); free(YWTrelo_h); free(YWTCrelo_h);
   free(WYTimhi_h); free(QWYTimhi_h); free(YWTimhi_h); free(YWTCimhi_h);
   free(WYTimlo_h); free(QWYTimlo_h); free(YWTimlo_h); free(YWTCimlo_h);
}
