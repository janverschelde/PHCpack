/* The file dbl4_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl4_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_baqr_kernels.h"
#include "quad_double_functions.h"
#include "dbl_baqr_flopcounts.h"

using namespace std;

/*
__global__ void dbl4_small_house
 ( double *x0hihi, double *x0lohi, double *x0hilo, double *x0lolo,
   double *x1hihi, double *x1lohi, double *x1hilo, double *x1lolo,
   int dim, int dimLog2,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
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
   ddg_div(shvhi[j],shvlo[j],prdhi[0],prdlo[0],&acchi,&acclo);
   vhi[j+1] = acchi;
   vlo[j+1] = acclo;
   if(j == 0) vhi[0] = 1.0;
   if(j == 0) vlo[0] = 0.0;
}

__global__ void cmplx4_small_house
 ( double *x0rehihi, double *x0relohi, double *x0rehilo, double *x0relolo,
   double *x0imhihi, double *x0imlohi, double *x0imhilo, double *x0imlolo,
   double *x1rehihi, double *x1relohi, double *x1rehilo, double *x1relolo,
   double *x1imhihi, double *x1imlohi, double *x1imhilo, double *x1imlolo,
   int dim, int dimLog2,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
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
   if(j == 0)
   {
      vrehi[0] = 1.0;
      vrelo[0] = 0.0;
      vimhi[0] = 0.0;
      vimlo[0] = 0.0;
   }
}

__global__ void dbl4_large_sum_of_squares
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *sumshihi, double *sumshilo, double *sumslohi, double *sumslolo,
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

void flopcount_dbl4_large_sum_of_squares
 ( int nblocks, int szt, int sztLog2,
   long long int *add, long long int *mul )
{
   *add += nblocks*szt*sztLog2;
   *mul += nblocks*szt;
}

__global__ void cmplx4_large_sum_of_squares
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int dim, int BS, int BSLog2 )
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

void flopcount_cmplx4_large_sum_of_squares
 ( int nblocks, int szt, int sztLog2,
   long long int *add, long long int *mul )
{
   *add += nblocks*szt*(1+sztLog2);
   *mul += 2*nblocks*szt;
}

__global__ void dbl4_sum_accumulator
 ( double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int nbsums, int nbsumsLog2,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo )
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

void flopcount_dbl4_sum_accumulator
 ( int nbt, int nbtLog2, long long int *add )
{
   *add += nbt*nbtLog2;
}

__global__ void dbl4_normalize
 ( int dim, int szt,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *v0hihi, double *v0lohi, double *v0hilo, double *v0lolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo )
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

void flopcount_dbl4_normalize ( int nblocks, int szt, long long int *div )
{
   *div += nblocks*szt;
}

__global__ void cmplx4_normalize
 ( int dim, int szt,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *inv0rehihi, double *inv0relohi,
   double *inv0rehilo, double *inv0relolo,
   double *inv0imhihi, double *inv0imlohi,
   double *inv0imhilo, double *inv0imlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo )
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

void flopcount_cmplx4_normalize
 ( int nblocks, int szt, long long int *add, long long int *mul )
{
   *add += 2*nblocks*szt;
   *mul += 4*nblocks*szt;
}

__global__ void dbl4_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi,
   double *betahilo, double *betalolo )
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

__global__ void cmplx4_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
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

__global__ void dbl4_small_betaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo )
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

__global__ void dbl4_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *RTdotvhihi, double *RTdotvlohi,
   double *RTdotvhilo, double *RTdotvlolo )
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

__global__ void cmplx4_RHdotv
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

__global__ void dbl4_sum_betaRTdotv
 ( int nrows,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *RTdotvhihi, double *RTdotvlohi,
   double *RTdotvhilo, double *RTdotvlolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo )
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

__global__ void cmplx4_sum_betaRHdotv
 ( int nrows,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *RTdotvrehihi, double *RTdotvrelohi,
   double *RTdotvrehilo, double *RTdotvrelolo,
   double *RTdotvimhihi, double *RTdotvimlohi,
   double *RTdotvimhilo, double *RTdotvimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
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

__global__ void dbl4_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo )
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

__global__ void cmplx4_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
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

__global__ void dbl4_beta_times_V
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo )
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

__global__ void cmplx4_beta_times_V
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo )
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

__global__ void dbl4_initialize_WYT
 ( int dim, int szt,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
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

__global__ void cmplx4_initialize_WYH
 ( int dim, int szt,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYTrehihi, double *WYTrelohi,
   double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi,
   double *WYTimhilo, double *WYTimlolo )
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

__global__ void dbl4_update_WYT
 ( int dim, int szt,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
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

__global__ void cmplx4_update_WYH
 ( int dim, int szt,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo )
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

__global__ void dbl4_beta_next_W
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
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

__global__ void cmplx4_beta_next_W
 ( int nrows, int szt, double *Bhi, double *Blo,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo )
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

__global__ void dbl4_small_WYT
 ( int nrows, int szt,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
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

__global__ void cmplx4_small_WYH
 ( int nrows, int szt,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *Yrehihi, double *Yrelohi, double *Yrehilo, double *Yrelolo,
   double *Yimhihi, double *Yimlohi, double *Yimhilo, double *Yimlolo,
   double *WYTrehihi, double *WYTrelohi,
   double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi,
   double *WYTimhilo, double *WYTimlolo )
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

__global__ void dbl4_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihi, double *Qlohi, double *Qhilo, double *Qlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo,
   double *QWYThihi, double *QWYTlohi, double *QWYThilo, double *QWYTlolo )
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

__global__ void cmplx4_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihi, double *Qrelohi, double *Qrehilo, double *Qrelolo,
   double *Qimhihi, double *Qimlohi, double *Qimhilo, double *Qimlolo,
   double *WYTrehihi, double *WYTrelohi, double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi, double *WYTimhilo, double *WYTimlolo,
   double *QWYTrehihi, double *QWYTrelohi,
   double *QWYTrehilo, double *QWYTrelolo,
   double *QWYTimhihi, double *QWYTimlohi,
   double *QWYTimhilo, double *QWYTimlolo )
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

__global__ void dbl4_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWThihi, double *YWTlohi, double *YWThilo, double *YWTlolo,
   double *Chihi, double *Clohi, double *Chilo, double *Clolo,
   double *YWTChihi, double *YWTClohi, double *YWTChilo, double *YWTClolo )
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

__global__ void cmplx4_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehihi, double *YWTrelohi, double *YWTrehilo, double *YWTrelolo,
   double *YWTimhihi, double *YWTimlohi, double *YWTimhilo, double *YWTimlolo,
   double *Crehihi, double *Crelohi, double *Crehilo, double *Crelolo,
   double *Cimhihi, double *Cimlohi, double *Cimhilo, double *Cimlolo,
   double *YWTCrehihi, double *YWTCrelohi,
   double *YWTCrehilo, double *YWTCrelolo,
   double *YWTCimhihi, double *YWTCimlohi,
   double *YWTCimhilo, double *YWTCimlolo )
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

__global__ void dbl4_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihi, double *Qlohi, double *Qhilo, double *Qlolo,
   double *QWYThihi, double *QWYTlohi, double *QWYThilo, double *QWYTlolo )
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

__global__ void cmplx4_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihi, double *Qrelohi, double *Qrehilo, double *Qrelolo,
   double *Qimhihi, double *Qimlohi, double *Qimhilo, double *Qimlolo,
   double *QWYTrehihi, double *QWYTrelohi,
   double *QWYTrehilo, double *QWYTrelolo,
   double *QWYTimhihi, double *QWYTimlohi,
   double *QWYTimhilo, double *QWYTimlolo )
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

__global__ void dbl4_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *YWTChihi, double *YWTClohi, double *YWTChilo, double *YWTClolo )
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

__global__ void cmplx4_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *YWTCrehihi, double *YWTCrelohi,
   double *YWTCrehilo, double *YWTCrelolo,
   double *YWTCimhihi, double *YWTCimlohi,
   double *YWTCimhilo, double *YWTCimlolo )
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

void GPU_dbl4_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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

void GPU_cmplx4_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahi_h[L],&betahi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalo_h[L],&betalo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
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

void GPU_dbl4_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *sumshihi_h, double *sumslohi_h,
   double *sumshilo_h, double *sumslolo_h,
   double *sumshihi_d, double *sumslohi_d,
   double *sumshilo_d, double *sumslolo_d,
   double *sigmahihi_h, double *sigmalohi_h, 
   double *sigmahilo_h, double *sigmalolo_h, 
   double *sigmahihi_d, double *sigmalohi_d,
   double *sigmahilo_d, double *sigmalolo_d,
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
   flopcount_dbl2_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

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
   flopcount_dbl2_sum_accumulator(nblocks,nblLog2,add);

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
      flopcount_dbl2_normalize(nblocks,szt,div);
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

void GPU_cmplx4_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *sumshihi_h, double *sumslohi_h,
   double *sumshilo_h, double *sumslolo_h,
   double *sumshihi_d, double *sumslohi_d,
   double *sumshilo_d, double *sumslolo_d,
   double *sigmahihi_h, double *sigmalohi_h,
   double *sigmahilo_h, double *sigmalolo_h,
   double *sigmahihi_d, double *sigmalohi_d,
   double *sigmahilo_d, double *sigmalolo_d,
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
   flopcount_cmplx2_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

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
   flopcount_dbl2_sum_accumulator(nblocks,nblLog2,add);

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
      flopcount_cmplx2_normalize(nblocks,szt,add,mul);
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

void GPU_dbl4_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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

void GPU_cmplx4_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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

void GPU_dbl4_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *RTdotvhihi_h, double *RTdotvlohi_h,
   double *RTdotvhilo_h, double *RTdotvlolo_h,
   double *RTdotvhihi_d, double *RTdotvlohi_d,
   double *RTdotvhilo_d, double *RTdotvlolo_d,
   double *whihi_h, double *wlohi_h, double *whilo_h, double *wlolo_h,
   double *whihi_d, double *wlohi_d, double *whilo_d, double *wlolo_d,
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

void GPU_cmplx4_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *RHdotvrehihi_h, double *RHdotvrelohi_h,
   double *RHdotvrehilo_h, double *RHdotvrelolo_h,
   double *RHdotvimhihi_h, double *RHdotvimlohi_h,
   double *RHdotvimhilo_h, double *RHdotvimlolo_h,
   double *RHdotvrehihi_d, double *RHdotvrelohi_d,
   double *RHdotvrehilo_d, double *RHdotvrelolo_d,
   double *RHdotvimhihi_d, double *RHdotvimlohi_d,
   double *RHdotvimhilo_d, double *RHdotvimlolo_d,
   double *wrehihi_h, double *wrelohi_h, double *wrehilo_h, double *wrelolo_h,
   double *wimhihi_h, double *wimlohi_h, double *wimhilo_h, double *wimlolo_h,
   double *wrehihi_d, double *wrelohi_d, double *wrehilo_d, double *wrelolo_d,
   double *wimhihi_d, double *wimlohi_d, double *wimhilo_d, double *wimlolo_d,
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

void GPU_dbl4_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhihi_h, double *Vlohi_h, double *Vhilo_h, double *Vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *Whihi_h, double *Wlohi_h, double *Whilo_h, double *Wlolo_h,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *WYThihi_h, double *WYTlohi_h, double *WYThilo_h, double *WYTlolo_h,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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

void GPU_cmplx4_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehihi_h, double *Vrelohi_h, double *Vrehilo_h, double *Vrelolo_h,
   double *Vimhihi_h, double *Vimlohi_h, double *Vimhilo_h, double *Vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *Wrehihi_h, double *Wrelohi_h, double *Wrehilo_h, double *Wrelolo_h,
   double *Wimhihi_h, double *Wimlohi_h, double *Wimhilo_h, double *Wimlolo_h,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *WYHrehihi_h, double *WYHrelohi_h,
   double *WYHrehilo_h, double *WYHrelolo_h,
   double *WYHimhihi_h, double *WYHimlohi_h,
   double *WYHimhilo_h, double *WYHimlolo_h,
   double *WYHrehihi_d, double *WYHrelohi_d,
   double *WYHrehilo_d, double *WYHrelolo_d,
   double *WYHimhihi_d, double *WYHimlohi_d,
   double *WYHimhilo_d, double *WYHimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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

void GPU_dbl4_small_WYT
 ( int nrows, int szt,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *Yhihi_d, double *Ylohi_d, double *Yhilo_d, double *Ylolo_d,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *WYThihi_h, double *WYTlohi_h, double *WYThilo_h, double *WYTlolo_h,
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

void GPU_cmplx4_small_WYH
 ( int nrows, int szt,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *Yrehihi_d, double *Yrelohi_d, double *Yrehilo_d, double *Yrelolo_d,
   double *Yimhihi_d, double *Yimlohi_d, double *Yimhilo_d, double *Yimlolo_d,
   double *WYTrehihi_d, double *WYTrelohi_d,
   double *WYTrehilo_d, double *WYTrelolo_d,
   double *WYTimhihi_d, double *WYTimlohi_d,
   double *WYTimhilo_d, double *WYTimlolo_d,
   double *WYTrehihi_h, double *WYTrelohi_h,
   double *WYTrehilo_h, double *WYTrelolo_h,
   double *WYTimhihi_h, double *WYTimlohi_h,
   double *WYTimhilo_h, double *WYTimlolo_h,
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

void GPU_dbl4_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhihi_d, double *Ylohi_d, double *Yhilo_d, double *Ylolo_d,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *YWThihi_d, double *YWTlohi_d, double *YWThilo_d, double *YWTlolo_d,
   double *YWThihi_h, double *YWTlohi_h, double *YWThilo_h, double *YWTlolo_h,
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

void GPU_cmplx4_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehihi_d, double *Yrelohi_d, double *Yrehilo_d, double *Yrelolo_d,
   double *Yimhihi_d, double *Yimlohi_d, double *Yimhilo_d, double *Yimlolo_d,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *YWTrehihi_d, double *YWTrelohi_d,
   double *YWTrehilo_d, double *YWTrelolo_d,
   double *YWTimhihi_d, double *YWTimlohi_d,
   double *YWTimhilo_d, double *YWTimlolo_d,
   double *YWTrehihi_h, double *YWTrelohi_h,
   double *YWTrehilo_h, double *YWTrelolo_h,
   double *YWTimhihi_h, double *YWTimlohi_h,
   double *YWTimhilo_h, double *YWTimlolo_h,
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

void GPU_dbl4_small_QWYT
 ( int dim, int szt, int idx,
   double *Qhihi_d, double *Qlohi_d, double *Qhilo_d, double *Qlolo_d,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *QWYThihi_d, double *QWYTlohi_d,
   double *QWYThilo_d, double *QWYTlolo_d,
   double *QWYThihi_h, double *QWYTlohi_h,
   double *QWYThilo_h, double *QWYTlolo_h,
   double *Qhihi_h, double *Qlohi_h, double *Qhilo_h, double *Qlolo_h,
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

void GPU_cmplx4_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehihi_d, double *Qrelohi_d, double *Qrehilo_d, double *Qrelolo_d,
   double *Qimhihi_d, double *Qimlohi_d, double *Qimhilo_d, double *Qimlolo_d,
   double *WYTrehihi_d, double *WYTrelohi_d,
   double *WYTrehilo_d, double *WYTrelolo_d,
   double *WYTimhihi_d, double *WYTimlohi_d,
   double *WYTimhilo_d, double *WYTimlolo_d,
   double *QWYTrehihi_d, double *QWYTrelohi_d,
   double *QWYTrehilo_d, double *QWYTrelolo_d,
   double *QWYTimhihi_d, double *QWYTimlohi_d,
   double *QWYTimhilo_d, double *QWYTimlolo_d,
   double *QWYTrehihi_h, double *QWYTrelohi_h,
   double *QWYTrehilo_h, double *QWYTrelolo_h,
   double *QWYTimhihi_h, double *QWYTimlohi_h,
   double *QWYTimhilo_h, double *QWYTimlolo_h,
   double *Qrehihi_h, double *Qrelohi_h, double *Qrehilo_h, double *Qrelolo_h,
   double *Qimhihi_h, double *Qimlohi_h, double *Qimhilo_h, double *Qimlolo_h,
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

void GPU_dbl4_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWThihi_d, double *YWTlohi_d, double *YWThilo_d, double *YWTlolo_d,
   double *Chihi_d, double *Clohi_d, double *Chilo_d, double *Clolo_d,
   double *YWTChihi_d, double *YWTClohi_d,
   double *YWTChilo_d, double *YWTClolo_d,
   double *YWTChihi_h, double *YWTClohi_h,
   double *YWTChilo_h, double *YWTClolo_h,
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

void GPU_cmplx4_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehihi_d, double *YWTrelohi_d,
   double *YWTrehilo_d, double *YWTrelolo_d,
   double *YWTimhihi_d, double *YWTimlohi_d,
   double *YWTimhilo_d, double *YWTimlolo_d,
   double *Crehihi_d, double *Crelohi_d, double *Crehilo_d, double *Crelolo_d,
   double *Cimhihi_d, double *Cimlohi_d, double *Cimhilo_d, double *Cimlolo_d,
   double *YWTCrehihi_d, double *YWTCrelohi_d,
   double *YWTCrehilo_d, double *YWTCrelolo_d,
   double *YWTCimhihi_d, double *YWTCimlohi_d,
   double *YWTCimhilo_d, double *YWTCimlolo_d,
   double *YWTCrehihi_h, double *YWTCrelohi_h,
   double *YWTCrehilo_h, double *YWTCrelolo_h,
   double *YWTCimhihi_h, double *YWTCimlohi_h,
   double *YWTCimhilo_h, double *YWTCimlolo_h,
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

void GPU_dbl4_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhihi_d, double *Qlohi_d, double *Qhilo_d, double *Qlolo_d,
   double *QWYThihi_d, double *QWYTlohi_d,
   double *QWYThilo_d, double *QWYTlolo_d,
   double *Qhihi_h, double *Qlohi_h, double *Qhilo_h, double *Qlolo_h,
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

void GPU_cmplx4_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehihi_d, double *Qrelohi_d, double *Qrehilo_d, double *Qrelolo_d,
   double *Qimhihi_d, double *Qimlohi_d, double *Qimhilo_d, double *Qimlolo_d,
   double *QWYTrehihi_d, double *QWYTrelohi_d,
   double *QWYTrehilo_d, double *QWYTrelolo_d,
   double *QWYTimhihi_d, double *QWYTimlohi_d,
   double *QWYTimhilo_d, double *QWYTimlolo_d,
   double *Qrehihi_h, double *Qrelohi_h, double *Qrehilo_h, double *Qrelolo_h,
   double *Qimhihi_h, double *Qimlohi_h, double *Qimhilo_h, double *Qimlolo_h,
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

void GPU_dbl4_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rhihi_d, double *Rlohi_d, double *Rhilo_d, double *Rlolo_d,
   double *YWTChihi_d, double *YWTClohi_d,
   double *YWTChilo_d, double *YWTClolo_d,
   double *Rhihi_h, double *Rlohi_h, double *Rhilo_h, double *Rlolo_h,
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

void GPU_cmplx4_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehihi_d, double *Rrelohi_d, double *Rrehilo_d, double *Rrelolo_d,
   double *Rimhihi_d, double *Rimlohi_d, double *Rimhilo_d, double *Rimlolo_d,
   double *YWTCrehihi_d, double *YWTCrelohi_d,
   double *YWTCrehilo_d, double *YWTCrelolo_d,
   double *YWTCimhihi_d, double *YWTCimlohi_d,
   double *YWTCimhilo_d, double *YWTCimlolo_d,
   double *Rrehihi_h, double *Rrelohi_h, double *Rrehilo_h, double *Rrelolo_h,
   double *Rimhihi_h, double *Rimlohi_h, double *Rimhilo_h, double *Rimlolo_h,
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
*/

void GPU_dbl4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;          // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Ahihi_h = new double[dim];    // A on the host
   double *Alohi_h = new double[dim]; 
   double *Ahilo_h = new double[dim];
   double *Alolo_h = new double[dim]; 
   double *Ahihi_d;                      // A on the device
   double *Alohi_d; 
   double *Ahilo_d; 
   double *Alolo_d; 
   double *Qhihi_h = new double[nrows2]; // Q on the host
   double *Qlohi_h = new double[nrows2]; 
   double *Qhilo_h = new double[nrows2]; 
   double *Qlolo_h = new double[nrows2]; 
   double *Qhihi_d;                      // Q on the device
   double *Qlohi_d;
   double *Qhilo_d;
   double *Qlolo_d;
   double *vhihi_h = new double[nrows];  // Householder vector
   double *vlohi_h = new double[nrows];
   double *vhilo_h = new double[nrows];
   double *vlolo_h = new double[nrows];
   double *betahihi_h = new double[szt]; //  beta on the host
   double *betalohi_h = new double[szt]; 
   double *betahilo_h = new double[szt]; 
   double *betalolo_h = new double[szt]; 
   double *betahihi_d;                   // beta on the device
   double *betalohi_d;
   double *betahilo_d;
   double *betalolo_d;
   double *Vhihi_h = new double[nrows*szt]; // V matrix
   double *Vlohi_h = new double[nrows*szt];
   double *Vhilo_h = new double[nrows*szt];
   double *Vlolo_h = new double[nrows*szt];
   double *Vhihi_d;                         // V on the device
   double *Vlohi_d;
   double *Vhilo_d;
   double *Vlolo_d;
   double *Whihi_h = new double[nrows*szt]; // W on the host
   double *Wlohi_h = new double[nrows*szt];
   double *Whilo_h = new double[nrows*szt];
   double *Wlolo_h = new double[nrows*szt];
   double *Whihi_d;                         // W on the device
   double *Wlohi_d;
   double *Whilo_d;
   double *Wlolo_d;
   double *WYThihi_h = new double[nrows2];  // W*Y^T 
   double *WYTlohi_h = new double[nrows2];
   double *WYThilo_h = new double[nrows2];
   double *WYTlolo_h = new double[nrows2];
   double *WYThihi_d;                       // WYT on the device
   double *WYTlohi_d;
   double *WYThilo_d;
   double *WYTlolo_d;
   double *YWThihi_h = new double[nrows2];  // Y*W^T
   double *YWTlohi_h = new double[nrows2];
   double *YWThilo_h = new double[nrows2];
   double *YWTlolo_h = new double[nrows2];
   double *YWThihi_d;                       // YWT on the device
   double *YWTlohi_d;
   double *YWThilo_d;
   double *YWTlolo_d;
   double *QWYThihi_h = new double[nrows2]; // Q*WY^T
   double *QWYTlohi_h = new double[nrows2];
   double *QWYThilo_h = new double[nrows2];
   double *QWYTlolo_h = new double[nrows2];
   double *QWYThihi_d;                      // QWYT on the device
   double *QWYTlohi_d;
   double *QWYThilo_d;
   double *QWYTlolo_d;
   double *YWTChihi_h = new double[dim];    // YWT*C on the host
   double *YWTClohi_h = new double[dim];
   double *YWTChilo_h = new double[dim];
   double *YWTClolo_h = new double[dim];
   double *YWTChihi_d;                      // YWTC on the device
   double *YWTClohi_d;
   double *YWTChilo_d;
   double *YWTClolo_d;
   double *RTdotvhihi_h = new double[nrows2]; // R^T dotted with v
   double *RTdotvlohi_h = new double[nrows2];
   double *RTdotvhilo_h = new double[nrows2];
   double *RTdotvlolo_h = new double[nrows2];
   double *RTdotvhihi_d;                      // RTdotv on the device
   double *RTdotvlohi_d;
   double *RTdotvhilo_d;
   double *RTdotvlolo_d;
   double *bRTvhihi_h = new double[nrows];  // beta*R^T*v
   double *bRTvlohi_h = new double[nrows];
   double *bRTvhilo_h = new double[nrows];
   double *bRTvlolo_h = new double[nrows];
   double *bRTvhihi_d;                      // bRTv on the device
   double *bRTvlohi_d;
   double *bRTvhilo_d;
   double *bRTvlolo_d;
   double *sumshihi_h = new double[nrows];  // subsums for large house
   double *sumslohi_h = new double[nrows];
   double *sumshilo_h = new double[nrows];
   double *sumslolo_h = new double[nrows];
   double *sumshihi_d;                      // sums on the device
   double *sumslohi_d;
   double *sumshilo_d;
   double *sumslolo_d;
   double sigmahihi_h,sigmalohi_h,sigmahilo_h,sigmalolo_h;
   double *sigmahihi_d;                     // sigma on the device
   double *sigmalohi_d;
   double *sigmahilo_d;
   double *sigmalolo_d;

   int ix = 0;                          // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Ahihi_h[ix]   = Ahihi[i][j];
         Alohi_h[ix]   = Alohi[i][j];
         Ahilo_h[ix]   = Ahilo[i][j];
         Alolo_h[ix++] = Alolo[i][j];
      }

   ix = 0;                              // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qhihi_h[ix]   = 1.0;
            Qlohi_h[ix]   = 0.0;
            Qhilo_h[ix]   = 0.0;
            Qlolo_h[ix++] = 0.0;
         }
         else
         {
            Qhihi_h[ix]   = 0.0;
            Qlohi_h[ix]   = 0.0;
            Qhilo_h[ix]   = 0.0;
            Qlolo_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Ahihi_d,sznum);
   cudaMalloc((void**)&Alohi_d,sznum);
   cudaMalloc((void**)&Ahilo_d,sznum);
   cudaMalloc((void**)&Alolo_d,sznum);
   cudaMemcpy(Ahihi_d,Ahihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alohi_d,Alohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahilo_d,Ahilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alolo_d,Alolo_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&betahihi_d,szbeta);
   cudaMalloc((void**)&betalohi_d,szbeta);
   cudaMalloc((void**)&betahilo_d,szbeta);
   cudaMalloc((void**)&betalolo_d,szbeta);

   for(int i=0; i<szt; i++)
   {
      betahihi_h[i] = 0.0;
      betalohi_h[i] = 0.0;
      betahilo_h[i] = 0.0;
      betalolo_h[i] = 0.0;
   }
   cudaMemcpy(betahihi_d,betahihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohi_d,betalohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilo_d,betahilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolo_d,betalolo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);  // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vhihi_d,szVandW + szpad); // padding only in allocation
   cudaMalloc((void**)&Vlohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vhilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vlolo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vhihi_h[ix] = 0.0; 
      Vlohi_h[ix] = 0.0; 
      Vhilo_h[ix] = 0.0; 
      Vlolo_h[ix++] = 0.0; 
   }
   Vhihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vhihi_d,Vhihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlohi_d,Vlohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhilo_d,Vhilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlolo_d,Vlolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Whihi_d,szVandW + szpad); // padding only in allocation
   cudaMalloc((void**)&Wlolo_d,szVandW + szpad); 
   cudaMalloc((void**)&Whilo_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlolo_d,szVandW + szpad); 

   cudaMalloc((void**)&RTdotvhihi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlohi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhilo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlolo_d,szVandW + szpad);
   cudaMalloc((void**)&bRTvhihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlolo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihi_d,szhouse);
   cudaMalloc((void**)&sumslohi_d,szhouse);
   cudaMalloc((void**)&sumshilo_d,szhouse);
   cudaMalloc((void**)&sumslolo_d,szhouse);
   cudaMalloc((void**)&sigmahihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalolo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYThihi_d,szWYT + szpad); // padding for W*Y^T product
   cudaMalloc((void**)&WYTlohi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThilo_d,szWYT + szpad); 
   cudaMalloc((void**)&WYTlolo_d,szWYT + szpad); 
   cudaMalloc((void**)&Qhihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qhilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qlolo_d,szWYT + szpad);
   cudaMemcpy(Qhihi_d,Qhihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlohi_d,Qlohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhilo_d,Qhilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlolo_d,Qlolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYThihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlolo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWThihi_d,szYWT + szpad); // padding for Y*W^T product
   cudaMalloc((void**)&YWTlohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWThilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTChihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTChilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTClolo_d,sznum + szpad);

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
            GPU_dbl4_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Ahihi_h,   Alohi_h,   Ahilo_h,   Alolo_h,
                   Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
                   vhihi_h,   vlohi_h,   vhilo_h,   vlolo_h,
                   Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
                betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            GPU_dbl4_small_leftRupdate
               (nrows,ncols,szt,colidx,k,L,
                   Ahihi_h,   Alohi_h,   Ahilo_h,   Alolo_h,
                   Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
                   Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
                betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                tileRlapms,addcnt,mulcnt,verbose);
         }
         else
         {
            GPU_dbl4_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                     Ahihi_h,     Alohi_h,     Ahilo_h,     Alolo_h,
                     Ahihi_d,     Alohi_d,     Ahilo_d,     Alolo_d,
                     vhihi_h,     vlohi_h,     vhilo_h,     vlolo_h,
                     Vhihi_d,     Vlohi_d,     Vhilo_d,     Vlolo_d,
                  betahihi_h,  betalohi_h,  betahilo_h,  betalolo_h,
                  betahihi_d,  betalohi_d,  betahilo_d,  betalolo_d,
                  sumshihi_h,  sumslohi_h,  sumshilo_h,  sumslolo_h,
                  sumshihi_d,  sumslohi_d,  sumshilo_d,  sumslolo_d,
                &sigmahihi_h,&sigmalohi_h,&sigmahilo_h,&sigmalolo_h,
                 sigmahihi_d, sigmalohi_d, sigmahilo_d, sigmalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            GPU_dbl4_medium_leftRupdate
               (nrows,ncols,szt,colidx,k,L,
                     Ahihi_h,     Alohi_h,     Ahilo_h,     Alolo_h,
                     Ahihi_d,     Alohi_d,     Ahilo_d,     Alolo_d,
                     Vhihi_d,     Vlohi_d,     Vhilo_d,     Vlolo_d,
                  betahihi_h,  betalohi_h,  betahilo_h,  betalolo_h,
                  betahihi_d,  betalohi_d,  betahilo_d,  betalolo_d,
                RTdotvhihi_h,RTdotvlohi_h,RTdotvhilo_h,RTdotvlolo_h,
                RTdotvhihi_d,RTdotvlohi_d,RTdotvhilo_d,RTdotvlolo_d,
                  bRTvhihi_h,  bRTvlohi_h,  bRTvhilo_h,  bRTvlolo_h,
                  bRTvhihi_d,  bRTvlohi_d,  bRTvhilo_d,  bRTvlolo_d,
                RTvlapms,tileRlapms,addcnt,mulcnt,verbose);
         }
      }
      GPU_dbl4_medium_VB_to_W
         (nrows,szt,szt,k,
             Vhihi_h,   Vlohi_h,   Vhilo_h,   Vlolo_h,
             Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
             Whihi_h,   Wlohi_h,   Whilo_h,   Wlolo_h,
             Whihi_d,   Wlohi_d,   Whilo_d,   Wlolo_d,
           WYThihi_h, WYTlohi_h, WYThilo_h, WYTlolo_h,
           WYThihi_d, WYTlohi_d, WYThilo_d, WYTlolo_d,
          betahihi_h,betalohi_h,betahilo_h,betalolo_h,
          betahihi_d,betalohi_d,betahilo_d,betalolo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);
/*
      GPU_dbl2_small_WYT
         (nrows-k*szt,szt,Whi_d,Wlo_d,Vhi_d,Vlo_d,WYThi_d,WYTlo_d,
          WYThi_h,WYTlo_h,WYTlapms,verbose);
 */
      GPU_dbl4_small_QWYT
         (nrows,szt,k,
             Qhihi_d,   Qlohi_d,   Qhilo_d,   Qlolo_d,
           WYThihi_d, WYTlohi_d, WYThilo_d, WYTlolo_d,
          QWYThihi_d,QWYTlohi_d,QWYThilo_d,QWYTlolo_d,
          QWYThihi_h,QWYTlohi_h,QWYThilo_h,QWYTlolo_h,
             Qhihi_h,   Qlohi_h,   Qhilo_h,   Qlolo_h,
          QWYTlapms,addcnt,mulcnt,verbose);

      GPU_dbl4_small_Qupdate
         (nrows,szt,k,
             Qhihi_d,   Qlohi_d,   Qhilo_d,   Qlolo_d,
          QWYThihi_d,QWYTlohi_d,QWYThilo_d,QWYTlolo_d,
             Qhihi_h,   Qlohi_h,   Qhilo_h,   Qlolo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                                           // update R
      {
         GPU_dbl4_small_YWT
            (nrows,szt,k,
               Vhihi_d,  Vlohi_d,  Vhilo_d,  Vlolo_d,
               Whihi_d,  Wlohi_d,  Whilo_d,  Wlolo_d,
             YWThihi_d,YWTlohi_d,YWThilo_d,YWTlolo_d,
             YWThihi_h,YWTlohi_h,YWThilo_h,YWTlolo_h,
             YWTlapms,addcnt,mulcnt,verbose);

         GPU_dbl4_small_YWTC
            (nrows,ncols,szt,k,
              YWThihi_d, YWTlohi_d, YWThilo_d, YWTlolo_d,
                Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
             YWTChihi_d,YWTClohi_d,YWTChilo_d,YWTClolo_d,
             YWTChihi_h,YWTClohi_h,YWTChilo_h,YWTClolo_h,
             YWTClapms,addcnt,mulcnt,verbose);

         GPU_dbl4_small_R_add_YWTC
            (nrows,ncols,szt,k,
                Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
             YWTChihi_d,YWTClohi_d,YWTChilo_d,YWTClolo_d,
                Ahihi_h,   Alohi_h,   Ahilo_h,   Alolo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qhihi_h,Qhihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlohi_h,Qlohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhilo_h,Qhilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlolo_h,Qlolo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qhihi[i][j] = Qhihi_h[ix];
         Qlohi[i][j] = Qlohi_h[ix];
         Qhilo[i][j] = Qhilo_h[ix];
         Qlolo[i][j] = Qlolo_h[ix++];
      }

   cudaMemcpy(Ahihi_h,Ahihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alohi_h,Alohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahilo_h,Ahilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alolo_h,Alolo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rhihi[i][j] = Ahihi_h[j*nrows+i];
         Rlohi[i][j] = Alohi_h[j*nrows+i];
         Rhilo[i][j] = Ahilo_h[j*nrows+i];
         Rlolo[i][j] = Alolo_h[j*nrows+i];
      }

   free(Ahihi_h); free(Alohi_h); free(Ahilo_h); free(Alolo_h);
   free(Qhihi_h); free(Qlohi_h); free(Qhilo_h); free(Qlolo_h); 
   free(vhihi_h); free(vlohi_h); free(vhilo_h); free(vlolo_h);
   free(Vhihi_h); free(Vlohi_h); free(Vhilo_h); free(Vlolo_h);
   free(Whihi_h); free(Wlohi_h); free(Whilo_h); free(Wlolo_h);
   free(sumshihi_h); free(sumslohi_h); free(sumshilo_h); free(sumslolo_h);

   free(RTdotvhihi_h); free(RTdotvlohi_h);
   free(RTdotvhilo_h); free(RTdotvlolo_h);
   free(bRTvhihi_h); free(bRTvlohi_h); free(bRTvhilo_h); free(bRTvlolo_h);
   free(WYThihi_h); free(QWYThihi_h); free(WYThilo_h); free(QWYThilo_h);
   free(YWThihi_h); free(YWTChihi_h); free(YWThilo_h); free(YWTChilo_h);
   free(WYTlohi_h); free(QWYTlohi_h); free(WYTlolo_h); free(QWYTlolo_h);
   free(YWTlohi_h); free(YWTClohi_h); free(YWTlolo_h); free(YWTClolo_h);
}


void GPU_cmplx4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;        // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Arehihi_h = new double[dim];  // the real parts of A
   double *Arelohi_h = new double[dim];
   double *Arehilo_h = new double[dim];
   double *Arelolo_h = new double[dim];
   double *Aimhihi_h = new double[dim];  // the imaginary parts of A
   double *Aimlohi_h = new double[dim]; 
   double *Aimhilo_h = new double[dim]; 
   double *Aimlolo_h = new double[dim]; 
   double *Arehihi_d;                    // Are on the device
   double *Arelohi_d;
   double *Arehilo_d;
   double *Arelolo_d;
   double *Aimhihi_d;                    // Aim on the device
   double *Aimlohi_d;
   double *Aimhilo_d;
   double *Aimlolo_d;
   double *Qrehihi_h = new double[nrows2]; // real parts of Q 
   double *Qrelohi_h = new double[nrows2];
   double *Qrehilo_h = new double[nrows2];
   double *Qrelolo_h = new double[nrows2];
   double *Qimhihi_h = new double[nrows2]; // imaginary parts of Q
   double *Qimlohi_h = new double[nrows2];
   double *Qimhilo_h = new double[nrows2];
   double *Qimlolo_h = new double[nrows2];
   double *Qrehihi_d;                      // Qre on the device
   double *Qrelohi_d;
   double *Qrehilo_d;
   double *Qrelolo_d;
   double *Qimhihi_d;                      // Qim on the device
   double *Qimlohi_d;
   double *Qimhilo_d;
   double *Qimlolo_d;
   double *vrehihi_h = new double[nrows];  // real parts of Householder v
   double *vrelohi_h = new double[nrows];
   double *vrehilo_h = new double[nrows];
   double *vrelolo_h = new double[nrows];
   double *vimhihi_h = new double[nrows];  // imaginary parts of Householder v
   double *vimlohi_h = new double[nrows]; 
   double *vimhilo_h = new double[nrows]; 
   double *vimlolo_h = new double[nrows]; 
   double *betahihi_h = new double[szt+1]; // beta
   double *betalohi_h = new double[szt+1]; 
   double *betahilo_h = new double[szt+1]; 
   double *betalolo_h = new double[szt+1]; 
   double *betahihi_d;                     // beta on the device
   double *betalohi_d;
   double *betahilo_d;
   double *betalolo_d;
   double *Vrehihi_h = new double[nrows*szt]; // real parts of V
   double *Vrelohi_h = new double[nrows*szt]; 
   double *Vrehilo_h = new double[nrows*szt]; 
   double *Vrelolo_h = new double[nrows*szt]; 
   double *Vimhihi_h = new double[nrows*szt]; // imaginary parts of V
   double *Vimlohi_h = new double[nrows*szt];
   double *Vimhilo_h = new double[nrows*szt];
   double *Vimlolo_h = new double[nrows*szt];
   double *Vrehihi_d;                         // Vre on device
   double *Vrelohi_d;
   double *Vrehilo_d;
   double *Vrelolo_d;
   double *Vimhihi_d;                         // Vim on device
   double *Vimlohi_d;
   double *Vimhilo_d;
   double *Vimlolo_d;
   double *Wrehihi_h = new double[nrows*szt]; // real parts of W
   double *Wrelohi_h = new double[nrows*szt];
   double *Wrehilo_h = new double[nrows*szt];
   double *Wrelolo_h = new double[nrows*szt];
   double *Wimhihi_h = new double[nrows*szt]; // imaginary parts of W
   double *Wimlohi_h = new double[nrows*szt];
   double *Wimhilo_h = new double[nrows*szt];
   double *Wimlolo_h = new double[nrows*szt];
   double *Wrehihi_d;                         // Wre on the device
   double *Wrelohi_d;
   double *Wrehilo_d;
   double *Wrelolo_d;
   double *Wimhihi_d;                         // Wim on the device
   double *Wimlohi_d;
   double *Wimhilo_d;
   double *Wimlolo_d;
   double *WYTrehihi_h = new double[nrows2];  // real parts of W*Y^T
   double *WYTrelohi_h = new double[nrows2];
   double *WYTrehilo_h = new double[nrows2];
   double *WYTrelolo_h = new double[nrows2];
   double *WYTimhihi_h = new double[nrows2];  // imaginary parts of W*Y^T
   double *WYTimlohi_h = new double[nrows2];
   double *WYTimhilo_h = new double[nrows2];
   double *WYTimlolo_h = new double[nrows2];
   double *WYTrehihi_d;                       // WYTre on the device 
   double *WYTrelohi_d;
   double *WYTrehilo_d;
   double *WYTrelolo_d;
   double *WYTimhihi_d;                       // WYTim on the device
   double *WYTimlohi_d;
   double *WYTimhilo_d;
   double *WYTimlolo_d;
   double *YWTrehihi_h = new double[nrows2];  // real parts of Y*W^T
   double *YWTrelohi_h = new double[nrows2];
   double *YWTrehilo_h = new double[nrows2];
   double *YWTrelolo_h = new double[nrows2];
   double *YWTimhihi_h = new double[nrows2];  // imaginary parts of Y*W^T
   double *YWTimlohi_h = new double[nrows2];
   double *YWTimhilo_h = new double[nrows2];
   double *YWTimlolo_h = new double[nrows2];
   double *YWTrehihi_d;                       // YWTre on the device
   double *YWTrelohi_d;
   double *YWTrehilo_d;
   double *YWTrelolo_d;
   double *YWTimhihi_d;                       // YWTim on the device
   double *YWTimlohi_d; 
   double *YWTimhilo_d; 
   double *YWTimlolo_d; 
   double *QWYTrehihi_h = new double[nrows2]; // real parts of Q*WY^T
   double *QWYTrelohi_h = new double[nrows2];
   double *QWYTrehilo_h = new double[nrows2];
   double *QWYTrelolo_h = new double[nrows2];
   double *QWYTimhihi_h = new double[nrows2]; // imaginary parts of Q*WY^T
   double *QWYTimlohi_h = new double[nrows2];
   double *QWYTimhilo_h = new double[nrows2];
   double *QWYTimlolo_h = new double[nrows2];
   double *QWYTrehihi_d;                      // QWYTre on the device
   double *QWYTrelohi_d;
   double *QWYTrehilo_d;
   double *QWYTrelolo_d;
   double *QWYTimhihi_d;                      // QWYTim on the device
   double *QWYTimlohi_d;
   double *QWYTimhilo_d;
   double *QWYTimlolo_d;
   double *YWTCrehihi_h = new double[dim];    // real parts of YWT*C
   double *YWTCrelohi_h = new double[dim]; 
   double *YWTCrehilo_h = new double[dim]; 
   double *YWTCrelolo_h = new double[dim]; 
   double *YWTCimhihi_h = new double[dim];    // imaginary parts of YWT*C
   double *YWTCimlohi_h = new double[dim];
   double *YWTCimhilo_h = new double[dim];
   double *YWTCimlolo_h = new double[dim];
   double *YWTCrehihi_d;                      // YWTCre on the device
   double *YWTCrelohi_d;
   double *YWTCrehilo_d;
   double *YWTCrelolo_d;
   double *YWTCimhihi_d;                      // YWTCim on the device
   double *YWTCimlohi_d;
   double *YWTCimhilo_d;
   double *YWTCimlolo_d;
   double *RHdotvrehihi_h = new double[nrows2]; // real R^H dotted with v
   double *RHdotvrelohi_h = new double[nrows2]; 
   double *RHdotvrehilo_h = new double[nrows2]; 
   double *RHdotvrelolo_h = new double[nrows2]; 
   double *RHdotvimhihi_h = new double[nrows2]; // imaginary R^H dotted with v
   double *RHdotvimlohi_h = new double[nrows2]; 
   double *RHdotvimhilo_h = new double[nrows2]; 
   double *RHdotvimlolo_h = new double[nrows2]; 
   double *RHdotvrehihi_d;                      // RHdotvre on the device
   double *RHdotvrelohi_d;
   double *RHdotvrehilo_d;
   double *RHdotvrelolo_d;
   double *RHdotvimhihi_d;                      // RHdotvim on the device
   double *RHdotvimlohi_d;
   double *RHdotvimhilo_d;
   double *RHdotvimlolo_d;
   double *bRHvrehihi_h = new double[nrows];  // real parts of beta*R^H*v
   double *bRHvrelohi_h = new double[nrows];
   double *bRHvrehilo_h = new double[nrows];
   double *bRHvrelolo_h = new double[nrows];
   double *bRHvimhihi_h = new double[nrows];  // imaginary parts of beta*R^H*v
   double *bRHvimlohi_h = new double[nrows];
   double *bRHvimhilo_h = new double[nrows];
   double *bRHvimlolo_h = new double[nrows];
   double *bRHvrehihi_d;                      // bRHvre on the device
   double *bRHvrelohi_d;
   double *bRHvrehilo_d;
   double *bRHvrelolo_d;
   double *bRHvimhihi_d;                      // bRHvim on the device
   double *bRHvimlohi_d; 
   double *bRHvimhilo_d; 
   double *bRHvimlolo_d; 
   double *sumshihi_h = new double[nrows];  // subsums for large house
   double *sumslohi_h = new double[nrows]; 
   double *sumshilo_h = new double[nrows]; 
   double *sumslolo_h = new double[nrows]; 
   double *sumshihi_d;                      // sums on the device
   double *sumslohi_d;
   double *sumshilo_d;
   double *sumslolo_d;
   double sigmahihi_h,sigmalohi_h,sigmahilo_h,sigmalolo_h;
   double *sigmahihi_d;                     // sigma on the device
   double *sigmalohi_d;
   double *sigmahilo_d;
   double *sigmalolo_d;

   int ix = 0;                            // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Arehihi_h[ix]   = Arehihi[i][j];
         Arelohi_h[ix]   = Arelohi[i][j];
         Arehilo_h[ix]   = Arehilo[i][j];
         Arelolo_h[ix]   = Arelolo[i][j];
         Aimhihi_h[ix]   = Aimhihi[i][j];
         Aimlohi_h[ix]   = Aimlohi[i][j];
         Aimhilo_h[ix]   = Aimhilo[i][j];
         Aimlolo_h[ix++] = Aimlolo[i][j];
      }

   ix = 0;                                // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qrehihi_h[ix]   = 1.0;
            Qrelohi_h[ix]   = 0.0;
            Qrehilo_h[ix]   = 0.0;
            Qrelolo_h[ix]   = 0.0;
            Qimhihi_h[ix]   = 0.0;
            Qimlohi_h[ix]   = 0.0;
            Qimhilo_h[ix]   = 0.0;
            Qimlolo_h[ix++] = 0.0;
         }
         else
         {
            Qrehihi_h[ix]   = 0.0;
            Qrelohi_h[ix]   = 0.0;
            Qrehilo_h[ix]   = 0.0;
            Qrelolo_h[ix]   = 0.0;
            Qimhihi_h[ix]   = 0.0;
            Qimlohi_h[ix]   = 0.0;
            Qimhilo_h[ix]   = 0.0;
            Qimlolo_h[ix++] = 0.0;
         }
         // cout << "Q[" << ix-1 << "] : "
         //      << Qre_h[ix-1] << "  " << Qim_h[ix-1] << endl;
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Arehihi_d,sznum);
   cudaMalloc((void**)&Arelohi_d,sznum);
   cudaMalloc((void**)&Arehilo_d,sznum);
   cudaMalloc((void**)&Arelolo_d,sznum);
   cudaMalloc((void**)&Aimhihi_d,sznum);
   cudaMalloc((void**)&Aimlohi_d,sznum);
   cudaMalloc((void**)&Aimhilo_d,sznum);
   cudaMalloc((void**)&Aimlolo_d,sznum);
   cudaMemcpy(Arehihi_d,Arehihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelohi_d,Arelohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arehilo_d,Arehilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelolo_d,Arelolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhihi_d,Aimhihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlohi_d,Aimlohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhilo_d,Aimhilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlolo_d,Aimlolo_h,sznum,cudaMemcpyHostToDevice);
   // allocate one extra beta for use in the cmplx2_normalize
   const size_t szbeta = (szt+1)*sizeof(double);
   cudaMalloc((void**)&betahihi_d,szbeta);
   cudaMalloc((void**)&betalohi_d,szbeta);
   cudaMalloc((void**)&betahilo_d,szbeta);
   cudaMalloc((void**)&betalolo_d,szbeta);

   for(int i=0; i<szt; i++)
   {
      betahihi_h[i] = 0.0;
      betalohi_h[i] = 0.0;
      betahilo_h[i] = 0.0;
      betalolo_h[i] = 0.0;
   }
   cudaMemcpy(betahihi_d,betahihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohi_d,betalohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilo_d,betahilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolo_d,betalolo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);    // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vrehihi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Vrelohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vrehilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vrelolo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlolo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vrehihi_h[ix]   = 0.0; 
      Vrelohi_h[ix]   = 0.0; 
      Vrehilo_h[ix]   = 0.0; 
      Vrelolo_h[ix]   = 0.0; 
      Vimhihi_h[ix]   = 0.0; 
      Vimlohi_h[ix]   = 0.0; 
      Vimhilo_h[ix]   = 0.0; 
      Vimlolo_h[ix++] = 0.0; 
   }
   Vrehihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vrehihi_d,Vrehihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelohi_d,Vrelohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrehilo_d,Vrehilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelolo_d,Vrelolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhihi_d,Vimhihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlohi_d,Vimlohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhilo_d,Vimhilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlolo_d,Vimlolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Wrehihi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Wrelohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wrehilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wrelolo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhihi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlolo_d,szVandW + szpad);

   cudaMalloc((void**)&RHdotvrehihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrehilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelolo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlolo_d,szVandW + szpad);
   cudaMalloc((void**)&bRHvrehihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrehilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelolo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlolo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihi_d,szhouse);
   cudaMalloc((void**)&sumslohi_d,szhouse);
   cudaMalloc((void**)&sumshilo_d,szhouse);
   cudaMalloc((void**)&sumslolo_d,szhouse);
   cudaMalloc((void**)&sigmahihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalolo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYTrehihi_d,szWYT + szpad); // padding for W*Y^T 
   cudaMalloc((void**)&WYTrelohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrehilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrelolo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhihi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlolo_d,szWYT + szpad);
   cudaMemcpy(Qrehihi_d,Qrehihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelohi_d,Qrelohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrehilo_d,Qrehilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelolo_d,Qrelolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhihi_d,Qimhihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlohi_d,Qimlohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhilo_d,Qimhilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlolo_d,Qimlolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYTrehihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrehilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelolo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlolo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWTrehihi_d,szYWT + szpad); // padding for Y*W^T
   cudaMalloc((void**)&YWTrelohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrehilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrelolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTCrehihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrehilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelolo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlolo_d,sznum + szpad);

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
            GPU_cmplx4_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                 Arehihi_h, Arelohi_h, Arehilo_h, Arelolo_h,
                 Aimhihi_h, Aimlohi_h, Aimhilo_h, Aimlolo_h,
                 Arehihi_d, Arelohi_d, Arehilo_d, Arelolo_d,
                 Aimhihi_d, Aimlohi_d, Aimhilo_d, Aimlolo_d,
                 vrehihi_h, vrelohi_h, vrehilo_h, vrelolo_h,
                 vimhihi_h, vimlohi_h, vimhilo_h, vimlolo_h,
                 Vrehihi_d, Vrelohi_d, Vrehilo_d, Vrelolo_d,
                 Vimhihi_d, Vimlohi_d, Vimhilo_d, Vimlolo_d,
                betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            GPU_cmplx4_small_leftRupdate
               (nrows,ncols,szt,colidx,k,L,
                 Arehihi_h, Arelohi_h, Arehilo_h, Arelolo_h,
                 Aimhihi_h, Aimlohi_h, Aimhilo_h, Aimlolo_h,
                 Arehihi_d, Arelohi_d, Arehilo_d, Arelolo_d,
                 Aimhihi_d, Aimlohi_d, Aimhilo_d, Aimlolo_d,
                 Vrehihi_d, Vrelohi_d, Vrehilo_d, Vrelolo_d,
                 Vimhihi_d, Vimlohi_d, Vimhilo_d, Vimlolo_d,
                betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                tileRlapms,addcnt,mulcnt,verbose);
         }
         else
         {
            GPU_cmplx4_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Arehihi_h,   Arelohi_h,   Arehilo_h,   Arelolo_h,
                   Aimhihi_h,   Aimlohi_h,   Aimhilo_h,   Aimlolo_h,
                   Arehihi_d,   Arelohi_d,   Arehilo_d,   Arelolo_d,
                   Aimhihi_d,   Aimlohi_d,   Aimhilo_d,   Aimlolo_d,
                   vrehihi_h,   vrelohi_h,   vrehilo_h,   vrelolo_h,
                   vimhihi_h,   vimlohi_h,   vimhilo_h,   vimlolo_h,
                   Vrehihi_d,   Vrelohi_d,   Vrehilo_d,   Vrelolo_d,
                   Vimhihi_d,   Vimlohi_d,   Vimhilo_d,   Vimlolo_d,
                  betahihi_h,  betalohi_h,  betahilo_h,  betalolo_h,
                  betahihi_d,  betalohi_d,  betahilo_d,  betalolo_d,
                  sumshihi_h,  sumslohi_h,  sumshilo_h,  sumslolo_h,
                  sumshihi_d,  sumslohi_d,  sumshilo_d,  sumslolo_d,
                &sigmahihi_h,&sigmalohi_h,&sigmahilo_h,&sigmalolo_h,
                 sigmahihi_d, sigmalohi_d, sigmahilo_d, sigmalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            GPU_cmplx4_medium_leftRupdate
               (nrows,ncols,szt,colidx,k,L,
                     Arehihi_h,     Arelohi_h,     Arehilo_h,     Arelolo_h,
                     Aimhihi_h,     Aimlohi_h,     Aimhilo_h,     Aimlolo_h,
                     Arehihi_d,     Arelohi_d,     Arehilo_d,     Arelolo_d,
                     Aimhihi_d,     Aimlohi_d,     Aimhilo_d,     Aimlolo_d,
                     Vrehihi_d,     Vrelohi_d,     Vrehilo_d,     Vrelolo_d,
                     Vimhihi_d,     Vimlohi_d,     Vimhilo_d,     Vimlolo_d,
                    betahihi_h,    betalohi_h,    betahilo_h,    betalolo_h,
                    betahihi_d,    betalohi_d,    betahilo_d,    betalolo_d,
                RHdotvrehihi_h,RHdotvrelohi_h,RHdotvrehilo_h,RHdotvrelolo_h,
                RHdotvimhihi_h,RHdotvimlohi_h,RHdotvimhilo_h,RHdotvimlolo_h,
                RHdotvrehihi_d,RHdotvrelohi_d,RHdotvrehilo_d,RHdotvrelolo_d,
                RHdotvimhihi_d,RHdotvimlohi_d,RHdotvimhilo_d,RHdotvimlolo_d,
                  bRHvrehihi_h,  bRHvrelohi_h,  bRHvrehilo_h,  bRHvrelolo_h,
                  bRHvimhihi_h,  bRHvimlohi_h,  bRHvimhilo_h,  bRHvimlolo_h,
                  bRHvrehihi_d,  bRHvrelohi_d,  bRHvrehilo_d,  bRHvrelolo_d,
                  bRHvimhihi_d,  bRHvimlohi_d,  bRHvimhilo_d,  bRHvimlolo_d,
                RHvlapms,tileRlapms,addcnt,mulcnt,verbose);
         }
      }
      GPU_cmplx4_medium_VB_to_W
         (nrows,szt,szt,k,
            Vrehihi_h,  Vrelohi_h,  Vrehilo_h,  Vrelolo_h,
            Vimhihi_h,  Vimlohi_h,  Vimhilo_h,  Vimlolo_h,
            Vrehihi_d,  Vrelohi_d,  Vrehilo_d,  Vrelolo_d,
            Vimhihi_d,  Vimlohi_d,  Vimhilo_d,  Vimlolo_d,
            Wrehihi_h,  Wrelohi_h,  Wrehilo_h,  Wrelolo_h,
            Wimhihi_h,  Wimlohi_h,  Wimhilo_h,  Wimlolo_h,
            Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
            Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
          WYTrehihi_h,WYTrelohi_h,WYTrehilo_h,WYTrelolo_h,
          WYTimhihi_h,WYTimlohi_h,WYTimhilo_h,WYTimlolo_h,
          WYTrehihi_d,WYTrelohi_d,WYTrehilo_d,WYTrelolo_d,
          WYTimhihi_d,WYTimlohi_d,WYTimhilo_d,WYTimlolo_d,
           betahihi_h, betalohi_h, betahilo_h, betalolo_h,
           betahihi_d, betalohi_d, betahilo_d, betalolo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);

/*
      GPU_cmplx2_small_WYH(nrows-k*szt,szt,
            Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
            Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
          WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d,
          WYTrehi_h,WYTrelo_h,WYTimhi_h,WYTimlo_h,WYHlapms,verbose);
 */

      GPU_cmplx4_small_QWYH(nrows,szt,k,
             Qrehihi_d,   Qrelohi_d,   Qrehilo_d,   Qrelolo_d,
             Qimhihi_d,   Qimlohi_d,   Qimhilo_d,   Qimlolo_d,
           WYTrehihi_d, WYTrelohi_d, WYTrehilo_d, WYTrelolo_d,
           WYTimhihi_d, WYTimlohi_d, WYTimhilo_d, WYTimlolo_d,
          QWYTrehihi_d,QWYTrelohi_d,QWYTrehilo_d,QWYTrelolo_d,
          QWYTimhihi_d,QWYTimlohi_d,QWYTimhilo_d,QWYTimlolo_d,
          QWYTrehihi_h,QWYTrelohi_h,QWYTrehilo_h,QWYTrelolo_h,
          QWYTimhihi_h,QWYTimlohi_h,QWYTimhilo_h,QWYTimlolo_h,
             Qrehihi_h,   Qrelohi_h,   Qrehilo_h,   Qrelolo_h,
             Qimhihi_h,   Qimlohi_h,   Qimhilo_h,   Qimlolo_h,
          QWYHlapms,addcnt,mulcnt,verbose);

      GPU_cmplx4_small_Qupdate
         (nrows,szt,k,
             Qrehihi_d,   Qrelohi_d,   Qrehilo_d,   Qrelolo_d,
             Qimhihi_d,   Qimlohi_d,   Qimhilo_d,   Qimlolo_d,
          QWYTrehihi_d,QWYTrelohi_d,QWYTrehilo_d,QWYTrelolo_d,
          QWYTimhihi_d,QWYTimlohi_d,QWYTimhilo_d,QWYTimlolo_d,
             Qrehihi_h,   Qrelohi_h,   Qrehilo_h,   Qrelolo_h,
             Qimhihi_h,   Qimlohi_h,   Qimhilo_h,   Qimlolo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                              // update R
      {
         GPU_cmplx4_small_YWH
            (nrows,szt,k,
               Vrehihi_d,  Vrelohi_d,  Vrehilo_d,  Vrelolo_d,
               Vimhihi_d,  Vimlohi_d,  Vimhilo_d,  Vimlolo_d,
               Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
               Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
             YWTrehihi_d,YWTrelohi_d,YWTrehilo_d,YWTrelolo_d,
             YWTimhihi_d,YWTimlohi_d,YWTimhilo_d,YWTimlolo_d,
             YWTrehihi_h,YWTrelohi_h,YWTrehilo_h,YWTrelolo_h,
             YWTimhihi_h,YWTimlohi_h,YWTimhilo_h,YWTimlolo_h,
             YWHlapms,addcnt,mulcnt,verbose);

         GPU_cmplx4_small_YWHC
            (nrows,ncols,szt,k,
              YWTrehihi_d, YWTrelohi_d, YWTrehilo_d, YWTrelolo_d,
              YWTimhihi_d, YWTimlohi_d, YWTimhilo_d, YWTimlolo_d,
                Arehihi_d,   Arelohi_d,   Arehilo_d,   Arelolo_d,
                Aimhihi_d,   Aimlohi_d,   Aimhilo_d,   Aimlolo_d,
             YWTCrehihi_d,YWTCrelohi_d,YWTCrehilo_d,YWTCrelolo_d,
             YWTCimhihi_d,YWTCimlohi_d,YWTCimhilo_d,YWTCimlolo_d,
             YWTCrehihi_h,YWTCrelohi_h,YWTCrehilo_h,YWTCrelolo_h,
             YWTCimhihi_h,YWTCimlohi_h,YWTCimhilo_h,YWTCimlolo_h,
             YWHClapms,addcnt,mulcnt,verbose);

         GPU_cmplx4_small_R_add_YWHC
            (nrows,ncols,szt,k,
                Arehihi_d,   Arelohi_d,   Arehilo_d,   Arelolo_d,
                Aimhihi_d,   Aimlohi_d,   Aimhilo_d,   Aimlolo_d,
             YWTCrehihi_d,YWTCrelohi_d,YWTCrehilo_d,YWTCrelolo_d,
             YWTCimhihi_d,YWTCimlohi_d,YWTCimhilo_d,YWTCimlolo_d,
                Arehihi_h,   Arelohi_h,   Arehilo_h,   Arelolo_h,
                Aimhihi_h,   Aimlohi_h,   Aimhilo_h,   Aimlolo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qrehihi_h,Qrehihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelohi_h,Qrelohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrehilo_h,Qrehilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelolo_h,Qrelolo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhihi_h,Qimhihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlohi_h,Qimlohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhilo_h,Qimhilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlolo_h,Qimlolo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qrehihi[i][j] = Qrehihi_h[ix];
         Qrelohi[i][j] = Qrelohi_h[ix];
         Qrehilo[i][j] = Qrehilo_h[ix];
         Qrelolo[i][j] = Qrelolo_h[ix];
         Qimhihi[i][j] = Qimhihi_h[ix];
         Qimlohi[i][j] = Qimlohi_h[ix];
         Qimhilo[i][j] = Qimhilo_h[ix];
         Qimlolo[i][j] = Qimlolo_h[ix++];
      }

   cudaMemcpy(Arehihi_h,Arehihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelohi_h,Arelohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arehilo_h,Arehilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelolo_h,Arelolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhihi_h,Aimhihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlohi_h,Aimlohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhilo_h,Aimhilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlolo_h,Aimlolo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rrehihi[i][j] = Arehihi_h[j*nrows+i];
         Rrelohi[i][j] = Arelohi_h[j*nrows+i];
         Rrehilo[i][j] = Arehilo_h[j*nrows+i];
         Rrelolo[i][j] = Arelolo_h[j*nrows+i];
         Rimhihi[i][j] = Aimhihi_h[j*nrows+i];
         Rimlohi[i][j] = Aimlohi_h[j*nrows+i];
         Rimhilo[i][j] = Aimhilo_h[j*nrows+i];
         Rimlolo[i][j] = Aimlolo_h[j*nrows+i];
      }

   free(Arehihi_h); free(Arelohi_h); free(Arehilo_h); free(Arelolo_h);
   free(Aimhihi_h); free(Aimlohi_h); free(Aimhilo_h); free(Aimlolo_h);
   free(Qrehihi_h); free(Qrelohi_h); free(Qrehilo_h); free(Qrelolo_h);
   free(Qimhihi_h); free(Qimlohi_h); free(Qimhilo_h); free(Qimlolo_h);
   free(vrehihi_h); free(vrelohi_h); free(vrehilo_h); free(vrelolo_h); 
   free(vimhihi_h); free(vimlohi_h); free(vimlohi_h); free(vimlolo_h);
   free(Vrehihi_h); free(Vrelohi_h); free(Vrehilo_h); free(Vrelolo_h);
   free(Vimhihi_h); free(Vimlohi_h); free(Vimhilo_h); free(Vimlolo_h);
   free(RHdotvrehihi_h); free(RHdotvrelohi_h);
   free(RHdotvrehilo_h); free(RHdotvrelolo_h);
   free(RHdotvimhihi_h); free(RHdotvimlohi_h);
   free(RHdotvimhilo_h); free(RHdotvimlolo_h);
   free(bRHvrehihi_h); free(bRHvrelohi_h);
   free(bRHvrehilo_h); free(bRHvrelolo_h);
   free(bRHvimhihi_h); free(bRHvimlohi_h);
   free(bRHvimhilo_h); free(bRHvimlolo_h);
   free(sumshihi_h); free(sumslohi_h);
   free(sumshilo_h); free(sumslolo_h);
   free(Wrehihi_h); free(Wrelohi_h); free(Wrehilo_h); free(Wrelolo_h);
   free(Wimhihi_h); free(Wimlohi_h); free(Wimhilo_h); free(Wimlolo_h);
   free(WYTrehihi_h); free(WYTrelohi_h); free(WYTrehilo_h); free(WYTrelolo_h);
   free(WYTimhihi_h); free(WYTimlohi_h); free(WYTimhilo_h); free(WYTimlolo_h);
   free(QWYTrehihi_h); free(QWYTrelohi_h);
   free(QWYTrehilo_h); free(QWYTrelolo_h);
   free(QWYTimhihi_h); free(QWYTimlohi_h);
   free(QWYTimhilo_h); free(QWYTimlolo_h);
   free(YWTrehihi_h); free(YWTrelohi_h); free(YWTrehilo_h); free(YWTrelolo_h);
   free(YWTimhihi_h); free(YWTimlohi_h); free(YWTimhilo_h); free(YWTimlolo_h);
   free(WYTrehihi_h); free(WYTrelohi_h); free(WYTrehilo_h); free(WYTrelolo_h);
   free(WYTimhihi_h); free(WYTimlohi_h); free(WYTimhilo_h); free(WYTimlolo_h);
   free(YWTCrehihi_h); free(YWTCrelohi_h);
   free(YWTCrehilo_h); free(YWTCrelolo_h);
   free(YWTCimhihi_h); free(YWTCimlohi_h);
   free(YWTCimhilo_h); free(YWTCimlolo_h);
}
