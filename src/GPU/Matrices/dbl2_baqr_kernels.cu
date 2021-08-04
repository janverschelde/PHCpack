/* The file dbl2_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl2_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_baqr_kernels.h"

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
   ddg_div(shvhi[j],shvlo[j],prdhi[0],prdlo[0],&acchi,&acclo);
   vhi[j+1] = acchi;
   vlo[j+1] = acclo;
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

__global__ void dbl2_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double *Vhi, double *Vlo, double *Whi, double *Wlo )
{
   const int tdx = threadIdx.x;        // index of thread in block
   double wrkhi,wrklo,pkhi,pklo,mypkhi,mypklo,zihi,zilo;
   int idx;

   __shared__ double shvhi[dd_shmemsize]; // one work vector
   __shared__ double shvlo[dd_shmemsize];
   __shared__ double shwhi[dd_shmemsize]; // the other work vector
   __shared__ double shwlo[dd_shmemsize];
   __shared__ double shphi[dd_shmemsize]; // to share Y^T*v
   __shared__ double shplo[dd_shmemsize]; 

   shvhi[tdx] = Vhi[tdx];
   shvlo[tdx] = Vlo[tdx];
   // wrk = -B[0]*shv[tdx];               // first column of W
   ddg_mul(Bhi[0],Blo[0],shvhi[tdx],shvlo[tdx],&wrkhi,&wrklo);
   ddg_minus(&wrkhi,&wrklo);
   Whi[tdx] = wrkhi;
   Wlo[tdx] = wrklo;

   for(int j=1; j<ncols; j++)          // compute column j of W
   {
      idx = j*nrows + tdx;
      shvhi[tdx] = Vhi[idx];           // j-th Householder vector
      shvlo[tdx] = Vlo[idx]; 

      for(int k=0; k<j; k++)
      {
         pkhi = 0.0;                   // k-th component of Y^T*v
         pklo = 0.0;
         idx = k*nrows + tdx;
         shwhi[tdx] = Vhi[idx];        // load V[k][i]
         shwlo[tdx] = Vlo[idx]; 
         // shp[tdx] = shw[tdx]*shv[tdx]; // V[k][i]*v[i]
         ddg_mul(shwhi[tdx],shwlo[tdx],shvhi[tdx],shvlo[tdx],
                 &shphi[tdx],&shplo[tdx]);

         __syncthreads();
         for(int i=0; i<nrows; i++) // pk = pk + shp[i];
            ddg_inc(&pkhi,&pklo,shphi[i],shplo[i]);

         __syncthreads();        // this synchronization is critical!
         if(tdx == k)
         {
            mypkhi = pkhi;
            mypklo = pklo;
         }
      }
      __syncthreads();
      shphi[tdx] = mypkhi;             // share p[k]
      shplo[tdx] = mypklo;
      __syncthreads();
      zihi = 0.0;                      // i-th component of W*p
      zilo = 0.0;
      for(int k=0; k<j; k++)
      {
         idx = k*nrows + tdx;
         shwhi[tdx] = Whi[idx];        // load W[k][i]
         shwlo[tdx] = Wlo[idx];
         // zi = zi + shw[tdx]*shp[k];
         ddg_mul(shwhi[tdx],shwlo[tdx],shphi[k],shplo[k],&wrkhi,&wrklo);
         ddg_inc(&zihi,&zilo,wrkhi,wrklo);
      }
      // zi = zi + shv[tdx];
      ddg_inc(&zihi,&zilo,shvhi[tdx],shvlo[tdx]);
      // wrk = -B[j]*zi;
      ddg_mul(Bhi[j],Blo[j],zihi,zilo,&wrkhi,&wrklo);
      ddg_minus(&wrkhi,&wrklo);
      idx = j*nrows + tdx;
      Whi[idx] = wrkhi;               // wrk is assigned to W[j][tdx]
      Wlo[idx] = wrklo;
      __syncthreads();
   }
}

__global__ void cmplx2_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double *Vrehi, double *Vrelo, double *Vimhi, double *Vimlo,
   double *Wrehi, double *Wrelo, double *Wimhi, double *Wimlo )
{
   const int tdx = threadIdx.x;        // index of thread in block
   double wrk_rehi,wrk_relo,wrk_imhi,wrk_imlo;
   double pk_rehi,pk_relo,pk_imhi,pk_imlo;
   double mypk_rehi,mypk_relo,mypk_imhi,mypk_imlo;
   double zi_rehi,zi_relo,zi_imhi,zi_imlo;
   int VWidx;

   __shared__ double shvrehi[cdd_shmemsize]; // one work vector
   __shared__ double shvrelo[cdd_shmemsize]; 
   __shared__ double shvimhi[cdd_shmemsize];
   __shared__ double shvimlo[cdd_shmemsize];
   __shared__ double shwrehi[cdd_shmemsize]; // the other work vector
   __shared__ double shwrelo[cdd_shmemsize];
   __shared__ double shwimhi[cdd_shmemsize];
   __shared__ double shwimlo[cdd_shmemsize];
   __shared__ double shprehi[cdd_shmemsize]; // to share Y^T*v
   __shared__ double shprelo[cdd_shmemsize];
   __shared__ double shpimhi[cdd_shmemsize];
   __shared__ double shpimlo[cdd_shmemsize];

   shvrehi[tdx] = Vrehi[tdx];
   shvrelo[tdx] = Vrelo[tdx];
   shvimhi[tdx] = Vimhi[tdx];
   shvimlo[tdx] = Vimlo[tdx];
   // wrk_re = -B[0]*shvre[tdx];            // first column of W
   ddg_mul(Bhi[0],Blo[0],shvrehi[tdx],shvrelo[tdx],&wrk_rehi,&wrk_relo);
   ddg_minus(&wrk_rehi,&wrk_relo);
   // wrk_im = -B[0]*shvim[tdx];
   ddg_mul(Bhi[0],Blo[0],shvimhi[tdx],shvimlo[tdx],&wrk_imhi,&wrk_imlo);
   ddg_minus(&wrk_imhi,&wrk_imlo);
   Wrehi[tdx] = wrk_rehi;
   Wrelo[tdx] = wrk_relo;
   Wimhi[tdx] = wrk_imhi;
   Wimlo[tdx] = wrk_imlo;

   for(int j=1; j<ncols; j++)          // compute column j of W
   {
      VWidx = j*nrows + tdx;
      shvrehi[tdx] = Vrehi[VWidx];     // j-th Householder vector
      shvrelo[tdx] = Vrelo[VWidx];
      shvimhi[tdx] = Vimhi[VWidx];
      shvimlo[tdx] = Vimlo[VWidx];

      for(int k=0; k<j; k++)
      {
         pk_rehi = 0.0;                // k-th component of Y^H*v
         pk_relo = 0.0;
         pk_imhi = 0.0;
         pk_imlo = 0.0;
         VWidx = k*nrows + tdx;
         shwrehi[tdx] = Vrehi[VWidx];  // load V[k][i]
         shwrelo[tdx] = Vrelo[VWidx];
         shwimhi[tdx] = Vimhi[VWidx];
         shwimlo[tdx] = Vimlo[VWidx];
         // shp[tdx] = shw[tdx]*shv[tdx]; V[k][i]*v[i], Hermitian transpose!
         // shpre[tdx] =   shwre[tdx]*shvre[tdx] + shwim[tdx]*shvim[tdx];
         ddg_mul(shwrehi[tdx],shwrelo[tdx],
                 shvrehi[tdx],shvrelo[tdx],&shprehi[tdx],&shprelo[tdx]);
         ddg_mul(shwimhi[tdx],shwimlo[tdx],
                 shvimhi[tdx],shvimlo[tdx],&wrk_imhi,&wrk_imlo);
         ddg_inc(&shprehi[tdx],&shprelo[tdx],wrk_imhi,wrk_imlo);
         // shpim[tdx] = - shwim[tdx]*shvre[tdx] + shwre[tdx]*shvim[tdx];
         ddg_mul(shwrehi[tdx],shwrelo[tdx],
                 shvimhi[tdx],shvimlo[tdx],&shpimhi[tdx],&shpimlo[tdx]);
         ddg_mul(shwimhi[tdx],shwimlo[tdx],
                 shvrehi[tdx],shvrelo[tdx],&wrk_imhi,&wrk_imlo);
         ddg_dec(&shpimhi[tdx],&shpimlo[tdx],wrk_imhi,wrk_imlo);
         __syncthreads();
         for(int i=0; i<nrows; i++)
         {
            // pk_re = pk_re + shpre[i];
            ddg_inc(&pk_rehi,&pk_relo,shprehi[i],shprelo[i]);
            // pk_im = pk_im + shpim[i];
            ddg_inc(&pk_imhi,&pk_imlo,shpimhi[i],shpimlo[i]);
         }
         __syncthreads();              // important synchronization
         if(tdx == k)
         {
            mypk_rehi = pk_rehi;
            mypk_relo = pk_relo;
            mypk_imhi = pk_imhi;
            mypk_imlo = pk_imlo;
         }
      }
      __syncthreads();
      shprehi[tdx] = mypk_rehi;        // share p[k]
      shprelo[tdx] = mypk_relo;
      shpimhi[tdx] = mypk_imhi;
      shpimlo[tdx] = mypk_imlo;
      __syncthreads();
      zi_rehi = 0.0;                   // i-th component of W*p
      zi_relo = 0.0;
      zi_imhi = 0.0;
      zi_imlo = 0.0;
      for(int k=0; k<j; k++)
      {
         VWidx = k*nrows + tdx;
         shwrehi[tdx] = Wrehi[VWidx];  // load W[k][i]
         shwrelo[tdx] = Wrelo[VWidx];
         shwimhi[tdx] = Wimhi[VWidx];
         shwimlo[tdx] = Wimlo[VWidx];
         // zi = zi + shw[tdx]*shp[k];
         // zi_re = zi_re + shwre[tdx]*shpre[k] - shwim[tdx]*shpim[k];
         ddg_mul(shwrehi[tdx],shwrelo[tdx],shprehi[k],shprelo[k],
                 &wrk_rehi,&wrk_relo);
         ddg_inc(&zi_rehi,&zi_relo,wrk_rehi,wrk_relo);
         ddg_mul(shwimhi[tdx],shwimlo[tdx],shpimhi[k],shpimlo[k],
                 &wrk_rehi,&wrk_relo);
         ddg_dec(&zi_rehi,&zi_relo,wrk_rehi,wrk_relo);
         // zi_im = zi_im + shwim[tdx]*shpre[k] + shwre[tdx]*shpim[k];
         ddg_mul(shwimhi[tdx],shwimlo[tdx],shprehi[k],shprelo[k],
                 &wrk_rehi,&wrk_relo);
         ddg_inc(&zi_imhi,&zi_imlo,wrk_rehi,wrk_relo);
         ddg_mul(shwrehi[tdx],shwrelo[tdx],shpimhi[k],shpimlo[k],
                 &wrk_rehi,&wrk_relo);
         ddg_inc(&zi_imhi,&zi_imlo,wrk_rehi,wrk_relo);
      }
      // zi_re = zi_re + shvre[tdx];
      ddg_inc(&zi_rehi,&zi_relo,shvrehi[tdx],shvrelo[tdx]);
      // zi_im = zi_im + shvim[tdx];
      ddg_inc(&zi_imhi,&zi_imlo,shvimhi[tdx],shvimlo[tdx]);
      // wrk_re = -B[j]*zi_re;
      ddg_mul(Bhi[j],Blo[j],zi_rehi,zi_relo,&wrk_rehi,&wrk_relo);
      ddg_minus(&wrk_rehi,&wrk_relo);
      // wrk_im = -B[j]*zi_im;
      ddg_mul(Bhi[j],Blo[j],zi_imhi,zi_imlo,&wrk_imhi,&wrk_imlo);
      ddg_minus(&wrk_imhi,&wrk_imlo);
      VWidx = j*nrows + tdx;
      Wrehi[VWidx] = wrk_rehi;        // wrk is assigned to W[j][tdx]
      Wrelo[VWidx] = wrk_relo;
      Wimhi[VWidx] = wrk_imhi;
      Wimlo[VWidx] = wrk_imlo;
      __syncthreads();
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

__global__ void cmplx2_small_WYT
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

__global__ void cmplx2_small_QWYT
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

__global__ void cmplx2_small_YWTC
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

__global__ void cmplx2_small_R_add_YWTC
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
   double *lapms, bool verbose )
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

void GPU_cmplx2_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehi_h, double *Arelo_h, double *Aimhi_h, double *Aimlo_h,
   double *Arehi_d, double *Arelo_d, double *Aimhi_d, double *Aimlo_d,
   double *vrehi_h, double *vrelo_h, double *vimhi_h, double *vimlo_h,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, bool verbose )
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

void GPU_dbl2_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, bool verbose )
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
   double *lapms, bool verbose )
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
   double *whi_h, double *wlo_h, double *whi_d, double *wlo_d,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix
   // total number of entries in R that will be modified
   const int sizenum = (nrows - colidx)*(endcol - colidx);
   const int nbrblocks = (int) ceil(sizenum/((double) szt));

   cudaEventRecord(start);           // 2nd argument: ncols -> endcol
   // changed second argument ncols into endcol
   // to avoid updating the next tile
   // dbl_medium_betaRTv<<<nbrblocks,szt>>>
   //   (nrows,endcol,szt,colidx,A_d,&V_d[L*nVrows+L],&beta_d[L],w_d);
   // number of threads must be ncols - colidx, not endcol - colidx
   dbl2_small_betaRTv<<<1,nrows-colidx>>> // nrows ...
     (nrows,endcol,szt,colidx,Ahi_d,Alo_d,
      &Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
      &betahi_d[L],&betalo_d[L],whi_d,wlo_d);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to update " << sizenum << " numbers ..." << endl;
      cout << "   nrows : " << nrows << "  endcol : " << endcol
           << "  szt : " << szt << "  colidx : " << colidx << endl;
   }
   dbl2_medium_subvbetaRTv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Ahi_d,Alo_d,&Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
       &betahi_d[L],&betalo_d[L],whi_d,wlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRTv = (endcol-colidx)*sizeof(double);

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

void GPU_dbl2_VB_to_W
 ( int nrows, int ncols, int szt,
   double *Vhi_h, double *Vlo_h, double *Vhi_d, double *Vlo_d,
   double *Whi_h, double *Wlo_h, double *Whi_d, double *Wlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   dbl2_VB_to_W<<<1,nrows>>>
      (nrows,ncols,betahi_d,betalo_d,Vhi_d,Vlo_d,Whi_d,Wlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = nrows*sizeof(double);
      const size_t szVandW = szt*szhouse;

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
         for(int i=0; i<nrows; i++) 
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
         for(int i=0; i<nrows; i++) 
         {
            cout << "W[" << i << "][" << j << "] : "
                 << Whi_h[ix] << "  " << Wlo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx2_VB_to_W
 ( int nrows, int ncols, int szt,
   double *Vrehi_h, double *Vrelo_h, double *Vimhi_h, double *Vimlo_h,
   double *Vrehi_d, double *Vrelo_d, double *Vimhi_d, double *Vimlo_d,
   double *Wrehi_h, double *Wrelo_h, double *Wimhi_h, double *Wimlo_h,
   double *Wrehi_d, double *Wrelo_d, double *Wimhi_d, double *Wimlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(verbose)
   {
      const size_t szhouse = nrows*sizeof(double);
      const size_t szVandW = szt*szhouse;

      cudaMemcpy(Vrehi_h,Vrehi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelo_h,Vrelo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhi_h,Vimhi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlo_h,Vimlo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<nrows; i++) 
         {
            cout << "V[" << i << "][" << j << "]re : "
                 << Vrehi_h[ix] << "  " << Vrelo_h[ix] << endl;
            cout << "V[" << i << "][" << j << "]im : "
                 << Vimhi_h[ix] << "  " << Vimlo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   cmplx2_VB_to_W<<<1,nrows>>>
      (nrows,ncols,betahi_d,betalo_d,Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
       Wrehi_d,Wrelo_d,Wimhi_d,Wimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = nrows*sizeof(double);
      const size_t szVandW = szt*szhouse;

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
         for(int i=0; i<nrows; i++) 
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
         for(int i=0; i<nrows; i++) 
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehi_h[ix] << "  " << Wrelo_h[ix] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhi_h[ix] << "  " << Wimlo_h[ix] << endl;
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

void GPU_cmplx2_small_WYT
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
   cmplx2_small_WYT<<<nbrblocks,szt>>>
      (nrows,szt,Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
                 Yrehi_d,  Yrelo_d,  Yimhi_d,  Yimlo_d,
               WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

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
   double *lapms, bool verbose )
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

void GPU_cmplx2_small_YWT
 ( int nrows, int szt, int idx,
   double *Yrehi_d, double *Yrelo_d, double *Yimhi_d, double *Yimlo_d,
   double *Wrehi_d, double *Wrelo_d, double *Wimhi_d, double *Wimlo_d,
   double *YWTrehi_d, double *YWTrelo_d, double *YWTimhi_d, double *YWTimlo_d,
   double *YWTrehi_h, double *YWTrelo_h, double *YWTimhi_h, double *YWTimlo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx2_small_WYT<<<nbrblocks,szt>>>
      (rowdim,szt,Yrehi_d,  Yrelo_d,  Yimhi_d,  Yimlo_d,
                  Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
                YWTrehi_d,YWTrelo_d,YWTimhi_d,YWTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

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
   double *lapms, bool verbose )
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

void GPU_cmplx2_small_QWYT
 ( int dim, int szt, int idx,
   double *Qrehi_d, double *Qrelo_d, double *Qimhi_d, double *Qimlo_d,
   double *WYTrehi_d, double *WYTrelo_d, double *WYTimhi_d, double *WYTimlo_d,
   double *QWYTrehi_d, double *QWYTrelo_d,
   double *QWYTimhi_d, double *QWYTimlo_d,
   double *QWYTrehi_h, double *QWYTrelo_h,
   double *QWYTimhi_h, double *QWYTimlo_h,
   double *Qrehi_h, double *Qrelo_h, double *Qimhi_h, double *Qimlo_h,
   double *lapms, bool verbose )
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
   cmplx2_small_QWYT<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
        WYTrehi_d, WYTrelo_d, WYTimhi_d, WYTimlo_d,
       QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

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
   double *YWTChi_h, double *YWTClo_h, double *lapms, bool verbose )
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

void GPU_cmplx2_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehi_d, double *YWTrelo_d, double *YWTimhi_d, double *YWTimlo_d,
   double *Crehi_d, double *Crelo_d, double *Cimhi_d, double *Cimlo_d,
   double *YWTCrehi_d, double *YWTCrelo_d,
   double *YWTCimhi_d, double *YWTCimlo_d,
   double *YWTCrehi_h, double *YWTCrelo_h,
   double *YWTCimhi_h, double *YWTCimlo_h, double *lapms, bool verbose )
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
   cmplx2_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWTrehi_d, YWTrelo_d, YWTimhi_d, YWTimlo_d,
          Crehi_d,   Crelo_d,   Cimhi_d,   Cimlo_d,
       YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

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
   double *Qhi_h, double *Qlo_h, double *lapms, bool verbose )
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
   double *lapms, bool verbose )
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
   double *lapms, bool verbose )
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

void GPU_cmplx2_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehi_d, double *Rrelo_d, double *Rimhi_d, double *Rimlo_d,
   double *YWTCrehi_d, double *YWTCrelo_d,
   double *YWTCimhi_d, double *YWTCimlo_d,
   double *Rrehi_h, double *Rrelo_h, double *Rimhi_h, double *Rimlo_h,
   double *lapms, bool verbose )
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
   cmplx2_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
          Rrehi_d,   Rrelo_d,   Rimhi_d,   Rimlo_d,
       YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

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
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, bool verbose )
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
   double *bRTvhi_h = new double[nrows];  // high doubles of beta*R^T*v
   double *bRTvlo_h = new double[nrows];  // low doubles of beta*R^T*v
   double *bRTvhi_d;                      // bRTvhi_h on the device
   double *bRTvlo_d;                      // bRTvlo_h on the device

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

   cudaMalloc((void**)&bRTvhi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlo_d,szhouse + szpad);

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

   *houselapms = 0.0;
   *tileRlapms = 0.0;
   *vb2Wlapms = 0.0;
   *WYTlapms = 0.0; *QWYTlapms = 0.0; *Qaddlapms = 0.0;
   *YWTlapms = 0.0; *YWTClapms = 0.0; *Raddlapms = 0.0;
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
         GPU_dbl2_small_house
            (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
             Ahi_h,Alo_h,Ahi_d,Alo_d,vhi_h,vlo_h,Vhi_d,Vlo_d,
             betahi_h,betalo_h,betahi_d,betalo_d,houselapms,verbose);
         if(nrows - colidx <= szt)
         {
            GPU_dbl2_small_leftRupdate
               (nrows,ncols,szt,colidx,k,L,Ahi_h,Alo_h,Ahi_d,Alo_d,
                Vhi_d,Vlo_d,betahi_h,betalo_h,betahi_d,betalo_d,
                tileRlapms,verbose);
         }
         else
         {
            GPU_dbl2_medium_leftRupdate
               (nrows,ncols,szt,colidx,k,L,Ahi_h,Alo_h,Ahi_d,Alo_d,
                Vhi_d,Vlo_d,betahi_h,betalo_h,betahi_d,betalo_d,
                bRTvhi_h,bRTvlo_h,bRTvhi_d,bRTvlo_d,tileRlapms,verbose);
         }
      }
      // changed nrows into nrows - k*szt and ncols into szt
      GPU_dbl2_VB_to_W
         (nrows-k*szt,szt,szt,Vhi_h,Vlo_h,Vhi_d,Vlo_d,Whi_h,Wlo_h,
          Whi_d,Wlo_d,betahi_h,betalo_h,betahi_d,betalo_d,vb2Wlapms,verbose);
      // update Q, WYT matrix has nrows - k*szt instead of nrows
      GPU_dbl2_small_WYT
         (nrows-k*szt,szt,Whi_d,Wlo_d,Vhi_d,Vlo_d,WYThi_d,WYTlo_d,
          WYThi_h,WYTlo_h,WYTlapms,verbose);
      GPU_dbl2_small_QWYT
         (nrows,szt,k,Qhi_d,Qlo_d,WYThi_d,WYTlo_d,QWYThi_d,QWYTlo_d,
          QWYThi_h,QWYTlo_h,Qhi_h,Qlo_h,QWYTlapms,verbose);
      GPU_dbl2_small_Qupdate
         (nrows,szt,k,Qhi_d,Qlo_d,QWYThi_d,QWYTlo_d,Qhi_h,Qlo_h,
          Qaddlapms,verbose);
      if(k < nbt-1)                                           // update R
      {
         GPU_dbl2_small_YWT
            (nrows,szt,k,Vhi_d,Vlo_d,Whi_d,Wlo_d,YWThi_d,YWTlo_d,
             YWThi_h,YWTlo_h,YWTlapms,verbose);
         GPU_dbl2_small_YWTC
            (nrows,ncols,szt,k,YWThi_d,YWTlo_d,Ahi_d,Alo_d,
             YWTChi_d,YWTClo_d,YWTChi_h,YWTClo_h,YWTClapms,verbose);
         GPU_dbl2_small_R_add_YWTC
            (nrows,ncols,szt,k,Ahi_d,Alo_d,YWTChi_d,YWTClo_d,
             Ahi_h,Alo_h,Raddlapms,verbose);
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
   free(Whi_h); free(Wlo_h); free(bRTvhi_h); free(bRTvlo_h);
   free(WYThi_h); free(QWYThi_h); free(YWThi_h); free(YWTChi_h);
   free(WYTlo_h); free(QWYTlo_h); free(YWTlo_h); free(YWTClo_h);
}

void GPU_cmplx2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, bool verbose )
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
   double *betahi_h = new double[szt];   // high doubles of beta
   double *betalo_h = new double[szt];   // low doubles of beta
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

   *houselapms = 0.0;
   *tileRlapms = 0.0;
   *vb2Wlapms = 0.0;
   *WYTlapms = 0.0; *QWYTlapms = 0.0; *Qaddlapms = 0.0;
   *YWTlapms = 0.0; *YWTClapms = 0.0; *Raddlapms = 0.0;
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
         GPU_cmplx2_small_house
            (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
             Arehi_h,Arelo_h,Aimhi_h,Aimlo_h,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
             vrehi_h,vrelo_h,vimhi_h,vimlo_h,Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
             betahi_h,betalo_h,betahi_d,betalo_d,houselapms,verbose);
         GPU_cmplx2_small_leftRupdate
            (nrows,ncols,szt,colidx,k,L,
             Arehi_h,Arelo_h,Aimhi_h,Aimlo_h,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
             Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
             betahi_h,betalo_h,betahi_d,betalo_d,tileRlapms,verbose);
      }
      GPU_cmplx2_VB_to_W(nrows-k*szt,szt,szt,
          Vrehi_h,Vrelo_h,Vimhi_h,Vimlo_h,Vrehi_d,Vrelo_d,Vimhi_d,Vimlo_d,
          Wrehi_h,Wrelo_h,Wimhi_h,Wimlo_h,Wrehi_d,Wrelo_d,Wimhi_d,Wimlo_d,
          betahi_h,betalo_h,betahi_d,betalo_d,vb2Wlapms,verbose);
      GPU_cmplx2_small_WYT(nrows-k*szt,szt,
            Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
            Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
          WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d,
          WYTrehi_h,WYTrelo_h,WYTimhi_h,WYTimlo_h,WYTlapms,verbose);
      GPU_cmplx2_small_QWYT(nrows,szt,k,
             Qrehi_d,   Qrelo_d,   Qimhi_d,   Qimlo_d,
           WYTrehi_d, WYTrelo_d, WYTimhi_d, WYTimlo_d,
          QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d,
          QWYTrehi_h,QWYTrelo_h,QWYTimhi_h,QWYTimlo_h,
             Qrehi_h,   Qrelo_h,   Qimhi_h,   Qimlo_h, QWYTlapms,verbose);
      GPU_cmplx2_small_Qupdate
         (nrows,szt,k,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
          QWYTrehi_d,QWYTrelo_d,QWYTimhi_d,QWYTimlo_d,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Qaddlapms,verbose);
      if(k < nbt-1)                              // update R
      {
         GPU_cmplx2_small_YWT
            (nrows,szt,k,
               Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
               Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
             YWTrehi_d,YWTrelo_d,YWTimhi_d,YWTimlo_d,
             YWTrehi_h,YWTrelo_h,YWTimhi_h,YWTimlo_h,YWTlapms,verbose);
         GPU_cmplx2_small_YWTC
            (nrows,ncols,szt,k,
              YWTrehi_d, YWTrelo_d, YWTimhi_d, YWTimlo_d,
                Arehi_d,   Arelo_d,   Aimhi_d,   Aimlo_d,
             YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d,
             YWTCrehi_h,YWTCrelo_h,YWTCimhi_h,YWTCimlo_h,YWTClapms,verbose);
         GPU_cmplx2_small_R_add_YWTC
            (nrows,ncols,szt,k,
                Arehi_d,   Arelo_d,   Aimhi_d,   Aimlo_d,
             YWTCrehi_d,YWTCrelo_d,YWTCimhi_d,YWTCimlo_d,
                Arehi_h,   Arelo_h,   Aimhi_h,   Aimlo_h,Raddlapms,verbose);
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
   free(Wrehi_h); free(Wimhi_h); free(Wrelo_h); free(Wimlo_h);
   free(WYTrehi_h); free(QWYTrehi_h); free(YWTrehi_h); free(YWTCrehi_h);
   free(WYTrelo_h); free(QWYTrelo_h); free(YWTrelo_h); free(YWTCrelo_h);
   free(WYTimhi_h); free(QWYTimhi_h); free(YWTimhi_h); free(YWTCimhi_h);
   free(WYTimlo_h); free(QWYTimlo_h); free(YWTimlo_h); free(YWTCimlo_h);
}
