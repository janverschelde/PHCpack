/* The file dbl_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "dbl_baqr_kernels.h"

using namespace std;

__global__ void dbl_small_house
 ( double *x0, double *x1, int dim, int dimLog2, double *v, double *beta )
{
   const int j = threadIdx.x;

   __shared__ double shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];

   bool stopflag = false;
   double mu,v0,v0p2;

   shv[j] = x1[j];              // reading of vector into shared memory
   prd[j] = shv[j]*shv[j];      // for the 2-norm computation

   v[j+1] = shv[j];             // copies x to v, in case beta is zero
   if(j == 0) v[0] = 1.0;

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {
      if(prd[0] == 0.0)                    // prd[0] is sigma of house
      {
         *beta = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   if(j == 0)                              // thread zero sets beta
   {
      mu = sqrt((*x0)*(*x0) + prd[0]);
      if(*x0 <= 0.0)
         v0 = *x0 - mu;
      else
         v0 = -prd[0]/(*x0 + mu);

      v0p2 = v0*v0;
      *beta = 2.0*v0p2/(prd[0] + v0p2);
      prd[0] = v0;                         // v0 needed for normalization
   }
   __syncthreads();
   shv[j] = shv[j]/prd[0];
   v[j+1] = shv[j];
   if(j == 0) v[0] = 1.0;
}

__global__ void cmplx_small_house
 ( double *x0re, double *x0im, double *x1re, double *x1im,
   int dim, int dimLog2, double *vre, double *vim, double *beta )
{
   const int j = threadIdx.x;

   __shared__ double shvre[d_shmemsize];
   __shared__ double shvim[d_shmemsize];
   __shared__ double prd[d_shmemsize];
   __shared__ double v0parts[2];

   bool stopflag = false;
   double mu,v0re,v0im,x0rad,sqrx0,sqrv0,inv0re,inv0im,zre,zim;

   shvre[j] = x1re[j];          // reading of vector into shared memory
   shvim[j] = x1im[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];

   vre[j+1] = shvre[j];         // copies x to v, in case beta is zero
   vim[j+1] = shvim[j];
   if(j == 0) vre[0] = 1.0;
   if(j == 0) vim[0] = 0.0;

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {
      if(prd[0] == 0.0)                    // prd[0] is sigma of house
      {
         *beta = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   if(j == 0)                              // thread zero sets beta
   {
      sqrx0 = (*x0re)*(*x0re) + (*x0im)*(*x0im);
      x0rad = sqrt(sqrx0);
      mu = sqrt(sqrx0 + prd[0]);

      if(x0rad == 0.0)
      {
         v0re = -mu;
         v0im = 0.0;
      }
      else
      {
         mu = mu/x0rad;
         v0re = (*x0re) - mu*(*x0re);
         v0im = (*x0im) - mu*(*x0im);
      }
      sqrv0 = v0re*v0re + v0im*v0im;
      *beta = 2.0*sqrv0/(prd[0] + sqrv0);

      prd[0] = sqrv0;                     // sqrv0 needed for normalization
      v0parts[0] = v0re;                  // share v0re with all threads
      v0parts[1] = v0im;                  // share v0im with all threads
   }
   __syncthreads();
   inv0re = v0parts[0]/prd[0];               // real part of 1/v[0]
   inv0im = -v0parts[1]/prd[0];              // imag part of 1/v[0]
   zre = shvre[j]*inv0re - shvim[j]*inv0im;  // real part of v[j]/v[0]
   zim = shvim[j]*inv0re + shvre[j]*inv0im;  // imag part of v[j]/v[0]
   vre[j+1] = zre;
   vim[j+1] = zim;
   if(j == 0) vre[0] = 1.0;
   if(j == 0) vim[0] = 0.0;
}

__global__ void dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *R, double *v, double *beta )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   int Rcolidx;
   double w,Rtdx;

   __shared__ double shv[d_shmemsize]; // slice of v

   shv[tdx] = v[tdx];
   __syncthreads();
   w = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rtdx = R[Roffset + i + tdx*nrows];
      w = w + Rtdx*shv[i];
   }
   w = (*beta)*w;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx = R[Rcolidx];
      Rtdx = Rtdx - shv[i]*w;
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k) R[Rcolidx] = Rtdx;
   }
}

__global__ void cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rre, double *Rim, double *vre, double *vim, double *beta )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   int Rcolidx;
   double w_re,w_im,Rtdx_re,Rtdx_im,acc;

   __shared__ double shvre[d_shmemsize]; // slice of vre
   __shared__ double shvim[d_shmemsize]; // slice of vim

   shvre[tdx] = vre[tdx];
   shvim[tdx] = vim[tdx];
   __syncthreads();
   w_re = 0.0;
   w_im = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx_re = Rre[Rcolidx];
      Rtdx_im = Rim[Rcolidx];
      // w = w + Rtdx*shv[i]; beware of the Hermitian transpose!
      w_re = w_re + Rtdx_re*shvre[i] + Rtdx_im*shvim[i];
      w_im = w_im - Rtdx_im*shvre[i] + Rtdx_re*shvim[i];
   }
   acc = *beta;
   w_re = acc*w_re;
   w_im = acc*w_im;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx_re = Rre[Rcolidx];
      Rtdx_im = Rim[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w; beware of the Hermitian transpose!
      Rtdx_re = Rtdx_re - (shvre[i]*w_re + shvim[i]*w_im);
      Rtdx_im = Rtdx_im - (shvim[i]*w_re - shvre[i]*w_im);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k) Rre[Rcolidx] = Rtdx_re;
      __syncthreads();
      if(tdx < ncols-k) Rim[Rcolidx] = Rtdx_im;
   }
}

__global__ void dbl_VB_to_W
 ( int nrows, int ncols, double *B, double *V, double *W )
{
   const int tdx = threadIdx.x;        // index of thread in block
   double wrk,pk,mypk,zi;

   __shared__ double shv[d_shmemsize]; // one work vector
   __shared__ double shw[d_shmemsize]; // the other work vector
   __shared__ double shp[d_shmemsize]; // to share Y^T*v

   shv[tdx] = V[tdx];
   wrk = -B[0]*shv[tdx];               // first column of W
   W[tdx] = wrk;

   for(int j=1; j<ncols; j++)          // compute column j of W
   {
      shv[tdx] = V[j*nrows + tdx];     // j-th Householder vector
      for(int k=0; k<j; k++)
      {
         pk = 0.0;                     // k-th component of Y^T*v
         shw[tdx] = V[k*nrows + tdx];  // load V[k][i]
         shp[tdx] = shw[tdx]*shv[tdx]; // V[k][i]*v[i]
         __syncthreads();
         for(int i=0; i<nrows; i++) pk = pk + shp[i];
         if(tdx == k) mypk = pk;
      }
      __syncthreads();
      shp[tdx] = mypk;                 // share p[k]
      __syncthreads();
      zi = 0.0;                        // i-th component of W*p
      for(int k=0; k<j; k++)
      {
         shw[tdx] = W[k*nrows + tdx];  // load W[k][i]
         zi = zi + shw[tdx]*shp[k];
      }
      zi = zi + shv[tdx];
      wrk = -B[j]*zi;
      W[j*nrows + tdx] = wrk;          // wrk is assigned to W[j][tdx]
      __syncthreads();
   }
}

__global__ void cmplx_VB_to_W
 ( int nrows, int ncols, double *B,
   double *Vre, double *Vim, double *Wre, double *Wim )
{
   const int tdx = threadIdx.x;        // index of thread in block
   double wrk_re,wrk_im,pk_re,pk_im,mypk_re,mypk_im,zi_re,zi_im;
   int VWidx;

   __shared__ double shvre[d_shmemsize]; // one work vector
   __shared__ double shvim[d_shmemsize];
   __shared__ double shwre[d_shmemsize]; // the other work vector
   __shared__ double shwim[d_shmemsize];
   __shared__ double shpre[d_shmemsize]; // to share Y^T*v
   __shared__ double shpim[d_shmemsize];

   shvre[tdx] = Vre[tdx];
   shvim[tdx] = Vim[tdx];
   wrk_re = -B[0]*shvre[tdx];            // first column of W
   wrk_im = -B[0]*shvim[tdx];
   Wre[tdx] = wrk_re;
   Wim[tdx] = wrk_im;

   for(int j=1; j<ncols; j++)          // compute column j of W
   {
      VWidx = j*nrows + tdx;
      shvre[tdx] = Vre[VWidx];         // j-th Householder vector
      shvim[tdx] = Vim[VWidx];

      for(int k=0; k<j; k++)
      {
         pk_re = 0.0;                  // k-th component of Y^H*v
         pk_im = 0.0;
         VWidx = k*nrows + tdx;
         shwre[tdx] = Vre[VWidx];      // load V[k][i]
         shwim[tdx] = Vim[VWidx];
         // shp[tdx] = shw[tdx]*shv[tdx]; V[k][i]*v[i], Hermitian transpose!
         shpre[tdx] =   shwre[tdx]*shvre[tdx] + shwim[tdx]*shvim[tdx];
         shpim[tdx] = - shwim[tdx]*shvre[tdx] + shwre[tdx]*shvim[tdx];
         __syncthreads();
         for(int i=0; i<nrows; i++)
         {
            pk_re = pk_re + shpre[i];
            pk_im = pk_im + shpim[i];
         }
         if(tdx == k) mypk_re = pk_re;
         if(tdx == k) mypk_im = pk_im;
      }
      __syncthreads();
      shpre[tdx] = mypk_re;            // share p[k]
      shpim[tdx] = mypk_im;
      __syncthreads();
      zi_re = 0.0;                     // i-th component of W*p
      zi_im = 0.0;
      for(int k=0; k<j; k++)
      {
         VWidx = k*nrows + tdx;
         shwre[tdx] = Wre[VWidx];      // load W[k][i]
         shwim[tdx] = Wim[VWidx];
         // zi = zi + shw[tdx]*shp[k];
         zi_re = zi_re + shvre[tdx]*shpre[k] - shvim[tdx]*shpim[k];
         zi_im = zi_im + shvim[tdx]*shpre[k] + shvre[tdx]*shpim[k];
      }
      zi_re = zi_re + shvre[tdx];
      zi_im = zi_im + shvim[tdx];
      wrk_re = -B[j]*zi_re;
      wrk_im = -B[j]*zi_im;
      VWidx = j*nrows + tdx;
      Wre[VWidx] = wrk_re;             // wrk is assigned to W[j][tdx]
      Wim[VWidx] = wrk_im;
      __syncthreads();
   }
}

__global__ void dbl_small_WYT
 ( int nrows, int szt, double *W, double *Y, double *WYT )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double result = 0.0;
   double a,b;

   for(int k=0; k<szt; k++)
   {
      a = W[k*nrows + row];   // if(nrows == szt) then row = bdx
      b = Y[k*nrows + col];   // if(nrows == szt) then col = tdx
      result = result + a*b;
   }
   __syncthreads();
   WYT[offset] = result;
}

__global__ void cmplx_small_WYT
 ( int nrows, int szt, double *Wre, double *Wim,
   double *Yre, double *Yim, double *WYTre, double *WYTim )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resultre = 0.0;
   double resultim = 0.0;
   double a_re,a_im,b_re,b_im;
   int Widx,Yidx;

   for(int k=0; k<szt; k++)
   {
      Widx = k*nrows + row;
      a_re = Wre[Widx];            // if(nrows == szt) then row = bdx
      a_im = Wim[Widx]; 
      Yidx = k*nrows + col;
      b_re = Yre[Yidx];            // if(nrows == szt) then col = tdx
      b_im = Yim[Yidx];
      // result = result + a*b;
      resultre = resultre + a_re*b_re - a_im*b_im;
      resultim = resultim + a_im*b_re + a_re*b_im;
   }
   __syncthreads();
   WYTre[offset] = resultre;
   WYTim[offset] = resultim;
}

__global__ void dbl_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Q, double *WYT, double *QWYT )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double result = 0.0;
   double a,b;

   for(int k=0; k<rowdim; k++)       // run over rowdim, not just szt
   {                                 // coloff shifts by col*row elements
      a = Q[row*dim + coloff + k];   // row = bdx, if dim == szt, coloff == 0
      b = WYT[k*rowdim + col];       // if(dim == szt) then col = tdx
      result = result + a*b;
   }
   __syncthreads();
   QWYT[offset] = result;            // no column offset in saving QWYT
}

__global__ void cmplx_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qre, double *Qim, double *WYTre, double *WYTim,
   double *QWYTre, double *QWYTim )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resultre = 0.0;
   double resultim = 0.0;
   double a_re,a_im,b_re,b_im;
   int Qidx,WYTidx;

   for(int k=0; k<rowdim; k++)          // run over rowdim, not just szt
   {                                    // coloff shifts by col*row elements
      Qidx = row*dim + coloff + k;
      a_re = Qre[Qidx];                 // row = bdx,
      a_im = Qim[Qidx];                 // if dim == szt, coloff == 0
      WYTidx = k*rowdim + col;
      b_re = WYTre[WYTidx];             // if(dim == szt) then col = tdx
      b_im = WYTim[WYTidx];
      // result = result + a*b;
      resultre = resultre + a_re*b_re - a_im*b_im;
      resultim = resultim + a_im*b_re + a_re*b_im;
   }
   __syncthreads();
   QWYTre[offset] = resultre;           // no column offset in saving QWYT
   QWYTim[offset] = resultim;
}

__global__ void dbl_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWT, double *C, double *YWTC )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C

   double result = 0.0;
   double a,b;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      a = YWT[row*rowdim + k];         // YWT is stored row by row
      b = C[colCoff0 + k];             // but C is stored column by column
      result = result + a*b;
   }
   __syncthreads();
   YWTC[(coloff + col)*nrows + (rowoff + row)] = result;
}

__global__ void cmplx_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWTre, double *YWTim,
   double *Cre, double *Cim, double *YWTCre, double *YWTCim )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C

   double resultre = 0.0;
   double resultim = 0.0;
   double a_re,a_im,b_re,b_im;
   int YWTidx,Cidx;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      YWTidx = row*rowdim + k;
      a_re = YWTre[YWTidx];           // YWT is stored row by row
      a_im = YWTim[YWTidx];
      Cidx = colCoff0 + k;
      b_re = Cre[Cidx];               // but C is stored column by column
      b_im = Cim[Cidx];
      // result = result + a*b;
      resultre = resultre + a_re*b_re - a_im*b_im;
      resultim = resultim + a_im*b_re + a_re*b_im;
   }
   __syncthreads();
   YWTCre[(coloff + col)*nrows + (rowoff + row)] = resultre;
   YWTCim[(coloff + col)*nrows + (rowoff + row)] = resultim;
}

__global__ void dbl_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff, double *Q, double *QWYT )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double a,b;

   a = Q[idx1];       // row = bdx, if dim == szt, coloff == 0
   b = QWYT[offset];  // if(dim == szt) then col = tdx
   a = a + b;

   __syncthreads();
   Q[idx1] = a;
}

__global__ void cmplx_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qre, double *Qim, double *QWYTre, double *QWYTim )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double a_re,a_im,b_re,b_im;

   a_re = Qre[idx1];       // row = bdx, if dim == szt, coloff == 0
   a_im = Qim[idx1];
   b_re = QWYTre[offset];  // if(dim == szt) then col = tdx
   b_im = QWYTim[offset];
   a_re = a_re + b_re;
   a_im = a_im + b_im;

   __syncthreads();
   Qre[idx1] = a_re;
   Qim[idx1] = a_im;
}

__global__ void dbl_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *R, double *YWTC )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double a,b;
   
   a = R[idx];
   b = YWTC[idx];
   a = a + b;
  
   __syncthreads();
   R[idx] = a;
}

__global__ void cmplx_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rre, double *Rim, double *YWTCre, double *YWTCim )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double a_re,a_im,b_re,b_im;
   
   a_re = Rre[idx];
   a_im = Rim[idx];
   b_re = YWTCre[idx];
   b_im = YWTCim[idx];
   a_re = a_re + b_re;
   a_im = a_im + b_im;
  
   __syncthreads();
   Rre[idx] = a_re;
   Rim[idx] = a_im;
}

void GPU_dbl_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *A_h, double *A_d,
   double *v_h, double *V_d, double *beta_h, double *beta_d,
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
      for(int i=0; i<L; i++) v_h[i] = 0.0; // insert zeros
      cudaMemcpy(&V_d[L*nVrows],v_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      beta_h[L] = 0.0; v_h[0] = 1.0;
      cudaMemcpy(&beta_d[L],&beta_h[L],sizeof(double),cudaMemcpyHostToDevice);
      cudaMemcpy(&V_d[L*nVrows+L],v_h,sizeof(double),cudaMemcpyHostToDevice);
   }
   else
   {
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl_small_house<<<1,nrows1>>>
         (&A_d[rowidx],&A_d[rowidx+1],nrows1,nrLog2,
          &V_d[L*nVrows+L],&beta_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&beta_h[L],&beta_d[L],sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(v_h,&V_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : " << beta_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : " << v_h[i] << endl;
   }
}

void GPU_cmplx_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *vre_h, double *vim_h, double *Vre_d, double *Vim_d,
   double *beta_h, double *beta_d, double *lapms, bool verbose )
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
         vre_h[i] = 0.0;
         vim_h[i] = 0.0;
      }
      cudaMemcpy(&Vre_d[L*nVrows],vre_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vim_d[L*nVrows],vim_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      beta_h[L] = 0.0; vre_h[0] = 1.0; vim_h[0] = 0.0;
      cudaMemcpy(&beta_d[L],&beta_h[L],sizeof(double),cudaMemcpyHostToDevice);
      cudaMemcpy(&Vre_d[L*nVrows+L],vre_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vim_d[L*nVrows+L],vim_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      cmplx_small_house<<<1,nrows1>>>
         (&Are_d[rowidx],&Aim_d[rowidx],&Are_d[rowidx+1],&Aim_d[rowidx+1],
          nrows1,nrLog2,&Vre_d[L*nVrows+L],&Vim_d[L*nVrows+L],&beta_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&beta_h[L],&beta_d[L],sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(vre_h,&Vre_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vim_h,&Vim_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : " << beta_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vre_h[i] << "  " << vim_h[i] << endl;
   }
}

void GPU_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
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
   dbl_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,A_d,&V_d[L*nVrows+L],&beta_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(A_h,A_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << A_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *Vre_d, double *Vim_d, double *beta_h, double *beta_d,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);
   cmplx_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,Are_d,Aim_d,
       &Vre_d[L*nVrows+L],&Vim_d[L*nVrows+L],&beta_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Are_h,Are_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aim_h,Aim_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Are_h[j*nrows+i] << "  "
                 << Aim_h[j*nrows+i] << endl;
   }
}

void GPU_dbl_VB_to_W
 ( int nrows, int ncols, int szt,
   double *V_h, double *V_d, double *W_h, double *W_d,
   double *beta_h, double *beta_d, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   dbl_VB_to_W<<<1,nrows>>>(nrows,ncols,beta_d,V_d,W_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = nrows*sizeof(double);
      const size_t szVandW = szt*szhouse;

      cudaMemcpy(beta_h,beta_d,szbeta,cudaMemcpyDeviceToHost);
      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : " << beta_h[j] << endl;
      cudaMemcpy(V_h,V_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<nrows; i++) 
            cout << "V[" << i << "][" << j << "] : " << V_h[ix++] << endl;
      cudaMemcpy(W_h,W_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<nrows; i++) 
            cout << "W[" << i << "][" << j << "] : " << W_h[ix++] << endl;
   }
}

void GPU_cmplx_VB_to_W
 ( int nrows, int ncols, int szt,
   double *Vre_h, double *Vim_h, double *Vre_d, double *Vim_d,
   double *Wre_h, double *Wim_h, double *Wre_d, double *Wim_d,
   double *beta_h, double *beta_d, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   cmplx_VB_to_W<<<1,nrows>>>(nrows,ncols,beta_d,Vre_d,Vim_d,Wre_d,Wim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = nrows*sizeof(double);
      const size_t szVandW = szt*szhouse;

      cudaMemcpy(beta_h,beta_d,szbeta,cudaMemcpyDeviceToHost);
      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : " << beta_h[j] << endl;

      cudaMemcpy(Vre_h,Vre_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vim_h,Vim_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<nrows; i++) 
            cout << "V[" << i << "][" << j << "] : "
                 << Vre_h[ix] << "  " << Vim_h[ix++] << endl;

      cudaMemcpy(Wre_h,Wre_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wim_h,Wim_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<nrows; i++) 
            cout << "W[" << i << "][" << j << "] : "
                 << Wre_h[ix] << "  " << Wim_h[ix++] << endl;
   }
}

void GPU_dbl_small_WYT
 ( int nrows, int szt, double *W_d, double *Y_d, double *WYT_d,
   double *WYT_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   dbl_small_WYT<<<nbrblocks,szt>>>(nrows,szt,W_d,Y_d,WYT_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYT_h,WYT_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYT_h[ix++] << endl;
   }
}

void GPU_cmplx_small_WYT
 ( int nrows, int szt, double *Wre_d, double *Wim_d,
   double *Yre_d, double *Yim_d, double *WYTre_d, double *WYTim_d,
   double *WYTre_h, double *WYTim_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   cmplx_small_WYT<<<nbrblocks,szt>>>
      (nrows,szt,Wre_d,Wim_d,Yre_d,Yim_d,WYTre_d,WYTim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYTre_h,WYTre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTim_h,WYTim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYTre_h[ix] << "  "
                 << WYTim_h[ix++] << endl;
   }
}

void GPU_dbl_small_YWT
 ( int nrows, int szt, int idx, double *Y_d, double *W_d, double *YWT_d,
   double *YWT_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl_small_WYT<<<nbrblocks,szt>>>(rowdim,szt,Y_d,W_d,YWT_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWT_h,YWT_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWT_h[ix++] << endl;
   }
}

void GPU_cmplx_small_YWT
 ( int nrows, int szt, int idx,
   double *Yre_d, double *Yim_d, double *Wre_d, double *Wim_d,
   double *YWTre_d, double *YWTim_d, double *YWTre_h, double *YWTim_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx_small_WYT<<<nbrblocks,szt>>>
      (rowdim,szt,Yre_d,Yim_d,Wre_d,Wim_d,YWTre_d,YWTim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWTre_h,YWTre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTim_h,YWTim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWTre_h[ix] << "  " << YWTim_h[ix++] << endl;
   }
}

void GPU_dbl_small_QWYT
 ( int dim, int szt, int idx, double *Q_d, double *WYT_d, double *QWYT_d,
   double *QWYT_h, double *Q_h, double *lapms, bool verbose )
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

      cudaMemcpy(Q_h,Q_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
            cout << "Q[" << i << "][" << j << "] : "
                 << Q_h[ix++] << endl;
   }

   cudaEventRecord(start);
   dbl_small_QWYT<<<nbrblocks,szt>>>(dim,rowdim,szt,coloff,Q_d,WYT_d,QWYT_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYT_h,QWYT_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYT_h[ix++] << endl;
   }
}

void GPU_cmplx_small_QWYT
 ( int dim, int szt, int idx, double *Qre_d, double *Qim_d,
   double *WYTre_d, double *WYTim_d, double *QWYTre_d, double *QWYTim_d,
   double *QWYTre_h, double *QWYTim_h, double *Qre_h, double *Qim_h,
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

      cudaMemcpy(Qre_h,Qre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qim_h,Qim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
            cout << "Q[" << i << "][" << j << "] : "
                 << Qre_h[ix] << "  " << Qim_h[ix++] << endl;
   }

   cudaEventRecord(start);
   cmplx_small_QWYT<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qre_d,Qim_d,WYTre_d,WYTim_d,
       QWYTre_d,QWYTim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYTre_h,QWYTre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTim_h,QWYTim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYTre_h[ix] << "  " << QWYTim_h[ix++] << endl;
   }
}

void GPU_dbl_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWT_d, double *C_d, double *YWTC_d, double *YWTC_h,
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

      double *C_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(C_h,C_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
            cout << "C_h[" << i << "][" << j << "] : "
                 << C_h[j*nrows+i] << endl;

      free(C_h);
   }

   cudaEventRecord(start);
   dbl_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,YWT_d,C_d,YWTC_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTC_h,YWTC_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "YWTC[" << i << "][" << j << "] : "
                 << YWTC_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTre_d, double *YWTim_d, double *Cre_d, double *Cim_d,
   double *YWTCre_d, double *YWTCim_d, double *YWTCre_h, double *YWTCim_h,
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

      double *Cre_h = new double[nrows*ncols];
      double *Cim_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Cre_h,Cre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cim_h,Cim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
            cout << "C_h[" << i << "][" << j << "] : "
                 << Cre_h[j*nrows+i] << "  "
                 << Cim_h[j*nrows+i] << endl;

      free(Cre_h); free(Cim_h);
   }
   cudaEventRecord(start);
   cmplx_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
       YWTre_d,YWTim_d,Cre_d,Cim_d,YWTCre_d,YWTCim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTCre_h,YWTCre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCim_h,YWTCim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "YWTC[" << i << "][" << j << "] : "
                 << YWTCre_h[j*nrows + i] << "  "
                 << YWTCim_h[j*nrows + i] << endl;
   }
}

void GPU_dbl_small_Qupdate
 ( int dim, int szt, int idx, double *Q_d, double *QWYT_d,
   double *Q_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl_small_Qupdate<<<nbrblocks,szt>>>(dim,rowdim,szt,coloff,Q_d,QWYT_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Q_h,Q_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
            cout << "Q[" << i << "][" << j << "] : "
                 << Q_h[ix++] << endl;
   }
}

void GPU_cmplx_small_Qupdate
 ( int dim, int szt, int idx, double *Qre_d, double *Qim_d,
   double *QWYTre_d, double *QWYTim_d, double *Qre_h, double *Qim_h,
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
   cmplx_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qre_d,Qim_d,QWYTre_d,QWYTim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qre_h,Qre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qim_h,Qim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
            cout << "Q[" << i << "][" << j << "] : "
                 << Qre_h[ix] << "  " << Qim_h[ix++] << endl;
   }
}

void GPU_dbl_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx, double *R_d, double *YWTC_d,
   double *R_h, double *lapms, bool verbose )
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
   dbl_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,R_d,YWTC_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(R_h,R_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << R_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rre_d, double *Rim_d, double *YWTCre_d, double *YWTCim_d,
   double *Rre_h, double *Rim_h, double *lapms, bool verbose )
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
   cmplx_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
       Rre_d,Rim_d,YWTCre_d,YWTCim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rre_h,Rre_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rim_h,Rim_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rre_h[j*nrows + i] << "  "
                 << Rim_h[j*nrows + i] << endl;
   }
}

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, bool verbose )
{
   const int dim = nrows*ncols;         // total number of doubles
   const int nrows2 = nrows*nrows;
   double *A_h = new double[dim];       // matrix A on the host
   double *A_d;                         // matrix on the device
   double *Q_h = new double[nrows2];    // orthogonal Q on host
   double *Q_d;                         // orthogonal Q on device
   double *v_h = new double[nrows];     // Householder vector on host
   double *beta_h = new double[szt];    // beta on the host
   double *beta_d;                      // beta on the device
   double *V_h = new double[nrows*szt]; // matrix of Householder vectors
   double *V_d;                         // Householder vectors on device
   double *W_h = new double[nrows*szt]; // the W matrix on the host
   double *W_d;                         // the W matrix on the device
   double *WYT_h = new double[nrows2];  // W*Y^T on the host
   double *WYT_d;                       // W*Y^T on the device
   double *YWT_h = new double[nrows2];  // Y*W^T on the host
   double *YWT_d;                       // Y*W^T on the device
   double *QWYT_h = new double[nrows2]; // Q*WY^T on the host
   double *QWYT_d;                      // Q*WY^T on the device
   double *YWTC_h = new double[dim];    // YWT*C on the host
   double *YWTC_d;                      // YWT*C on the device

   int ix = 0;                          // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++) A_h[ix++] = A[i][j];

   ix = 0;                              // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
            Q_h[ix++] = 1.0;
         else
            Q_h[ix++] = 0.0;
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&A_d,sznum);
   cudaMemcpy(A_d,A_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&beta_d,szbeta);
   for(int i=0; i<szt; i++) beta_h[i] = 0.0;
   cudaMemcpy(beta_d,beta_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);  // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&V_d,szVandW + szpad); // padding only in allocation
   ix = 0;
   for(int i=0; i<nrows*szt; i++) V_h[ix++] = 0.0; 
   V_h[--ix] = 1.0; // initialize last vector for square tiles
   cudaMemcpy(V_d,V_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&W_d,szVandW + szpad); // padding only in allocation

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYT_d,szWYT + szpad); // padding for W*Y^T product
   cudaMalloc((void**)&Q_d,szWYT);
   cudaMemcpy(Q_d,Q_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYT_d,szWYT);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWT_d,szYWT + szpad); // padding for Y*W^T product
   cudaMalloc((void**)&YWTC_d,sznum + szpad);

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
         GPU_dbl_small_house
            (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
             A_h,A_d,v_h,V_d,beta_h,beta_d,houselapms,verbose);
         GPU_dbl_small_leftRupdate
            (nrows,ncols,szt,colidx,k,L,A_h,A_d,V_d,beta_h,beta_d,
             tileRlapms,verbose);
      }
      // changed nrows into nrows - k*szt and ncols into szt
      GPU_dbl_VB_to_W
         (nrows-k*szt,szt,szt,V_h,V_d,W_h,W_d,beta_h,beta_d,vb2Wlapms,verbose);
      // update Q, WYT matrix has nrows - k*szt instead of nrows
      GPU_dbl_small_WYT(nrows-k*szt,szt,W_d,V_d,WYT_d,WYT_h,WYTlapms,verbose);
      GPU_dbl_small_QWYT
         (nrows,szt,k,Q_d,WYT_d,QWYT_d,QWYT_h,Q_h,QWYTlapms,verbose);
      GPU_dbl_small_Qupdate
         (nrows,szt,k,Q_d,QWYT_d,Q_h,Qaddlapms,verbose);
      if(k < nbt-1)                                           // update R
      {
         GPU_dbl_small_YWT(nrows,szt,k,V_d,W_d,YWT_d,YWT_h,YWTlapms,verbose);
         GPU_dbl_small_YWTC
            (nrows,ncols,szt,k,YWT_d,A_d,YWTC_d,YWTC_h,YWTClapms,verbose);
         GPU_dbl_small_R_add_YWTC
            (nrows,ncols,szt,k,A_d,YWTC_d,A_h,Raddlapms,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Q_h,Q_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++) Q[i][j] = Q_h[ix++];

   cudaMemcpy(A_h,A_d,sznum,cudaMemcpyDeviceToHost);
   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
         R[i][j] = A_h[j*nrows+i];

   free(A_h); free(Q_h); free(v_h); free(V_h); free(W_h);
   free(WYT_h); free(QWYT_h); free(YWT_h); free(YWTC_h);
}

void GPU_cmplx_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, bool verbose )
{
   const int dim = nrows*ncols;           // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Are_h = new double[dim];       // real parts of A on the host
   double *Aim_h = new double[dim];       // imaginary parts of A on the host
   double *Are_d;                         // Are on the device
   double *Aim_d;                         // Aim on the device
   double *Qre_h = new double[nrows2];    // real parts of Q on host
   double *Qim_h = new double[nrows2];    // imaginary parts of Q on host
   double *Qre_d;                         // Qre on device
   double *Qim_d;                         // Qim on device
   double *vre_h = new double[nrows];     // real parts of Householder vector
   double *vim_h = new double[nrows];     // imaginary parts on host
   double *beta_h = new double[szt];      // beta on the host
   double *beta_d;                        // beta on the device
   double *Vre_h = new double[nrows*szt]; // real parts of Householder vectors
   double *Vim_h = new double[nrows*szt]; // imaginary parts
   double *Vre_d;                         // Vre on device
   double *Vim_d;                         // Vim on device
   double *Wre_h = new double[nrows*szt]; // real parts of the W matrix
   double *Wim_h = new double[nrows*szt]; // imaginary parts of the W matrix
   double *Wre_d;                         // Wre on the device
   double *Wim_d;                         // Wim on the device
   double *WYTre_h = new double[nrows2];  // real parts of W*Y^T on the host
   double *WYTim_h = new double[nrows2];  // imaginary parts of W*Y^T
   double *WYTre_d;                       // WYTre on the device 
   double *WYTim_d;                       // WYTim on the device
   double *YWTre_h = new double[nrows2];  // real parts of Y*W^T on the host
   double *YWTim_h = new double[nrows2];  // imginary parts of Y*W^T
   double *YWTre_d;                       // YWTre on the device
   double *YWTim_d;                       // YWTim on the device
   double *QWYTre_h = new double[nrows2]; // real parts of Q*WY^T on the host
   double *QWYTim_h = new double[nrows2]; // imaginary parts of Q*WY^T
   double *QWYTre_d;                      // QWYTre on the device
   double *QWYTim_d;                      // QWYTim on the device
   double *YWTCre_h = new double[dim];    // real parts of YWT*C on the host
   double *YWTCim_h = new double[dim];    // imaginary parts of YWT*C
   double *YWTCre_d;                      // YWTCre on the device
   double *YWTCim_d;                      // YWTCim on the device

   int ix = 0;                            // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Are_h[ix]   = Are[i][j];
         Aim_h[ix++] = Aim[i][j];
      }

   ix = 0;                                // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qre_h[ix]   = 1.0;
            Qim_h[ix++] = 0.0;
         }
         else
         {
            Qre_h[ix]   = 0.0;
            Qim_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Are_d,sznum);
   cudaMalloc((void**)&Aim_d,sznum);
   cudaMemcpy(Are_d,Are_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aim_d,Aim_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&beta_d,szbeta);
   for(int i=0; i<szt; i++) beta_h[i] = 0.0;
   cudaMemcpy(beta_d,beta_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);    // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vre_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Vim_d,szVandW + szpad); // padding added
   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vre_h[ix]   = 0.0; 
      Vim_h[ix++] = 0.0; 
   }
   Vre_h[--ix] = 1.0; // initialize last vector for square tiles
   cudaMemcpy(Vre_d,Vre_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vim_d,Vim_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Wre_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Wim_d,szVandW + szpad); // padding added

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYTre_d,szWYT + szpad); // padding for W*Y^T 
   cudaMalloc((void**)&WYTim_d,szWYT + szpad);
   cudaMalloc((void**)&Qre_d,szWYT);
   cudaMalloc((void**)&Qim_d,szWYT);
   cudaMemcpy(Qre_d,Qre_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qim_d,Qim_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYTre_d,szWYT);
   cudaMalloc((void**)&QWYTim_d,szWYT);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWTre_d,szYWT + szpad); // padding for Y*W^T
   cudaMalloc((void**)&YWTim_d,szYWT + szpad);
   cudaMalloc((void**)&YWTCre_d,sznum + szpad);
   cudaMalloc((void**)&YWTCim_d,sznum + szpad);

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
         GPU_cmplx_small_house
            (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
             Are_h,Aim_h,Are_d,Aim_d,vre_h,vim_h,Vre_d,Vim_d,
             beta_h,beta_d,houselapms,verbose);
         GPU_cmplx_small_leftRupdate
            (nrows,ncols,szt,colidx,k,L,Are_h,Aim_h,Are_d,Aim_d,
             Vre_d,Vim_d,beta_h,beta_d,tileRlapms,verbose);
      }
      GPU_cmplx_VB_to_W
         (nrows-k*szt,szt,szt,Vre_h,Vim_h,Vre_d,Vim_d,Wre_h,Wim_h,
          Wre_d,Wim_d,beta_h,beta_d,vb2Wlapms,verbose);
      GPU_cmplx_small_WYT
         (nrows-k*szt,szt,Wre_d,Wim_d,Vre_d,Vim_d,WYTre_d,WYTim_d,
          WYTre_h,WYTim_h,WYTlapms,verbose);
      GPU_cmplx_small_QWYT
         (nrows,szt,k,Qre_d,Qim_d,WYTre_d,WYTim_d,QWYTre_d,QWYTim_d,
          QWYTre_h,QWYTim_h,Qre_h,Qim_h,QWYTlapms,verbose);
      GPU_cmplx_small_Qupdate
         (nrows,szt,k,Qre_d,Qim_d,QWYTre_d,QWYTim_d,Qre_h,Qim_h,
          Qaddlapms,verbose);
      if(k < nbt-1)                              // update R
      {
         GPU_cmplx_small_YWT
            (nrows,szt,k,Vre_d,Vim_d,Wre_d,Wim_d,YWTre_d,YWTim_d,
             YWTre_h,YWTim_h,YWTlapms,verbose);
         GPU_cmplx_small_YWTC
            (nrows,ncols,szt,k,YWTre_d,YWTim_d,Are_d,Aim_d,YWTCre_d,
             YWTCim_d,YWTCre_h,YWTCim_h,YWTClapms,verbose);
         GPU_cmplx_small_R_add_YWTC
            (nrows,ncols,szt,k,Are_d,Aim_d,YWTCre_d,YWTCim_d,Are_h,Aim_h,
             Raddlapms,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qre_h,Qre_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qim_h,Qim_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qre[i][j] = Qre_h[ix];
         Qim[i][j] = Qim_h[ix++];
      }

   cudaMemcpy(Are_h,Are_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aim_h,Aim_d,sznum,cudaMemcpyDeviceToHost);
   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rre[i][j] = Are_h[j*nrows+i];
         Rim[i][j] = Aim_h[j*nrows+i];
      }

   free(Are_h); free(Aim_h); free(Qre_h); free(Qim_h);
   free(vre_h); free(vim_h); free(Vre_h); free(Vim_h);
   free(Wre_h); free(Wim_h);
   free(WYTre_h); free(QWYTre_h); free(YWTre_h); free(YWTCre_h);
   free(WYTim_h); free(QWYTim_h); free(YWTim_h); free(YWTCim_h);
}
