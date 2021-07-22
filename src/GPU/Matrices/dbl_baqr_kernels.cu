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

__global__ void dbl_factors_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *R, double *v, double *beta )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int idx = bdx*szt + tdx;        // global thread index
   const int Roffset = (k+bdx)*szt + k;
   int Rcolidx;
   double w,Rtdx;

   __shared__ double shv[d_shmemsize]; // slice of v

   if(idx < nrows - k)   // nrows - k threads in all blocks work
   {
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
         R[Rcolidx] = Rtdx;
      }
   }
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
      if(tdx < ncols-k) R[Rcolidx] = Rtdx;
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

__global__ void dbl_small_QWYT
 ( int dim, int szt, int coloff, double *Q, double *WYT, double *QWYT )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / dim;
   const int col = offset % dim;       // thread 0 computes QWYT[row][col]

   double result = 0.0;
   double a,b;

   for(int k=0; k<dim; k++)          // run over dim, not just szt
   {                                 // coloff shifts by col*row elements
      a = Q[k*dim + row*(1+coloff)]; // row = bdx, if dim == szt, coloff == 0
      b = WYT[k*dim + col];          // if(dim == szt) then col = tdx
      result = result + a*b;
   }
   __syncthreads();
   QWYT[offset + row*coloff] = result;
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
      a = YWT[row*ncols + k];          // YWT is stored row by row
      b = C[colCoff0 + k];             // but C is stored column by column
      result = result + a*b;
   }
   __syncthreads();
   YWTC[(coloff + col)*nrows + (rowoff + row)] = result;
}

__global__ void dbl_small_Qupdate
 ( int dim, int szt, int coloff, double *Q, double *QWYT )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / dim;

   double result = 0.0;
   double a,b;

   a = Q[row*coloff + offset];     // row = bdx, if dim == szt, coloff == 0
   b = QWYT[row*coloff + offset];  // if(dim == szt) then col = tdx
   result = a + b;

   __syncthreads();
   Q[row*coloff + offset] = result;
}

void GPU_dbl_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *x0_d, double *A_h, double *A_d,
   double *v_h, double *V_d, double *beta_h, double *beta_d,
   double *lapms, bool verbose )
{
   int nrLog2 = ceil(log2((double) nrows1));
   int rowidx = colidx*(nrows+1);       // start of number in A_h

   if(verbose)
   {
      cout << "nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  nbt : " << nbt << endl;
      cout << "k : " << k 
           << "  L : " << L
           << "  nrows1 : " << nrows1
           << "  colidx : " << colidx
           << "  rowidx : " << rowidx << endl;
   }

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaMemcpy(x0_d,&A_h[rowidx],sizeof(double),cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl_small_house<<<1,nrows1>>>
      (x0_d,&A_d[rowidx+1],nrows1,nrLog2,&V_d[L*nrows+L],&beta_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
 
   if(verbose)
   {
      const size_t szhouse = nrows*sizeof(double);

      cudaMemcpy(&beta_h[L],&beta_d[L],sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(v_h,&V_d[L*nrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << L << "] : " << beta_h[L] << endl;
      for(int i=0; i<nrows; i++)
         cout << "v[" << i << "] : " << v_h[i] << endl;
   }
}

void GPU_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);           // 2nd argument: ncols -> szt
   dbl_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,szt,szt,colidx,A_d,&V_d[L*nrows+L],&beta_d[L]);
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

void GPU_dbl_small_WYT
 ( int nrows, int szt, double *W_d, double *Y_d, double *WYT_d,
   double *WYT_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

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

void GPU_dbl_small_YWT
 ( int nrows, int szt, double *Y_d, double *W_d, double *YWT_d,
   double *YWT_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   dbl_small_WYT<<<nbrblocks,szt>>>(nrows,szt,Y_d,W_d,YWT_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(YWT_h,YWT_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWT_h[ix++] << endl;
   }
}

void GPU_dbl_small_QWYT
 ( int dim, int szt, int idx, double *Q_d, double *WYT_d, double *QWYT_d,
   double *QWYT_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl_small_QWYT<<<nbrblocks,szt>>>(dim,szt,coloff,Q_d,WYT_d,QWYT_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(QWYT_h,QWYT_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYT_h[ix++] << endl;
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
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
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
   dbl_small_Qupdate<<<nbrblocks,szt>>>(dim,szt,coloff,Q_d,QWYT_d);
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

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms,
   double *walltimesec, bool verbose )
{
   const int dim = nrows*ncols;         // total number of doubles
   const int nrows2 = nrows*nrows;
   double *A_h = new double[dim];       // matrix A on the host
   double *A_d;                         // matrix on the device
   double *v_h = new double[nrows];     // Householder vector on host
   double *x0_d;                        // first element for house on device
   double *beta_h = new double[szt];    // beta on the host
   double *beta_d;                      // beta on the device
   double *V_h = new double[nrows*szt]; // matrix of Householder vectors
   double *V_d;                         // Householder vectors on device
   double *W_h = new double[nrows*szt]; // the W matrix on the host
   double *W_d;                         // the W matrix 
   double *WYT_h = new double[nrows2];  // W*Y^T on host
   double *WYT_d;                       // W*Y^T on device
   double *YWT_h = new double[nrows2];  // Y*W^T on host
   double *YWT_d;                       // Y*W^T on device
   double *Q_h = new double[nrows2];    // orthogonal Q on host
   double *Q_d;                         // orthogonal Q on device
   double *QWYT_h = new double[nrows2]; // Q*WY^T on host
   double *QWYT_d;                      // Q*WY^T on device
   double *YWTC_h = new double[dim];    // YWT*C on host
   double *YWTC_d;                      // YWT*C on device

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
   cudaMalloc((void**)&x0_d,sizeof(double));

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
   *WYTlapms = 0.0;
   *YWTlapms = 0.0;
   *QWYTlapms = 0.0;
   *YWTClapms = 0.0;
   *Qaddlapms = 0.0;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      int colidx,nrows1;

      for(int L=0; L<szt; L++)  // L runs over the columns in one block
      {
         colidx = k*szt + L;              // index of the current column
         nrows1 = nrows - colidx - 1;     // #rows in Householder vector - 1
         if(nrows1 > 0)
         {
            GPU_dbl_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                x0_d,A_h,A_d,v_h,V_d,beta_h,beta_d,houselapms,verbose);

            GPU_dbl_small_leftRupdate
               (nrows,ncols,szt,colidx,L,A_h,A_d,V_d,beta_h,beta_d,
                tileRlapms,verbose);
         }
      }
      GPU_dbl_VB_to_W
         (nrows,ncols,szt,V_h,V_d,W_h,W_d,beta_h,beta_d,vb2Wlapms,verbose);
      // update Q
      GPU_dbl_small_WYT(nrows,szt,W_d,V_d,WYT_d,WYT_h,WYTlapms,verbose);
      GPU_dbl_small_QWYT
        (nrows,szt,k,Q_d,WYT_d,QWYT_d,QWYT_h,QWYTlapms,verbose);
      GPU_dbl_small_Qupdate
        (nrows,szt,k,Q_d,QWYT_d,Q_h,Qaddlapms,verbose);
      if(k < nbt-1) // update R
      {
         GPU_dbl_small_YWT(nrows,szt,V_d,W_d,YWT_d,YWT_h,YWTlapms,verbose);
         GPU_dbl_small_YWTC
            (nrows,ncols,szt,k,YWT_d,A_d,YWTC_d,YWTC_h,YWTClapms,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Q_h,Q_d,szWYT,cudaMemcpyHostToDevice);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++) Q[i][j] = Q_h[ix++];

   cudaMemcpy(A_h,A_d,sznum,cudaMemcpyHostToDevice);
   ix = 0;                                           // copy columns of R
   for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++) R[i][j] = A_h[ix++];

   free(A_h); free(v_h); free(V_h); free(W_h);
   free(WYT_h); free(QWYT_h); free(YWT_h); free(YWTC_h);
}
