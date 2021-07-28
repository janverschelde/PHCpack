/* The file dbl2_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl2_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "dbl2_baqr_kernels.h"

using namespace std;

__global__ void dbl2_small_house
 ( double *x0hi, double *x0lo, double *x1hi, double *x1lo,
   int dim, int dimLog2,
   double *vhi, double *vlo, double *betahi, double *betalo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_small_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *Rhi, double *Rlo,
   double *vhi, double *vlo, double *betahi, double *betalo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double *Vhi, double *Vlo, double *Whi, double *Wlo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_small_WYT
 ( int nrows, int szt, double *Whi, double *Wlo, double *Yhi, double *Ylo,
   double *WYThi, double *WYTlo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhi, double *Qlo, double *WYThi, double *WYTlo,
   double *QWYThi, double *QWYTlo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWThi, double *YWTlo,
   double *Chi, double *Clo, double *YWTChi, double *YWTClo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_small_Qupdate
 ( int dim, int szt, int coloff,
   double *Qhi, double *Qlo, double *QWYThi, double *QWYTlo )
{
   const int j = threadIdx.x;
}

__global__ void dbl2_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhi, double *Rlo, double *YWTChi, double *YWTClo )
{
   const int j = threadIdx.x;
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
      cudaMemcpy(Ahi_h,Alo_d,sznum,cudaMemcpyDeviceToHost);
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

void GPU_dbl2_small_YWTC
 ( int nrows, int ncols, int szt, int idx, double *YWThi_d, double *YWTlo_d,
   double *Chi_d, double *Clo_d, double *YWTChi_d, double *YWTClo_d,
   double *YWTChi_h, double *YWTClo_h, double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
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
      (dim,szt,coloff,Qhi_d,Qlo_d,QWYThi_d,QWYTlo_d);
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

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYThi_d,szWYT + szpad); // padding for W*Y^T product
   cudaMalloc((void**)&WYTlo_d,szWYT + szpad); 
   cudaMalloc((void**)&Qhi_d,szWYT);
   cudaMalloc((void**)&Qlo_d,szWYT);
   cudaMemcpy(Qhi_d,Qhi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlo_d,Qlo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYThi_d,szWYT);
   cudaMalloc((void**)&QWYTlo_d,szWYT);

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
         GPU_dbl2_small_leftRupdate
            (nrows,ncols,szt,colidx,k,L,Ahi_h,Alo_h,Ahi_d,Alo_d,
             Vhi_d,Vlo_d,betahi_h,betalo_d,betahi_d,betalo_d,
             tileRlapms,verbose);
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
         (nrows,szt,k,Qhi_d,Qlo_d,WYThi_d,WYThi_d,QWYThi_d,QWYTlo_d,
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
   free(Whi_h); free(Wlo_h);
   free(WYThi_h); free(QWYThi_h); free(YWThi_h); free(YWTChi_h);
   free(WYTlo_h); free(QWYTlo_h); free(YWTlo_h); free(YWTClo_h);
}
