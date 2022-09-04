// The file dbl4_bals_kernels.cu defines the functions with prototypes in
// the file dbl4_bals_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_baqr_kernels.h"
#include "dbl4_tabs_kernels.h"
#include "dbl4_bals_kernels.h"

using namespace std;

__global__ void dbl4_bals_tail
 ( int ncols, int szt,
   double *Ahihi, double *Alohi, double *Ahilo, double *Alolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajhihi;               // register for Ahihi[idx][j]
   double Ajlohi;               // register for Alohi[idx][j]
   double Ajhilo;               // register for Ahilo[idx][j]
   double Ajlolo;               // register for Alolo[idx][j]
   double xjhihi;               // register for xhihi[j]
   double xjlohi;               // register for xlohi[j]
   double xjhilo;               // register for xhilo[j]
   double xjlolo;               // register for xlolo[j]
   double bihihi = bhihi[idx];  // register for bhihi[idx]
   double bilohi = blohi[idx];  // register for blohi[idx]
   double bihilo = bhilo[idx];  // register for bhilo[idx]
   double bilolo = blolo[idx];  // register for blolo[idx]
   double acchihi,acclohi,acchilo,acclolo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajhihi = Ahihi[offset+j];
      Ajlohi = Alohi[offset+j];
      Ajhilo = Ahilo[offset+j];
      Ajlolo = Alolo[offset+j];
      xjhihi = xhihi[j];
      xjlohi = xlohi[j];
      xjhilo = xhilo[j];
      xjlolo = xlolo[j];
      // bi = bi - Aj*xj;
      qdg_mul(Ajhihi,Ajlohi,Ajhilo,Ajlolo,
              xjhihi,xjlohi,xjhilo,xjlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&bihihi,&bilohi,&bihilo,&bilolo,
              acchihi,acclohi,acchilo,acclolo);
   }
   bhihi[idx] = bihihi;
   blohi[idx] = bilohi;
   bhilo[idx] = bihilo;
   blolo[idx] = bilolo;
}

__global__ void dbl4_bals_qtb
 ( int ncols, int szt,
   double *Qthihi, double *Qtlohi, double *Qthilo, double *Qtlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *rhihi, double *rlohi, double *rhilo, double *rlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjhihi;           // register for Q^Thihi[idx][j]
   double Qjlohi;           // register for Q^Tlohi[idx][j]
   double Qjhilo;           // register for Q^Thilo[idx][j]
   double Qjlolo;           // register for Q^Tlolo[idx][j]
   double bjhihi;           // register for bhihi[j]
   double bjlohi;           // register for blohi[j]
   double bjhilo;           // register for bhilo[j]
   double bjlolo;           // register for blolo[j]
   double rihihi = 0.0;     // register for result, rhihi[idx]
   double rilohi = 0.0;     // register for result, rlohi[idx]
   double rihilo = 0.0;     // register for result, rhilo[idx]
   double rilolo = 0.0;     // register for result, rlolo[idx]
   double acchihi,acclohi,acchilo,acclolo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjhihi = Qthihi[offset+j];
      Qjlohi = Qtlohi[offset+j];
      Qjhilo = Qthilo[offset+j];
      Qjlolo = Qtlolo[offset+j];
      bjhihi = bhihi[j];
      bjlohi = blohi[j];
      bjhilo = bhilo[j];
      bjlolo = blolo[j];
      // ri = ri + Qj*bj;
      qdg_mul(Qjhihi,Qjlohi,Qjhilo,Qjlolo,
              bjhihi,bjlohi,bjhilo,bjlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&rihihi,&rilohi,&rihilo,&rilolo,
              acchihi,acclohi,acchilo,acclolo);
   }
   rhihi[idx] = rihihi;
   rlohi[idx] = rilohi;
   rhilo[idx] = rihilo;
   rlolo[idx] = rilolo;
}

void GPU_dbl4_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo, 
   double *xhihi, double *xlohi, double *xhilo, double *xlolo, bool verbose )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;

   if(verbose) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_dbl4_blocked_houseqr
      (nrows,ncols,szt,nbt,Ahihi,Alohi,Ahilo,Alolo,
       Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if(verbose)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(verbose)
   {
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihi[i][j] << "  " << Rlohi[i][j] << endl << "  "
                 << Rhilo[i][j] << "  " << Rlolo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << bhihi[i] << "  " << blohi[i] << endl << "  "
              << bhilo[i] << "  " << blolo[i] << endl;
   }
   double **workRhihi = new double*[nrows]; // work around because
   double **workRlohi = new double*[nrows]; // solver modifies R ...
   double **workRhilo = new double*[nrows]; 
   double **workRlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRhihi[i] = new double[ncols];
      workRlohi[i] = new double[ncols];
      workRhilo[i] = new double[ncols];
      workRlolo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRhihi[i][j] = Rhihi[i][j];
         workRlohi[i][j] = Rlohi[i][j];
         workRhilo[i][j] = Rhilo[i][j];
         workRlolo[i][j] = Rlolo[i][j];
      }
   }
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   GPU_dbl4_upper_tiled_solver
      (ncols,szt,nbt,workRhihi,workRlohi,workRhilo,workRlolo,
       bhihi,blohi,bhilo,blolo,xhihi,xlohi,xhilo,xlolo,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&addover,&bsmulcnt,&mulover,&bsdivcnt,&divover);

   if(verbose)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihi[i][j] << "  " << Rlohi[i][j] << endl << "  "
                 << Rhilo[i][j] << "  " << Rlolo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
               << bhihi[i] << "  " << blohi[i] << endl << "  "
               << bhilo[i] << "  " << blolo[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRhihi[i]); free(workRlohi[i]);
      free(workRhilo[i]); free(workRlolo[i]);
   }
   free(workRhihi); free(workRlohi);
   free(workRhilo); free(workRlolo);
}

void GPU_dbl4_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   bool verbose )
{
   if(verbose)
   {
      cout << "GPU_dbl4_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshihi[k][i] << "  " << rhslohi[k][i] << "  "
                 << rhshilo[k][i] << "  " << rhslolo[k][i] << endl;
      }
   }
   double *bhihi_d;
   double *blohi_d;
   double *bhilo_d;
   double *blolo_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&bhihi_d,szrhs);
   cudaMalloc((void**)&blohi_d,szrhs);
   cudaMalloc((void**)&bhilo_d,szrhs);
   cudaMalloc((void**)&blolo_d,szrhs);

   double *xhihi_d;
   double *xlohi_d;
   double *xhilo_d;
   double *xlolo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xhihi_d,szsol);
   cudaMalloc((void**)&xlohi_d,szsol);
   cudaMalloc((void**)&xhilo_d,szsol);
   cudaMalloc((void**)&xlolo_d,szsol);
   cudaMemcpy(xhihi_d,solhihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlohi_d,sollohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xhilo_d,solhilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlolo_d,sollolo[stage-1],szsol,cudaMemcpyHostToDevice);

   double *Ahihi_d;
   double *Alohi_d;
   double *Ahilo_d;
   double *Alolo_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Ahihi_d,szmat);
   cudaMalloc((void**)&Alohi_d,szmat);
   cudaMalloc((void**)&Ahilo_d,szmat);
   cudaMalloc((void**)&Alolo_d,szmat);

   double *Ahihi_h = new double[szmat];
   double *Alohi_h = new double[szmat];
   double *Ahilo_h = new double[szmat];
   double *Alolo_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(verbose)
         cout << "GPU_dbl4_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Ahihi_h[idx]   = mathihi[k-stage+1][i][j];
            Alohi_h[idx]   = matlohi[k-stage+1][i][j];
            Ahilo_h[idx]   = mathilo[k-stage+1][i][j];
            Alolo_h[idx++] = matlolo[k-stage+1][i][j];
         }

      cudaMemcpy(bhihi_d,rhshihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blohi_d,rhslohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bhilo_d,rhshilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blolo_d,rhslolo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahihi_d,Ahihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alohi_d,Alohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahilo_d,Ahilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alolo_d,Alolo_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      dbl4_bals_tail<<<nbt,szt>>>
          (ncols,szt,Ahihi_d,Alohi_d,Ahilo_d,Alolo_d,
           xhihi_d,xlohi_d,xhilo_d,xlolo_d,bhihi_d,blohi_d,bhilo_d,blolo_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhshihi[k],bhihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslohi[k],blohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhshilo[k],bhilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslolo[k],blolo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Ahihi_h); free(Alohi_h); free(Ahilo_h); free(Alolo_h);

   if(verbose)
   {
      cout << "GPU_dbl4_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshihi[k][i] << "  " << rhslohi[k][i] << "  "
                 << rhshilo[k][i] << "  " << rhslolo[k][i] << endl;
      }
   }
}

void GPU_dbl4_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo, bool verbose )
{
   double *bhihi_d;
   double *blohi_d;
   double *bhilo_d;
   double *blolo_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&bhihi_d,szrhs);
   cudaMalloc((void**)&blohi_d,szrhs);
   cudaMalloc((void**)&bhilo_d,szrhs);
   cudaMalloc((void**)&blolo_d,szrhs);
   cudaMemcpy(bhihi_d,bhihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blohi_d,blohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bhilo_d,bhilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blolo_d,blolo,szrhs,cudaMemcpyHostToDevice);

   double *rhihi_d;
   double *rlohi_d;
   double *rhilo_d;
   double *rlolo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rhihi_d,szsol);
   cudaMalloc((void**)&rlohi_d,szsol);
   cudaMalloc((void**)&rhilo_d,szsol);
   cudaMalloc((void**)&rlolo_d,szsol);

   double *Qthihi_d;
   double *Qtlohi_d;
   double *Qthilo_d;
   double *Qtlolo_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&Qthihi_d,szmat);
   cudaMalloc((void**)&Qtlohi_d,szmat);
   cudaMalloc((void**)&Qthilo_d,szmat);
   cudaMalloc((void**)&Qtlolo_d,szmat);

   double *Qthihi_h = new double[szmat];
   double *Qtlohi_h = new double[szmat];
   double *Qthilo_h = new double[szmat];
   double *Qtlolo_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         Qthihi_h[idx]   = Qhihi[j][i];
         Qtlohi_h[idx]   = Qlohi[j][i];
         Qthilo_h[idx]   = Qhilo[j][i];
         Qtlolo_h[idx++] = Qlolo[j][i];
      }

   cudaMemcpy(Qthihi_d,Qthihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlohi_d,Qtlohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qthilo_d,Qthilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlolo_d,Qtlolo_h,szmat,cudaMemcpyHostToDevice);

   dbl4_bals_qtb<<<nbt,szt>>>
      (ncols,szt,Qthihi_d,Qtlohi_d,Qthilo_d,Qtlolo_d,
       bhihi_d,blohi_d,bhilo_d,blolo_d,rhihi_d,rlohi_d,rhilo_d,rlolo_d);

   cudaMemcpy(bhihi,rhihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blohi,rlohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhilo,rhilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blolo,rlolo_d,szrhs,cudaMemcpyDeviceToHost);

   free(Qthihi_h); free(Qtlohi_h); free(Qthilo_h); free(Qtlolo_h);
}

void GPU_dbl4_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   const bool bvrb = (vrblvl > 0);

   double **Ahihi = new double*[nrows];
   double **Alohi = new double*[nrows];
   double **Ahilo = new double*[nrows];
   double **Alolo = new double*[nrows];
   double *bhihi = new double[nrows];
   double *blohi = new double[nrows];
   double *bhilo = new double[nrows];
   double *blolo = new double[nrows];
   double *xhihi = new double[ncols];
   double *xlohi = new double[ncols];
   double *xhilo = new double[ncols];
   double *xlolo = new double[ncols];

   double **workRhihi = new double*[nrows]; // GPU upper solver changes R
   double **workRlohi = new double*[nrows];
   double **workRhilo = new double*[nrows];
   double **workRlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRhihi[i] = new double[ncols];
      workRlohi[i] = new double[ncols];
      workRhilo[i] = new double[ncols];
      workRlolo[i] = new double[ncols];
   }
   if(vrblvl)
   {
      cout << "GPU_dbl4_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshihi[k][i] << "  " << rhslohi[k][i] << "  "
                 << rhshilo[k][i] << "  " << rhslolo[k][i] << endl;
      }
   }

   for(int i=0; i<nrows; i++)
   {
      Ahihi[i] = new double[ncols];
      Alohi[i] = new double[ncols];
      Ahilo[i] = new double[ncols];
      Alolo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         Ahihi[i][j] = mathihi[0][i][j];
         Alohi[i][j] = matlohi[0][i][j];
         Ahilo[i][j] = mathilo[0][i][j];
         Alolo[i][j] = matlolo[0][i][j];
      }
      bhihi[i] = rhshihi[0][i];
      blohi[i] = rhslohi[0][i];
      bhilo[i] = rhshilo[0][i];
      blolo[i] = rhslolo[0][i];

      for(int j=0; j<ncols; j++)
      {
         Rhihi[i][j] = mathihi[0][i][j];
         Rlohi[i][j] = matlohi[0][i][j];
         Rhilo[i][j] = mathilo[0][i][j];
         Rlolo[i][j] = matlolo[0][i][j];
      }
   }
   GPU_dbl4_bals_head
      (nrows,ncols,szt,nbt,Ahihi,Alohi,Ahilo,Alolo,
       Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
       bhihi,blohi,bhilo,blolo,xhihi,xlohi,xhilo,xlolo,bvrb);

   for(int j=0; j<ncols; j++)
   {
      solhihi[0][j] = xhihi[j];
      sollohi[0][j] = xlohi[j];
      solhilo[0][j] = xhilo[j];
      sollolo[0][j] = xlolo[j];
   }
   for(int stage=1; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      GPU_dbl4_bals_tail
         (nrows,ncols,szt,nbt,degp1,stage,
          mathihi,matlohi,mathilo,matlolo,rhshihi,rhslohi,rhshilo,rhslolo,
          solhihi,sollohi,solhilo,sollolo,bvrb);

      if(vrblvl > 0)
      {
         cout << "blocks of rhs before assignment :" << endl;
         for(int k=0; k<degp1; k++)
         {
            for(int i=0; i<nrows; i++)
               cout << "rhs[" << k << "][" << i << "] : "
                    << rhshihi[k][i] << "  " << rhslohi[k][i] << "  "
                    << rhshilo[k][i] << "  " << rhslolo[k][i] << endl;
         }
      }

      for(int i=0; i<nrows; i++) 
      {
         cout << "assigning component " << i
              << ", stage = " << stage << endl;
         bhihi[i] = rhshihi[stage][i];
         blohi[i] = rhslohi[stage][i];
         bhilo[i] = rhshilo[stage][i];
         blolo[i] = rhslolo[stage][i];
         cout << "b[" << i << "] : "
              << bhihi[i] << "  " << blohi[i] << endl << "  "
              << bhilo[i] << "  " << blolo[i] << endl;
      }
      double bstimelapsed_d;
      double elapsedms,invlapsed,mullapsed,sublapsed;
      long long int bsaddcnt = 0;
      long long int bsmulcnt = 0;
      long long int bsdivcnt = 0;

      if(vrblvl > 0)
         cout << "-> GPU multiplies rhs with Q^T ..." << endl;

      GPU_dbl4_bals_qtb
         (ncols,szt,nbt,Qhihi,Qlohi,Qhilo,Qlolo,
                        bhihi,blohi,bhilo,blolo,bvrb);

      if(vrblvl > 0)
      {
         for(int i=0; i<nrows; i++)
            cout << "Qtb[" << i << "] : "
                 << bhihi[i] << "  " << blohi[i] << "  "
                 << bhilo[i] << "  " << blolo[i] << endl;
      }
      if(vrblvl > 0)
      {
         cout << "-> GPU solves an upper triangular system ..." << endl;
 
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
               cout << "R[" << i << "][" << j << "] : "
                    << Rhihi[i][j] << "  " << Rlohi[i][j] << "  "
                    << Rhilo[i][j] << "  " << Rlolo[i][j] << endl;
      }
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            workRhihi[i][j] = Rhihi[i][j];
            workRlohi[i][j] = Rlohi[i][j];
            workRhilo[i][j] = Rhilo[i][j];
            workRlolo[i][j] = Rlolo[i][j];
         }

      double addover = 0.0;
      double mulover = 0.0;
      double divover = 0.0;

      GPU_dbl4_upper_tiled_solver
         (ncols,szt,nbt,workRhihi,workRlohi,workRhilo,workRlolo,
          bhihi,blohi,bhilo,blolo,xhihi,xlohi,xhilo,xlolo,
          &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
          &bsaddcnt,&addover,&bsmulcnt,&mulover,&bsdivcnt,&divover);

     if(vrblvl > 0)
        for(int i=0; i<ncols; i++)
           cout << "x[" << i << "] : "
                << xhihi[i] << "  " << xlohi[i] << "  "
                << xhilo[i] << "  " << xlolo[i] << endl;

      for(int j=0; j<ncols; j++)
      {
         solhihi[stage][j] = xhihi[j];
         sollohi[stage][j] = xlohi[j];
         solhilo[stage][j] = xhilo[j];
         sollolo[stage][j] = xlolo[j];
      }
   }

   for(int i=0; i<nrows; i++)
   {
      free(Ahihi[i]); free(Alohi[i]); free(Ahilo[i]); free(Alolo[i]);
      free(workRhihi[i]); free(workRlohi[i]);
      free(workRhilo[i]); free(workRlolo[i]);
   }
   free(Ahihi); free(bhihi); free(xhihi); free(workRhihi);
   free(Alohi); free(blohi); free(xlohi); free(workRlohi);
   free(Ahilo); free(bhilo); free(xhilo); free(workRhilo);
   free(Alolo); free(blolo); free(xlolo); free(workRlolo);
}
