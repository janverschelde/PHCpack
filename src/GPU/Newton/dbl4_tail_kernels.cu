// The file dbl4_tail_kernels.cu defines the functions with prototypes in
// the file dbl4_tail_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

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

__global__ void cmplx4_bals_tail
 ( int ncols, int szt,
   double *Arehihi, double *Arelohi, double *Arehilo, double *Arelolo,
   double *Aimhihi, double *Aimlohi, double *Aimhilo, double *Aimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo, 
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajrehihi;             // register for Arehihi[idx][j]
   double Ajrelohi;             // register for Arelohi[idx][j]
   double Ajrehilo;             // register for Arehilo[idx][j]
   double Ajrelolo;             // register for Arelolo[idx][j]
   double Ajimhihi;             // register for Aimhihi[idx][j]
   double Ajimlohi;             // register for Aimlohi[idx][j]
   double Ajimhilo;             // register for Aimhilo[idx][j]
   double Ajimlolo;             // register for Aimlolo[idx][j]
   double xjrehihi;             // register for xrehihi[j]
   double xjrelohi;             // register for xrelohi[j]
   double xjrehilo;             // register for xrehilo[j]
   double xjrelolo;             // register for xrelolo[j]
   double xjimhihi;             // register for ximhihi[j]
   double xjimlohi;             // register for ximlohi[j]
   double xjimhilo;             // register for ximhilo[j]
   double xjimlolo;             // register for ximlolo[j]
   double birehihi = brehihi[idx];  // register for brehihi[idx]
   double birelohi = brelohi[idx];  // register for brelohi[idx]
   double birehilo = brehilo[idx];  // register for brehilo[idx]
   double birelolo = brelolo[idx];  // register for brelolo[idx]
   double biimhihi = bimhihi[idx];  // register for bimhihi[idx]
   double biimlohi = bimlohi[idx];  // register for bimlohi[idx]
   double biimhilo = bimhilo[idx];  // register for bimhilo[idx]
   double biimlolo = bimlolo[idx];  // register for bimlolo[idx]
   double acchihi,acclohi,acchilo,acclolo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajrehihi = Arehihi[offset+j];
      Ajrelohi = Arelohi[offset+j];
      Ajrehilo = Arehilo[offset+j];
      Ajrelolo = Arelolo[offset+j];
      Ajimhihi = Aimhihi[offset+j];
      Ajimlohi = Aimlohi[offset+j];
      Ajimhilo = Aimhilo[offset+j];
      Ajimlolo = Aimlolo[offset+j];
      xjrehihi = xrehihi[j];
      xjrelohi = xrelohi[j];
      xjrehilo = xrehilo[j];
      xjrelolo = xrelolo[j];
      xjimhihi = ximhihi[j];
      xjimlohi = ximlohi[j];
      xjimhilo = ximhilo[j];
      xjimlolo = ximlolo[j];
      // bi = bi - Aj*xj;
      // zre = Ajre*xjre - Ajim*xjim;
      // bire = bire - zre;
      qdg_mul(Ajrehihi,Ajrelohi,Ajrehilo,Ajrelolo,
              xjrehihi,xjrelohi,xjrehilo,xjrelolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&birehihi,&birelohi,&birehilo,&birelolo,
              acchihi,acclohi,acchilo,acclolo);
      qdg_mul(Ajimhihi,Ajimlohi,Ajimhilo,Ajimlolo,
              xjimhihi,xjimlohi,xjimhilo,xjimlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&birehihi,&birelohi,&birehilo,&birelolo,
              acchihi,acclohi,acchilo,acclolo);
      // zim = Ajre*xjim + Ajim*xjre;
      // biim = biim - zim;
      qdg_mul(Ajrehihi,Ajrelohi,Ajrehilo,Ajrelolo,
              xjimhihi,xjimlohi,xjimhilo,xjimlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&biimhihi,&biimlohi,&biimhilo,&biimlolo,
              acchihi,acclohi,acchilo,acclolo);
      qdg_mul(Ajimhihi,Ajimlohi,Ajimhilo,Ajimlolo,
              xjrehihi,xjrelohi,xjrehilo,xjrelolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&biimhihi,&biimlohi,&biimhilo,&biimlolo,
              acchihi,acclohi,acchilo,acclolo);
   }
   brehihi[idx] = birehihi;
   brelohi[idx] = birelohi;
   brehilo[idx] = birehilo;
   brelolo[idx] = birelolo;
   bimhihi[idx] = biimhihi;
   bimlohi[idx] = biimlohi;
   bimhilo[idx] = biimhilo;
   bimlolo[idx] = biimlolo;
}

void write_dbl4_balsflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for b = b - A*X: " << lapsms
        << " milliseconds." << endl;

   long long int flopcnt;
   if(ctype == 0)
      flopcnt = 89*ncols*ncols + 336*ncols*ncols;
      // as many + as * in one inner product
   else
      flopcnt = 4*89*ncols*ncols + 4*336*ncols*ncols;
      // for complex *: 2 ops for +, 6 for *, which is 8 in total

   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 4*ncols*ncols;
   else
      bytecnt = 8*ncols*ncols;

   cout << "    Total number of bytes : " << bytecnt << endl;

   double intensity = ((double) flopcnt)/bytecnt;
   cout << "     Arithmetic intensity : "
        << scientific << setprecision(3) << intensity
        << " #flops/#bytes" << endl;

   double kernflops = 1000.0*((double) flopcnt)/lapsms;
   // double wallflops = ((double) flopcnt)/timelapsed;
   const int gigacnt = pow(2.0,30);

   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
/*
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
 */
}

void GPU_dbl4_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double *totupdlapsedms, int vrblvl )
{
   if(vrblvl > 1)
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
      if(vrblvl > 1)
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

      if(vrblvl > 1)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cudaEvent_t start,stop;       // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl4_bals_tail<<<nbt,szt>>>
          (ncols,szt,Ahihi_d,Alohi_d,Ahilo_d,Alolo_d,
           xhihi_d,xlohi_d,xhilo_d,xlolo_d,bhihi_d,blohi_d,bhilo_d,blolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);

      *totupdlapsedms += milliseconds;

      if(vrblvl > 0) write_dbl4_balsflops(0,ncols,milliseconds);
      
      if(vrblvl > 1)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhshihi[k],bhihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslohi[k],blohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhshilo[k],bhilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslolo[k],blolo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Ahihi_h); free(Alohi_h); free(Ahilo_h); free(Alolo_h);

   if(vrblvl > 1)
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
   cudaFree(bhihi_d); cudaFree(blohi_d);
   cudaFree(bhilo_d); cudaFree(blolo_d);
   cudaFree(xhihi_d); cudaFree(xlohi_d);
   cudaFree(xhilo_d); cudaFree(xlolo_d);
   cudaFree(Ahihi_d); cudaFree(Alohi_d);
   cudaFree(Ahilo_d); cudaFree(Alolo_d);
}

void GPU_cmplx4_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double *totupdlapsedms, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "GPU_cmplx4_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehihi[k][i] << "  " << rhsrelohi[k][i] << endl << "  "
                 << rhsrehilo[k][i] << "  " << rhsrelolo[k][i] << endl << "  "
                 << rhsimhihi[k][i] << "  " << rhsimlohi[k][i] << endl << "  "
                 << rhsimhilo[k][i] << "  " << rhsimlolo[k][i] << endl;
      }
   }
   double *brehihi_d;
   double *brelohi_d;
   double *brehilo_d;
   double *brelolo_d;
   double *bimhihi_d;
   double *bimlohi_d;
   double *bimhilo_d;
   double *bimlolo_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&brehihi_d,szrhs);
   cudaMalloc((void**)&brelohi_d,szrhs);
   cudaMalloc((void**)&brehilo_d,szrhs);
   cudaMalloc((void**)&brelolo_d,szrhs);
   cudaMalloc((void**)&bimhihi_d,szrhs);
   cudaMalloc((void**)&bimlohi_d,szrhs);
   cudaMalloc((void**)&bimhilo_d,szrhs);
   cudaMalloc((void**)&bimlolo_d,szrhs);

   double *xrehihi_d;
   double *xrelohi_d;
   double *xrehilo_d;
   double *xrelolo_d;
   double *ximhihi_d;
   double *ximlohi_d;
   double *ximhilo_d;
   double *ximlolo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xrehihi_d,szsol);
   cudaMalloc((void**)&xrelohi_d,szsol);
   cudaMalloc((void**)&xrehilo_d,szsol);
   cudaMalloc((void**)&xrelolo_d,szsol);
   cudaMalloc((void**)&ximhihi_d,szsol);
   cudaMalloc((void**)&ximlohi_d,szsol);
   cudaMalloc((void**)&ximhilo_d,szsol);
   cudaMalloc((void**)&ximlolo_d,szsol);
   cudaMemcpy(xrehihi_d,solrehihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelohi_d,solrelohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehilo_d,solrehilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelolo_d,solrelolo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhihi_d,solimhihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlohi_d,solimlohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhilo_d,solimhilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlolo_d,solimlolo[stage-1],szsol,cudaMemcpyHostToDevice);

   double *Arehihi_d;
   double *Arelohi_d;
   double *Arehilo_d;
   double *Arelolo_d;
   double *Aimhihi_d;
   double *Aimlohi_d;
   double *Aimhilo_d;
   double *Aimlolo_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Arehihi_d,szmat);
   cudaMalloc((void**)&Arelohi_d,szmat);
   cudaMalloc((void**)&Arehilo_d,szmat);
   cudaMalloc((void**)&Arelolo_d,szmat);
   cudaMalloc((void**)&Aimhihi_d,szmat);
   cudaMalloc((void**)&Aimlohi_d,szmat);
   cudaMalloc((void**)&Aimhilo_d,szmat);
   cudaMalloc((void**)&Aimlolo_d,szmat);

   double *Arehihi_h = new double[szmat];
   double *Arelohi_h = new double[szmat];
   double *Arehilo_h = new double[szmat];
   double *Arelolo_h = new double[szmat];
   double *Aimhihi_h = new double[szmat];
   double *Aimlohi_h = new double[szmat];
   double *Aimhilo_h = new double[szmat];
   double *Aimlolo_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(vrblvl > 1)
         cout << "GPU_cmplx4_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Arehihi_h[idx]   = matrehihi[k-stage+1][i][j];
            Arelohi_h[idx]   = matrelohi[k-stage+1][i][j];
            Arehilo_h[idx]   = matrehilo[k-stage+1][i][j];
            Arelolo_h[idx]   = matrelolo[k-stage+1][i][j];
            Aimhihi_h[idx]   = matimhihi[k-stage+1][i][j];
            Aimlohi_h[idx]   = matimlohi[k-stage+1][i][j];
            Aimhilo_h[idx]   = matimhilo[k-stage+1][i][j];
            Aimlolo_h[idx++] = matimlolo[k-stage+1][i][j];
         }
      
      cudaMemcpy(brehihi_d,rhsrehihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelohi_d,rhsrelohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brehilo_d,rhsrehilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelolo_d,rhsrelolo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhihi_d,rhsimhihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlohi_d,rhsimlohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhilo_d,rhsimhilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlolo_d,rhsimlolo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehihi_d,Arehihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelohi_d,Arelohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehilo_d,Arehilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelolo_d,Arelolo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhihi_d,Aimhihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlohi_d,Aimlohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhilo_d,Aimhilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlolo_d,Aimlolo_h,szmat,cudaMemcpyHostToDevice);

      if(vrblvl > 1)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cudaEvent_t start,stop;       // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      cmplx4_bals_tail<<<nbt,szt>>>
         (ncols,szt,Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
                    Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
          xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
          ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
          brehihi_d,brelohi_d,brehilo_d,brelolo_d,
          bimhihi_d,bimlohi_d,bimhilo_d,bimlolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);

      *totupdlapsedms += milliseconds;

      if(vrblvl > 0) write_dbl4_balsflops(1,ncols,milliseconds);
      
      if(vrblvl > 1)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhsrehihi[k],brehihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelohi[k],brelohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrehilo[k],brehilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelolo[k],brelolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhihi[k],bimhihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlohi[k],bimlohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhilo[k],bimhilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlolo[k],bimlolo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Arehihi_h); free(Aimhihi_h);
   free(Arelohi_h); free(Aimlohi_h);
   free(Arehilo_h); free(Aimhilo_h);
   free(Arelolo_h); free(Aimlolo_h);

   if(vrblvl > 1)
   {
      cout << "GPU_cmplx4_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehihi[k][i] << "  " << rhsrelohi[k][i] << endl << "  "
                 << rhsrehilo[k][i] << "  " << rhsrelolo[k][i] << endl << "  "
                 << rhsimhihi[k][i] << "  " << rhsimlohi[k][i] << endl << "  "
                 << rhsimhilo[k][i] << "  " << rhsimlolo[k][i] << endl;
      }
   }
   cudaFree(brehihi_d); cudaFree(brelohi_d);
   cudaFree(brehilo_d); cudaFree(brelolo_d);
   cudaFree(bimhihi_d); cudaFree(bimlohi_d);
   cudaFree(bimhilo_d); cudaFree(bimlolo_d);
   cudaFree(xrehihi_d); cudaFree(xrelohi_d);
   cudaFree(xrehilo_d); cudaFree(xrelolo_d);
   cudaFree(ximhihi_d); cudaFree(ximlohi_d);
   cudaFree(ximhilo_d); cudaFree(ximlolo_d);
   cudaFree(Arehihi_d); cudaFree(Arelohi_d);
   cudaFree(Arehilo_d); cudaFree(Arelolo_d);
   cudaFree(Aimhihi_d); cudaFree(Aimlohi_d);
   cudaFree(Aimhilo_d); cudaFree(Aimlolo_d);
}

void GPU_dbl4_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   double *lapms, long long int *add, long long int *mul, int vrblvl )
{
   double *rhihi_d;
   double *rlohi_d;
   double *rhilo_d;
   double *rlolo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhihi_d,szrhs);
   cudaMalloc((void**)&rlohi_d,szrhs);
   cudaMalloc((void**)&rhilo_d,szrhs);
   cudaMalloc((void**)&rlolo_d,szrhs);

   double *xhihi_d;
   double *xlohi_d;
   double *xhilo_d;
   double *xlolo_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&xhihi_d,szsol);
   cudaMalloc((void**)&xlohi_d,szsol);
   cudaMalloc((void**)&xhilo_d,szsol);
   cudaMalloc((void**)&xlolo_d,szsol);

   double *Ahihi_d;
   double *Alohi_d;
   double *Ahilo_d;
   double *Alolo_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&Ahihi_d,szmat);
   cudaMalloc((void**)&Alohi_d,szmat);
   cudaMalloc((void**)&Ahilo_d,szmat);
   cudaMalloc((void**)&Alolo_d,szmat);

   double *Ahihi_h = new double[dim*dim];
   double *Alohi_h = new double[dim*dim];
   double *Ahilo_h = new double[dim*dim];
   double *Alolo_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   if(vrblvl > 0)
      cout << "GPU_dbl4_linear_residue for deg+1 : " << degp1 << endl;

   for(int i=tailidx; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rhihi_d,rhshihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlohi_d,rhslohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rhilo_d,rhshilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlolo_d,rhslolo[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=(i-tailidx); j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Ahihi_h[idx]   = mathihi[j][i1][j1];
               Alohi_h[idx]   = matlohi[j][i1][j1];
               Ahilo_h[idx]   = mathilo[j][i1][j1];
               Alolo_h[idx++] = matlolo[j][i1][j1];
            }
      
         cudaMemcpy(Ahihi_d,Ahihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alohi_d,Alohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahilo_d,Ahilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alolo_d,Alolo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(xhihi_d,solhihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlohi_d,sollohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xhilo_d,solhilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlolo_d,sollolo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahihi_d,Ahihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alohi_d,Alohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahilo_d,Ahilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alolo_d,Alolo_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 1)
            cout << "GPU_dbl4_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         dbl4_bals_tail<<<nbt,szt>>>
            (dim,szt,Ahihi_d,Alohi_d,Ahilo_d,Alolo_d,
                     xhihi_d,xlohi_d,xhilo_d,xlolo_d,
                     rhihi_d,rlohi_d,rhilo_d,rlolo_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *lapms += milliseconds;
         flopcount_dbl_bals_tail(dim,add,mul);

         if(vrblvl > 0) write_dbl4_balsflops(0,dim,milliseconds);
      }
      cudaMemcpy(resvechihi[i],rhihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclohi[i],rlohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvechilo[i],rhilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclolo[i],rlolo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<degp1; i++) 
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << solhihi[i][j] << "  " << sollohi[i][j] << endl;
            cout << solhilo[i][j] << "  " << sollolo[i][j] << endl;
         }
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << resvechihi[i][j] << "  " << resveclohi[i][j] << endl;
            cout << resvechilo[i][j] << "  " << resveclolo[i][j] << endl;
         }
      }
   }
   *resmaxhihi = 0.0;
   *resmaxlohi = 0.0;
   *resmaxhilo = 0.0;
   *resmaxlolo = 0.0;
 
   for(int i=tailidx; i<degp1; i++) 
   {
      double *rihihi = resvechihi[i];
      double *rilohi = resveclohi[i];
      double *rihilo = resvechilo[i];
      double *rilolo = resveclolo[i];

      for(int j=0; j<dim; j++)
         if(abs(rihihi[j]) > *resmaxhihi)
         {
            *resmaxhihi = abs(rihihi[j]);
            *resmaxlohi = abs(rilohi[j]);
            *resmaxhilo = abs(rihilo[j]);
            *resmaxlolo = abs(rilolo[j]);
         }
   }
   free(Ahihi_h); free(Alohi_h); free(Ahilo_h); free(Alolo_h);

   cudaFree(rhihi_d); cudaFree(rlohi_d);
   cudaFree(rhilo_d); cudaFree(rlolo_d);
   cudaFree(xhihi_d); cudaFree(xlohi_d);
   cudaFree(xhilo_d); cudaFree(xlolo_d);
   cudaFree(Ahihi_d); cudaFree(Alohi_d);
   cudaFree(Ahilo_d); cudaFree(Alolo_d);
}

void GPU_cmplx4_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi, 
   double **rhsimhilo, double **rhsimlolo, 
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   double *lapms, long long int *add, long long int *mul, int vrblvl )
{
   double *rrehihi_d;
   double *rrelohi_d;
   double *rrehilo_d;
   double *rrelolo_d;
   double *rimhihi_d;
   double *rimlohi_d;
   double *rimhilo_d;
   double *rimlolo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rrehihi_d,szrhs);
   cudaMalloc((void**)&rrelohi_d,szrhs);
   cudaMalloc((void**)&rrehilo_d,szrhs);
   cudaMalloc((void**)&rrelolo_d,szrhs);
   cudaMalloc((void**)&rimhihi_d,szrhs);
   cudaMalloc((void**)&rimlohi_d,szrhs);
   cudaMalloc((void**)&rimhilo_d,szrhs);
   cudaMalloc((void**)&rimlolo_d,szrhs);

   double *xrehihi_d;
   double *xrelohi_d;
   double *xrehilo_d;
   double *xrelolo_d;
   double *ximhihi_d;
   double *ximlohi_d;
   double *ximhilo_d;
   double *ximlolo_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&xrehihi_d,szsol);
   cudaMalloc((void**)&xrelohi_d,szsol);
   cudaMalloc((void**)&xrehilo_d,szsol);
   cudaMalloc((void**)&xrelolo_d,szsol);
   cudaMalloc((void**)&ximhihi_d,szsol);
   cudaMalloc((void**)&ximlohi_d,szsol);
   cudaMalloc((void**)&ximhilo_d,szsol);
   cudaMalloc((void**)&ximlolo_d,szsol);

   double *Arehihi_d;
   double *Arelohi_d;
   double *Arehilo_d;
   double *Arelolo_d;
   double *Aimhihi_d;
   double *Aimlohi_d;
   double *Aimhilo_d;
   double *Aimlolo_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&Arehihi_d,szmat);
   cudaMalloc((void**)&Arelohi_d,szmat);
   cudaMalloc((void**)&Arehilo_d,szmat);
   cudaMalloc((void**)&Arelolo_d,szmat);
   cudaMalloc((void**)&Aimhihi_d,szmat);
   cudaMalloc((void**)&Aimlohi_d,szmat);
   cudaMalloc((void**)&Aimhilo_d,szmat);
   cudaMalloc((void**)&Aimlolo_d,szmat);

   double *Arehihi_h = new double[dim*dim];
   double *Arelohi_h = new double[dim*dim];
   double *Arehilo_h = new double[dim*dim];
   double *Arelolo_h = new double[dim*dim];
   double *Aimhihi_h = new double[dim*dim];
   double *Aimlohi_h = new double[dim*dim];
   double *Aimhilo_h = new double[dim*dim];
   double *Aimlolo_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   if(vrblvl > 0)
      cout << "GPU_cmplx4_linear_residue for deg+1 : " << degp1 << endl;

   for(int i=tailidx; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rrehihi_d,rhsrehihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelohi_d,rhsrelohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrehilo_d,rhsrehilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelolo_d,rhsrelolo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhihi_d,rhsimhihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlohi_d,rhsimlohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhilo_d,rhsimhilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlolo_d,rhsimlolo[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=(i-tailidx); j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Arehihi_h[idx]   = matrehihi[j][i1][j1];
               Arelohi_h[idx]   = matrelohi[j][i1][j1];
               Arehilo_h[idx]   = matrehilo[j][i1][j1];
               Arelolo_h[idx]   = matrelolo[j][i1][j1];
               Aimhihi_h[idx]   = matimhihi[j][i1][j1];
               Aimlohi_h[idx]   = matimlohi[j][i1][j1];
               Aimhilo_h[idx]   = matimhilo[j][i1][j1];
               Aimlolo_h[idx++] = matimlolo[j][i1][j1];
            }
      
         cudaMemcpy(Arehihi_d,Arehihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arelohi_d,Arelohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arehilo_d,Arehilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arelolo_d,Arelolo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimhihi_d,Aimhihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimlohi_d,Aimlohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimhilo_d,Aimhilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimlolo_d,Aimlolo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(xrehihi_d,solrehihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelohi_d,solrelohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrehilo_d,solrehilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelolo_d,solrelolo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhihi_d,solimhihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlohi_d,solimlohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhilo_d,solimhilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlolo_d,solimlolo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(Arehihi_d,Arehihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arelohi_d,Arelohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arehilo_d,Arehilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arelolo_d,Arelolo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimhihi_d,Aimhihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimlohi_d,Aimlohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimhilo_d,Aimhilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimlolo_d,Aimlolo_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 1)
            cout << "GPU_cmplx4_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         cmplx4_bals_tail<<<nbt,szt>>>
            (dim,szt,Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
                     Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
                     xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                     ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
                     rrehihi_d,rrelohi_d,rrehilo_d,rrelolo_d,
                     rimhihi_d,rimlohi_d,rimhilo_d,rimlolo_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *lapms += milliseconds;
         flopcount_cmplx_bals_tail(dim,add,mul);

         if(vrblvl > 0) write_dbl4_balsflops(1,dim,milliseconds);
      }
      cudaMemcpy(resvecrehihi[i],rrehihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelohi[i],rrelohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrehilo[i],rrehilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelolo[i],rrelolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhihi[i],rimhihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlohi[i],rimlohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhilo[i],rimhilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlolo[i],rimlolo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=tailidx; i<degp1; i++)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solrehihi[i][j] << "  " << solrelohi[i][j] << endl << "  "
                 << solrehilo[i][j] << "  " << solrelolo[i][j] << endl << "  "
                 << solimhihi[i][j] << "  " << solimlohi[i][j] << endl << "  "
                 << solimhilo[i][j] << "  " << solimlolo[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << resvecrehihi[i][j] << "  "
                 << resvecrelohi[i][j] << endl << "  "
                 << resvecrehilo[i][j] << "  "
                 << resvecrelolo[i][j] << endl << "  "
                 << resvecimhihi[i][j] << "  "
                 << resvecimlohi[i][j] << endl << "  "
                 << resvecimhilo[i][j] << "  "
                 << resvecimlolo[i][j] << endl;
      }
   }
   *resmaxhihi = 0.0; *resmaxlohi = 0.0;
   *resmaxhilo = 0.0; *resmaxlolo = 0.0;

   for(int i=tailidx; i<degp1; i++)
   {
      double *rirehihi = resvecrehihi[i];
      double *rirelohi = resvecrelohi[i];
      double *rirehilo = resvecrehilo[i];
      double *rirelolo = resvecrelolo[i];
      double *riimhihi = resvecimhihi[i];
      double *riimlohi = resvecimlohi[i];
      double *riimhilo = resvecimhilo[i];
      double *riimlolo = resvecimlolo[i];

      for(int j=0; j<dim; j++)
         if(abs(rirehihi[j]) + abs(riimhihi[j]) > *resmaxhihi)
         {
            *resmaxhihi = abs(rirehihi[j]) + abs(riimhihi[j]);
            *resmaxlohi = abs(rirelohi[j]) + abs(riimlohi[j]);
            *resmaxhilo = abs(rirehilo[j]) + abs(riimhilo[j]);
            *resmaxlolo = abs(rirelolo[j]) + abs(riimlolo[j]);
         }
   }
   free(Arehihi_h); free(Arelohi_h); free(Arehilo_h); free(Arelolo_h);
   free(Aimhihi_h); free(Aimlohi_h); free(Aimhilo_h); free(Aimlolo_h);

   cudaFree(rrehihi_d); cudaFree(rrelohi_d);
   cudaFree(rrehilo_d); cudaFree(rrelolo_d);
   cudaFree(rimhihi_d); cudaFree(rimlohi_d);
   cudaFree(rimhilo_d); cudaFree(rimlolo_d);
   cudaFree(xrehihi_d); cudaFree(xrelohi_d);
   cudaFree(xrehilo_d); cudaFree(xrelolo_d);
   cudaFree(ximhihi_d); cudaFree(ximlohi_d);
   cudaFree(ximhilo_d); cudaFree(ximlolo_d);
   cudaFree(Arehihi_d); cudaFree(Arelohi_d);
   cudaFree(Arehilo_d); cudaFree(Arelolo_d);
   cudaFree(Aimhihi_d); cudaFree(Aimlohi_d);
   cudaFree(Aimhilo_d); cudaFree(Aimlolo_d);
}
