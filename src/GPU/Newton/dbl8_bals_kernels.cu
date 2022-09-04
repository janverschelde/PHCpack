// The file dbl8_bals_kernels.cu defines the functions with prototypes in
// the file dbl8_bals_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#endif
#include "dbl8_baqr_kernels.h"
#include "dbl8_tabs_kernels.h"
#include "dbl8_bals_kernels.h"

using namespace std;

__global__ void dbl8_bals_tail
 ( int ncols, int szt,
   double *Ahihihi, double *Alohihi, double *Ahilohi, double *Alolohi,
   double *Ahihilo, double *Alohilo, double *Ahilolo, double *Alololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajhihihi;                 // register for Ahihihi[idx][j]
   double Ajlohihi;                 // register for Alohihi[idx][j]
   double Ajhilohi;                 // register for Ahilohi[idx][j]
   double Ajlolohi;                 // register for Alolohi[idx][j]
   double Ajhihilo;                 // register for Ahihilo[idx][j]
   double Ajlohilo;                 // register for Alohilo[idx][j]
   double Ajhilolo;                 // register for Ahilolo[idx][j]
   double Ajlololo;                 // register for Alololo[idx][j]
   double xjhihihi;                 // register for xhihihi[j]
   double xjlohihi;                 // register for xlohihi[j]
   double xjhilohi;                 // register for xhilohi[j]
   double xjlolohi;                 // register for xlolohi[j]
   double xjhihilo;                 // register for xhihilo[j]
   double xjlohilo;                 // register for xlohilo[j]
   double xjhilolo;                 // register for xhilolo[j]
   double xjlololo;                 // register for xlololo[j]
   double bihihihi = bhihihi[idx];  // register for bhihihi[idx]
   double bilohihi = blohihi[idx];  // register for blohihi[idx]
   double bihilohi = bhilohi[idx];  // register for bhilohi[idx]
   double bilolohi = blolohi[idx];  // register for blolohi[idx]
   double bihihilo = bhihilo[idx];  // register for bhihilo[idx]
   double bilohilo = blohilo[idx];  // register for blohilo[idx]
   double bihilolo = bhilolo[idx];  // register for bhilolo[idx]
   double bilololo = blololo[idx];  // register for blololo[idx]
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajhihihi = Ahihihi[offset+j];
      Ajlohihi = Alohihi[offset+j];
      Ajhilohi = Ahilohi[offset+j];
      Ajlolohi = Alolohi[offset+j];
      Ajhihilo = Ahihilo[offset+j];
      Ajlohilo = Alohilo[offset+j];
      Ajhilolo = Ahilolo[offset+j];
      Ajlololo = Alololo[offset+j];
      xjhihihi = xhihihi[j];
      xjlohihi = xlohihi[j];
      xjhilohi = xhilohi[j];
      xjlolohi = xlolohi[j];
      xjhihilo = xhihilo[j];
      xjlohilo = xlohilo[j];
      xjhilolo = xhilolo[j];
      xjlololo = xlololo[j];
      // bi = bi - Aj*xj;
      odg_mul(Ajhihihi,Ajlohihi,Ajhilohi,Ajlolohi,
              Ajhihilo,Ajlohilo,Ajhilolo,Ajlololo,
              xjhihihi,xjlohihi,xjhilohi,xjlolohi,
              xjhihilo,xjlohilo,xjhilolo,xjlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&bihihihi,&bilohihi,&bihilohi,&bilolohi,
              &bihihilo,&bilohilo,&bihilolo,&bilololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
   }
   bhihihi[idx] = bihihihi;
   blohihi[idx] = bilohihi;
   bhilohi[idx] = bihilohi;
   blolohi[idx] = bilolohi;
   bhihilo[idx] = bihihilo;
   blohilo[idx] = bilohilo;
   bhilolo[idx] = bihilolo;
   blololo[idx] = bilololo;
}

__global__ void dbl8_bals_qtb
 ( int ncols, int szt,
   double *Qthihihi, double *Qtlohihi, double *Qthilohi, double *Qtlolohi,
   double *Qthihilo, double *Qtlohilo, double *Qthilolo, double *Qtlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *rhihihi, double *rlohihi, double *rhilohi, double *rlolohi,
   double *rhihilo, double *rlohilo, double *rhilolo, double *rlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjhihihi;           // register for Q^Thihihi[idx][j]
   double Qjlohihi;           // register for Q^Tlohihi[idx][j]
   double Qjhilohi;           // register for Q^Thilohi[idx][j]
   double Qjlolohi;           // register for Q^Tlolohi[idx][j]
   double Qjhihilo;           // register for Q^Thihilo[idx][j]
   double Qjlohilo;           // register for Q^Tlohilo[idx][j]
   double Qjhilolo;           // register for Q^Thilolo[idx][j]
   double Qjlololo;           // register for Q^Tlololo[idx][j]
   double bjhihihi;           // register for bhihihi[j]
   double bjlohihi;           // register for blohihi[j]
   double bjhilohi;           // register for bhilohi[j]
   double bjlolohi;           // register for blolohi[j]
   double bjhihilo;           // register for bhihilo[j]
   double bjlohilo;           // register for blohilo[j]
   double bjhilolo;           // register for bhilolo[j]
   double bjlololo;           // register for blololo[j]
   double rihihihi = 0.0;     // register for result, rhihihi[idx]
   double rilohihi = 0.0;     // register for result, rlohihi[idx]
   double rihilohi = 0.0;     // register for result, rhilohi[idx]
   double rilolohi = 0.0;     // register for result, rlolohi[idx]
   double rihihilo = 0.0;     // register for result, rhihilo[idx]
   double rilohilo = 0.0;     // register for result, rlohilo[idx]
   double rihilolo = 0.0;     // register for result, rhilolo[idx]
   double rilololo = 0.0;     // register for result, rlololo[idx]
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjhihihi = Qthihihi[offset+j];
      Qjlohihi = Qtlohihi[offset+j];
      Qjhilohi = Qthilohi[offset+j];
      Qjlolohi = Qtlolohi[offset+j];
      Qjhihilo = Qthihilo[offset+j];
      Qjlohilo = Qtlohilo[offset+j];
      Qjhilolo = Qthilolo[offset+j];
      Qjlololo = Qtlololo[offset+j];
      bjhihihi = bhihihi[j];
      bjlohihi = blohihi[j];
      bjhilohi = bhilohi[j];
      bjlolohi = blolohi[j];
      bjhihilo = bhihilo[j];
      bjlohilo = blohilo[j];
      bjhilolo = bhilolo[j];
      bjlololo = blololo[j];
      // ri = ri + Qj*bj;
      odg_mul(Qjhihihi,Qjlohihi,Qjhilohi,Qjlolohi,
              Qjhihilo,Qjlohilo,Qjhilolo,Qjlololo,
              bjhihihi,bjlohihi,bjhilohi,bjlolohi,
              bjhihilo,bjlohilo,bjhilolo,bjlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&rihihihi,&rilohihi,&rihilohi,&rilolohi,
              &rihihilo,&rilohilo,&rihilolo,&rilololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
   }
   rhihihi[idx] = rihihihi;
   rlohihi[idx] = rilohihi;
   rhilohi[idx] = rihilohi;
   rlolohi[idx] = rilolohi;
   rhihilo[idx] = rihihilo;
   rlohilo[idx] = rilohilo;
   rhilolo[idx] = rihilolo;
   rlololo[idx] = rilololo;
}

void GPU_dbl8_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   bool verbose )
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

   GPU_dbl8_blocked_houseqr
      (nrows,ncols,szt,nbt,
       Ahihihi,Alohihi,Ahilohi,Alolohi,Ahihilo,Alohilo,Ahilolo,Alololo,
       Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
       Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
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
                 << Rhihihi[i][j] << "  " << Rlohihi[i][j] << endl << "  "
                 << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl << "  "
                 << Rhihilo[i][j] << "  " << Rlohilo[i][j] << endl << "  "
                 << Rhilolo[i][j] << "  " << Rlololo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << bhihihi[i] << "  " << blohihi[i] << endl << "  "
              << bhilohi[i] << "  " << blolohi[i] << endl << "  "
              << bhihilo[i] << "  " << blohilo[i] << endl << "  "
              << bhilolo[i] << "  " << blololo[i] << endl;
   }
   double **workRhihihi = new double*[nrows]; // work around because
   double **workRlohihi = new double*[nrows]; // solver modifies R ...
   double **workRhilohi = new double*[nrows]; 
   double **workRlolohi = new double*[nrows];
   double **workRhihilo = new double*[nrows];
   double **workRlohilo = new double*[nrows];
   double **workRhilolo = new double*[nrows]; 
   double **workRlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRhihihi[i] = new double[ncols];
      workRlohihi[i] = new double[ncols];
      workRhilohi[i] = new double[ncols];
      workRlolohi[i] = new double[ncols];
      workRhihilo[i] = new double[ncols];
      workRlohilo[i] = new double[ncols];
      workRhilolo[i] = new double[ncols];
      workRlololo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRhihihi[i][j] = Rhihihi[i][j];
         workRlohihi[i][j] = Rlohihi[i][j];
         workRhilohi[i][j] = Rhilohi[i][j];
         workRlolohi[i][j] = Rlolohi[i][j];
         workRhihilo[i][j] = Rhihilo[i][j];
         workRlohilo[i][j] = Rlohilo[i][j];
         workRhilolo[i][j] = Rhilolo[i][j];
         workRlololo[i][j] = Rlololo[i][j];
      }
   }
   GPU_dbl8_upper_tiled_solver
      (ncols,szt,nbt,
       workRhihihi,workRlohihi,workRhilohi,workRlolohi,
       workRhihilo,workRlohilo,workRhilolo,workRlololo,
       bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
       xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);

   if(verbose)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;

      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihihi[i][j] << "  " << Rlohihi[i][j] << endl << "  "
                 << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl << "  "
                 << Rhihilo[i][j] << "  " << Rlohilo[i][j] << endl << "  "
                 << Rhilolo[i][j] << "  " << Rlololo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
               << bhihihi[i] << "  " << blohihi[i] << endl << "  "
               << bhilohi[i] << "  " << blolohi[i] << endl << "  "
               << bhihilo[i] << "  " << blohilo[i] << endl << "  "
               << bhilolo[i] << "  " << blololo[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRhihihi[i]); free(workRlohihi[i]);
      free(workRhilohi[i]); free(workRlolohi[i]);
      free(workRhihilo[i]); free(workRlohilo[i]);
      free(workRhilolo[i]); free(workRlololo[i]);
   }
   free(workRhihihi); free(workRlohihi);
   free(workRhilohi); free(workRlolohi);
   free(workRhihilo); free(workRlohilo);
   free(workRhilolo); free(workRlololo);
}

void GPU_dbl8_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo, bool verbose )
{
   if(verbose)
   {
      cout << "GPU_dbl8_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshihihi[k][i] << "  " << rhslohihi[k][i] << "  "
                 << rhshilohi[k][i] << "  " << rhslolohi[k][i] << endl
                 << "  "
                 << rhshihilo[k][i] << "  " << rhslohilo[k][i] << "  "
                 << rhshilolo[k][i] << "  " << rhslololo[k][i] << endl;
      }
   }
   double *bhihihi_d;
   double *blohihi_d;
   double *bhilohi_d;
   double *blolohi_d;
   double *bhihilo_d;
   double *blohilo_d;
   double *bhilolo_d;
   double *blololo_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&bhihihi_d,szrhs);
   cudaMalloc((void**)&blohihi_d,szrhs);
   cudaMalloc((void**)&bhilohi_d,szrhs);
   cudaMalloc((void**)&blolohi_d,szrhs);
   cudaMalloc((void**)&bhihilo_d,szrhs);
   cudaMalloc((void**)&blohilo_d,szrhs);
   cudaMalloc((void**)&bhilolo_d,szrhs);
   cudaMalloc((void**)&blololo_d,szrhs);

   double *xhihihi_d;
   double *xlohihi_d;
   double *xhilohi_d;
   double *xlolohi_d;
   double *xhihilo_d;
   double *xlohilo_d;
   double *xhilolo_d;
   double *xlololo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xhihihi_d,szsol);
   cudaMalloc((void**)&xlohihi_d,szsol);
   cudaMalloc((void**)&xhilohi_d,szsol);
   cudaMalloc((void**)&xlolohi_d,szsol);
   cudaMalloc((void**)&xhihilo_d,szsol);
   cudaMalloc((void**)&xlohilo_d,szsol);
   cudaMalloc((void**)&xhilolo_d,szsol);
   cudaMalloc((void**)&xlololo_d,szsol);
   cudaMemcpy(xhihihi_d,solhihihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlohihi_d,sollohihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xhilohi_d,solhilohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlolohi_d,sollolohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xhihilo_d,solhihilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlohilo_d,sollohilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xhilolo_d,solhilolo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlololo_d,sollololo[stage-1],szsol,cudaMemcpyHostToDevice);

   double *Ahihihi_d;
   double *Alohihi_d;
   double *Ahilohi_d;
   double *Alolohi_d;
   double *Ahihilo_d;
   double *Alohilo_d;
   double *Ahilolo_d;
   double *Alololo_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Ahihihi_d,szmat);
   cudaMalloc((void**)&Alohihi_d,szmat);
   cudaMalloc((void**)&Ahilohi_d,szmat);
   cudaMalloc((void**)&Alolohi_d,szmat);
   cudaMalloc((void**)&Ahihilo_d,szmat);
   cudaMalloc((void**)&Alohilo_d,szmat);
   cudaMalloc((void**)&Ahilolo_d,szmat);
   cudaMalloc((void**)&Alololo_d,szmat);

   double *Ahihihi_h = new double[szmat];
   double *Alohihi_h = new double[szmat];
   double *Ahilohi_h = new double[szmat];
   double *Alolohi_h = new double[szmat];
   double *Ahihilo_h = new double[szmat];
   double *Alohilo_h = new double[szmat];
   double *Ahilolo_h = new double[szmat];
   double *Alololo_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(verbose)
         cout << "GPU_dbl8_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Ahihihi_h[idx]   = mathihihi[k-stage+1][i][j];
            Alohihi_h[idx]   = matlohihi[k-stage+1][i][j];
            Ahilohi_h[idx]   = mathilohi[k-stage+1][i][j];
            Alolohi_h[idx]   = matlolohi[k-stage+1][i][j];
            Ahihilo_h[idx]   = mathihilo[k-stage+1][i][j];
            Alohilo_h[idx]   = matlohilo[k-stage+1][i][j];
            Ahilolo_h[idx]   = mathilolo[k-stage+1][i][j];
            Alololo_h[idx++] = matlololo[k-stage+1][i][j];
         }

      cudaMemcpy(bhihihi_d,rhshihihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blohihi_d,rhslohihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bhilohi_d,rhshilohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blolohi_d,rhslolohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bhihilo_d,rhshihilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blohilo_d,rhslohilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bhilolo_d,rhshilolo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blololo_d,rhslololo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahihihi_d,Ahihihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alohihi_d,Alohihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahilohi_d,Ahilohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alolohi_d,Alolohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahihilo_d,Ahihilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alohilo_d,Alohilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahilolo_d,Ahilolo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alololo_d,Alololo_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      dbl8_bals_tail<<<nbt,szt>>>
          (ncols,szt,
           Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
           Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
           xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
           xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d,
           bhihihi_d,blohihi_d,bhilohi_d,blolohi_d,
           bhihilo_d,blohilo_d,bhilolo_d,blololo_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhshihihi[k],bhihihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslohihi[k],blohihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhshilohi[k],bhilohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslolohi[k],blolohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhshihilo[k],bhihilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslohilo[k],blohilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhshilolo[k],bhilolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslololo[k],blololo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Ahihihi_h); free(Alohihi_h); free(Ahilohi_h); free(Alolohi_h);
   free(Ahihilo_h); free(Alohilo_h); free(Ahilolo_h); free(Alololo_h);

   if(verbose)
   {
      cout << "GPU_dbl8_bals_tail copied blocks of rhs :" << endl;

      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshihihi[k][i] << "  " << rhslohihi[k][i] << "  "
                 << rhshilohi[k][i] << "  " << rhslolohi[k][i] << endl
                 << "  "
                 << rhshihilo[k][i] << "  " << rhslohilo[k][i] << "  "
                 << rhshilolo[k][i] << "  " << rhslololo[k][i] << endl;
      }
   }
}

void GPU_dbl8_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   bool verbose )
{
   double *bhihihi_d;
   double *blohihi_d;
   double *bhilohi_d;
   double *blolohi_d;
   double *bhihilo_d;
   double *blohilo_d;
   double *bhilolo_d;
   double *blololo_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&bhihihi_d,szrhs);
   cudaMalloc((void**)&blohihi_d,szrhs);
   cudaMalloc((void**)&bhilohi_d,szrhs);
   cudaMalloc((void**)&blolohi_d,szrhs);
   cudaMalloc((void**)&bhihilo_d,szrhs);
   cudaMalloc((void**)&blohilo_d,szrhs);
   cudaMalloc((void**)&bhilolo_d,szrhs);
   cudaMalloc((void**)&blololo_d,szrhs);
   cudaMemcpy(bhihihi_d,bhihihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blohihi_d,blohihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bhilohi_d,bhilohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blolohi_d,blolohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bhihilo_d,bhihilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blohilo_d,blohilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bhilolo_d,bhilolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blololo_d,blololo,szrhs,cudaMemcpyHostToDevice);

   double *rhihihi_d;
   double *rlohihi_d;
   double *rhilohi_d;
   double *rlolohi_d;
   double *rhihilo_d;
   double *rlohilo_d;
   double *rhilolo_d;
   double *rlololo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rhihihi_d,szsol);
   cudaMalloc((void**)&rlohihi_d,szsol);
   cudaMalloc((void**)&rhilohi_d,szsol);
   cudaMalloc((void**)&rlolohi_d,szsol);
   cudaMalloc((void**)&rhihilo_d,szsol);
   cudaMalloc((void**)&rlohilo_d,szsol);
   cudaMalloc((void**)&rhilolo_d,szsol);
   cudaMalloc((void**)&rlololo_d,szsol);

   double *Qthihihi_d;
   double *Qtlohihi_d;
   double *Qthilohi_d;
   double *Qtlolohi_d;
   double *Qthihilo_d;
   double *Qtlohilo_d;
   double *Qthilolo_d;
   double *Qtlololo_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&Qthihihi_d,szmat);
   cudaMalloc((void**)&Qtlohihi_d,szmat);
   cudaMalloc((void**)&Qthilohi_d,szmat);
   cudaMalloc((void**)&Qtlolohi_d,szmat);
   cudaMalloc((void**)&Qthihilo_d,szmat);
   cudaMalloc((void**)&Qtlohilo_d,szmat);
   cudaMalloc((void**)&Qthilolo_d,szmat);
   cudaMalloc((void**)&Qtlololo_d,szmat);

   double *Qthihihi_h = new double[szmat];
   double *Qtlohihi_h = new double[szmat];
   double *Qthilohi_h = new double[szmat];
   double *Qtlolohi_h = new double[szmat];
   double *Qthihilo_h = new double[szmat];
   double *Qtlohilo_h = new double[szmat];
   double *Qthilolo_h = new double[szmat];
   double *Qtlololo_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         Qthihihi_h[idx]   = Qhihihi[j][i];
         Qtlohihi_h[idx]   = Qlohihi[j][i];
         Qthilohi_h[idx]   = Qhilohi[j][i];
         Qtlolohi_h[idx]   = Qlolohi[j][i];
         Qthihilo_h[idx]   = Qhihilo[j][i];
         Qtlohilo_h[idx]   = Qlohilo[j][i];
         Qthilolo_h[idx]   = Qhilolo[j][i];
         Qtlololo_h[idx++] = Qlololo[j][i];
      }

   cudaMemcpy(Qthihihi_d,Qthihihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlohihi_d,Qtlohihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qthilohi_d,Qthilohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlolohi_d,Qtlolohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qthihilo_d,Qthihilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlohilo_d,Qtlohilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qthilolo_d,Qthilolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlololo_d,Qtlololo_h,szmat,cudaMemcpyHostToDevice);

   dbl8_bals_qtb<<<nbt,szt>>>
      (ncols,szt,
       Qthihihi_d,Qtlohihi_d,Qthilohi_d,Qtlolohi_d,
       Qthihilo_d,Qtlohilo_d,Qthilolo_d,Qtlololo_d,
       bhihihi_d,blohihi_d,bhilohi_d,blolohi_d,
       bhihilo_d,blohilo_d,bhilolo_d,blololo_d,
       rhihihi_d,rlohihi_d,rhilohi_d,rlolohi_d,
       rhihilo_d,rlohilo_d,rhilolo_d,rlololo_d);

   cudaMemcpy(bhihihi,rhihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blohihi,rlohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhilohi,rhilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blolohi,rlolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhihilo,rhihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blohilo,rlohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhilolo,rhilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blololo,rlololo_d,szrhs,cudaMemcpyDeviceToHost);

   free(Qthihihi_h); free(Qtlohihi_h); free(Qthilohi_h); free(Qtlolohi_h);
   free(Qthihilo_h); free(Qtlohilo_h); free(Qthilolo_h); free(Qtlololo_h);
}

void GPU_dbl8_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   const bool bvrb = (vrblvl > 0);

   double **Ahihihi = new double*[nrows];
   double **Alohihi = new double*[nrows];
   double **Ahilohi = new double*[nrows];
   double **Alolohi = new double*[nrows];
   double **Ahihilo = new double*[nrows];
   double **Alohilo = new double*[nrows];
   double **Ahilolo = new double*[nrows];
   double **Alololo = new double*[nrows];
   double *bhihihi = new double[nrows];
   double *blohihi = new double[nrows];
   double *bhilohi = new double[nrows];
   double *blolohi = new double[nrows];
   double *bhihilo = new double[nrows];
   double *blohilo = new double[nrows];
   double *bhilolo = new double[nrows];
   double *blololo = new double[nrows];
   double *xhihihi = new double[ncols];
   double *xlohihi = new double[ncols];
   double *xhilohi = new double[ncols];
   double *xlolohi = new double[ncols];
   double *xhihilo = new double[ncols];
   double *xlohilo = new double[ncols];
   double *xhilolo = new double[ncols];
   double *xlololo = new double[ncols];

   double **workRhihihi = new double*[nrows]; // GPU upper solver changes R
   double **workRlohihi = new double*[nrows];
   double **workRhilohi = new double*[nrows];
   double **workRlolohi = new double*[nrows];
   double **workRhihilo = new double*[nrows];
   double **workRlohilo = new double*[nrows];
   double **workRhilolo = new double*[nrows];
   double **workRlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRhihihi[i] = new double[ncols];
      workRlohihi[i] = new double[ncols];
      workRhilohi[i] = new double[ncols];
      workRlolohi[i] = new double[ncols];
      workRhihilo[i] = new double[ncols];
      workRlohilo[i] = new double[ncols];
      workRhilolo[i] = new double[ncols];
      workRlololo[i] = new double[ncols];
   }
   if(vrblvl)
   {
      cout << "GPU_dbl8_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshihihi[k][i] << "  " << rhslohihi[k][i] << "  "
                 << rhshilohi[k][i] << "  " << rhslolohi[k][i] << endl
                 << "  "
                 << rhshihilo[k][i] << "  " << rhslohilo[k][i] << "  "
                 << rhshilolo[k][i] << "  " << rhslololo[k][i] << endl;
      }
   }

   for(int i=0; i<nrows; i++)
   {
      Ahihihi[i] = new double[ncols];
      Alohihi[i] = new double[ncols];
      Ahilohi[i] = new double[ncols];
      Alolohi[i] = new double[ncols];
      Ahihilo[i] = new double[ncols];
      Alohilo[i] = new double[ncols];
      Ahilolo[i] = new double[ncols];
      Alololo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         Ahihihi[i][j] = mathihihi[0][i][j];
         Alohihi[i][j] = matlohihi[0][i][j];
         Ahilohi[i][j] = mathilohi[0][i][j];
         Alolohi[i][j] = matlolohi[0][i][j];
         Ahihilo[i][j] = mathihilo[0][i][j];
         Alohilo[i][j] = matlohilo[0][i][j];
         Ahilolo[i][j] = mathilolo[0][i][j];
         Alololo[i][j] = matlololo[0][i][j];
      }
      bhihihi[i] = rhshihihi[0][i];
      blohihi[i] = rhslohihi[0][i];
      bhilohi[i] = rhshilohi[0][i];
      blolohi[i] = rhslolohi[0][i];
      bhihilo[i] = rhshihilo[0][i];
      blohilo[i] = rhslohilo[0][i];
      bhilolo[i] = rhshilolo[0][i];
      blololo[i] = rhslololo[0][i];

      for(int j=0; j<ncols; j++)
      {
         Rhihihi[i][j] = mathihihi[0][i][j];
         Rlohihi[i][j] = matlohihi[0][i][j];
         Rhilohi[i][j] = mathilohi[0][i][j];
         Rlolohi[i][j] = matlolohi[0][i][j];
         Rhihilo[i][j] = mathihilo[0][i][j];
         Rlohilo[i][j] = matlohilo[0][i][j];
         Rhilolo[i][j] = mathilolo[0][i][j];
         Rlololo[i][j] = matlololo[0][i][j];
      }
   }
   GPU_dbl8_bals_head
      (nrows,ncols,szt,nbt,
       Ahihihi,Alohihi,Ahilohi,Alolohi,Ahihilo,Alohilo,Ahilolo,Alololo,
       Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
       Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
       bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
       xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,bvrb);

   for(int j=0; j<ncols; j++)
   {
      solhihihi[0][j] = xhihihi[j];
      sollohihi[0][j] = xlohihi[j];
      solhilohi[0][j] = xhilohi[j];
      sollolohi[0][j] = xlolohi[j];
      solhihilo[0][j] = xhihilo[j];
      sollohilo[0][j] = xlohilo[j];
      solhilolo[0][j] = xhilolo[j];
      sollololo[0][j] = xlololo[j];
   }
   for(int stage=1; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      GPU_dbl8_bals_tail
         (nrows,ncols,szt,nbt,degp1,stage,
          mathihihi,matlohihi,mathilohi,matlolohi,
          mathihilo,matlohilo,mathilolo,matlololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,bvrb);

      if(vrblvl > 0)
      {
         cout << "blocks of rhs before assignment :" << endl;
         for(int k=0; k<degp1; k++)
         {
            for(int i=0; i<nrows; i++)
               cout << "rhs[" << k << "][" << i << "] : "
                    << rhshihihi[k][i] << "  " << rhslohihi[k][i] << "  "
                    << rhshilohi[k][i] << "  " << rhslolohi[k][i] << endl
                    << "  "
                    << rhshihilo[k][i] << "  " << rhslohilo[k][i] << "  "
                    << rhshilolo[k][i] << "  " << rhslololo[k][i] << endl;
         }
      }

      for(int i=0; i<nrows; i++) 
      {
         cout << "assigning component " << i
              << ", stage = " << stage << endl;

         bhihihi[i] = rhshihihi[stage][i];
         blohihi[i] = rhslohihi[stage][i];
         bhilohi[i] = rhshilohi[stage][i];
         blolohi[i] = rhslolohi[stage][i];
         bhihilo[i] = rhshihilo[stage][i];
         blohilo[i] = rhslohilo[stage][i];
         bhilolo[i] = rhshilolo[stage][i];
         blololo[i] = rhslololo[stage][i];

         cout << "b[" << i << "] : "
              << bhihihi[i] << "  " << blohihi[i] << endl << "  "
              << bhilohi[i] << "  " << blolohi[i] << endl << "  "
              << bhihilo[i] << "  " << blohilo[i] << endl << "  "
              << bhilolo[i] << "  " << blololo[i] << endl;
      }
      double bstimelapsed_d;
      double elapsedms,invlapsed,mullapsed,sublapsed;
      long long int bsaddcnt = 0;
      long long int bsmulcnt = 0;
      long long int bsdivcnt = 0;

      if(vrblvl > 0)
         cout << "-> GPU multiplies rhs with Q^T ..." << endl;

      GPU_dbl8_bals_qtb
         (ncols,szt,nbt,
          Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
          bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
          bvrb);

      if(vrblvl > 0)
      {
         for(int i=0; i<nrows; i++)
            cout << "Qtb[" << i << "] : "
                 << bhihihi[i] << "  " << blohihi[i] << "  "
                 << bhilohi[i] << "  " << blolohi[i] << endl << "  "
                 << bhihilo[i] << "  " << blohilo[i] << "  "
                 << bhilolo[i] << "  " << blololo[i] << endl;
      }
      if(vrblvl > 0)
      {
         cout << "-> GPU solves an upper triangular system ..." << endl;
 
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
               cout << "R[" << i << "][" << j << "] : "
                    << Rhihihi[i][j] << "  " << Rlohihi[i][j] << "  "
                    << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl << "  "
                    << Rhihilo[i][j] << "  " << Rlohilo[i][j] << "  "
                    << Rhilolo[i][j] << "  " << Rlololo[i][j] << endl;
      }
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            workRhihihi[i][j] = Rhihihi[i][j];
            workRlohihi[i][j] = Rlohihi[i][j];
            workRhilohi[i][j] = Rhilohi[i][j];
            workRlolohi[i][j] = Rlolohi[i][j];
            workRhihilo[i][j] = Rhihilo[i][j];
            workRlohilo[i][j] = Rlohilo[i][j];
            workRhilolo[i][j] = Rhilolo[i][j];
            workRlololo[i][j] = Rlololo[i][j];
         }

      GPU_dbl8_upper_tiled_solver
         (ncols,szt,nbt,
          workRhihihi,workRlohihi,workRhilohi,workRlolohi,
          workRhihilo,workRlohilo,workRhilolo,workRlololo,
          bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
          xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
          &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

     if(vrblvl > 0)
        for(int i=0; i<ncols; i++)
           cout << "x[" << i << "] : "
                << xhihihi[i] << "  " << xlohihi[i] << "  "
                << xhilohi[i] << "  " << xlolohi[i] << endl << "  "
                << xhihilo[i] << "  " << xlohilo[i] << "  "
                << xhilolo[i] << "  " << xlololo[i] << endl;

      for(int j=0; j<ncols; j++)
      {
         solhihihi[stage][j] = xhihihi[j];
         sollohihi[stage][j] = xlohihi[j];
         solhilohi[stage][j] = xhilohi[j];
         sollolohi[stage][j] = xlolohi[j];
         solhihilo[stage][j] = xhihilo[j];
         sollohilo[stage][j] = xlohilo[j];
         solhilolo[stage][j] = xhilolo[j];
         sollololo[stage][j] = xlololo[j];
      }
   }

   for(int i=0; i<nrows; i++)
   {
      free(Ahihihi[i]); free(Alohihi[i]); free(Ahilohi[i]); free(Alolohi[i]);
      free(Ahihilo[i]); free(Alohilo[i]); free(Ahilolo[i]); free(Alololo[i]);
      free(workRhihihi[i]); free(workRlohihi[i]);
      free(workRhilohi[i]); free(workRlolohi[i]);
      free(workRhihilo[i]); free(workRlohilo[i]);
      free(workRhilolo[i]); free(workRlololo[i]);
   }
   free(Ahihihi); free(bhihihi); free(xhihihi); free(workRhihihi);
   free(Alohihi); free(blohihi); free(xlohihi); free(workRlohihi);
   free(Ahilohi); free(bhilohi); free(xhilohi); free(workRhilohi);
   free(Alolohi); free(blolohi); free(xlolohi); free(workRlolohi);
   free(Ahihilo); free(bhihilo); free(xhihilo); free(workRhihilo);
   free(Alohilo); free(blohilo); free(xlohilo); free(workRlohilo);
   free(Ahilolo); free(bhilolo); free(xhilolo); free(workRhilolo);
   free(Alololo); free(blololo); free(xlololo); free(workRlololo);
}
