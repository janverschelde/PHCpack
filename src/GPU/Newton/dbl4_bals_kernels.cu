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
#include "dbl4_tail_kernels.h"
#include "dbl4_bals_kernels.h"
#include "write_dbl4_bstimeflops.h"
#include "write_dbl4_qrtimeflops.h"
#include "dbl_onenorms_host.h"

using namespace std;

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

__global__ void cmplx4_bals_qhb
 ( int ncols, int szt,
   double *QHrehihi, double *QHrelohi, double *QHrehilo, double *QHrelolo, 
   double *QHimhihi, double *QHimlohi, double *QHimhilo, double *QHimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *rrehihi, double *rrelohi, double *rrehilo, double *rrelolo,
   double *rimhihi, double *rimlohi, double *rimhilo, double *rimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjrehihi;         // registers for real part of Q^H[idx][j]
   double Qjrelohi;
   double Qjrehilo;
   double Qjrelolo;
   double Qjimhihi;         // registers for imaginary part of Q^H[idx][j]
   double Qjimlohi;
   double Qjimhilo;
   double Qjimlolo;
   double bjrehihi;         // register for brehihi[j]
   double bjrelohi;         // register for brelohi[j]
   double bjrehilo;         // register for brehilo[j]
   double bjrelolo;         // register for brelolo[j]
   double bjimhihi;         // register for bimhihi[j]
   double bjimlohi;         // register for bimlohi[j]
   double bjimhilo;         // register for bimhilo[j]
   double bjimlolo;         // register for bimlolo[j]
   double rirehihi = 0.0;   // register for result, rrehihi[idx]
   double rirelohi = 0.0;   // register for result, rrelohi[idx]
   double rirehilo = 0.0;   // register for result, rrehilo[idx]
   double rirelolo = 0.0;   // register for result, rrelolo[idx]
   double riimhihi = 0.0;   // register for result, rimhihi[idx]
   double riimlohi = 0.0;   // register for result, rimlohi[idx]
   double riimhilo = 0.0;   // register for result, rimhilo[idx]
   double riimlolo = 0.0;   // register for result, rimlolo[idx]
   double acchihi,acclohi,acchilo,acclolo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjrehihi = QHrehihi[offset+j];
      Qjrelohi = QHrelohi[offset+j];
      Qjrehilo = QHrehilo[offset+j];
      Qjrelolo = QHrelolo[offset+j];
      Qjimhihi = QHimhihi[offset+j];
      Qjimlohi = QHimlohi[offset+j];
      Qjimhilo = QHimhilo[offset+j];
      Qjimlolo = QHimlolo[offset+j];
      bjrehihi = brehihi[j];
      bjrelohi = brelohi[j];
      bjrehilo = brehilo[j];
      bjrelolo = brelolo[j];
      bjimhihi = bimhihi[j];
      bjimlohi = bimlohi[j];
      bjimhilo = bimhilo[j];
      bjimlolo = bimlolo[j];
      // ri = ri + Qj*bj;
      // zre = Qjre*bjre - Qjim*bjim;
      // rire = rire + zre;
      qdg_mul(Qjrehihi,Qjrelohi,Qjrehilo,Qjrelolo,
              bjrehihi,bjrelohi,bjrehilo,bjrelolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&rirehihi,&rirelohi,&rirehilo,&rirelolo,
              acchihi,acclohi,acchilo,acclolo);
      qdg_mul(Qjimhihi,Qjimlohi,Qjimhilo,Qjimlolo,
              bjimhihi,bjimlohi,bjimhilo,bjimlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&rirehihi,&rirelohi,&rirehilo,&rirelolo,
              acchihi,acclohi,acchilo,acclolo);
      // zim = Qjre*bjim + Qjim*bjre;
      // riim = riim + zim;
      qdg_mul(Qjrehihi,Qjrelohi,Qjrehilo,Qjrelolo,
              bjimhihi,bjimlohi,bjimhilo,bjimlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&riimhihi,&riimlohi,&riimhilo,&riimlolo,
              acchihi,acclohi,acchilo,acclolo);
      qdg_mul(Qjimhihi,Qjimlohi,Qjimhilo,Qjimlolo,
              bjrehihi,bjrelohi,bjrehilo,bjrelolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&riimhihi,&riimlohi,&riimhilo,&riimlolo,
              acchihi,acclohi,acchilo,acclolo);
   }
   rrehihi[idx] = rirehihi; rrelohi[idx] = rirelohi;
   rrehilo[idx] = rirehilo; rrelolo[idx] = rirelolo;
   rimhihi[idx] = riimhihi; rimlohi[idx] = riimlohi;
   rimhilo[idx] = riimhilo; rimlolo[idx] = riimlolo;
}

void GPU_dbl4_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo, 
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;
   bool verbose = (vrblvl > 1);

   if(vrblvl > 0) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_dbl4_blocked_houseqr
      (nrows,ncols,szt,nbt,Ahihi,Alohi,Ahilo,Alolo,
       Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl4_qrtimeflops
         (0,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^T ..." << endl;

   GPU_dbl4_bals_qtb
      (ncols,szt,nbt,Qhihi,Qlohi,Qhilo,Qlolo,
                     bhihi,blohi,bhilo,blolo,totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "Qtb[" << i << "] : "
              << bhihi[i] << "  " << blohi[i] << endl << "  "
              << bhilo[i] << "  " << blolo[i] << endl;
   }
   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if(vrblvl > 0)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(vrblvl > 1)
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

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl4_bstimeflops
         (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,addover,bsmulcnt,mulover,bsdivcnt,divover);

   if(vrblvl > 1)
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

void GPU_cmplx4_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo, 
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;
   bool verbose = (vrblvl > 1);

   if(vrblvl > 0) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_cmplx4_blocked_houseqr
      (nrows,ncols,szt,nbt,
       Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
       Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
       Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl4_qrtimeflops
         (1,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0)
      cout << "-> GPU multiplies rhs with Q^H ..." << endl;

   GPU_cmplx4_bals_qhb
      (ncols,szt,nbt,
       Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
       brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
       totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "QHb[" << i << "] : "
              << brehihi[i] << "  " << brelohi[i] << endl << "  "
              << brehilo[i] << "  " << brelolo[i] << endl << "  "
              << bimhihi[i] << "  " << bimlohi[i] << endl << "  "
              << bimhilo[i] << "  " << bimlolo[i] << endl;
   }
   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   if(vrblvl > 0)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rrehihi[i][j] << "  " << Rrelohi[i][j] << endl << "  "
                 << Rrehilo[i][j] << "  " << Rrelolo[i][j] << endl << "  "
                 << Rimhihi[i][j] << "  " << Rimlohi[i][j] << endl << "  "
                 << Rimhilo[i][j] << "  " << Rimlolo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << brehihi[i] << "  " << brelohi[i] << endl << "  "
              << brehilo[i] << "  " << brelolo[i] << endl << "  "
              << bimhihi[i] << "  " << bimlohi[i] << endl << "  "
              << bimhilo[i] << "  " << bimlolo[i] << endl;
   }
   double **workRrehihi = new double*[nrows]; // work around ...
   double **workRrelohi = new double*[nrows];
   double **workRrehilo = new double*[nrows]; 
   double **workRrelolo = new double*[nrows];
   double **workRimhihi = new double*[nrows];
   double **workRimlohi = new double*[nrows];
   double **workRimhilo = new double*[nrows];
   double **workRimlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRrehihi[i] = new double[ncols]; workRrelohi[i] = new double[ncols];
      workRrehilo[i] = new double[ncols]; workRrelolo[i] = new double[ncols];
      workRimhihi[i] = new double[ncols]; workRimlohi[i] = new double[ncols];
      workRimhilo[i] = new double[ncols]; workRimlolo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRrehihi[i][j] = Rrehihi[i][j]; workRrelohi[i][j] = Rrelohi[i][j];
         workRrehilo[i][j] = Rrehilo[i][j]; workRrelolo[i][j] = Rrelolo[i][j];
         workRimhihi[i][j] = Rimhihi[i][j]; workRimlohi[i][j] = Rimlohi[i][j];
         workRimhilo[i][j] = Rimhilo[i][j]; workRimlolo[i][j] = Rimlolo[i][j];
      }
   }
   GPU_cmplx4_upper_tiled_solver
      (ncols,szt,nbt,
       workRrehihi,workRrelohi,workRrehilo,workRrelolo,
       workRimhihi,workRimlohi,workRimhilo,workRimlolo,
       brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
       xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&addover,&bsmulcnt,&mulover,&bsdivcnt,&divover);

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl4_bstimeflops
         (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,addover,bsmulcnt,mulover,bsdivcnt,divover);

   if(vrblvl > 1)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rrehihi[i][j] << "  " << Rrelohi[i][j] << endl << "  "
                 << Rrehilo[i][j] << "  " << Rrelolo[i][j] << endl << "  "
                 << Rimhihi[i][j] << "  " << Rimlohi[i][j] << endl << "  "
                 << Rimhilo[i][j] << "  " << Rimlolo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << brehihi[i] << "  " << brelohi[i] << endl << "  "
              << brehilo[i] << "  " << brelolo[i] << endl << "  "
              << bimhihi[i] << "  " << bimlohi[i] << endl << "  "
              << bimhilo[i] << "  " << bimlolo[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRrehihi[i]); free(workRimhihi[i]);
      free(workRrehilo[i]); free(workRimhilo[i]);
      free(workRrelohi[i]); free(workRimlohi[i]);
      free(workRrelolo[i]); free(workRimlolo[i]);
   }
   free(workRrehihi); free(workRimhihi);
   free(workRrehilo); free(workRimhilo);
   free(workRrelohi); free(workRimlohi);
   free(workRrelolo); free(workRimlolo);
}

void write_dbl4_qtbflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for Q^T*b : " << lapsms << " milliseconds." << endl;

   long long int flopcnt;
   const long long int longncols2 = ncols*ncols; // to avoid overflow
   if(ctype == 0)
      flopcnt = 89*longncols2 + 336*longncols2;
      // as many + as * in one inner product
   else
      flopcnt = 4*89*longncols2 + 4*336*longncols2;
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

void GPU_dbl4_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *totqtblapsedms, int vrblvl )
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

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   dbl4_bals_qtb<<<nbt,szt>>>
      (ncols,szt,Qthihi_d,Qtlohi_d,Qthilo_d,Qtlolo_d,
       bhihi_d,blohi_d,bhilo_d,blolo_d,rhihi_d,rlohi_d,rhilo_d,rlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(bhihi,rhihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blohi,rlohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhilo,rhilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blolo,rlolo_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl4_qtbflops(0,ncols,milliseconds);

   free(Qthihi_h); free(Qtlohi_h); free(Qthilo_h); free(Qtlolo_h);

   cudaFree(bhihi_d); cudaFree(blohi_d); cudaFree(bhilo_d); cudaFree(blolo_d);
   cudaFree(rhihi_d); cudaFree(rlohi_d); cudaFree(rhilo_d); cudaFree(rlolo_d);
   cudaFree(Qthihi_d); cudaFree(Qtlohi_d);
   cudaFree(Qthilo_d); cudaFree(Qtlolo_d);
}

void GPU_cmplx4_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *totqtblapsedms, int vrblvl )
{
   double *brehihi_d;
   double *brelohi_d;
   double *brehilo_d;
   double *brelolo_d;
   double *bimhihi_d;
   double *bimlohi_d;
   double *bimhilo_d;
   double *bimlolo_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&brehihi_d,szrhs);
   cudaMalloc((void**)&brelohi_d,szrhs);
   cudaMalloc((void**)&brehilo_d,szrhs);
   cudaMalloc((void**)&brelolo_d,szrhs);
   cudaMalloc((void**)&bimhihi_d,szrhs);
   cudaMalloc((void**)&bimlohi_d,szrhs);
   cudaMalloc((void**)&bimhilo_d,szrhs);
   cudaMalloc((void**)&bimlolo_d,szrhs);
   cudaMemcpy(brehihi_d,brehihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelohi_d,brelohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brehilo_d,brehilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelolo_d,brelolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhihi_d,bimhihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlohi_d,bimlohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhilo_d,bimhilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlolo_d,bimlolo,szrhs,cudaMemcpyHostToDevice);

   double *rrehihi_d;
   double *rrelohi_d;
   double *rrehilo_d;
   double *rrelolo_d;
   double *rimhihi_d;
   double *rimlohi_d;
   double *rimhilo_d;
   double *rimlolo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rrehihi_d,szsol);
   cudaMalloc((void**)&rrelohi_d,szsol);
   cudaMalloc((void**)&rrehilo_d,szsol);
   cudaMalloc((void**)&rrelolo_d,szsol);
   cudaMalloc((void**)&rimhihi_d,szsol);
   cudaMalloc((void**)&rimlohi_d,szsol);
   cudaMalloc((void**)&rimhilo_d,szsol);
   cudaMalloc((void**)&rimlolo_d,szsol);

   double *QHrehihi_d;
   double *QHrelohi_d;
   double *QHrehilo_d;
   double *QHrelolo_d;
   double *QHimhihi_d;
   double *QHimlohi_d;
   double *QHimhilo_d;
   double *QHimlolo_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&QHrehihi_d,szmat);
   cudaMalloc((void**)&QHrelohi_d,szmat);
   cudaMalloc((void**)&QHrehilo_d,szmat);
   cudaMalloc((void**)&QHrelolo_d,szmat);
   cudaMalloc((void**)&QHimhihi_d,szmat);
   cudaMalloc((void**)&QHimlohi_d,szmat);
   cudaMalloc((void**)&QHimhilo_d,szmat);
   cudaMalloc((void**)&QHimlolo_d,szmat);

   double *QHrehihi_h = new double[szmat];
   double *QHrelohi_h = new double[szmat];
   double *QHrehilo_h = new double[szmat];
   double *QHrelolo_h = new double[szmat];
   double *QHimhihi_h = new double[szmat];
   double *QHimlohi_h = new double[szmat];
   double *QHimhilo_h = new double[szmat];
   double *QHimlolo_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         QHrehihi_h[idx]   = Qrehihi[j][i];
         QHrelohi_h[idx]   = Qrelohi[j][i];
         QHrehilo_h[idx]   = Qrehilo[j][i];
         QHrelolo_h[idx]   = Qrelolo[j][i];
         QHimhihi_h[idx]   = - Qimhihi[j][i]; // Hermitian transpose !
         QHimlohi_h[idx]   = - Qimlohi[j][i]; // Hermitian transpose !
         QHimhilo_h[idx]   = - Qimhilo[j][i]; // Hermitian transpose !
         QHimlolo_h[idx++] = - Qimlolo[j][i]; // Hermitian transpose !
      }

   cudaMemcpy(QHrehihi_d,QHrehihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelohi_d,QHrelohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrehilo_d,QHrehilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelolo_d,QHrelolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhihi_d,QHimhihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlohi_d,QHimlohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhilo_d,QHimhilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlolo_d,QHimlolo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   cmplx4_bals_qhb<<<nbt,szt>>>
      (ncols,szt,QHrehihi_d,QHrelohi_d,QHrehilo_d,QHrelolo_d,
                 QHimhihi_d,QHimlohi_d,QHimhilo_d,QHimlolo_d,
       brehihi_d,brelohi_d,brehilo_d,brelolo_d,
       bimhihi_d,bimlohi_d,bimhilo_d,bimlolo_d,
       rrehihi_d,rrelohi_d,rrehilo_d,rrelolo_d,
       rimhihi_d,rimlohi_d,rimhilo_d,rimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(brehihi,rrehihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelohi,rrelohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brehilo,rrehilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelolo,rrelolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhihi,rimhihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlohi,rimlohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhilo,rimhilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlolo,rimlolo_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl4_qtbflops(1,ncols,milliseconds);

   cudaFree(brehihi_d); cudaFree(brelohi_d);
   cudaFree(brehilo_d); cudaFree(brelolo_d);
   cudaFree(bimhihi_d); cudaFree(bimlohi_d);
   cudaFree(bimhilo_d); cudaFree(bimlolo_d);
   cudaFree(rrehihi_d); cudaFree(rrelohi_d);
   cudaFree(rrehilo_d); cudaFree(rrelolo_d);
   cudaFree(rimhihi_d); cudaFree(rimlohi_d);
   cudaFree(rimhilo_d); cudaFree(rimlolo_d);
   cudaFree(QHrehihi_d); cudaFree(QHrelohi_d);
   cudaFree(QHrehilo_d); cudaFree(QHrelolo_d);
   cudaFree(QHimhihi_d); cudaFree(QHimlohi_d);
   cudaFree(QHimhilo_d); cudaFree(QHimlolo_d);

   free(QHrehihi_h); free(QHimhihi_h);
   free(QHrelohi_h); free(QHimlohi_h);
   free(QHrehilo_h); free(QHimhilo_h);
   free(QHrelolo_h); free(QHimlolo_h);
}

void GPU_dbl4_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped substitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

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
   if(vrblvl > 1)
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
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "-> skipping GPU_dbl4_bals_head ..." << endl;
   }
   else
   {
      CPU_dbl_onenorm(nrows,rhshihi[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm < 1.0e-56)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "-> no skipping GPU_dbl4_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "-> skip call to GPU_dbl4_bals_head ..." << endl;

            *noqr = true;

            for(int j=0; j<ncols; j++)
            {
               solhihi[0][j] = 0.0; sollohi[0][j] = 0.0;
               solhilo[0][j] = 0.0; sollolo[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "-> calling GPU_dbl4_bals_head ..." << endl;

         double **Ahihi = new double*[nrows];
         double **Alohi = new double*[nrows];
         double **Ahilo = new double*[nrows];
         double **Alolo = new double*[nrows];

         for(int i=0; i<nrows; i++)
         {
            Ahihi[i] = new double[ncols];
            Alohi[i] = new double[ncols];
            Ahilo[i] = new double[ncols];
            Alolo[i] = new double[ncols];

            for(int j=0; j<ncols; j++)
            {
               Ahihi[i][j] = mathihi[0][i][j]; Alohi[i][j] = matlohi[0][i][j];
               Ahilo[i][j] = mathilo[0][i][j]; Alolo[i][j] = matlolo[0][i][j];
            }
            bhihi[i] = rhshihi[0][i]; blohi[i] = rhslohi[0][i];
            bhilo[i] = rhshilo[0][i]; blolo[i] = rhslolo[0][i];

            for(int j=0; j<ncols; j++)
            {
               Rhihi[i][j] = mathihi[0][i][j]; Rlohi[i][j] = matlohi[0][i][j];
               Rhilo[i][j] = mathilo[0][i][j]; Rlolo[i][j] = matlolo[0][i][j];
            }
         }
         GPU_dbl4_bals_head
            (nrows,ncols,szt,nbt,Ahihi,Alohi,Ahilo,Alolo,
             Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
             bhihi,blohi,bhilo,blolo,xhihi,xlohi,xhilo,xlolo,
             totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;
   
         for(int j=0; j<ncols; j++)
         {
            solhihi[0][j] = xhihi[j];
            sollohi[0][j] = xlohi[j];
            solhilo[0][j] = xhilo[j];
            sollolo[0][j] = xlolo[j];
         }

         for(int i=0; i<nrows; i++)
         {
            free(Ahihi[i]); free(Alohi[i]); free(Ahilo[i]); free(Alolo[i]);
         }
         free(Ahihi); free(Alohi); free(Ahilo); free(Alolo);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xshihi = solhihi[stage-1];   // solution to do the update with
      CPU_dbl_onenorm(dim,xshihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-56)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << stage-1 << "] ..." << endl;

         GPU_dbl4_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,
             mathihi,matlohi,mathilo,matlolo,rhshihi,rhslohi,rhshilo,rhslolo,
             solhihi,sollohi,solhilo,sollolo,totupdlapsedms,vrblvl);

         if(vrblvl > 1)
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
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         bhihi[i] = rhshihi[stage][i];
         blohi[i] = rhslohi[stage][i];
         bhilo[i] = rhshilo[stage][i];
         blolo[i] = rhslolo[stage][i];
         // cout << "b[" << i << "] : "
         //      << bhihi[i] << "  " << blohi[i] << endl << "  "
         //      << bhilo[i] << "  " << blolo[i] << endl;
      }
      CPU_dbl_onenorm(nrows,bhihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-56) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xhihi[i] = 0.0; xlohi[i] = 0.0;
            xhilo[i] = 0.0; xlolo[i] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "-> run backsubstitution for x[" << stage << "] ..."
                 << endl;

         if(firstbs)
         {
            *newtail = stage;
            firstbs = false;
         }

         double bstimelapsed_d;
         double elapsedms,invlapsed,mullapsed,sublapsed;
         long long int bsaddcnt = 0;
         long long int bsmulcnt = 0;
         long long int bsdivcnt = 0;

         if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^T ..." << endl;

         GPU_dbl4_bals_qtb
            (ncols,szt,nbt,Qhihi,Qlohi,Qhilo,Qlolo,
                           bhihi,blohi,bhilo,blolo,totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "Qtb[" << i << "] : "
                    << bhihi[i] << "  " << blohi[i] << "  "
                    << bhilo[i] << "  " << blolo[i] << endl;
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
 
            if(vrblvl > 1)
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

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl4_bstimeflops
               (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,addover,bsmulcnt,mulover,
                bsdivcnt,divover);

         if(vrblvl > 1)
            for(int i=0; i<ncols; i++)
               cout << "x[" << i << "] : "
                    << xhihi[i] << "  " << xlohi[i] << "  "
                    << xhilo[i] << "  " << xlolo[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solhihi[stage][j] = xhihi[j];
         sollohi[stage][j] = xlohi[j];
         solhilo[stage][j] = xhilo[j];
         sollolo[stage][j] = xlolo[j];
      }
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;

   for(int i=0; i<nrows; i++)
   {
      free(workRhihi[i]); free(workRlohi[i]);
      free(workRhilo[i]); free(workRlolo[i]);
   }
   free(bhihi); free(xhihi); free(workRhihi);
   free(blohi); free(xlohi); free(workRlohi);
   free(bhilo); free(xhilo); free(workRhilo);
   free(blolo); free(xlolo); free(workRlolo);
}

void GPU_cmplx4_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo, 
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi, 
   double **solimhilo, double **solimlolo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   double *brehihi = new double[nrows];
   double *brelohi = new double[nrows];
   double *brehilo = new double[nrows];
   double *brelolo = new double[nrows];
   double *bimhihi = new double[nrows];
   double *bimlohi = new double[nrows];
   double *bimhilo = new double[nrows];
   double *bimlolo = new double[nrows];
   double *xrehihi = new double[ncols];
   double *xrelohi = new double[ncols];
   double *xrehilo = new double[ncols];
   double *xrelolo = new double[ncols];
   double *ximhihi = new double[ncols];
   double *ximlohi = new double[ncols];
   double *ximhilo = new double[ncols];
   double *ximlolo = new double[ncols];

   double **workRrehihi = new double*[nrows]; // GPU upper solver changes R
   double **workRrelohi = new double*[nrows];
   double **workRrehilo = new double*[nrows];
   double **workRrelolo = new double*[nrows];
   double **workRimhihi = new double*[nrows]; 
   double **workRimlohi = new double*[nrows]; 
   double **workRimhilo = new double*[nrows]; 
   double **workRimlolo = new double*[nrows]; 

   for(int i=0; i<nrows; i++)
   {
      workRrehihi[i] = new double[ncols];
      workRrelohi[i] = new double[ncols];
      workRrehilo[i] = new double[ncols];
      workRrelolo[i] = new double[ncols];
      workRimhihi[i] = new double[ncols];
      workRimlohi[i] = new double[ncols];
      workRimhilo[i] = new double[ncols];
      workRimlolo[i] = new double[ncols];
   }
   if(vrblvl > 1)
   {
      cout << "GPU_cmplx4_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehihi[k][i] << "  " << rhsrelohi[k][i] << endl
                 << "  "
                 << rhsrehilo[k][i] << "  " << rhsrelolo[k][i] << endl
                 << "  "
                 << rhsimhihi[k][i] << "  " << rhsimlohi[k][i] << endl
                 << "  "
                 << rhsimhilo[k][i] << "  " << rhsimlolo[k][i] << endl;
         }
   }
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "-> skipping GPU_cmplx4_bals_head ..." << endl;
   }
   else
   {
      CPU_cmplx_onenorm(nrows,rhsrehihi[0],rhsimhihi[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm < 1.0e-56)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "-> no skipping GPU_cmplx4_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "-> skip call to GPU_cmplx4_bals_head ..." << endl;

            *noqr = true;

            for(int j=0; j<ncols; j++)
            {
               solrehihi[0][j] = 0.0; solrelohi[0][j] = 0.0;
               solrehilo[0][j] = 0.0; solrelolo[0][j] = 0.0;
               solimhihi[0][j] = 0.0; solimlohi[0][j] = 0.0;
               solimhilo[0][j] = 0.0; solimlolo[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "-> calling GPU_cmplx4_bals_head ..." << endl;

         double **Arehihi = new double*[nrows];
         double **Arelohi = new double*[nrows];
         double **Arehilo = new double*[nrows];
         double **Arelolo = new double*[nrows];
         double **Aimhihi = new double*[nrows];
         double **Aimlohi = new double*[nrows];
         double **Aimhilo = new double*[nrows];
         double **Aimlolo = new double*[nrows];

         for(int i=0; i<nrows; i++)
         {
            Arehihi[i] = new double[ncols]; Arelohi[i] = new double[ncols];
            Arehilo[i] = new double[ncols]; Arelolo[i] = new double[ncols];
            Aimhihi[i] = new double[ncols]; Aimlohi[i] = new double[ncols];
            Aimhilo[i] = new double[ncols]; Aimlolo[i] = new double[ncols];

            for(int j=0; j<ncols; j++)
            {
               Arehihi[i][j] = matrehihi[0][i][j];
               Arelohi[i][j] = matrelohi[0][i][j];
               Arehilo[i][j] = matrehilo[0][i][j];
               Arelolo[i][j] = matrelolo[0][i][j];
               Aimhihi[i][j] = matimhihi[0][i][j];
               Aimlohi[i][j] = matimlohi[0][i][j];
               Aimhilo[i][j] = matimhilo[0][i][j];
               Aimlolo[i][j] = matimlolo[0][i][j];
            }
            brehihi[i] = rhsrehihi[0][i]; brelohi[i] = rhsrelohi[0][i];
            brehilo[i] = rhsrehilo[0][i]; brelolo[i] = rhsrelolo[0][i];
            bimhihi[i] = rhsimhihi[0][i]; bimlohi[i] = rhsimlohi[0][i];
            bimhilo[i] = rhsimhilo[0][i]; bimlolo[i] = rhsimlolo[0][i];

            for(int j=0; j<ncols; j++)
            {
               Rrehihi[i][j] = matrehihi[0][i][j];
               Rrelohi[i][j] = matrelohi[0][i][j];
               Rrehilo[i][j] = matrehilo[0][i][j];
               Rrelolo[i][j] = matrelolo[0][i][j];
               Rimhihi[i][j] = matimhihi[0][i][j];
               Rimlohi[i][j] = matimlohi[0][i][j];
               Rimhilo[i][j] = matimhilo[0][i][j];
               Rimlolo[i][j] = matimlolo[0][i][j];
            }
         }
         GPU_cmplx4_bals_head
            (nrows,ncols,szt,nbt,
             Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
             Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
             Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
             brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
             xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
             totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         for(int j=0; j<ncols; j++)
         {
            solrehihi[0][j] = xrehihi[j];
            solrelohi[0][j] = xrelohi[j];
            solrehilo[0][j] = xrehilo[j];
            solrelolo[0][j] = xrelolo[j];
            solimhihi[0][j] = ximhihi[j];
            solimlohi[0][j] = ximlohi[j];
            solimhilo[0][j] = ximhilo[j];
            solimlolo[0][j] = ximlolo[j];
         }
         for(int i=0; i<nrows; i++)
         {
            free(Arehihi[i]); free(Arelohi[i]);
            free(Arehilo[i]); free(Arelolo[i]);
            free(Aimhihi[i]); free(Aimlohi[i]);
            free(Aimhilo[i]); free(Aimlolo[i]);
         }
         free(Arehihi); free(Arelohi); free(Arehilo); free(Arelolo);
         free(Aimhihi); free(Aimlohi); free(Aimhilo); free(Aimlolo);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xrs = solrehihi[stage-1];  // solution to do the update with
      double *xis = solimhihi[stage-1];

      CPU_cmplx_onenorm(dim,xrs,xis,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-56)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << stage-1 << "] ..." << endl;

         GPU_cmplx4_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,
             matrehihi,matrelohi,matrehilo,matrelolo,
             matimhihi,matimlohi,matimhilo,matimlolo,
             rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
             rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
             solrehihi,solrelohi,solrehilo,solrelolo,
             solimhihi,solimlohi,solimhilo,solimlolo,totupdlapsedms,vrblvl);

         if(vrblvl > 1)
         {
            cout << "blocks of rhs before assignment :" << endl;
            for(int k=0; k<degp1; k++)
            {
               for(int i=0; i<nrows; i++)
                  cout << "rhs[" << k << "][" << i << "] : "
                       << rhsrehihi[k][i] << "  "
                       << rhsrelohi[k][i] << endl << "  " 
                       << rhsrehilo[k][i] << "  "
                       << rhsrelolo[k][i] << endl << "  " 
                       << rhsimhihi[k][i] << "  "
                       << rhsimlohi[k][i] << endl << "  "
                       << rhsimhilo[k][i] << "  "
                       << rhsimlolo[k][i] << endl;
            }
         }
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         brehihi[i] = rhsrehihi[stage][i];
         brelohi[i] = rhsrelohi[stage][i];
         brehilo[i] = rhsrehilo[stage][i];
         brelolo[i] = rhsrelolo[stage][i];
         bimhihi[i] = rhsimhihi[stage][i];
         bimlohi[i] = rhsimlohi[stage][i];
         bimhilo[i] = rhsimhilo[stage][i];
         bimlolo[i] = rhsimlolo[stage][i];
         // cout << "b[" << i << "] : "
         //      << brehihi[i] << "  " << brelohi[i] << endl << "  "
         //      << brehilo[i] << "  " << brelolo[i] << endl << "  "
         //      << bimhihi[i] << "  " << bimlohi[i] << endl << "  "
         //      << bimhilo[i] << "  " << bimlolo[i] << endl;
      }
      CPU_cmplx_onenorm(dim,brehihi,bimhihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-56) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitutions for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xrehihi[i] = 0.0; xrelohi[i] = 0.0;
            xrehilo[i] = 0.0; xrelolo[i] = 0.0;
            ximhihi[i] = 0.0; ximlohi[i] = 0.0;
            ximhilo[i] = 0.0; ximlolo[i] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "-> run backsubstitutions for x[" << stage << "] ..."
                 << endl;

         if(firstbs)
         {
            *newtail = stage;
            firstbs = false;
         }
         double bstimelapsed_d;
         double elapsedms,invlapsed,mullapsed,sublapsed;
         long long int bsaddcnt = 0;
         long long int bsmulcnt = 0;
         long long int bsdivcnt = 0;
         double addover = 0.0;
         double mulover = 0.0;
         double divover = 0.0;

         if(vrblvl > 0)
            cout << "-> GPU multiplies rhs with Q^H ..." << endl;

         GPU_cmplx4_bals_qhb
            (ncols,szt,nbt,
             Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
             brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
             totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "QHb[" << i << "] : "
                    << brehihi[i] << "  " << brelohi[i] << endl << "  "
                    << brehilo[i] << "  " << brelolo[i] << endl << "  "
                    << bimhihi[i] << "  " << bimlohi[i] << endl << "  "
                    << bimhilo[i] << "  " << bimlolo[i] << endl;
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
 
            if(vrblvl > 1)
               for(int i=0; i<nrows; i++)
                  for(int j=0; j<ncols; j++)
                     cout << "R[" << i << "][" << j << "] : "
                          << Rrehihi[i][j] << "  " << Rrelohi[i][j] << endl
                          << "  "
                          << Rrehilo[i][j] << "  " << Rrelolo[i][j] << endl
                          << "  "
                          << Rimhihi[i][j] << "  " << Rimlohi[i][j] << endl
                          << "  "
                          << Rimhilo[i][j] << "  " << Rimlolo[i][j] << endl;
         }
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
            {
               workRrehihi[i][j] = Rrehihi[i][j];
               workRrelohi[i][j] = Rrelohi[i][j];
               workRrehilo[i][j] = Rrehilo[i][j];
               workRrelolo[i][j] = Rrelolo[i][j];
               workRimhihi[i][j] = Rimhihi[i][j];
               workRimlohi[i][j] = Rimlohi[i][j];
               workRimhilo[i][j] = Rimhilo[i][j];
               workRimlolo[i][j] = Rimlolo[i][j];
            }

         GPU_cmplx4_upper_tiled_solver
            (ncols,szt,nbt,
             workRrehihi,workRrelohi,workRrehilo,workRrelolo,
             workRimhihi,workRimlohi,workRimhilo,workRimlolo,
             brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
             xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
             &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
             &bsaddcnt,&addover,&bsmulcnt,&mulover,&bsdivcnt,&divover);

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
             write_dbl4_bstimeflops
                (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,
                 bstimelapsed_d,bsaddcnt,addover,bsmulcnt,mulover,
                 bsdivcnt,divover);

         if(vrblvl > 1)
            for(int i=0; i<ncols; i++)
               cout << "x[" << i << "] : "
                    << xrehihi[i] << "  " << xrelohi[i] << endl << "  "
                    << xrehilo[i] << "  " << xrelolo[i] << endl << "  "
                    << ximhihi[i] << "  " << ximlohi[i] << endl << "  "
                    << ximhilo[i] << "  " << ximlolo[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solrehihi[stage][j] = xrehihi[j]; solrelohi[stage][j] = xrelohi[j];
         solrehilo[stage][j] = xrehilo[j]; solrelolo[stage][j] = xrelolo[j];
         solimhihi[stage][j] = ximhihi[j]; solimlohi[stage][j] = ximlohi[j];
         solimhilo[stage][j] = ximhilo[j]; solimlolo[stage][j] = ximlolo[j];
      }
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;

   for(int i=0; i<nrows; i++)
   {
      free(workRrehihi[i]); free(workRrelohi[i]);
      free(workRrehilo[i]); free(workRrelolo[i]);
      free(workRimhihi[i]); free(workRimlohi[i]);
      free(workRimhilo[i]); free(workRimlolo[i]);
   }
   free(brehihi); free(xrehihi); free(workRrehihi);
   free(brelohi); free(xrelohi); free(workRrelohi);
   free(brehilo); free(xrehilo); free(workRrehilo);
   free(brelolo); free(xrelolo); free(workRrelolo);
   free(bimhihi); free(ximhihi); free(workRimhihi);
   free(bimlohi); free(ximlohi); free(workRimlohi);
   free(bimhilo); free(ximhilo); free(workRimhilo);
   free(bimlolo); free(ximlolo); free(workRimlolo);
}
