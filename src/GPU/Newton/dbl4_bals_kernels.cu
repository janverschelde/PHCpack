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
      ddg_dec(&biimhihi,&biimlohi,acchilo,acclolo);
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

   GPU_cmplx4_blocked_houseqr
      (nrows,ncols,szt,nbt,
       Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
       Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
       Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   if(verbose)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(verbose)
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

   if(verbose)
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
   double **solimhilo, double **solimlolo, bool verbose )
{
   if(verbose)
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
      if(verbose)
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

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cmplx4_bals_tail<<<nbt,szt>>>
         (ncols,szt,Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
                    Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
          xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
          ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
          brehihi_d,brelohi_d,brehilo_d,brelolo_d,
          bimhihi_d,bimlohi_d,bimhilo_d,bimlolo_d);
      
      if(verbose)
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

   if(verbose)
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

void GPU_cmplx4_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   bool verbose )
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

   cmplx4_bals_qhb<<<nbt,szt>>>
      (ncols,szt,QHrehihi_d,QHrelohi_d,QHrehilo_d,QHrelolo_d,
                 QHimhihi_d,QHimlohi_d,QHimhilo_d,QHimlolo_d,
       brehihi_d,brelohi_d,brehilo_d,brelolo_d,
       bimhihi_d,bimlohi_d,bimhilo_d,bimlolo_d,
       rrehihi_d,rrelohi_d,rrehilo_d,rrelolo_d,
       rimhihi_d,rimlohi_d,rimhilo_d,rimlolo_d);

   cudaMemcpy(brehihi,rrehihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelohi,rrelohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brehilo,rrehilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelolo,rrelolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhihi,rimhihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlohi,rimlohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhilo,rimhilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlolo,rimlolo_d,szrhs,cudaMemcpyDeviceToHost);

   free(QHrehihi_h); free(QHimhihi_h);
   free(QHrelohi_h); free(QHimlohi_h);
   free(QHrehilo_h); free(QHimhilo_h);
   free(QHrelolo_h); free(QHimlolo_h);
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

void GPU_cmplx4_bals_solve
 ( int dim, int degp1, int szt, int nbt,
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
   double **solimhilo, double **solimlolo, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   const bool bvrb = (vrblvl > 0);

   double **Arehihi = new double*[nrows];
   double **Arelohi = new double*[nrows];
   double **Arehilo = new double*[nrows];
   double **Arelolo = new double*[nrows];
   double **Aimhihi = new double*[nrows];
   double **Aimlohi = new double*[nrows];
   double **Aimhilo = new double*[nrows];
   double **Aimlolo = new double*[nrows];

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
   if(vrblvl)
   {
      cout << "GPU_cmplx4_bals_solve blocks of rhs :" << endl;
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
       xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,bvrb);

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
   for(int stage=1; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      GPU_cmplx4_bals_tail
         (nrows,ncols,szt,nbt,degp1,stage,
          matrehihi,matrelohi,matrehilo,matrelolo,
          matimhihi,matimlohi,matimhilo,matimlolo,
          rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
          rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
          solrehihi,solrelohi,solrehilo,solrelolo,
          solimhihi,solimlohi,solimhilo,solimlolo,bvrb);

      if(vrblvl > 0)
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

      for(int i=0; i<nrows; i++) 
      {
         cout << "assigning component " << i
              << ", stage = " << stage << endl;
         brehihi[i] = rhsrehihi[stage][i];
         brelohi[i] = rhsrelohi[stage][i];
         brehilo[i] = rhsrehilo[stage][i];
         brelolo[i] = rhsrelolo[stage][i];
         bimhihi[i] = rhsimhihi[stage][i];
         bimlohi[i] = rhsimlohi[stage][i];
         bimhilo[i] = rhsimhilo[stage][i];
         bimlolo[i] = rhsimlolo[stage][i];
         cout << "b[" << i << "] : "
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
         cout << "-> GPU multiplies rhs with Q^H ..." << endl;

      GPU_cmplx4_bals_qhb
         (ncols,szt,nbt,
          Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
          brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
          bvrb);

      if(vrblvl > 0)
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
 
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
               cout << "R[" << i << "][" << j << "] : "
                    << Rrehihi[i][j] << "  " << Rrelohi[i][j] << endl << "  "
                    << Rrehilo[i][j] << "  " << Rrelolo[i][j] << endl << "  "
                    << Rimhihi[i][j] << "  " << Rimlohi[i][j] << endl << "  "
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

     if(vrblvl > 0)
        for(int i=0; i<ncols; i++)
           cout << "x[" << i << "] : "
                << xrehihi[i] << "  " << xrelohi[i] << endl << "  "
                << xrehilo[i] << "  " << xrelolo[i] << endl << "  "
                << ximhihi[i] << "  " << ximlohi[i] << endl << "  "
                << ximhilo[i] << "  " << ximlolo[i] << endl;

      for(int j=0; j<ncols; j++)
      {
         solrehihi[stage][j] = xrehihi[j]; solrelohi[stage][j] = xrelohi[j];
         solrehilo[stage][j] = xrehilo[j]; solrelolo[stage][j] = xrelolo[j];
         solimhihi[stage][j] = ximhihi[j]; solimlohi[stage][j] = ximlohi[j];
         solimhilo[stage][j] = ximhilo[j]; solimlolo[stage][j] = ximlolo[j];
      }
   }

   for(int i=0; i<nrows; i++)
   {
      free(Arehihi[i]); free(workRrehihi[i]);
      free(Arelohi[i]); free(workRrelohi[i]);
      free(Arehilo[i]); free(workRrehihi[i]);
      free(Arelolo[i]); free(workRrelohi[i]);
      free(Aimhihi[i]); free(workRimhihi[i]);
      free(Aimlohi[i]); free(workRimlohi[i]);
      free(Aimhilo[i]); free(workRimhilo[i]);
      free(Aimlolo[i]); free(workRimlolo[i]);
   }
   free(Arehihi); free(brehihi); free(xrehihi); free(workRrehihi);
   free(Arelohi); free(brelohi); free(xrelohi); free(workRrelohi);
   free(Arehilo); free(brehilo); free(xrehilo); free(workRrehilo);
   free(Arelolo); free(brelolo); free(xrelolo); free(workRrelolo);
   free(Aimhihi); free(bimhihi); free(ximhihi); free(workRimhihi);
   free(Aimlohi); free(bimlohi); free(ximlohi); free(workRimlohi);
   free(Aimhilo); free(bimhilo); free(ximhilo); free(workRimhilo);
   free(Aimlolo); free(bimlolo); free(ximlolo); free(workRimlolo);
}
