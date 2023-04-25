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
#include "dbl8_tail_kernels.h"
#include "dbl8_bals_kernels.h"
#include "write_dbl8_bstimeflops.h"
#include "write_dbl8_qrtimeflops.h"
#include "dbl_onenorms_host.h"

using namespace std;

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

__global__ void cmplx8_bals_qhb
 ( int ncols, int szt,
   double *QHrehihihi, double *QHrelohihi,
   double *QHrehilohi, double *QHrelolohi, 
   double *QHrehihilo, double *QHrelohilo,
   double *QHrehilolo, double *QHrelololo, 
   double *QHimhihihi, double *QHimlohihi,
   double *QHimhilohi, double *QHimlolohi,
   double *QHimhihilo, double *QHimlohilo,
   double *QHimhilolo, double *QHimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *rrehihihi, double *rrelohihi, double *rrehilohi, double *rrelolohi,
   double *rrehihilo, double *rrelohilo, double *rrehilolo, double *rrelololo,
   double *rimhihihi, double *rimlohihi, double *rimhilohi, double *rimlolohi,
   double *rimhihilo, double *rimlohilo, double *rimhilolo, double *rimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjrehihihi;         // registers for real part of Q^H[idx][j]
   double Qjrelohihi;
   double Qjrehilohi;
   double Qjrelolohi;
   double Qjrehihilo;
   double Qjrelohilo;
   double Qjrehilolo;
   double Qjrelololo;
   double Qjimhihihi;         // registers for imaginary part of Q^H[idx][j]
   double Qjimlohihi;
   double Qjimhilohi;
   double Qjimlolohi;
   double Qjimhihilo; 
   double Qjimlohilo;
   double Qjimhilolo;
   double Qjimlololo;
   double bjrehihihi;         // register for brehihihi[j]
   double bjrelohihi;         // register for brelohihi[j]
   double bjrehilohi;         // register for brehilohi[j]
   double bjrelolohi;         // register for brelolohi[j]
   double bjrehihilo;         // register for brehihilo[j]
   double bjrelohilo;         // register for brelohilo[j]
   double bjrehilolo;         // register for brehilolo[j]
   double bjrelololo;         // register for brelololo[j]
   double bjimhihihi;         // register for bimhihihi[j]
   double bjimlohihi;         // register for bimlohihi[j]
   double bjimhilohi;         // register for bimhilohi[j]
   double bjimlolohi;         // register for bimlolohi[j]
   double bjimhihilo;         // register for bimhihilo[j]
   double bjimlohilo;         // register for bimlohilo[j]
   double bjimhilolo;         // register for bimhilolo[j]
   double bjimlololo;         // register for bimlololo[j]
   double rirehihihi = 0.0;   // register for result, rrehihihi[idx]
   double rirelohihi = 0.0;   // register for result, rrelohihi[idx]
   double rirehilohi = 0.0;   // register for result, rrehilohi[idx]
   double rirelolohi = 0.0;   // register for result, rrelolohi[idx]
   double rirehihilo = 0.0;   // register for result, rrehihilo[idx]
   double rirelohilo = 0.0;   // register for result, rrelohilo[idx]
   double rirehilolo = 0.0;   // register for result, rrehilolo[idx]
   double rirelololo = 0.0;   // register for result, rrelololo[idx]
   double riimhihihi = 0.0;   // register for result, rimhihihi[idx]
   double riimlohihi = 0.0;   // register for result, rimlohihi[idx]
   double riimhilohi = 0.0;   // register for result, rimhilohi[idx]
   double riimlolohi = 0.0;   // register for result, rimlolohi[idx]
   double riimhihilo = 0.0;   // register for result, rimhihilo[idx]
   double riimlohilo = 0.0;   // register for result, rimlohilo[idx]
   double riimhilolo = 0.0;   // register for result, rimhilolo[idx]
   double riimlololo = 0.0;   // register for result, rimlololo[idx]
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjrehihihi = QHrehihihi[offset+j];
      Qjrelohihi = QHrelohihi[offset+j];
      Qjrehilohi = QHrehilohi[offset+j];
      Qjrelolohi = QHrelolohi[offset+j];
      Qjrehihilo = QHrehihilo[offset+j];
      Qjrelohilo = QHrelohilo[offset+j];
      Qjrehilolo = QHrehilolo[offset+j];
      Qjrelololo = QHrelololo[offset+j];
      Qjimhihihi = QHimhihihi[offset+j];
      Qjimlohihi = QHimlohihi[offset+j];
      Qjimhilohi = QHimhilohi[offset+j];
      Qjimlolohi = QHimlolohi[offset+j];
      Qjimhihilo = QHimhihilo[offset+j];
      Qjimlohilo = QHimlohilo[offset+j];
      Qjimhilolo = QHimhilolo[offset+j];
      Qjimlololo = QHimlololo[offset+j];
      bjrehihihi = brehihihi[j];
      bjrelohihi = brelohihi[j];
      bjrehilohi = brehilohi[j];
      bjrelolohi = brelolohi[j];
      bjrehihilo = brehihilo[j];
      bjrelohilo = brelohilo[j];
      bjrehilolo = brehilolo[j];
      bjrelololo = brelololo[j];
      bjimhihihi = bimhihihi[j];
      bjimlohihi = bimlohihi[j];
      bjimhilohi = bimhilohi[j];
      bjimlolohi = bimlolohi[j];
      bjimhihilo = bimhihilo[j];
      bjimlohilo = bimlohilo[j];
      bjimhilolo = bimhilolo[j];
      bjimlololo = bimlololo[j];
      // ri = ri + Qj*bj;
      // zre = Qjre*bjre - Qjim*bjim;
      // rire = rire + zre;
      odg_mul(Qjrehihihi,Qjrelohihi,Qjrehilohi,Qjrelolohi,
              Qjrehihilo,Qjrelohilo,Qjrehilolo,Qjrelololo,
              bjrehihihi,bjrelohihi,bjrehilohi,bjrelolohi,
              bjrehihilo,bjrelohilo,bjrehilolo,bjrelololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&rirehihihi,&rirelohihi,&rirehilohi,&rirelolohi,
              &rirehihilo,&rirelohilo,&rirehilolo,&rirelololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      odg_mul(Qjimhihihi,Qjimlohihi,Qjimhilohi,Qjimlolohi,
              Qjimhihilo,Qjimlohilo,Qjimhilolo,Qjimlololo,
              bjimhihihi,bjimlohihi,bjimhilohi,bjimlolohi,
              bjimhihilo,bjimlohilo,bjimhilolo,bjimlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&rirehihihi,&rirelohihi,&rirehilohi,&rirelolohi,
              &rirehihilo,&rirelohilo,&rirehilolo,&rirelololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      // zim = Qjre*bjim + Qjim*bjre;
      // riim = riim + zim;
      odg_mul(Qjrehihihi,Qjrelohihi,Qjrehilohi,Qjrelolohi,
              Qjrehihilo,Qjrelohilo,Qjrehilolo,Qjrelololo,
              bjimhihihi,bjimlohihi,bjimhilohi,bjimlolohi,
              bjimhihilo,bjimlohilo,bjimhilolo,bjimlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&riimhihihi,&riimlohihi,&riimhilohi,&riimlolohi,
              &riimhihilo,&riimlohilo,&riimhilolo,&riimlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      odg_mul(Qjimhihihi,Qjimlohihi,Qjimhilohi,Qjimlolohi,
              Qjimhihilo,Qjimlohilo,Qjimhilolo,Qjimlololo,
              bjrehihihi,bjrelohihi,bjrehilohi,bjrelolohi,
              bjrehihilo,bjrelohilo,bjrehilolo,bjrelololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&riimhihihi,&riimlohihi,&riimhilohi,&riimlolohi,
              &riimhihilo,&riimlohilo,&riimhilolo,&riimlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
   }
   rrehihihi[idx] = rirehihihi; rrelohihi[idx] = rirelohihi;
   rrehilohi[idx] = rirehilohi; rrelolohi[idx] = rirelolohi;
   rrehihilo[idx] = rirehihilo; rrelohilo[idx] = rirelohilo;
   rrehilolo[idx] = rirehilolo; rrelololo[idx] = rirelololo;
   rimhihihi[idx] = riimhihihi; rimlohihi[idx] = riimlohihi;
   rimhilohi[idx] = riimhilohi; rimlolohi[idx] = riimlolohi;
   rimhihilo[idx] = riimhihilo; rimlohilo[idx] = riimlohilo;
   rimhilolo[idx] = riimhilolo; rimlololo[idx] = riimlololo;
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

   GPU_dbl8_blocked_houseqr
      (nrows,ncols,szt,nbt,
       Ahihihi,Alohihi,Ahilohi,Alolohi,Ahihilo,Alohilo,Ahilolo,Alololo,
       Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
       Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl8_qrtimeflops
         (0,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^T ..." << endl;

   GPU_dbl8_bals_qtb
      (ncols,szt,nbt,
       Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
       bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
       totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "Qtb[" << i << "] : "
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
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(vrblvl > 1)
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

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl8_bstimeflops
         (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,bsmulcnt,bsdivcnt);

   if(vrblvl > 1)
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

void GPU_cmplx8_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi, 
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo, 
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
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

   GPU_cmplx8_blocked_houseqr
      (nrows,ncols,szt,nbt,
       Arehihihi,Arelohihi,Arehilohi,Arelolohi,
       Arehihilo,Arelohilo,Arehilolo,Arelololo,
       Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
       Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
       Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
       Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
       Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
       Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
       Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
       Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
       Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
       Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl8_qrtimeflops
         (1,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^H ..." << endl;

   GPU_cmplx8_bals_qhb
      (ncols,szt,nbt,
       Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
       Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
       Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
       Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
       brehihihi,brelohihi,brehilohi,brelolohi,
       brehihilo,brelohilo,brehilolo,brelololo,
       bimhihihi,bimlohihi,bimhilohi,bimlolohi,
       bimhihilo,bimlohilo,bimhilolo,bimlololo,totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "QHb[" << i << "] : "
              << brehihihi[i] << "  " << brelohihi[i] << endl << "  "
              << brehilohi[i] << "  " << brelolohi[i] << endl << "  "
              << brehihilo[i] << "  " << brelohilo[i] << endl << "  "
              << brehilolo[i] << "  " << brelololo[i] << endl << "  "
              << bimhihihi[i] << "  " << bimlohihi[i] << endl << "  "
              << bimhilohi[i] << "  " << bimlolohi[i] << endl << "  "
              << bimhihilo[i] << "  " << bimlohilo[i] << endl << "  "
              << bimhilolo[i] << "  " << bimlololo[i] << endl;
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
                 << Rrehihihi[i][j] << "  " << Rrelohihi[i][j] << endl
                 << "  "
                 << Rrehilohi[i][j] << "  " << Rrelolohi[i][j] << endl
                 << "  "
                 << Rrehihilo[i][j] << "  " << Rrelohilo[i][j] << endl
                 << "  "
                 << Rrehilolo[i][j] << "  " << Rrelololo[i][j] << endl
                 << "  "
                 << Rimhihihi[i][j] << "  " << Rimlohihi[i][j] << endl
                 << "  "
                 << Rimhilohi[i][j] << "  " << Rimlolohi[i][j] << endl
                 << "  "
                 << Rimhihilo[i][j] << "  " << Rimlohilo[i][j] << endl
                 << "  "
                 << Rimhilolo[i][j] << "  " << Rimlololo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << brehihihi[i] << "  " << brelohihi[i] << endl
              << "  "
              << brehilohi[i] << "  " << brelolohi[i] << endl
              << "  "
              << brehihilo[i] << "  " << brelohilo[i] << endl
              << "  "
              << brehilolo[i] << "  " << brelololo[i] << endl
              << "  "
              << bimhihihi[i] << "  " << bimlohihi[i] << endl
              << "  "
              << bimhilohi[i] << "  " << bimlolohi[i] << endl
              << "  "
              << bimhihilo[i] << "  " << bimlohilo[i] << endl
              << "  "
              << bimhilolo[i] << "  " << bimlololo[i] << endl;
   }
   double **workRrehihihi = new double*[nrows]; // work around ...
   double **workRrelohihi = new double*[nrows];
   double **workRrehilohi = new double*[nrows]; 
   double **workRrelolohi = new double*[nrows];
   double **workRrehihilo = new double*[nrows];
   double **workRrelohilo = new double*[nrows];
   double **workRrehilolo = new double*[nrows]; 
   double **workRrelololo = new double*[nrows];
   double **workRimhihihi = new double*[nrows];
   double **workRimlohihi = new double*[nrows];
   double **workRimhilohi = new double*[nrows];
   double **workRimlolohi = new double*[nrows];
   double **workRimhihilo = new double*[nrows];
   double **workRimlohilo = new double*[nrows];
   double **workRimhilolo = new double*[nrows];
   double **workRimlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRrehihihi[i] = new double[ncols];
      workRrelohihi[i] = new double[ncols];
      workRrehilohi[i] = new double[ncols];
      workRrelolohi[i] = new double[ncols];
      workRrehihilo[i] = new double[ncols];
      workRrelohilo[i] = new double[ncols];
      workRrehilolo[i] = new double[ncols];
      workRrelololo[i] = new double[ncols];
      workRimhihihi[i] = new double[ncols];
      workRimlohihi[i] = new double[ncols];
      workRimhilohi[i] = new double[ncols];
      workRimlolohi[i] = new double[ncols];
      workRimhihilo[i] = new double[ncols];
      workRimlohilo[i] = new double[ncols];
      workRimhilolo[i] = new double[ncols];
      workRimlololo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRrehihihi[i][j] = Rrehihihi[i][j];
         workRrelohihi[i][j] = Rrelohihi[i][j];
         workRrehilohi[i][j] = Rrehilohi[i][j];
         workRrelolohi[i][j] = Rrelolohi[i][j];
         workRrehihilo[i][j] = Rrehihilo[i][j];
         workRrelohilo[i][j] = Rrelohilo[i][j];
         workRrehilolo[i][j] = Rrehilolo[i][j];
         workRrelololo[i][j] = Rrelololo[i][j];
         workRimhihihi[i][j] = Rimhihihi[i][j];
         workRimlohihi[i][j] = Rimlohihi[i][j];
         workRimhilohi[i][j] = Rimhilohi[i][j];
         workRimlolohi[i][j] = Rimlolohi[i][j];
         workRimhihilo[i][j] = Rimhihilo[i][j];
         workRimlohilo[i][j] = Rimlohilo[i][j];
         workRimhilolo[i][j] = Rimhilolo[i][j];
         workRimlololo[i][j] = Rimlololo[i][j];
      }
   }
   GPU_cmplx8_upper_tiled_solver
      (ncols,szt,nbt,
       workRrehihihi,workRrelohihi,workRrehilohi,workRrelolohi,
       workRrehihilo,workRrelohilo,workRrehilolo,workRrelololo,
       workRimhihihi,workRimlohihi,workRimhilohi,workRimlolohi,
       workRimhihilo,workRimlohilo,workRimhilolo,workRimlololo,
       brehihihi,brelohihi,brehilohi,brelolohi,
       brehihilo,brelohilo,brehilolo,brelololo,
       bimhihihi,bimlohihi,bimhilohi,bimlolohi,
       bimhihilo,bimlohilo,bimhilolo,bimlololo,
       xrehihihi,xrelohihi,xrehilohi,xrelolohi,
       xrehihilo,xrelohilo,xrehilolo,xrelololo,
       ximhihihi,ximlohihi,ximhilohi,ximlolohi,
       ximhihilo,ximlohilo,ximhilolo,ximlololo,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl8_bstimeflops
         (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,bsmulcnt,bsdivcnt);

   if(vrblvl > 1)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rrehihihi[i][j] << "  " << Rrelohihi[i][j] << endl
                 << "  "
                 << Rrehilohi[i][j] << "  " << Rrelolohi[i][j] << endl
                 << "  "
                 << Rrehihilo[i][j] << "  " << Rrelohilo[i][j] << endl
                 << "  "
                 << Rrehilolo[i][j] << "  " << Rrelololo[i][j] << endl
                 << "  "
                 << Rimhihihi[i][j] << "  " << Rimlohihi[i][j] << endl
                 << "  "
                 << Rimhilohi[i][j] << "  " << Rimlolohi[i][j] << endl
                 << "  "
                 << Rimhihilo[i][j] << "  " << Rimlohilo[i][j] << endl
                 << "  "
                 << Rimhilolo[i][j] << "  " << Rimlololo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << brehihihi[i] << "  " << brelohihi[i] << endl
              << "  "
              << brehilohi[i] << "  " << brelolohi[i] << endl
              << "  "
              << brehihilo[i] << "  " << brelohilo[i] << endl
              << "  "
              << brehilolo[i] << "  " << brelololo[i] << endl
              << "  "
              << bimhihihi[i] << "  " << bimlohihi[i] << endl
              << "  "
              << bimhilohi[i] << "  " << bimlolohi[i] << endl
              << "  "
              << bimhihilo[i] << "  " << bimlohilo[i] << endl
              << "  "
              << bimhilolo[i] << "  " << bimlololo[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRrehihihi[i]); free(workRimhihihi[i]);
      free(workRrehilohi[i]); free(workRimhilohi[i]);
      free(workRrelohihi[i]); free(workRimlohihi[i]);
      free(workRrelolohi[i]); free(workRimlolohi[i]);
      free(workRrehihilo[i]); free(workRimhihilo[i]);
      free(workRrehilolo[i]); free(workRimhilolo[i]);
      free(workRrelohilo[i]); free(workRimlohilo[i]);
      free(workRrelololo[i]); free(workRimlololo[i]);
   }
   free(workRrehihihi); free(workRimhihihi);
   free(workRrehilohi); free(workRimhilohi);
   free(workRrelohihi); free(workRimlohihi);
   free(workRrelolohi); free(workRimlolohi);
   free(workRrehihilo); free(workRimhihilo);
   free(workRrehilolo); free(workRimhilolo);
   free(workRrelohilo); free(workRimlohilo);
   free(workRrelololo); free(workRimlololo);
}

void write_dbl8_qtbflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for Q^T*b : " << lapsms << " milliseconds." << endl;

   long long int flopcnt;
   const long long int longncols2 = ncols*ncols; // to avoid overflow
   if(ctype == 0)
      flopcnt = 270*longncols2 + 1742*longncols2;
      // as many + as * in one inner product
   else
      flopcnt = 4*270*longncols2 + 4*1742*longncols2;
      // for complex *: 2 ops for +, 6 for *, which is 8 in total

   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 8*ncols*ncols;
   else
      bytecnt = 16*ncols*ncols;

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

void GPU_dbl8_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *totqtblapsedms, int vrblvl )
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

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   dbl8_bals_qtb<<<nbt,szt>>>
      (ncols,szt,
       Qthihihi_d,Qtlohihi_d,Qthilohi_d,Qtlolohi_d,
       Qthihilo_d,Qtlohilo_d,Qthilolo_d,Qtlololo_d,
       bhihihi_d,blohihi_d,bhilohi_d,blolohi_d,
       bhihilo_d,blohilo_d,bhilolo_d,blololo_d,
       rhihihi_d,rlohihi_d,rhilohi_d,rlolohi_d,
       rhihilo_d,rlohilo_d,rhilolo_d,rlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(bhihihi,rhihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blohihi,rlohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhilohi,rhilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blolohi,rlolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhihilo,rhihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blohilo,rlohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bhilolo,rhilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blololo,rlololo_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl8_qtbflops(0,ncols,milliseconds);

   free(Qthihihi_h); free(Qtlohihi_h); free(Qthilohi_h); free(Qtlolohi_h);
   free(Qthihilo_h); free(Qtlohilo_h); free(Qthilolo_h); free(Qtlololo_h);

   cudaFree(bhihihi_d); cudaFree(blohihi_d);
   cudaFree(bhilohi_d); cudaFree(blolohi_d);
   cudaFree(bhihilo_d); cudaFree(blohilo_d);
   cudaFree(bhilolo_d); cudaFree(blololo_d);
   cudaFree(rhihihi_d); cudaFree(rlohihi_d);
   cudaFree(rhilohi_d); cudaFree(rlolohi_d);
   cudaFree(rhihilo_d); cudaFree(rlohilo_d);
   cudaFree(rhilolo_d); cudaFree(rlololo_d);
   cudaFree(Qthihihi_d); cudaFree(Qtlohihi_d);
   cudaFree(Qthilohi_d); cudaFree(Qtlolohi_d);
   cudaFree(Qthihilo_d); cudaFree(Qtlohilo_d);
   cudaFree(Qthilolo_d); cudaFree(Qtlololo_d);
}

void GPU_cmplx8_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *totqtblapsedms, int vrblvl )
{
   double *brehihihi_d;
   double *brelohihi_d;
   double *brehilohi_d;
   double *brelolohi_d;
   double *brehihilo_d;
   double *brelohilo_d;
   double *brehilolo_d;
   double *brelololo_d;
   double *bimhihihi_d;
   double *bimlohihi_d;
   double *bimhilohi_d;
   double *bimlolohi_d;
   double *bimhihilo_d;
   double *bimlohilo_d;
   double *bimhilolo_d;
   double *bimlololo_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&brehihihi_d,szrhs);
   cudaMalloc((void**)&brelohihi_d,szrhs);
   cudaMalloc((void**)&brehilohi_d,szrhs);
   cudaMalloc((void**)&brelolohi_d,szrhs);
   cudaMalloc((void**)&brehihilo_d,szrhs);
   cudaMalloc((void**)&brelohilo_d,szrhs);
   cudaMalloc((void**)&brehilolo_d,szrhs);
   cudaMalloc((void**)&brelololo_d,szrhs);
   cudaMalloc((void**)&bimhihihi_d,szrhs);
   cudaMalloc((void**)&bimlohihi_d,szrhs);
   cudaMalloc((void**)&bimhilohi_d,szrhs);
   cudaMalloc((void**)&bimlolohi_d,szrhs);
   cudaMalloc((void**)&bimhihilo_d,szrhs);
   cudaMalloc((void**)&bimlohilo_d,szrhs);
   cudaMalloc((void**)&bimhilolo_d,szrhs);
   cudaMalloc((void**)&bimlololo_d,szrhs);
   cudaMemcpy(brehihihi_d,brehihihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelohihi_d,brelohihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brehilohi_d,brehilohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelolohi_d,brelolohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brehihilo_d,brehihilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelohilo_d,brelohilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brehilolo_d,brehilolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelololo_d,brelololo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhihihi_d,bimhihihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlohihi_d,bimlohihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhilohi_d,bimhilohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlolohi_d,bimlolohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhihilo_d,bimhihilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlohilo_d,bimlohilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhilolo_d,bimhilolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlololo_d,bimlololo,szrhs,cudaMemcpyHostToDevice);

   double *rrehihihi_d;
   double *rrelohihi_d;
   double *rrehilohi_d;
   double *rrelolohi_d;
   double *rrehihilo_d;
   double *rrelohilo_d;
   double *rrehilolo_d;
   double *rrelololo_d;
   double *rimhihihi_d;
   double *rimlohihi_d;
   double *rimhilohi_d;
   double *rimlolohi_d;
   double *rimhihilo_d;
   double *rimlohilo_d;
   double *rimhilolo_d;
   double *rimlololo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rrehihihi_d,szsol);
   cudaMalloc((void**)&rrelohihi_d,szsol);
   cudaMalloc((void**)&rrehilohi_d,szsol);
   cudaMalloc((void**)&rrelolohi_d,szsol);
   cudaMalloc((void**)&rrehihilo_d,szsol);
   cudaMalloc((void**)&rrelohilo_d,szsol);
   cudaMalloc((void**)&rrehilolo_d,szsol);
   cudaMalloc((void**)&rrelololo_d,szsol);
   cudaMalloc((void**)&rimhihihi_d,szsol);
   cudaMalloc((void**)&rimlohihi_d,szsol);
   cudaMalloc((void**)&rimhilohi_d,szsol);
   cudaMalloc((void**)&rimlolohi_d,szsol);
   cudaMalloc((void**)&rimhihilo_d,szsol);
   cudaMalloc((void**)&rimlohilo_d,szsol);
   cudaMalloc((void**)&rimhilolo_d,szsol);
   cudaMalloc((void**)&rimlololo_d,szsol);

   double *QHrehihihi_d;
   double *QHrelohihi_d;
   double *QHrehilohi_d;
   double *QHrelolohi_d;
   double *QHrehihilo_d;
   double *QHrelohilo_d;
   double *QHrehilolo_d;
   double *QHrelololo_d;
   double *QHimhihihi_d;
   double *QHimlohihi_d;
   double *QHimhilohi_d;
   double *QHimlolohi_d;
   double *QHimhihilo_d;
   double *QHimlohilo_d;
   double *QHimhilolo_d;
   double *QHimlololo_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&QHrehihihi_d,szmat);
   cudaMalloc((void**)&QHrelohihi_d,szmat);
   cudaMalloc((void**)&QHrehilohi_d,szmat);
   cudaMalloc((void**)&QHrelolohi_d,szmat);
   cudaMalloc((void**)&QHrehihilo_d,szmat);
   cudaMalloc((void**)&QHrelohilo_d,szmat);
   cudaMalloc((void**)&QHrehilolo_d,szmat);
   cudaMalloc((void**)&QHrelololo_d,szmat);
   cudaMalloc((void**)&QHimhihihi_d,szmat);
   cudaMalloc((void**)&QHimlohihi_d,szmat);
   cudaMalloc((void**)&QHimhilohi_d,szmat);
   cudaMalloc((void**)&QHimlolohi_d,szmat);
   cudaMalloc((void**)&QHimhihilo_d,szmat);
   cudaMalloc((void**)&QHimlohilo_d,szmat);
   cudaMalloc((void**)&QHimhilolo_d,szmat);
   cudaMalloc((void**)&QHimlololo_d,szmat);

   double *QHrehihihi_h = new double[szmat];
   double *QHrelohihi_h = new double[szmat];
   double *QHrehilohi_h = new double[szmat];
   double *QHrelolohi_h = new double[szmat];
   double *QHrehihilo_h = new double[szmat];
   double *QHrelohilo_h = new double[szmat];
   double *QHrehilolo_h = new double[szmat];
   double *QHrelololo_h = new double[szmat];
   double *QHimhihihi_h = new double[szmat];
   double *QHimlohihi_h = new double[szmat];
   double *QHimhilohi_h = new double[szmat];
   double *QHimlolohi_h = new double[szmat];
   double *QHimhihilo_h = new double[szmat];
   double *QHimlohilo_h = new double[szmat];
   double *QHimhilolo_h = new double[szmat];
   double *QHimlololo_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         QHrehihihi_h[idx]   = Qrehihihi[j][i];
         QHrelohihi_h[idx]   = Qrelohihi[j][i];
         QHrehilohi_h[idx]   = Qrehilohi[j][i];
         QHrelolohi_h[idx]   = Qrelolohi[j][i];
         QHrehihilo_h[idx]   = Qrehihilo[j][i];
         QHrelohilo_h[idx]   = Qrelohilo[j][i];
         QHrehilolo_h[idx]   = Qrehilolo[j][i];
         QHrelololo_h[idx]   = Qrelololo[j][i];
         QHimhihihi_h[idx]   = - Qimhihihi[j][i]; // Hermitian transpose !
         QHimlohihi_h[idx]   = - Qimlohihi[j][i]; // Hermitian transpose !
         QHimhilohi_h[idx]   = - Qimhilohi[j][i]; // Hermitian transpose !
         QHimlolohi_h[idx]   = - Qimlolohi[j][i]; // Hermitian transpose !
         QHimhihilo_h[idx]   = - Qimhihilo[j][i]; // Hermitian transpose !
         QHimlohilo_h[idx]   = - Qimlohilo[j][i]; // Hermitian transpose !
         QHimhilolo_h[idx]   = - Qimhilolo[j][i]; // Hermitian transpose !
         QHimlololo_h[idx++] = - Qimlololo[j][i]; // Hermitian transpose !
      }

   cudaMemcpy(QHrehihihi_d,QHrehihihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelohihi_d,QHrelohihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrehilohi_d,QHrehilohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelolohi_d,QHrelolohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrehihilo_d,QHrehihilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelohilo_d,QHrelohilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrehilolo_d,QHrehilolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelololo_d,QHrelololo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhihihi_d,QHimhihihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlohihi_d,QHimlohihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhilohi_d,QHimhilohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlolohi_d,QHimlolohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhihilo_d,QHimhihilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlohilo_d,QHimlohilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhilolo_d,QHimhilolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlololo_d,QHimlololo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   cmplx8_bals_qhb<<<nbt,szt>>>
      (ncols,szt,QHrehihihi_d,QHrelohihi_d,QHrehilohi_d,QHrelolohi_d,
                 QHrehihilo_d,QHrelohilo_d,QHrehilolo_d,QHrelololo_d,
                 QHimhihihi_d,QHimlohihi_d,QHimhilohi_d,QHimlolohi_d,
                 QHimhihilo_d,QHimlohilo_d,QHimhilolo_d,QHimlololo_d,
       brehihihi_d,brelohihi_d,brehilohi_d,brelolohi_d,
       brehihilo_d,brelohilo_d,brehilolo_d,brelololo_d,
       bimhihihi_d,bimlohihi_d,bimhilohi_d,bimlolohi_d,
       bimhihilo_d,bimlohilo_d,bimhilolo_d,bimlololo_d,
       rrehihihi_d,rrelohihi_d,rrehilohi_d,rrelolohi_d,
       rrehihilo_d,rrelohilo_d,rrehilolo_d,rrelololo_d,
       rimhihihi_d,rimlohihi_d,rimhilohi_d,rimlolohi_d,
       rimhihilo_d,rimlohilo_d,rimhilolo_d,rimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(brehihihi,rrehihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelohihi,rrelohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brehilohi,rrehilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelolohi,rrelolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brehihilo,rrehihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelohilo,rrelohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brehilolo,rrehilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelololo,rrelololo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhihihi,rimhihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlohihi,rimlohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhilohi,rimhilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlolohi,rimlolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhihilo,rimhihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlohilo,rimlohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhilolo,rimhilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlololo,rimlololo_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl8_qtbflops(1,ncols,milliseconds);

   cudaFree(brehihihi_d); cudaFree(brelohihi_d);
   cudaFree(brehilohi_d); cudaFree(brelolohi_d);
   cudaFree(brehihilo_d); cudaFree(brelohilo_d);
   cudaFree(brehilolo_d); cudaFree(brelololo_d);
   cudaFree(bimhihihi_d); cudaFree(bimlohihi_d);
   cudaFree(bimhilohi_d); cudaFree(bimlolohi_d);
   cudaFree(bimhihilo_d); cudaFree(bimlohilo_d);
   cudaFree(bimhilolo_d); cudaFree(bimlololo_d);
   cudaFree(rrehihihi_d); cudaFree(rrelohihi_d);
   cudaFree(rrehilohi_d); cudaFree(rrelolohi_d);
   cudaFree(rrehihilo_d); cudaFree(rrelohilo_d);
   cudaFree(rrehilolo_d); cudaFree(rrelololo_d);
   cudaFree(rimhihihi_d); cudaFree(rimlohihi_d);
   cudaFree(rimhilohi_d); cudaFree(rimlolohi_d);
   cudaFree(rimhihilo_d); cudaFree(rimlohilo_d);
   cudaFree(rimhilolo_d); cudaFree(rimlololo_d);
   cudaFree(QHrehihihi_d); cudaFree(QHrelohihi_d);
   cudaFree(QHrehilohi_d); cudaFree(QHrelolohi_d);
   cudaFree(QHrehihilo_d); cudaFree(QHrelohilo_d);
   cudaFree(QHrehilolo_d); cudaFree(QHrelololo_d);
   cudaFree(QHimhihihi_d); cudaFree(QHimlohihi_d);
   cudaFree(QHimhilohi_d); cudaFree(QHimlolohi_d);
   cudaFree(QHimhihilo_d); cudaFree(QHimlohilo_d);
   cudaFree(QHimhilolo_d); cudaFree(QHimlololo_d);

   free(QHrehihihi_h); free(QHimhihihi_h);
   free(QHrelohihi_h); free(QHimlohihi_h);
   free(QHrehilohi_h); free(QHimhilohi_h);
   free(QHrelolohi_h); free(QHimlolohi_h);
   free(QHrehihilo_h); free(QHimhihilo_h);
   free(QHrelohilo_h); free(QHimlohilo_h);
   free(QHrehilolo_h); free(QHimhilolo_h);
   free(QHrelololo_h); free(QHimlololo_h);
}

void GPU_dbl8_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
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
   double **solhilolo, double **sollololo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped substitutions
   const double prevnorm = 1.0e+120;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

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
   if(vrblvl > 1)
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
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping GPU_dbl8_bals_head ..." << endl;
   }
   else
   {
      CPU_dbl_onenorm(nrows,rhshihihi[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm < 1.0e-120)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "no skipping GPU_dbl8_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "skip call to GPU_dbl8_bals_head ..." << endl;

            *noqr = true;

            for(int j=0; j<ncols; j++)
            {
               solhihihi[0][j] = 0.0; sollohihi[0][j] = 0.0;
               solhilohi[0][j] = 0.0; sollolohi[0][j] = 0.0;
               solhihilo[0][j] = 0.0; sollohilo[0][j] = 0.0;
               solhilolo[0][j] = 0.0; sollololo[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "calling GPU_dbl8_bals_head ..." << endl;

         double **Ahihihi = new double*[nrows];
         double **Alohihi = new double*[nrows];
         double **Ahilohi = new double*[nrows];
         double **Alolohi = new double*[nrows];
         double **Ahihilo = new double*[nrows];
         double **Alohilo = new double*[nrows];
         double **Ahilolo = new double*[nrows];
         double **Alololo = new double*[nrows];

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
             xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
             totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         for(int j=0; j<ncols; j++)
         {
            solhihihi[0][j] = xhihihi[j]; sollohihi[0][j] = xlohihi[j];
            solhilohi[0][j] = xhilohi[j]; sollolohi[0][j] = xlolohi[j];
            solhihilo[0][j] = xhihilo[j]; sollohilo[0][j] = xlohilo[j];
            solhilolo[0][j] = xhilolo[j]; sollololo[0][j] = xlololo[j];
         }
         for(int i=0; i<nrows; i++)
         {
            free(Ahihihi[i]); free(Alohihi[i]);
            free(Ahilohi[i]); free(Alolohi[i]);
            free(Ahihilo[i]); free(Alohilo[i]);
            free(Ahilolo[i]); free(Alololo[i]);
         }
         free(Ahihihi); free(Alohihi); free(Ahilohi); free(Alolohi);
         free(Ahihilo); free(Alohilo); free(Ahilolo); free(Alololo);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xshihihi = solhihihi[stage-1]; // solution used in update
      CPU_dbl_onenorm(dim,xshihihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-120)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << stage-1 << "] ..." << endl;

         GPU_dbl8_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,
             mathihihi,matlohihi,mathilohi,matlolohi,
             mathihilo,matlohilo,mathilolo,matlololo,
             rhshihihi,rhslohihi,rhshilohi,rhslolohi,
             rhshihilo,rhslohilo,rhshilolo,rhslololo,
             solhihihi,sollohihi,solhilohi,sollolohi,
             solhihilo,sollohilo,solhilolo,sollololo,totupdlapsedms,vrblvl);

         if(vrblvl > 1)
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
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;

         bhihihi[i] = rhshihihi[stage][i]; blohihi[i] = rhslohihi[stage][i];
         bhilohi[i] = rhshilohi[stage][i]; blolohi[i] = rhslolohi[stage][i];
         bhihilo[i] = rhshihilo[stage][i]; blohilo[i] = rhslohilo[stage][i];
         bhilolo[i] = rhshilolo[stage][i]; blololo[i] = rhslololo[stage][i];

         // cout << "b[" << i << "] : "
         //      << bhihihi[i] << "  " << blohihi[i] << endl << "  "
         //      << bhilohi[i] << "  " << blolohi[i] << endl << "  "
         //      << bhihilo[i] << "  " << blohilo[i] << endl << "  "
         //      << bhilolo[i] << "  " << blololo[i] << endl;
      }
      CPU_dbl_onenorm(nrows,bhihihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-120) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xhihihi[i] = 0.0; xlohihi[i] = 0.0;
            xhilohi[i] = 0.0; xlolohi[i] = 0.0;
            xhihilo[i] = 0.0; xlohilo[i] = 0.0;
            xhilolo[i] = 0.0; xlololo[i] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+120; // nrm*1.0e+8;

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

         if(vrblvl > 0)
            cout << "-> GPU multiplies rhs with Q^T ..." << endl;

         GPU_dbl8_bals_qtb
            (ncols,szt,nbt,
             Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
             bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
             totqtblapsedms,vrblvl);

         if(vrblvl > 1)
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
 
            if(vrblvl > 1)
               for(int i=0; i<nrows; i++)
                  for(int j=0; j<ncols; j++)
                     cout << "R[" << i << "][" << j << "] : "
                          << Rhihihi[i][j] << "  " << Rlohihi[i][j] << "  "
                          << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl
                          << "  "
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

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl8_bstimeflops
               (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,bsmulcnt,bsdivcnt);

         if(vrblvl > 1)
            for(int i=0; i<ncols; i++)
               cout << "x[" << i << "] : "
                    << xhihihi[i] << "  " << xlohihi[i] << "  "
                    << xhilohi[i] << "  " << xlolohi[i] << endl << "  "
                    << xhihilo[i] << "  " << xlohilo[i] << "  "
                    << xhilolo[i] << "  " << xlololo[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solhihihi[stage][j] = xhihihi[j]; sollohihi[stage][j] = xlohihi[j];
         solhilohi[stage][j] = xhilohi[j]; sollolohi[stage][j] = xlolohi[j];
         solhihilo[stage][j] = xhihilo[j]; sollohilo[stage][j] = xlohilo[j];
         solhilolo[stage][j] = xhilolo[j]; sollololo[stage][j] = xlololo[j];
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
      free(workRhihihi[i]); free(workRlohihi[i]);
      free(workRhilohi[i]); free(workRlolohi[i]);
      free(workRhihilo[i]); free(workRlohilo[i]);
      free(workRhilolo[i]); free(workRlololo[i]);
   }
   free(bhihihi); free(xhihihi); free(workRhihihi);
   free(blohihi); free(xlohihi); free(workRlohihi);
   free(bhilohi); free(xhilohi); free(workRhilohi);
   free(blolohi); free(xlolohi); free(workRlolohi);
   free(bhihilo); free(xhihilo); free(workRhihilo);
   free(blohilo); free(xlohilo); free(workRlohilo);
   free(bhilolo); free(xhilolo); free(workRhilolo);
   free(blololo); free(xlololo); free(workRlololo);
}

void GPU_cmplx8_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi, 
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo, 
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi, 
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo, 
   double **solimhilolo, double **solimlololo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped substitutions
   const double prevnorm = 1.0e+120;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   double *brehihihi = new double[nrows];
   double *brelohihi = new double[nrows];
   double *brehilohi = new double[nrows];
   double *brelolohi = new double[nrows];
   double *brehihilo = new double[nrows];
   double *brelohilo = new double[nrows];
   double *brehilolo = new double[nrows];
   double *brelololo = new double[nrows];
   double *bimhihihi = new double[nrows];
   double *bimlohihi = new double[nrows];
   double *bimhilohi = new double[nrows];
   double *bimlolohi = new double[nrows];
   double *bimhihilo = new double[nrows];
   double *bimlohilo = new double[nrows];
   double *bimhilolo = new double[nrows];
   double *bimlololo = new double[nrows];
   double *xrehihihi = new double[ncols];
   double *xrelohihi = new double[ncols];
   double *xrehilohi = new double[ncols];
   double *xrelolohi = new double[ncols];
   double *xrehihilo = new double[ncols];
   double *xrelohilo = new double[ncols];
   double *xrehilolo = new double[ncols];
   double *xrelololo = new double[ncols];
   double *ximhihihi = new double[ncols];
   double *ximlohihi = new double[ncols];
   double *ximhilohi = new double[ncols];
   double *ximlolohi = new double[ncols];
   double *ximhihilo = new double[ncols];
   double *ximlohilo = new double[ncols];
   double *ximhilolo = new double[ncols];
   double *ximlololo = new double[ncols];

   double **workRrehihihi = new double*[nrows]; // GPU upper solver changes R
   double **workRrelohihi = new double*[nrows];
   double **workRrehilohi = new double*[nrows];
   double **workRrelolohi = new double*[nrows];
   double **workRrehihilo = new double*[nrows];
   double **workRrelohilo = new double*[nrows];
   double **workRrehilolo = new double*[nrows];
   double **workRrelololo = new double*[nrows];
   double **workRimhihihi = new double*[nrows]; 
   double **workRimlohihi = new double*[nrows]; 
   double **workRimhilohi = new double*[nrows]; 
   double **workRimlolohi = new double*[nrows]; 
   double **workRimhihilo = new double*[nrows]; 
   double **workRimlohilo = new double*[nrows]; 
   double **workRimhilolo = new double*[nrows]; 
   double **workRimlololo = new double*[nrows]; 

   for(int i=0; i<nrows; i++)
   {
      workRrehihihi[i] = new double[ncols];
      workRrelohihi[i] = new double[ncols];
      workRrehilohi[i] = new double[ncols];
      workRrelolohi[i] = new double[ncols];
      workRrehihilo[i] = new double[ncols];
      workRrelohilo[i] = new double[ncols];
      workRrehilolo[i] = new double[ncols];
      workRrelololo[i] = new double[ncols];
      workRimhihihi[i] = new double[ncols];
      workRimlohihi[i] = new double[ncols];
      workRimhilohi[i] = new double[ncols];
      workRimlolohi[i] = new double[ncols];
      workRimhihilo[i] = new double[ncols];
      workRimlohilo[i] = new double[ncols];
      workRimhilolo[i] = new double[ncols];
      workRimlololo[i] = new double[ncols];
   }
   if(vrblvl > 1)
   {
      cout << "GPU_cmplx8_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehihihi[k][i] << "  " << rhsrelohihi[k][i] << endl
                 << "  "
                 << rhsrehilohi[k][i] << "  " << rhsrelolohi[k][i] << endl
                 << "  "
                 << rhsrehihilo[k][i] << "  " << rhsrelohilo[k][i] << endl
                 << "  "
                 << rhsrehilolo[k][i] << "  " << rhsrelololo[k][i] << endl
                 << "  "
                 << rhsimhihihi[k][i] << "  " << rhsimlohihi[k][i] << endl
                 << "  "
                 << rhsimhilohi[k][i] << "  " << rhsimlolohi[k][i] << endl
                 << "  "
                 << rhsimhihilo[k][i] << "  " << rhsimlohilo[k][i] << endl
                 << "  "
                 << rhsimhilolo[k][i] << "  " << rhsimlololo[k][i] << endl;
      }
   }
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "-> calling GPU_cmplx8_bals_head ..." << endl;
   }
   else
   {
      CPU_cmplx_onenorm(nrows,rhsrehihihi[0],rhsimhihihi[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm < 1.0e-120)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "-> no skipping GPU_cmplx8_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "-> skip call to GPU_cmplx8_bals_head ..." << endl;

            *noqr = true;

            for(int j=0; j<ncols; j++)
            {
               solrehihihi[0][j] = 0.0; solrelohihi[0][j] = 0.0;
               solrehilohi[0][j] = 0.0; solrelolohi[0][j] = 0.0;
               solrehihilo[0][j] = 0.0; solrelohilo[0][j] = 0.0;
               solrehilolo[0][j] = 0.0; solrelololo[0][j] = 0.0;
               solimhihihi[0][j] = 0.0; solimlohihi[0][j] = 0.0;
               solimhilohi[0][j] = 0.0; solimlolohi[0][j] = 0.0;
               solimhihilo[0][j] = 0.0; solimlohilo[0][j] = 0.0;
               solimhilolo[0][j] = 0.0; solimlololo[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "-> calling GPU_cmplx8_bals_head ..." << endl;

         double **Arehihihi = new double*[nrows];
         double **Arelohihi = new double*[nrows];
         double **Arehilohi = new double*[nrows];
         double **Arelolohi = new double*[nrows];
         double **Arehihilo = new double*[nrows];
         double **Arelohilo = new double*[nrows];
         double **Arehilolo = new double*[nrows];
         double **Arelololo = new double*[nrows];
         double **Aimhihihi = new double*[nrows];
         double **Aimlohihi = new double*[nrows];
         double **Aimhilohi = new double*[nrows];
         double **Aimlolohi = new double*[nrows];
         double **Aimhihilo = new double*[nrows];
         double **Aimlohilo = new double*[nrows];
         double **Aimhilolo = new double*[nrows];
         double **Aimlololo = new double*[nrows];

         for(int i=0; i<nrows; i++)
         {
            Arehihihi[i] = new double[ncols];
            Arelohihi[i] = new double[ncols];
            Arehilohi[i] = new double[ncols];
            Arelolohi[i] = new double[ncols];
            Arehihilo[i] = new double[ncols];
            Arelohilo[i] = new double[ncols];
            Arehilolo[i] = new double[ncols];
            Arelololo[i] = new double[ncols];
            Aimhihihi[i] = new double[ncols];
            Aimlohihi[i] = new double[ncols];
            Aimhilohi[i] = new double[ncols];
            Aimlolohi[i] = new double[ncols];
            Aimhihilo[i] = new double[ncols];
            Aimlohilo[i] = new double[ncols];
            Aimhilolo[i] = new double[ncols];
            Aimlololo[i] = new double[ncols];

            for(int j=0; j<ncols; j++)
            {
               Arehihihi[i][j] = matrehihihi[0][i][j];
               Arelohihi[i][j] = matrelohihi[0][i][j];
               Arehilohi[i][j] = matrehilohi[0][i][j];
               Arelolohi[i][j] = matrelolohi[0][i][j];
               Arehihilo[i][j] = matrehihilo[0][i][j];
               Arelohilo[i][j] = matrelohilo[0][i][j];
               Arehilolo[i][j] = matrehilolo[0][i][j];
               Arelololo[i][j] = matrelololo[0][i][j];
               Aimhihihi[i][j] = matimhihihi[0][i][j];
               Aimlohihi[i][j] = matimlohihi[0][i][j];
               Aimhilohi[i][j] = matimhilohi[0][i][j];
               Aimlolohi[i][j] = matimlolohi[0][i][j];
               Aimhihilo[i][j] = matimhihilo[0][i][j];
               Aimlohilo[i][j] = matimlohilo[0][i][j];
               Aimhilolo[i][j] = matimhilolo[0][i][j];
               Aimlololo[i][j] = matimlololo[0][i][j];
            }
            brehihihi[i] = rhsrehihihi[0][i]; brelohihi[i] = rhsrelohihi[0][i];
            brehilohi[i] = rhsrehilohi[0][i]; brelolohi[i] = rhsrelolohi[0][i];
            brehihilo[i] = rhsrehihilo[0][i]; brelohilo[i] = rhsrelohilo[0][i];
            brehilolo[i] = rhsrehilolo[0][i]; brelololo[i] = rhsrelololo[0][i];
            bimhihihi[i] = rhsimhihihi[0][i]; bimlohihi[i] = rhsimlohihi[0][i];
            bimhilohi[i] = rhsimhilohi[0][i]; bimlolohi[i] = rhsimlolohi[0][i];
            bimhihilo[i] = rhsimhihilo[0][i]; bimlohilo[i] = rhsimlohilo[0][i];
            bimhilolo[i] = rhsimhilolo[0][i]; bimlololo[i] = rhsimlololo[0][i];

            for(int j=0; j<ncols; j++)
            {
               Rrehihihi[i][j] = matrehihihi[0][i][j];
               Rrelohihi[i][j] = matrelohihi[0][i][j];
               Rrehilohi[i][j] = matrehilohi[0][i][j];
               Rrelolohi[i][j] = matrelolohi[0][i][j];
               Rrehihilo[i][j] = matrehihilo[0][i][j];
               Rrelohilo[i][j] = matrelohilo[0][i][j];
               Rrehilolo[i][j] = matrehilolo[0][i][j];
               Rrelololo[i][j] = matrelololo[0][i][j];
               Rimhihihi[i][j] = matimhihihi[0][i][j];
               Rimlohihi[i][j] = matimlohihi[0][i][j];
               Rimhilohi[i][j] = matimhilohi[0][i][j];
               Rimlolohi[i][j] = matimlolohi[0][i][j];
               Rimhihilo[i][j] = matimhihilo[0][i][j];
               Rimlohilo[i][j] = matimlohilo[0][i][j];
               Rimhilolo[i][j] = matimhilolo[0][i][j];
               Rimlololo[i][j] = matimlololo[0][i][j];
            }
         }
         GPU_cmplx8_bals_head
            (nrows,ncols,szt,nbt,
             Arehihihi,Arelohihi,Arehilohi,Arelolohi,
             Arehihilo,Arelohilo,Arehilolo,Arelololo,
             Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
             Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
             Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
             Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
             Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
             Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
             Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
             Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
             Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
             Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
             brehihihi,brelohihi,brehilohi,brelolohi,
             brehihilo,brelohilo,brehilolo,brelololo,
             bimhihihi,bimlohihi,bimhilohi,bimlolohi,
             bimhihilo,bimlohilo,bimhilolo,bimlololo,
             xrehihihi,xrelohihi,xrehilohi,xrelolohi,
             xrehihilo,xrelohilo,xrehilolo,xrelololo,
             ximhihihi,ximlohihi,ximhilohi,ximlolohi,
             ximhihilo,ximlohilo,ximhilolo,ximlololo,
             totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         for(int j=0; j<ncols; j++)
         {
            solrehihihi[0][j] = xrehihihi[j];
            solrelohihi[0][j] = xrelohihi[j];
            solrehilohi[0][j] = xrehilohi[j];
            solrelolohi[0][j] = xrelolohi[j];
            solrehihilo[0][j] = xrehihilo[j];
            solrelohilo[0][j] = xrelohilo[j];
            solrehilolo[0][j] = xrehilolo[j];
            solrelololo[0][j] = xrelololo[j];
            solimhihihi[0][j] = ximhihihi[j];
            solimlohihi[0][j] = ximlohihi[j];
            solimhilohi[0][j] = ximhilohi[j];
            solimlolohi[0][j] = ximlolohi[j];
            solimhihilo[0][j] = ximhihilo[j];
            solimlohilo[0][j] = ximlohilo[j];
            solimhilolo[0][j] = ximhilolo[j];
            solimlololo[0][j] = ximlololo[j];
         }
         for(int i=0; i<nrows; i++)
         {
            free(Arehihihi[i]); free(Arelohihi[i]);
            free(Arehilohi[i]); free(Arelolohi[i]);
            free(Arehihilo[i]); free(Arelohilo[i]);
            free(Arehilolo[i]); free(Arelololo[i]);
            free(Aimhihihi[i]); free(Aimlohihi[i]);
            free(Aimhilohi[i]); free(Aimlolohi[i]);
            free(Aimhihilo[i]); free(Aimlohilo[i]);
            free(Aimhilolo[i]); free(Aimlololo[i]);
         }
         free(Arehihihi); free(Arelohihi); free(Arehilohi); free(Arelolohi);
         free(Arehihilo); free(Arelohilo); free(Arehilolo); free(Arelololo);
         free(Aimhihihi); free(Aimlohihi); free(Aimhilohi); free(Aimlolohi);
         free(Aimhihilo); free(Aimlohilo); free(Aimhilolo); free(Aimlololo);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xrshihihi = solrehihihi[stage-1]; // solution used in update
      double *xishihihi = solimhihihi[stage-1];
      CPU_cmplx_onenorm(dim,xrshihihi,xishihihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-120)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << stage-1 << "] ..." << endl;

         GPU_cmplx8_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,
             matrehihihi,matrelohihi,matrehilohi,matrelolohi,
             matrehihilo,matrelohilo,matrehilolo,matrelololo,
             matimhihihi,matimlohihi,matimhilohi,matimlolohi,
             matimhihilo,matimlohilo,matimhilolo,matimlololo,
             rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
             rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
             rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
             rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
             solrehihihi,solrelohihi,solrehilohi,solrelolohi,
             solrehihilo,solrelohilo,solrehilolo,solrelololo,
             solimhihihi,solimlohihi,solimhilohi,solimlolohi,
             solimhihilo,solimlohilo,solimhilolo,solimlololo,
             totupdlapsedms,vrblvl);

         if(vrblvl > 1)
         {
            cout << "blocks of rhs before assignment :" << endl;
            for(int k=0; k<degp1; k++)
            {
               for(int i=0; i<nrows; i++)
                  cout << "rhs[" << k << "][" << i << "] : "
                       << rhsrehihihi[k][i] << "  "
                       << rhsrelohihi[k][i] << endl << "  " 
                       << rhsrehilohi[k][i] << "  "
                       << rhsrelolohi[k][i] << endl << "  " 
                       << rhsrehihilo[k][i] << "  "
                       << rhsrelohilo[k][i] << endl << "  " 
                       << rhsrehilolo[k][i] << "  "
                       << rhsrelololo[k][i] << endl << "  " 
                       << rhsimhihihi[k][i] << "  "
                       << rhsimlohihi[k][i] << endl << "  "
                       << rhsimhilohi[k][i] << "  "
                       << rhsimlolohi[k][i] << endl << "  "
                       << rhsimhihilo[k][i] << "  "
                       << rhsimlohilo[k][i] << endl << "  "
                       << rhsimhilolo[k][i] << "  "
                       << rhsimlololo[k][i] << endl;
            }
         }
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         brehihihi[i] = rhsrehihihi[stage][i];
         brelohihi[i] = rhsrelohihi[stage][i];
         brehilohi[i] = rhsrehilohi[stage][i];
         brelolohi[i] = rhsrelolohi[stage][i];
         brehihilo[i] = rhsrehihilo[stage][i];
         brelohilo[i] = rhsrelohilo[stage][i];
         brehilolo[i] = rhsrehilolo[stage][i];
         brelololo[i] = rhsrelololo[stage][i];
         bimhihihi[i] = rhsimhihihi[stage][i];
         bimlohihi[i] = rhsimlohihi[stage][i];
         bimhilohi[i] = rhsimhilohi[stage][i];
         bimlolohi[i] = rhsimlolohi[stage][i];
         bimhihilo[i] = rhsimhihilo[stage][i];
         bimlohilo[i] = rhsimlohilo[stage][i];
         bimhilolo[i] = rhsimhilolo[stage][i];
         bimlololo[i] = rhsimlololo[stage][i];
/*
         cout << "b[" << i << "] : "
              << brehihihi[i] << "  " << brelohihi[i] << endl << "  "
              << brehilohi[i] << "  " << brelolohi[i] << endl << "  "
              << brehihilo[i] << "  " << brelohilo[i] << endl << "  "
              << brehilolo[i] << "  " << brelololo[i] << endl << "  "
              << bimhihihi[i] << "  " << bimlohihi[i] << endl << "  "
              << bimhilohi[i] << "  " << bimlolohi[i] << endl << "  "
              << bimhihilo[i] << "  " << bimlohilo[i] << endl << "  "
              << bimhilolo[i] << "  " << bimlololo[i] << endl;
 */
      }
      CPU_cmplx_onenorm(dim,brehihihi,bimhihihi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-120) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitutions for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xrehihihi[i] = 0.0; xrelohihi[i] = 0.0;
            xrehilohi[i] = 0.0; xrelolohi[i] = 0.0;
            xrehihilo[i] = 0.0; xrelohilo[i] = 0.0;
            xrehilolo[i] = 0.0; xrelololo[i] = 0.0;
            ximhihihi[i] = 0.0; ximlohihi[i] = 0.0;
            ximhilohi[i] = 0.0; ximlolohi[i] = 0.0;
            ximhihilo[i] = 0.0; ximlohilo[i] = 0.0;
            ximhilolo[i] = 0.0; ximlololo[i] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+120; // nrm*1.0e+8;

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

         if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^H ..." << endl;

         GPU_cmplx8_bals_qhb
            (ncols,szt,nbt,
             Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
             Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
             Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
             Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
             brehihihi,brelohihi,brehilohi,brelolohi,
             brehihilo,brelohilo,brehilolo,brelololo,
             bimhihihi,bimlohihi,bimhilohi,bimlolohi,
             bimhihilo,bimlohilo,bimhilolo,bimlololo,totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "QHb[" << i << "] : "
                    << brehihihi[i] << "  " << brelohihi[i] << endl << "  "
                    << brehilohi[i] << "  " << brelolohi[i] << endl << "  "
                    << brehihilo[i] << "  " << brelohilo[i] << endl << "  "
                    << brehilolo[i] << "  " << brelololo[i] << endl << "  "
                    << bimhihihi[i] << "  " << bimlohihi[i] << endl << "  "
                    << bimhilohi[i] << "  " << bimlolohi[i] << endl << "  "
                    << bimhihilo[i] << "  " << bimlohilo[i] << endl << "  "
                    << bimhilolo[i] << "  " << bimlololo[i] << endl << "  ";
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
 
            if(vrblvl > 1)
               for(int i=0; i<nrows; i++)
                  for(int j=0; j<ncols; j++)
                     cout << "R[" << i << "][" << j << "] : "
                          << Rrehihihi[i][j] << "  " << Rrelohihi[i][j]
                          << endl << "  "
                          << Rrehilohi[i][j] << "  " << Rrelolohi[i][j]
                          << endl << "  "
                          << Rrehihilo[i][j] << "  " << Rrelohilo[i][j]
                          << endl << "  "
                          << Rrehilolo[i][j] << "  " << Rrelololo[i][j]
                          << endl << "  "
                          << Rimhihihi[i][j] << "  " << Rimlohihi[i][j]
                          << endl << "  "
                          << Rimhilohi[i][j] << "  " << Rimlolohi[i][j]
                          << endl << "  "
                          << Rimhihilo[i][j] << "  " << Rimlohilo[i][j]
                          << endl << "  "
                          << Rimhilolo[i][j] << "  " << Rimlololo[i][j]
                          << endl;
         }
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
            {
               workRrehihihi[i][j] = Rrehihihi[i][j];
               workRrelohihi[i][j] = Rrelohihi[i][j];
               workRrehilohi[i][j] = Rrehilohi[i][j];
               workRrelolohi[i][j] = Rrelolohi[i][j];
               workRrehihilo[i][j] = Rrehihilo[i][j];
               workRrelohilo[i][j] = Rrelohilo[i][j];
               workRrehilolo[i][j] = Rrehilolo[i][j];
               workRrelololo[i][j] = Rrelololo[i][j];
               workRimhihihi[i][j] = Rimhihihi[i][j];
               workRimlohihi[i][j] = Rimlohihi[i][j];
               workRimhilohi[i][j] = Rimhilohi[i][j];
               workRimlolohi[i][j] = Rimlolohi[i][j];
               workRimhihilo[i][j] = Rimhihilo[i][j];
               workRimlohilo[i][j] = Rimlohilo[i][j];
               workRimhilolo[i][j] = Rimhilolo[i][j];
               workRimlololo[i][j] = Rimlololo[i][j];
            }

         GPU_cmplx8_upper_tiled_solver
            (ncols,szt,nbt,
             workRrehihihi,workRrelohihi,workRrehilohi,workRrelolohi,
             workRrehihilo,workRrelohilo,workRrehilolo,workRrelololo,
             workRimhihihi,workRimlohihi,workRimhilohi,workRimlolohi,
             workRimhihilo,workRimlohilo,workRimhilolo,workRimlololo,
             brehihihi,brelohihi,brehilohi,brelolohi,
             brehihilo,brelohilo,brehilolo,brelololo,
             bimhihihi,bimlohihi,bimhilohi,bimlolohi,
             bimhihilo,bimlohilo,bimhilolo,bimlololo,
             xrehihihi,xrelohihi,xrehilohi,xrelolohi,
             xrehihilo,xrelohilo,xrehilolo,xrelololo,
             ximhihihi,ximlohihi,ximhilohi,ximlolohi,
             ximhihilo,ximlohilo,ximhilolo,ximlololo,
             &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
             &bsaddcnt,&bsmulcnt,&bsdivcnt);

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl8_bstimeflops
               (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,bsmulcnt,bsdivcnt);

        if(vrblvl > 1)
           for(int i=0; i<ncols; i++)
              cout << "x[" << i << "] : "
                   << xrehihihi[i] << "  " << xrelohihi[i] << endl << "  "
                   << xrehilohi[i] << "  " << xrelolohi[i] << endl << "  "
                   << xrehihilo[i] << "  " << xrelohilo[i] << endl << "  "
                   << xrehilolo[i] << "  " << xrelololo[i] << endl << "  "
                   << ximhihihi[i] << "  " << ximlohihi[i] << endl << "  "
                   << ximhilohi[i] << "  " << ximlolohi[i] << endl << "  "
                   << ximhihilo[i] << "  " << ximlohilo[i] << endl << "  "
                   << ximhilolo[i] << "  " << ximlololo[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solrehihihi[stage][j] = xrehihihi[j];
         solrelohihi[stage][j] = xrelohihi[j];
         solrehilohi[stage][j] = xrehilohi[j];
         solrelolohi[stage][j] = xrelolohi[j];
         solrehihilo[stage][j] = xrehihilo[j];
         solrelohilo[stage][j] = xrelohilo[j];
         solrehilolo[stage][j] = xrehilolo[j];
         solrelololo[stage][j] = xrelololo[j];
         solimhihihi[stage][j] = ximhihihi[j];
         solimlohihi[stage][j] = ximlohihi[j];
         solimhilohi[stage][j] = ximhilohi[j];
         solimlolohi[stage][j] = ximlolohi[j];
         solimhihilo[stage][j] = ximhihilo[j];
         solimlohilo[stage][j] = ximlohilo[j];
         solimhilolo[stage][j] = ximhilolo[j];
         solimlololo[stage][j] = ximlololo[j];
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
      free(workRrehihihi[i]); free(workRrelohihi[i]);
      free(workRrehilohi[i]); free(workRrelolohi[i]);
      free(workRrehihilo[i]); free(workRrelohilo[i]);
      free(workRrehilolo[i]); free(workRrelololo[i]);
      free(workRimhihihi[i]); free(workRimlohihi[i]);
      free(workRimhilohi[i]); free(workRimlolohi[i]);
      free(workRimhihilo[i]); free(workRimlohilo[i]);
      free(workRimhilolo[i]); free(workRimlololo[i]);
   }
   free(brehihihi); free(xrehihihi); free(workRrehihihi);
   free(brelohihi); free(xrelohihi); free(workRrelohihi);
   free(brehilohi); free(xrehilohi); free(workRrehilohi);
   free(brelolohi); free(xrelolohi); free(workRrelolohi);
   free(brehihilo); free(xrehihilo); free(workRrehihilo);
   free(brelohilo); free(xrelohilo); free(workRrelohilo);
   free(brehilolo); free(xrehilolo); free(workRrehilolo);
   free(brelololo); free(xrelololo); free(workRrelololo);
   free(bimhihihi); free(ximhihihi); free(workRimhihihi);
   free(bimlohihi); free(ximlohihi); free(workRimlohihi);
   free(bimhilohi); free(ximhilohi); free(workRimhilohi);
   free(bimlolohi); free(ximlolohi); free(workRimlolohi);
   free(bimhihilo); free(ximhihilo); free(workRimhihilo);
   free(bimlohilo); free(ximlohilo); free(workRimlohilo);
   free(bimhilolo); free(ximhilolo); free(workRimhilolo);
   free(bimlololo); free(ximlololo); free(workRimlololo);
}
