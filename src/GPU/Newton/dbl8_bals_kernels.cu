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

__global__ void cmplx8_bals_tail
 ( int ncols, int szt,
   double *Arehihihi, double *Arelohihi, double *Arehilohi, double *Arelolohi,
   double *Arehihilo, double *Arelohilo, double *Arehilolo, double *Arelololo,
   double *Aimhihihi, double *Aimlohihi, double *Aimhilohi, double *Aimlolohi,
   double *Aimhihilo, double *Aimlohilo, double *Aimhilolo, double *Aimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi, 
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo, 
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajrehihihi;             // register for Arehihihi[idx][j]
   double Ajrelohihi;             // register for Arelohihi[idx][j]
   double Ajrehilohi;             // register for Arehilohi[idx][j]
   double Ajrelolohi;             // register for Arelolohi[idx][j]
   double Ajrehihilo;             // register for Arehihilo[idx][j]
   double Ajrelohilo;             // register for Arelohilo[idx][j]
   double Ajrehilolo;             // register for Arehilolo[idx][j]
   double Ajrelololo;             // register for Arelololo[idx][j]
   double Ajimhihihi;             // register for Aimhihihi[idx][j]
   double Ajimlohihi;             // register for Aimlohihi[idx][j]
   double Ajimhilohi;             // register for Aimhilohi[idx][j]
   double Ajimlolohi;             // register for Aimlolohi[idx][j]
   double Ajimhihilo;             // register for Aimhihilo[idx][j]
   double Ajimlohilo;             // register for Aimlohilo[idx][j]
   double Ajimhilolo;             // register for Aimhilolo[idx][j]
   double Ajimlololo;             // register for Aimlololo[idx][j]
   double xjrehihihi;             // register for xrehihihi[j]
   double xjrelohihi;             // register for xrelohihi[j]
   double xjrehilohi;             // register for xrehilohi[j]
   double xjrelolohi;             // register for xrelolohi[j]
   double xjrehihilo;             // register for xrehihilo[j]
   double xjrelohilo;             // register for xrelohilo[j]
   double xjrehilolo;             // register for xrehilolo[j]
   double xjrelololo;             // register for xrelololo[j]
   double xjimhihihi;             // register for ximhihihi[j]
   double xjimlohihi;             // register for ximlohihi[j]
   double xjimhilohi;             // register for ximhilohi[j]
   double xjimlolohi;             // register for ximlolohi[j]
   double xjimhihilo;             // register for ximhihilo[j]
   double xjimlohilo;             // register for ximlohilo[j]
   double xjimhilolo;             // register for ximhilolo[j]
   double xjimlololo;             // register for ximlololo[j]
   double birehihihi = brehihihi[idx];  // register for brehihihi[idx]
   double birelohihi = brelohihi[idx];  // register for brelohihi[idx]
   double birehilohi = brehilohi[idx];  // register for brehilohi[idx]
   double birelolohi = brelolohi[idx];  // register for brelolohi[idx]
   double birehihilo = brehihilo[idx];  // register for brehihilo[idx]
   double birelohilo = brelohilo[idx];  // register for brelohilo[idx]
   double birehilolo = brehilolo[idx];  // register for brehilolo[idx]
   double birelololo = brelololo[idx];  // register for brelololo[idx]
   double biimhihihi = bimhihihi[idx];  // register for bimhihihi[idx]
   double biimlohihi = bimlohihi[idx];  // register for bimlohihi[idx]
   double biimhilohi = bimhilohi[idx];  // register for bimhilohi[idx]
   double biimlolohi = bimlolohi[idx];  // register for bimlolohi[idx]
   double biimhihilo = bimhihilo[idx];  // register for bimhihilo[idx]
   double biimlohilo = bimlohilo[idx];  // register for bimlohilo[idx]
   double biimhilolo = bimhilolo[idx];  // register for bimhilolo[idx]
   double biimlololo = bimlololo[idx];  // register for bimlololo[idx]
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajrehihihi = Arehihihi[offset+j];
      Ajrelohihi = Arelohihi[offset+j];
      Ajrehilohi = Arehilohi[offset+j];
      Ajrelolohi = Arelolohi[offset+j];
      Ajrehihilo = Arehihilo[offset+j];
      Ajrelohilo = Arelohilo[offset+j];
      Ajrehilolo = Arehilolo[offset+j];
      Ajrelololo = Arelololo[offset+j];
      Ajimhihihi = Aimhihihi[offset+j];
      Ajimlohihi = Aimlohihi[offset+j];
      Ajimhilohi = Aimhilohi[offset+j];
      Ajimlolohi = Aimlolohi[offset+j];
      Ajimhihilo = Aimhihilo[offset+j];
      Ajimlohilo = Aimlohilo[offset+j];
      Ajimhilolo = Aimhilolo[offset+j];
      Ajimlololo = Aimlololo[offset+j];
      xjrehihihi = xrehihihi[j];
      xjrelohihi = xrelohihi[j];
      xjrehilohi = xrehilohi[j];
      xjrelolohi = xrelolohi[j];
      xjrehihilo = xrehihilo[j];
      xjrelohilo = xrelohilo[j];
      xjrehilolo = xrehilolo[j];
      xjrelololo = xrelololo[j];
      xjimhihihi = ximhihihi[j];
      xjimlohihi = ximlohihi[j];
      xjimhilohi = ximhilohi[j];
      xjimlolohi = ximlolohi[j];
      xjimhihilo = ximhihilo[j];
      xjimlohilo = ximlohilo[j];
      xjimhilolo = ximhilolo[j];
      xjimlololo = ximlololo[j];
      // bi = bi - Aj*xj;
      // zre = Ajre*xjre - Ajim*xjim;
      // bire = bire - zre;
      odg_mul(Ajrehihihi,Ajrelohihi,Ajrehilohi,Ajrelolohi,
              Ajrehihilo,Ajrelohilo,Ajrehilolo,Ajrelololo,
              xjrehihihi,xjrelohihi,xjrehilohi,xjrelolohi,
              xjrehihilo,xjrelohilo,xjrehilolo,xjrelololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&birehihihi,&birelohihi,&birehilohi,&birelolohi,
              &birehihilo,&birelohilo,&birehilolo,&birelololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      odg_mul(Ajimhihihi,Ajimlohihi,Ajimhilohi,Ajimlolohi,
              Ajimhihilo,Ajimlohilo,Ajimhilolo,Ajimlololo,
              xjimhihihi,xjimlohihi,xjimhilohi,xjimlolohi,
              xjimhihilo,xjimlohilo,xjimhilolo,xjimlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&birehihihi,&birelohihi,&birehilohi,&birelolohi,
              &birehihilo,&birelohilo,&birehilolo,&birelololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      // zim = Ajre*xjim + Ajim*xjre;
      // biim = biim - zim;
      odg_mul(Ajrehihihi,Ajrelohihi,Ajrehilohi,Ajrelolohi,
              Ajrehihilo,Ajrelohilo,Ajrehilolo,Ajrelololo,
              xjimhihihi,xjimlohihi,xjimhilohi,xjimlolohi,
              xjimhihilo,xjimlohilo,xjimhilolo,xjimlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&biimhihihi,&biimlohihi,&biimhilohi,&biimlolohi,
              &biimhihilo,&biimlohilo,&biimhilolo,&biimlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      odg_mul(Ajimhihihi,Ajimlohihi,Ajimhilohi,Ajimlolohi,
              Ajimhihilo,Ajimlohilo,Ajimhilolo,Ajimlololo,
              xjrehihihi,xjrelohihi,xjrehilohi,xjrelolohi,
              xjrehihilo,xjrelohilo,xjrehilolo,xjrelololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&biimhihihi,&biimlohihi,&biimhilohi,&biimlolohi,
              &biimhihilo,&biimlohilo,&biimhilolo,&biimlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
   }
   brehihihi[idx] = birehihihi;
   brelohihi[idx] = birelohihi;
   brehilohi[idx] = birehilohi;
   brelolohi[idx] = birelolohi;
   brehihilo[idx] = birehihilo;
   brelohilo[idx] = birelohilo;
   brehilolo[idx] = birehilolo;
   brelololo[idx] = birelololo;
   bimhihihi[idx] = biimhihihi;
   bimlohihi[idx] = biimlohihi;
   bimhilohi[idx] = biimhilohi;
   bimlolohi[idx] = biimlolohi;
   bimhihilo[idx] = biimhihilo;
   bimlohilo[idx] = biimlohilo;
   bimhilolo[idx] = biimhilolo;
   bimlololo[idx] = biimlololo;
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

   if(verbose) cout << "-> GPU multiplies rhs with Q^T ..." << endl;

   GPU_dbl8_bals_qtb
      (ncols,szt,nbt,
       Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
       bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
       verbose);

   if(verbose)
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

   if(verbose) cout << "-> GPU multiplies rhs with Q^H ..." << endl;

   GPU_cmplx8_bals_qhb
      (ncols,szt,nbt,
       Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
       Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
       Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
       Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
       brehihihi,brelohihi,brehilohi,brelolohi,
       brehihilo,brelohilo,brehilolo,brelololo,
       bimhihihi,bimlohihi,bimhilohi,bimlolohi,
       bimhihilo,bimlohilo,bimhilolo,bimlololo,verbose);

   if(verbose)
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

   if(verbose)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(verbose)
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

   if(verbose)
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

void GPU_cmplx8_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
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
   double **solimhilolo, double **solimlololo, bool verbose )
{
   if(verbose)
   {
      cout << "GPU_cmplx8_bals_tail input blocks of rhs :" << endl;
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
   const size_t szrhs = nrows*sizeof(double);
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

   double *xrehihihi_d;
   double *xrelohihi_d;
   double *xrehilohi_d;
   double *xrelolohi_d;
   double *xrehihilo_d;
   double *xrelohilo_d;
   double *xrehilolo_d;
   double *xrelololo_d;
   double *ximhihihi_d;
   double *ximlohihi_d;
   double *ximhilohi_d;
   double *ximlolohi_d;
   double *ximhihilo_d;
   double *ximlohilo_d;
   double *ximhilolo_d;
   double *ximlololo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xrehihihi_d,szsol);
   cudaMalloc((void**)&xrelohihi_d,szsol);
   cudaMalloc((void**)&xrehilohi_d,szsol);
   cudaMalloc((void**)&xrelolohi_d,szsol);
   cudaMalloc((void**)&xrehihilo_d,szsol);
   cudaMalloc((void**)&xrelohilo_d,szsol);
   cudaMalloc((void**)&xrehilolo_d,szsol);
   cudaMalloc((void**)&xrelololo_d,szsol);
   cudaMalloc((void**)&ximhihihi_d,szsol);
   cudaMalloc((void**)&ximlohihi_d,szsol);
   cudaMalloc((void**)&ximhilohi_d,szsol);
   cudaMalloc((void**)&ximlolohi_d,szsol);
   cudaMalloc((void**)&ximhihilo_d,szsol);
   cudaMalloc((void**)&ximlohilo_d,szsol);
   cudaMalloc((void**)&ximhilolo_d,szsol);
   cudaMalloc((void**)&ximlololo_d,szsol);
   cudaMemcpy(xrehihihi_d,solrehihihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelohihi_d,solrelohihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehilohi_d,solrehilohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelolohi_d,solrelolohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehihilo_d,solrehihilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelohilo_d,solrelohilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehilolo_d,solrehilolo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelololo_d,solrelololo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhihihi_d,solimhihihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlohihi_d,solimlohihi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhilohi_d,solimhilohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlolohi_d,solimlolohi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhihilo_d,solimhihilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlohilo_d,solimlohilo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhilolo_d,solimhilolo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlololo_d,solimlololo[stage-1],szsol,cudaMemcpyHostToDevice);

   double *Arehihihi_d;
   double *Arelohihi_d;
   double *Arehilohi_d;
   double *Arelolohi_d;
   double *Arehihilo_d;
   double *Arelohilo_d;
   double *Arehilolo_d;
   double *Arelololo_d;
   double *Aimhihihi_d;
   double *Aimlohihi_d;
   double *Aimhilohi_d;
   double *Aimlolohi_d;
   double *Aimhihilo_d;
   double *Aimlohilo_d;
   double *Aimhilolo_d;
   double *Aimlololo_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Arehihihi_d,szmat);
   cudaMalloc((void**)&Arelohihi_d,szmat);
   cudaMalloc((void**)&Arehilohi_d,szmat);
   cudaMalloc((void**)&Arelolohi_d,szmat);
   cudaMalloc((void**)&Arehihilo_d,szmat);
   cudaMalloc((void**)&Arelohilo_d,szmat);
   cudaMalloc((void**)&Arehilolo_d,szmat);
   cudaMalloc((void**)&Arelololo_d,szmat);
   cudaMalloc((void**)&Aimhihihi_d,szmat);
   cudaMalloc((void**)&Aimlohihi_d,szmat);
   cudaMalloc((void**)&Aimhilohi_d,szmat);
   cudaMalloc((void**)&Aimlolohi_d,szmat);
   cudaMalloc((void**)&Aimhihilo_d,szmat);
   cudaMalloc((void**)&Aimlohilo_d,szmat);
   cudaMalloc((void**)&Aimhilolo_d,szmat);
   cudaMalloc((void**)&Aimlololo_d,szmat);

   double *Arehihihi_h = new double[szmat];
   double *Arelohihi_h = new double[szmat];
   double *Arehilohi_h = new double[szmat];
   double *Arelolohi_h = new double[szmat];
   double *Arehihilo_h = new double[szmat];
   double *Arelohilo_h = new double[szmat];
   double *Arehilolo_h = new double[szmat];
   double *Arelololo_h = new double[szmat];
   double *Aimhihihi_h = new double[szmat];
   double *Aimlohihi_h = new double[szmat];
   double *Aimhilohi_h = new double[szmat];
   double *Aimlolohi_h = new double[szmat];
   double *Aimhihilo_h = new double[szmat];
   double *Aimlohilo_h = new double[szmat];
   double *Aimhilolo_h = new double[szmat];
   double *Aimlololo_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(verbose)
         cout << "GPU_cmplx8_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Arehihihi_h[idx]   = matrehihihi[k-stage+1][i][j];
            Arelohihi_h[idx]   = matrelohihi[k-stage+1][i][j];
            Arehilohi_h[idx]   = matrehilohi[k-stage+1][i][j];
            Arelolohi_h[idx]   = matrelolohi[k-stage+1][i][j];
            Arehihilo_h[idx]   = matrehihilo[k-stage+1][i][j];
            Arelohilo_h[idx]   = matrelohilo[k-stage+1][i][j];
            Arehilolo_h[idx]   = matrehilolo[k-stage+1][i][j];
            Arelololo_h[idx]   = matrelololo[k-stage+1][i][j];
            Aimhihihi_h[idx]   = matimhihihi[k-stage+1][i][j];
            Aimlohihi_h[idx]   = matimlohihi[k-stage+1][i][j];
            Aimhilohi_h[idx]   = matimhilohi[k-stage+1][i][j];
            Aimlolohi_h[idx]   = matimlolohi[k-stage+1][i][j];
            Aimhihilo_h[idx]   = matimhihilo[k-stage+1][i][j];
            Aimlohilo_h[idx]   = matimlohilo[k-stage+1][i][j];
            Aimhilolo_h[idx]   = matimhilolo[k-stage+1][i][j];
            Aimlololo_h[idx++] = matimlololo[k-stage+1][i][j];
         }
      
      cudaMemcpy(brehihihi_d,rhsrehihihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelohihi_d,rhsrelohihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brehilohi_d,rhsrehilohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelolohi_d,rhsrelolohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brehihilo_d,rhsrehihilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelohilo_d,rhsrelohilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brehilolo_d,rhsrehilolo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelololo_d,rhsrelololo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhihihi_d,rhsimhihihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlohihi_d,rhsimlohihi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhilohi_d,rhsimhilohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlolohi_d,rhsimlolohi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhihilo_d,rhsimhihilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlohilo_d,rhsimlohilo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhilolo_d,rhsimhilolo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlololo_d,rhsimlololo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehihihi_d,Arehihihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelohihi_d,Arelohihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehilohi_d,Arehilohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelolohi_d,Arelolohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehihilo_d,Arehihilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelohilo_d,Arelohilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehilolo_d,Arehilolo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelololo_d,Arelololo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhihihi_d,Aimhihihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlohihi_d,Aimlohihi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhilohi_d,Aimhilohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlolohi_d,Aimlolohi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhihilo_d,Aimhihilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlohilo_d,Aimlohilo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhilolo_d,Aimhilolo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlololo_d,Aimlololo_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cmplx8_bals_tail<<<nbt,szt>>>
         (ncols,szt,Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
                    Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
                    Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
                    Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
          xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
          xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
          ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
          ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d,
          brehihihi_d,brelohihi_d,brehilohi_d,brelolohi_d,
          brehihilo_d,brelohilo_d,brehilolo_d,brelololo_d,
          bimhihihi_d,bimlohihi_d,bimhilohi_d,bimlolohi_d,
          bimhihilo_d,bimlohilo_d,bimhilolo_d,bimlololo_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhsrehihihi[k],brehihihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelohihi[k],brelohihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrehilohi[k],brehilohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelolohi[k],brelolohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrehihilo[k],brehihilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelohilo[k],brelohilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrehilolo[k],brehilolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelololo[k],brelololo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhihihi[k],bimhihihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlohihi[k],bimlohihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhilohi[k],bimhilohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlolohi[k],bimlolohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhihilo[k],bimhihilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlohilo[k],bimlohilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhilolo[k],bimhilolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlololo[k],bimlololo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Arehihihi_h); free(Aimhihihi_h);
   free(Arelohihi_h); free(Aimlohihi_h);
   free(Arehilohi_h); free(Aimhilohi_h);
   free(Arelolohi_h); free(Aimlolohi_h);
   free(Arehihilo_h); free(Aimhihilo_h);
   free(Arelohilo_h); free(Aimlohilo_h);
   free(Arehilolo_h); free(Aimhilolo_h);
   free(Arelololo_h); free(Aimlololo_h);

   if(verbose)
   {
      cout << "GPU_cmplx8_bals_tail copied blocks of rhs :" << endl;
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
   bool verbose )
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
   const bool bvrb = (vrblvl > 1);

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
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;

         bhihihi[i] = rhshihihi[stage][i];
         blohihi[i] = rhslohihi[stage][i];
         bhilohi[i] = rhshilohi[stage][i];
         blolohi[i] = rhslolohi[stage][i];
         bhihilo[i] = rhshihilo[stage][i];
         blohilo[i] = rhslohilo[stage][i];
         bhilolo[i] = rhshilolo[stage][i];
         blololo[i] = rhslololo[stage][i];

         // cout << "b[" << i << "] : "
         //      << bhihihi[i] << "  " << blohihi[i] << endl << "  "
         //      << bhilohi[i] << "  " << blolohi[i] << endl << "  "
         //      << bhihilo[i] << "  " << blohilo[i] << endl << "  "
         //      << bhilolo[i] << "  " << blololo[i] << endl;
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

     if(vrblvl > 1)
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

void GPU_cmplx8_bals_solve
 ( int dim, int degp1, int szt, int nbt,
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
   double **solimhilolo, double **solimlololo, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   const bool bvrb = (vrblvl > 1);

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
       ximhihilo,ximlohilo,ximhilolo,ximlololo,bvrb);

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
   for(int stage=1; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

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
          solimhihilo,solimlohilo,solimhilolo,solimlololo,bvrb);

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
          bimhihilo,bimlohilo,bimhilolo,bimlololo,bvrb);

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
   for(int i=0; i<nrows; i++)
   {
      free(Arehihihi[i]); free(workRrehihihi[i]);
      free(Arelohihi[i]); free(workRrelohihi[i]);
      free(Arehilohi[i]); free(workRrehilohi[i]);
      free(Arelolohi[i]); free(workRrelolohi[i]);
      free(Arehihilo[i]); free(workRrehihilo[i]);
      free(Arelohilo[i]); free(workRrelohilo[i]);
      free(Arehilolo[i]); free(workRrehilolo[i]);
      free(Arelololo[i]); free(workRrelololo[i]);
      free(Aimhihihi[i]); free(workRimhihihi[i]);
      free(Aimlohihi[i]); free(workRimlohihi[i]);
      free(Aimhilohi[i]); free(workRimhilohi[i]);
      free(Aimlolohi[i]); free(workRimlolohi[i]);
      free(Aimhihilo[i]); free(workRimhihilo[i]);
      free(Aimlohilo[i]); free(workRimlohilo[i]);
      free(Aimhilolo[i]); free(workRimhilolo[i]);
      free(Aimlololo[i]); free(workRimlololo[i]);
   }
   free(Arehihihi); free(brehihihi); free(xrehihihi); free(workRrehihihi);
   free(Arelohihi); free(brelohihi); free(xrelohihi); free(workRrelohihi);
   free(Arehilohi); free(brehilohi); free(xrehilohi); free(workRrehilohi);
   free(Arelolohi); free(brelolohi); free(xrelolohi); free(workRrelolohi);
   free(Arehihilo); free(brehihilo); free(xrehihilo); free(workRrehihilo);
   free(Arelohilo); free(brelohilo); free(xrelohilo); free(workRrelohilo);
   free(Arehilolo); free(brehilolo); free(xrehilolo); free(workRrehilolo);
   free(Arelololo); free(brelololo); free(xrelololo); free(workRrelololo);
   free(Aimhihihi); free(bimhihihi); free(ximhihihi); free(workRimhihihi);
   free(Aimlohihi); free(bimlohihi); free(ximlohihi); free(workRimlohihi);
   free(Aimhilohi); free(bimhilohi); free(ximhilohi); free(workRimhilohi);
   free(Aimlolohi); free(bimlolohi); free(ximlolohi); free(workRimlolohi);
   free(Aimhihilo); free(bimhihilo); free(ximhihilo); free(workRimhihilo);
   free(Aimlohilo); free(bimlohilo); free(ximlohilo); free(workRimlohilo);
   free(Aimhilolo); free(bimhilolo); free(ximhilolo); free(workRimhilolo);
   free(Aimlololo); free(bimlololo); free(ximlololo); free(workRimlololo);
}
