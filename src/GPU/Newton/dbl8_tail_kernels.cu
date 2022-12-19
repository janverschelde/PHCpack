// The file dbl8_tail_kernels.cu defines the functions with prototypes in
// the file dbl8_tail_kernels.h.

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
#include "dbl8_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

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

void write_dbl8_balsflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for b = b - A*x : " << lapsms
        << " milliseconds." << endl;

   long long int flopcnt;
   if(ctype == 0)
      flopcnt = 270*ncols*ncols + 1742*ncols*ncols;
      // as many + as * in one inner product
   else
      flopcnt = 4*270*ncols*ncols + 4*1742*ncols*ncols;
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
   double **solhilolo, double **sollololo,
   double *totupdlapsedms, int vrblvl )
{
   if(vrblvl > 1)
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
      if(vrblvl > 1)
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

      if(vrblvl > 1)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cudaEvent_t start,stop;       // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl8_bals_tail<<<nbt,szt>>>
          (ncols,szt,
           Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
           Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
           xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
           xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d,
           bhihihi_d,blohihi_d,bhilohi_d,blolohi_d,
           bhihilo_d,blohilo_d,bhilolo_d,blololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);

      *totupdlapsedms += milliseconds;

      if(vrblvl > 0) write_dbl8_balsflops(0,ncols,milliseconds);
      
      if(vrblvl > 1)
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

   if(vrblvl > 1)
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
   cudaFree(bhihihi_d); cudaFree(blohihi_d);
   cudaFree(bhilohi_d); cudaFree(blolohi_d);
   cudaFree(bhihilo_d); cudaFree(blohilo_d);
   cudaFree(bhilolo_d); cudaFree(blololo_d);
   cudaFree(xhihihi_d); cudaFree(xlohihi_d);
   cudaFree(xhilohi_d); cudaFree(xlolohi_d);
   cudaFree(xhihilo_d); cudaFree(xlohilo_d);
   cudaFree(xhilolo_d); cudaFree(xlololo_d);
   cudaFree(Ahihihi_d); cudaFree(Alohihi_d);
   cudaFree(Ahilohi_d); cudaFree(Alolohi_d);
   cudaFree(Ahihilo_d); cudaFree(Alohilo_d);
   cudaFree(Ahilolo_d); cudaFree(Alololo_d);
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
   double **solimhilolo, double **solimlololo,
   double *totupdlapsedms, int vrblvl )
{
   if(vrblvl > 1)
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
      if(vrblvl > 1)
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

      if(vrblvl > 1)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cudaEvent_t start,stop;       // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
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
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);

      *totupdlapsedms += milliseconds;

      if(vrblvl > 0) write_dbl8_balsflops(1,ncols,milliseconds);
      
      if(vrblvl > 1)
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

   if(vrblvl > 1)
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
   cudaFree(brehihihi_d); cudaFree(brelohihi_d);
   cudaFree(brehilohi_d); cudaFree(brelolohi_d);
   cudaFree(brehihilo_d); cudaFree(brelohilo_d);
   cudaFree(brehilolo_d); cudaFree(brelololo_d);
   cudaFree(bimhihihi_d); cudaFree(bimlohihi_d);
   cudaFree(bimhilohi_d); cudaFree(bimlolohi_d);
   cudaFree(bimhihilo_d); cudaFree(bimlohilo_d);
   cudaFree(bimhilolo_d); cudaFree(bimlololo_d);
   cudaFree(xrehihihi_d); cudaFree(xrelohihi_d);
   cudaFree(xrehilohi_d); cudaFree(xrelolohi_d);
   cudaFree(xrehihilo_d); cudaFree(xrelohilo_d);
   cudaFree(xrehilolo_d); cudaFree(xrelololo_d);
   cudaFree(ximhihihi_d); cudaFree(ximlohihi_d);
   cudaFree(ximhilohi_d); cudaFree(ximlolohi_d);
   cudaFree(ximhihilo_d); cudaFree(ximlohilo_d);
   cudaFree(ximhilolo_d); cudaFree(ximlololo_d);
   cudaFree(Arehihihi_d); cudaFree(Arelohihi_d);
   cudaFree(Arehilohi_d); cudaFree(Arelolohi_d);
   cudaFree(Arehihilo_d); cudaFree(Arelohilo_d);
   cudaFree(Arehilolo_d); cudaFree(Arelololo_d);
   cudaFree(Aimhihihi_d); cudaFree(Aimlohihi_d);
   cudaFree(Aimhilohi_d); cudaFree(Aimlolohi_d);
   cudaFree(Aimhihilo_d); cudaFree(Aimlohilo_d);
   cudaFree(Aimhilolo_d); cudaFree(Aimlololo_d);
}

void GPU_dbl8_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
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
   double **solhilolo, double **sollololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo,
   double *totreslapsedms, long long int *add, long long int *mul,
   int vrblvl )
{
   double *rhihihi_d;
   double *rlohihi_d;
   double *rhilohi_d;
   double *rlolohi_d;
   double *rhihilo_d;
   double *rlohilo_d;
   double *rhilolo_d;
   double *rlololo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhihihi_d,szrhs);
   cudaMalloc((void**)&rlohihi_d,szrhs);
   cudaMalloc((void**)&rhilohi_d,szrhs);
   cudaMalloc((void**)&rlolohi_d,szrhs);
   cudaMalloc((void**)&rhihilo_d,szrhs);
   cudaMalloc((void**)&rlohilo_d,szrhs);
   cudaMalloc((void**)&rhilolo_d,szrhs);
   cudaMalloc((void**)&rlololo_d,szrhs);

   double *xhihihi_d;
   double *xlohihi_d;
   double *xhilohi_d;
   double *xlolohi_d;
   double *xhihilo_d;
   double *xlohilo_d;
   double *xhilolo_d;
   double *xlololo_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&xhihihi_d,szsol);
   cudaMalloc((void**)&xlohihi_d,szsol);
   cudaMalloc((void**)&xhilohi_d,szsol);
   cudaMalloc((void**)&xlolohi_d,szsol);
   cudaMalloc((void**)&xhihilo_d,szsol);
   cudaMalloc((void**)&xlohilo_d,szsol);
   cudaMalloc((void**)&xhilolo_d,szsol);
   cudaMalloc((void**)&xlololo_d,szsol);

   double *Ahihihi_d;
   double *Alohihi_d;
   double *Ahilohi_d;
   double *Alolohi_d;
   double *Ahihilo_d;
   double *Alohilo_d;
   double *Ahilolo_d;
   double *Alololo_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&Ahihihi_d,szmat);
   cudaMalloc((void**)&Alohihi_d,szmat);
   cudaMalloc((void**)&Ahilohi_d,szmat);
   cudaMalloc((void**)&Alolohi_d,szmat);
   cudaMalloc((void**)&Ahihilo_d,szmat);
   cudaMalloc((void**)&Alohilo_d,szmat);
   cudaMalloc((void**)&Ahilolo_d,szmat);
   cudaMalloc((void**)&Alololo_d,szmat);

   double *Ahihihi_h = new double[dim*dim];
   double *Alohihi_h = new double[dim*dim];
   double *Ahilohi_h = new double[dim*dim];
   double *Alolohi_h = new double[dim*dim];
   double *Ahihilo_h = new double[dim*dim];
   double *Alohilo_h = new double[dim*dim];
   double *Ahilolo_h = new double[dim*dim];
   double *Alololo_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   if(vrblvl > 0)
      cout << "GPU_dbl8_linear_residue for deg+1 : " << degp1 << endl;

   for(int i=tailidx; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rhihihi_d,rhshihihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlohihi_d,rhslohihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rhilohi_d,rhshilohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlolohi_d,rhslolohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rhihilo_d,rhshihilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlohilo_d,rhslohilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rhilolo_d,rhshilolo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlololo_d,rhslololo[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=(i-tailidx); j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Ahihihi_h[idx]   = mathihihi[j][i1][j1];
               Alohihi_h[idx]   = matlohihi[j][i1][j1];
               Ahilohi_h[idx]   = mathilohi[j][i1][j1];
               Alolohi_h[idx]   = matlolohi[j][i1][j1];
               Ahihilo_h[idx]   = mathihilo[j][i1][j1];
               Alohilo_h[idx]   = matlohilo[j][i1][j1];
               Ahilolo_h[idx]   = mathilolo[j][i1][j1];
               Alololo_h[idx++] = matlololo[j][i1][j1];
            }
      
         cudaMemcpy(Ahihihi_d,Ahihihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alohihi_d,Alohihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahilohi_d,Ahilohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alolohi_d,Alolohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahihilo_d,Ahihilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alohilo_d,Alohilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahilolo_d,Ahilolo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alololo_d,Alololo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(xhihihi_d,solhihihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlohihi_d,sollohihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xhilohi_d,solhilohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlolohi_d,sollolohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xhihilo_d,solhihilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlohilo_d,sollohilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xhilolo_d,solhilolo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlololo_d,sollololo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahihihi_d,Ahihihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alohihi_d,Alohihi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahilohi_d,Ahilohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alolohi_d,Alolohi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahihilo_d,Ahihilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alohilo_d,Alohilo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahilolo_d,Ahilolo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alololo_d,Alololo_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 1)
            cout << "GPU_dbl8_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         dbl8_bals_tail<<<nbt,szt>>>
            (dim,szt,Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
                     Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
                     xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
                     xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d,
                     rhihihi_d,rlohihi_d,rhilohi_d,rlolohi_d,
                     rhihilo_d,rlohilo_d,rhilolo_d,rlololo_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *totreslapsedms += milliseconds;
         flopcount_dbl_bals_tail(dim,add,mul);

         if(vrblvl > 0) write_dbl8_balsflops(0,dim,milliseconds);
      }
      cudaMemcpy(resvechihihi[i],rhihihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclohihi[i],rlohihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvechilohi[i],rhilohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclolohi[i],rlolohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvechihilo[i],rhihilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclohilo[i],rlohilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvechilolo[i],rhilolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclololo[i],rlololo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=tailidx; i<degp1; i++)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << solhihihi[i][j] << "  " << sollohihi[i][j] << endl;
            cout << solhilohi[i][j] << "  " << sollolohi[i][j] << endl;
            cout << solhihilo[i][j] << "  " << sollohilo[i][j] << endl;
            cout << solhilolo[i][j] << "  " << sollololo[i][j] << endl;
         }
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << resvechihihi[i][j] << "  " << resveclohihi[i][j] << endl;
            cout << resvechilohi[i][j] << "  " << resveclolohi[i][j] << endl;
            cout << resvechihilo[i][j] << "  " << resveclohilo[i][j] << endl;
            cout << resvechilolo[i][j] << "  " << resveclololo[i][j] << endl;
         }
      }
   }
   *resmaxhihihi = 0.0; *resmaxlohihi = 0.0;
   *resmaxhilohi = 0.0; *resmaxlolohi = 0.0;
   *resmaxhihilo = 0.0; *resmaxlohilo = 0.0;
   *resmaxhilolo = 0.0; *resmaxlololo = 0.0;
   
   for(int i=tailidx; i<degp1; i++)
   {
      double *rihihihi = resvechihihi[i];
      double *rilohihi = resveclohihi[i];
      double *rihilohi = resvechilohi[i];
      double *rilolohi = resveclolohi[i];
      double *rihihilo = resvechihilo[i];
      double *rilohilo = resveclohilo[i];
      double *rihilolo = resvechilolo[i];
      double *rilololo = resveclololo[i];

      for(int j=0; j<dim; j++)
         if(abs(rihihihi[j]) > *resmaxhihihi)
         {
            *resmaxhihihi = abs(rihihihi[j]);
            *resmaxlohihi = abs(rilohihi[j]);
            *resmaxhilohi = abs(rihilohi[j]);
            *resmaxlolohi = abs(rilolohi[j]);
            *resmaxhihilo = abs(rihihilo[j]);
            *resmaxlohilo = abs(rilohilo[j]);
            *resmaxhilolo = abs(rihilolo[j]);
            *resmaxlololo = abs(rilololo[j]);
         }
   }
   free(Ahihihi_h); free(Alohihi_h); free(Ahilohi_h); free(Alolohi_h);
   free(Ahihilo_h); free(Alohilo_h); free(Ahilolo_h); free(Alololo_h);

   cudaFree(rhihihi_d); cudaFree(rlohihi_d);
   cudaFree(rhilohi_d); cudaFree(rlolohi_d);
   cudaFree(rhihilo_d); cudaFree(rlohilo_d);
   cudaFree(rhilolo_d); cudaFree(rlololo_d);
   cudaFree(xhihihi_d); cudaFree(xlohihi_d);
   cudaFree(xhilohi_d); cudaFree(xlolohi_d);
   cudaFree(xhihilo_d); cudaFree(xlohilo_d);
   cudaFree(xhilolo_d); cudaFree(xlololo_d);
   cudaFree(Ahihihi_d); cudaFree(Alohihi_d);
   cudaFree(Ahilohi_d); cudaFree(Alolohi_d);
   cudaFree(Ahihilo_d); cudaFree(Alohilo_d);
   cudaFree(Ahilolo_d); cudaFree(Alololo_d);
}

void GPU_cmplx8_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
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
   double **solimhilolo, double **solimlololo,
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo,
   double *totreslapsedms, long long int *add, long long int *mul,
   int vrblvl )
{
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
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rrehihihi_d,szrhs);
   cudaMalloc((void**)&rrelohihi_d,szrhs);
   cudaMalloc((void**)&rrehilohi_d,szrhs);
   cudaMalloc((void**)&rrelolohi_d,szrhs);
   cudaMalloc((void**)&rrehihilo_d,szrhs);
   cudaMalloc((void**)&rrelohilo_d,szrhs);
   cudaMalloc((void**)&rrehilolo_d,szrhs);
   cudaMalloc((void**)&rrelololo_d,szrhs);
   cudaMalloc((void**)&rimhihihi_d,szrhs);
   cudaMalloc((void**)&rimlohihi_d,szrhs);
   cudaMalloc((void**)&rimhilohi_d,szrhs);
   cudaMalloc((void**)&rimlolohi_d,szrhs);
   cudaMalloc((void**)&rimhihilo_d,szrhs);
   cudaMalloc((void**)&rimlohilo_d,szrhs);
   cudaMalloc((void**)&rimhilolo_d,szrhs);
   cudaMalloc((void**)&rimlololo_d,szrhs);

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
   const size_t szsol = dim*sizeof(double);
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
   const size_t szmat = dim*dim*sizeof(double);
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

   double *Arehihihi_h = new double[dim*dim];
   double *Arelohihi_h = new double[dim*dim];
   double *Arehilohi_h = new double[dim*dim];
   double *Arelolohi_h = new double[dim*dim];
   double *Arehihilo_h = new double[dim*dim];
   double *Arelohilo_h = new double[dim*dim];
   double *Arehilolo_h = new double[dim*dim];
   double *Arelololo_h = new double[dim*dim];
   double *Aimhihihi_h = new double[dim*dim];
   double *Aimlohihi_h = new double[dim*dim];
   double *Aimhilohi_h = new double[dim*dim];
   double *Aimlolohi_h = new double[dim*dim];
   double *Aimhihilo_h = new double[dim*dim];
   double *Aimlohilo_h = new double[dim*dim];
   double *Aimhilolo_h = new double[dim*dim];
   double *Aimlololo_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   if(vrblvl > 0)
      cout << "GPU_cmplx8_linear_residue for deg+1 : " << degp1 << endl;

   for(int i=tailidx; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rrehihihi_d,rhsrehihihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelohihi_d,rhsrelohihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrehilohi_d,rhsrehilohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelolohi_d,rhsrelolohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrehihilo_d,rhsrehihilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelohilo_d,rhsrelohilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrehilolo_d,rhsrehilolo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelololo_d,rhsrelololo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhihihi_d,rhsimhihihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlohihi_d,rhsimlohihi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhilohi_d,rhsimhilohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlolohi_d,rhsimlolohi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhihilo_d,rhsimhihilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlohilo_d,rhsimlohilo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhilolo_d,rhsimhilolo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlololo_d,rhsimlololo[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=(i-tailidx); j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Arehihihi_h[idx]   = matrehihihi[j][i1][j1];
               Arelohihi_h[idx]   = matrelohihi[j][i1][j1];
               Arehilohi_h[idx]   = matrehilohi[j][i1][j1];
               Arelolohi_h[idx]   = matrelolohi[j][i1][j1];
               Arehihilo_h[idx]   = matrehihilo[j][i1][j1];
               Arelohilo_h[idx]   = matrelohilo[j][i1][j1];
               Arehilolo_h[idx]   = matrehilolo[j][i1][j1];
               Arelololo_h[idx]   = matrelololo[j][i1][j1];
               Aimhihihi_h[idx]   = matimhihihi[j][i1][j1];
               Aimlohihi_h[idx]   = matimlohihi[j][i1][j1];
               Aimhilohi_h[idx]   = matimhilohi[j][i1][j1];
               Aimlolohi_h[idx]   = matimlolohi[j][i1][j1];
               Aimhihilo_h[idx]   = matimhihilo[j][i1][j1];
               Aimlohilo_h[idx]   = matimlohilo[j][i1][j1];
               Aimhilolo_h[idx]   = matimhilolo[j][i1][j1];
               Aimlololo_h[idx++] = matimlololo[j][i1][j1];
            }
      
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
         cudaMemcpy(xrehihihi_d,solrehihihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelohihi_d,solrelohihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrehilohi_d,solrehilohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelolohi_d,solrelolohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrehihilo_d,solrehihilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelohilo_d,solrelohilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrehilolo_d,solrehilolo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelololo_d,solrelololo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhihihi_d,solimhihihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlohihi_d,solimlohihi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhilohi_d,solimhilohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlolohi_d,solimlolohi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhihilo_d,solimhihilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlohilo_d,solimlohilo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhilolo_d,solimhilolo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlololo_d,solimlololo[i-j],szsol,cudaMemcpyHostToDevice);
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

         if(vrblvl > 1)
            cout << "GPU_cmplx8_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         cmplx8_bals_tail<<<nbt,szt>>>
            (dim,szt,Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
                     Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
                     Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
                     Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
                     xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
                     xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
                     ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
                     ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d,
                     rrehihihi_d,rrelohihi_d,rrehilohi_d,rrelolohi_d,
                     rrehihilo_d,rrelohilo_d,rrehilolo_d,rrelololo_d,
                     rimhihihi_d,rimlohihi_d,rimhilohi_d,rimlolohi_d,
                     rimhihilo_d,rimlohilo_d,rimhilolo_d,rimlololo_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *totreslapsedms += milliseconds;
         flopcount_cmplx_bals_tail(dim,add,mul);

         if(vrblvl > 0) write_dbl8_balsflops(1,dim,milliseconds);
      }
      cudaMemcpy(resvecrehihihi[i],rrehihihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelohihi[i],rrelohihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrehilohi[i],rrehilohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelolohi[i],rrelolohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrehihilo[i],rrehihilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelohilo[i],rrelohilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrehilolo[i],rrehilolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelololo[i],rrelololo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhihihi[i],rimhihihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlohihi[i],rimlohihi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhilohi[i],rimhilohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlolohi[i],rimlolohi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhihilo[i],rimhihilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlohilo[i],rimlohilo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhilolo[i],rimhilolo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlololo[i],rimlololo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=tailidx; i<degp1; i++)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solrehihihi[i][j] << "  " << solrelohihi[i][j] << endl
                 << "  "
                 << solrehilohi[i][j] << "  " << solrelolohi[i][j] << endl
                 << "  "
                 << solrehihilo[i][j] << "  " << solrelohilo[i][j] << endl
                 << "  "
                 << solrehilolo[i][j] << "  " << solrelololo[i][j] << endl
                 << "  "
                 << solimhihihi[i][j] << "  " << solimlohihi[i][j] << endl
                 << "  "
                 << solimhilohi[i][j] << "  " << solimlolohi[i][j] << endl
                 << "  "
                 << solimhihilo[i][j] << "  " << solimlohilo[i][j] << endl
                 << "  "
                 << solimhilolo[i][j] << "  " << solimlololo[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << resvecrehihihi[i][j] << "  "
                 << resvecrelohihi[i][j] << endl << "  "
                 << resvecrehilohi[i][j] << "  "
                 << resvecrelolohi[i][j] << endl << "  "
                 << resvecrehihilo[i][j] << "  "
                 << resvecrelohilo[i][j] << endl << "  "
                 << resvecrehilolo[i][j] << "  "
                 << resvecrelololo[i][j] << endl << "  "
                 << resvecimhihihi[i][j] << "  "
                 << resvecimlohihi[i][j] << endl << "  "
                 << resvecimhilohi[i][j] << "  "
                 << resvecimlolohi[i][j] << endl << "  "
                 << resvecimhihilo[i][j] << "  "
                 << resvecimlohilo[i][j] << endl << "  "
                 << resvecimhilolo[i][j] << "  "
                 << resvecimlololo[i][j] << endl;
      }
   }
   *resmaxhihihi = 0.0; *resmaxlohihi = 0.0;
   *resmaxhilohi = 0.0; *resmaxlolohi = 0.0;
   *resmaxhihilo = 0.0; *resmaxlohilo = 0.0;
   *resmaxhilolo = 0.0; *resmaxlololo = 0.0;

   for(int i=tailidx; i<degp1; i++)
   {
      double *rirehihihi = resvecrehihihi[i];
      double *rirelohihi = resvecrelohihi[i];
      double *rirehilohi = resvecrehilohi[i];
      double *rirelolohi = resvecrelolohi[i];
      double *rirehihilo = resvecrehihilo[i];
      double *rirelohilo = resvecrelohilo[i];
      double *rirehilolo = resvecrehilolo[i];
      double *rirelololo = resvecrelololo[i];
      double *riimhihihi = resvecimhihihi[i];
      double *riimlohihi = resvecimlohihi[i];
      double *riimhilohi = resvecimhilohi[i];
      double *riimlolohi = resvecimlolohi[i];
      double *riimhihilo = resvecimhihilo[i];
      double *riimlohilo = resvecimlohilo[i];
      double *riimhilolo = resvecimhilolo[i];
      double *riimlololo = resvecimlololo[i];

      for(int j=0; j<dim; j++)
         if(abs(rirehihihi[j]) + abs(riimhihihi[j]) > *resmaxhihihi)
         {
            *resmaxhihihi = abs(rirehihihi[j]) + abs(riimhihihi[j]);
            *resmaxlohihi = abs(rirelohihi[j]) + abs(riimlohihi[j]);
            *resmaxhilohi = abs(rirehilohi[j]) + abs(riimhilohi[j]);
            *resmaxlolohi = abs(rirelolohi[j]) + abs(riimlolohi[j]);
            *resmaxhihilo = abs(rirehihilo[j]) + abs(riimhihilo[j]);
            *resmaxlohilo = abs(rirelohilo[j]) + abs(riimlohilo[j]);
            *resmaxhilolo = abs(rirehilolo[j]) + abs(riimhilolo[j]);
            *resmaxlololo = abs(rirelololo[j]) + abs(riimlololo[j]);
         }
   }
   free(Arehihihi_h); free(Arelohihi_h); free(Arehilohi_h); free(Arelolohi_h);
   free(Arehihilo_h); free(Arelohilo_h); free(Arehilolo_h); free(Arelololo_h);
   free(Aimhihihi_h); free(Aimlohihi_h); free(Aimhilohi_h); free(Aimlolohi_h);
   free(Aimhihilo_h); free(Aimlohilo_h); free(Aimhilolo_h); free(Aimlololo_h);

   cudaFree(rrehihihi_d); cudaFree(rrelohihi_d);
   cudaFree(rrehilohi_d); cudaFree(rrelolohi_d);
   cudaFree(rrehihilo_d); cudaFree(rrelohilo_d);
   cudaFree(rrehilolo_d); cudaFree(rrelololo_d);
   cudaFree(rimhihihi_d); cudaFree(rimlohihi_d);
   cudaFree(rimhilohi_d); cudaFree(rimlolohi_d);
   cudaFree(rimhihilo_d); cudaFree(rimlohilo_d);
   cudaFree(rimhilolo_d); cudaFree(rimlololo_d);
   cudaFree(xrehihihi_d); cudaFree(xrelohihi_d);
   cudaFree(xrehilohi_d); cudaFree(xrelolohi_d);
   cudaFree(xrehihilo_d); cudaFree(xrelohilo_d);
   cudaFree(xrehilolo_d); cudaFree(xrelololo_d);
   cudaFree(ximhihihi_d); cudaFree(ximlohihi_d);
   cudaFree(ximhilohi_d); cudaFree(ximlolohi_d);
   cudaFree(ximhihilo_d); cudaFree(ximlohilo_d);
   cudaFree(ximhilolo_d); cudaFree(ximlololo_d);
   cudaFree(Arehihihi_d); cudaFree(Arelohihi_d);
   cudaFree(Arehilohi_d); cudaFree(Arelolohi_d);
   cudaFree(Arehihilo_d); cudaFree(Arelohilo_d);
   cudaFree(Arehilolo_d); cudaFree(Arelololo_d);
   cudaFree(Aimhihihi_d); cudaFree(Aimlohihi_d);
   cudaFree(Aimhilohi_d); cudaFree(Aimlolohi_d);
   cudaFree(Aimhihilo_d); cudaFree(Aimlohilo_d);
   cudaFree(Aimhilolo_d); cudaFree(Aimlololo_d);
}
