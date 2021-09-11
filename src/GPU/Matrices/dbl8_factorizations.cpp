/* The file dbl8_factorizations.cpp defines functions specified in
 * the file dbl8_factorizations.h. */

#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"

// #include <iostream>
// using namespace std;

void CPU_dbl8_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Bhihihi, double **Blohihi, double **Bhilohi, double **Blolohi,
   double **Bhihilo, double **Blohilo, double **Bhilolo, double **Blololo,
   double **Chihihi, double **Clohihi, double **Chilohi, double **Clolohi,
   double **Chihilo, double **Clohilo, double **Chilolo, double **Clololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Chihihi[i][j] = 0.0; Clohihi[i][j] = 0.0;
         Chilohi[i][j] = 0.0; Clolohi[i][j] = 0.0;
         Chihilo[i][j] = 0.0; Clohilo[i][j] = 0.0;
         Chilolo[i][j] = 0.0; Clololo[i][j] = 0.0;

         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            odf_mul(Ahihihi[i][k],Alohihi[i][k],Ahilohi[i][k],Alolohi[i][k],
                    Ahihilo[i][k],Alohilo[i][k],Ahilolo[i][k],Alololo[i][k],
                    Bhihihi[k][j],Blohihi[k][j],Bhilohi[k][j],Blolohi[k][j],
                    Bhihilo[k][j],Blohilo[k][j],Bhilolo[k][j],Blololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&Chihihi[i][j],&Clohihi[i][j],&Chilohi[i][j],&Clolohi[i][j],
                    &Chihilo[i][j],&Clohilo[i][j],&Chilolo[i][j],&Clololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }
}

void CPU_cmplx8_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Brehihihi, double **Brelohihi,
   double **Brehilohi, double **Brelolohi,
   double **Brehihilo, double **Brelohilo,
   double **Brehilolo, double **Brelololo,
   double **Bimhihihi, double **Bimlohihi,
   double **Bimhilohi, double **Bimlolohi,
   double **Bimhihilo, double **Bimlohilo,
   double **Bimhilolo, double **Bimlololo,
   double **Crehihihi, double **Crelohihi,
   double **Crehilohi, double **Crelolohi,
   double **Crehihilo, double **Crelohilo,
   double **Crehilolo, double **Crelololo,
   double **Cimhihihi, double **Cimlohihi,
   double **Cimhilohi, double **Cimlolohi,
   double **Cimhihilo, double **Cimlohilo,
   double **Cimhilolo, double **Cimlololo )
{
   double zrehihihi,zrelohihi,zrehilohi,zrelolohi;
   double zrehihilo,zrelohilo,zrehilolo,zrelololo;
   double zimhihihi,zimlohihi,zimhilohi,zimlolohi;
   double zimhihilo,zimlohilo,zimhilolo,zimlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Crehihihi[i][j] = 0.0; Crelohihi[i][j] = 0.0;
         Crehilohi[i][j] = 0.0; Crelolohi[i][j] = 0.0;
         Crehihilo[i][j] = 0.0; Crelohilo[i][j] = 0.0;
         Crehilolo[i][j] = 0.0; Crelololo[i][j] = 0.0;
         Cimhihihi[i][j] = 0.0; Cimlohihi[i][j] = 0.0;
         Cimhilohi[i][j] = 0.0; Cimlolohi[i][j] = 0.0;
         Cimhihilo[i][j] = 0.0; Cimlohilo[i][j] = 0.0;
         Cimhilolo[i][j] = 0.0; Cimlololo[i][j] = 0.0;

         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            // zre = Are[i][k]*Bre[k][j] - Aim[i][k]*Bim[k][j];
            odf_mul(Arehihihi[i][k],Arelohihi[i][k],
                    Arehilohi[i][k],Arelolohi[i][k],
                    Arehihilo[i][k],Arelohilo[i][k],
                    Arehilolo[i][k],Arelololo[i][k],
                    Brehihihi[k][j],Brelohihi[k][j],
                    Brehilohi[k][j],Brelolohi[k][j],
                    Brehihilo[k][j],Brelohilo[k][j],
                    Brehilolo[k][j],Brelololo[k][j],
                    &zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                    &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo);
            odf_mul(Aimhihihi[i][k],Aimlohihi[i][k],
                    Aimhilohi[i][k],Aimlolohi[i][k],
                    Aimhihilo[i][k],Aimlohilo[i][k],
                    Aimhilolo[i][k],Aimlololo[i][k],
                    Bimhihihi[k][j],Bimlohihi[k][j],
                    Bimhilohi[k][j],Bimlolohi[k][j],
                    Bimhihilo[k][j],Bimlohilo[k][j],
                    Bimhilolo[k][j],Bimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_dec(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                    &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            // zim = Aim[i][k]*Bre[k][j] + Are[i][k]*Bim[k][j];
            odf_mul(Aimhihihi[i][k],Aimlohihi[i][k],
                    Aimhilohi[i][k],Aimlolohi[i][k],
                    Aimhihilo[i][k],Aimlohilo[i][k],
                    Aimhilolo[i][k],Aimlololo[i][k],
                    Brehihihi[k][j],Brelohihi[k][j],
                    Brehilohi[k][j],Brelolohi[k][j],
                    Brehihilo[k][j],Brelohilo[k][j],
                    Brehilolo[k][j],Brelololo[k][j],
                    &zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                    &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo);
            odf_mul(Arehihihi[i][k],Arelohihi[i][k],
                    Arehilohi[i][k],Arelolohi[i][k],
                    Arehihilo[i][k],Arelohilo[i][k],
                    Arehilolo[i][k],Arelololo[i][k],
                    Bimhihihi[k][j],Bimlohihi[k][j],
                    Bimhilohi[k][j],Bimlolohi[k][j],
                    Bimhihilo[k][j],Bimlohilo[k][j],
                    Bimhilolo[k][j],Bimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                    &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            // Cre[i][j] = Cre[i][j] + zre;
            odf_inc(&Crehihihi[i][j],&Crelohihi[i][j],
                    &Crehilohi[i][j],&Crelolohi[i][j],
                    &Crehihilo[i][j],&Crelohilo[i][j],
                    &Crehilolo[i][j],&Crelololo[i][j],
                    zrehihihi,zrelohihi,zrehilohi,zrelolohi,
                    zrehihilo,zrelohilo,zrehilolo,zrelololo);
            // Cim[i][j] = Cim[i][j] + zim;
            odf_inc(&Cimhihihi[i][j],&Cimlohihi[i][j],
                    &Cimhilohi[i][j],&Cimlolohi[i][j],
                    &Cimhihilo[i][j],&Cimlohilo[i][j],
                    &Cimhilolo[i][j],&Cimlololo[i][j],
                    zimhihihi,zimlohihi,zimhilohi,zimlolohi,
                    zimhihilo,zimlohilo,zimhilolo,zimlololo);
         }
      }
}

void CPU_dbl8_factors_forward
 ( int dim,
   double **Lhihihi, double **Llohihi, double **Lhilohi, double **Llolohi,
   double **Lhihilo, double **Llohilo, double **Lhilolo, double **Llololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   xhihihi[0] = bhihihi[0]; xlohihi[0] = blohihi[0];
   xhilohi[0] = bhilohi[0]; xlolohi[0] = blolohi[0];
   xhihilo[0] = bhihilo[0]; xlohilo[0] = blohilo[0];
   xhilolo[0] = bhilolo[0]; xlololo[0] = blololo[0];

   for(int i=1; i<dim; i++)
   {
      xhihihi[i] = bhihihi[i]; xlohihi[i] = blohihi[i];
      xhilohi[i] = bhilohi[i]; xlolohi[i] = blolohi[i];
      xhihilo[i] = bhihilo[i]; xlohilo[i] = blohilo[i];
      xhilolo[i] = bhilolo[i]; xlololo[i] = blololo[i];

      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         odf_mul(Lhihihi[i][j],Llohihi[i][j],Lhilohi[i][j],Llolohi[i][j],
                 Lhihilo[i][j],Llohilo[i][j],Lhilolo[i][j],Llololo[i][j],
                 xhihihi[j],   xlohihi[j],   xhilohi[j],   xlolohi[j],
                 xhihilo[j],   xlohilo[j],   xhilolo[j],   xlololo[j],
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odf_dec(&xhihihi[i],&xlohihi[i],&xhilohi[i],&xlolohi[i],
                 &xhihilo[i],&xlohilo[i],&xhilolo[i],&xlololo[i],
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
      }
   }
}

void CPU_cmplx8_factors_forward
 ( int dim,
   double **Lrehihihi, double **Lrelohihi,
   double **Lrehilohi, double **Lrelolohi,
   double **Lrehihilo, double **Lrelohilo,
   double **Lrehilolo, double **Lrelololo,
   double **Limhihihi, double **Limlohihi,
   double **Limhilohi, double **Limlolohi,
   double **Limhihilo, double **Limlohilo,
   double **Limhilolo, double **Limlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo )
{
   double acc1rehihihi,acc1relohihi,acc1rehilohi,acc1relolohi;
   double acc1rehihilo,acc1relohilo,acc1rehilolo,acc1relololo;
   double acc1imhihihi,acc1imlohihi,acc1imhilohi,acc1imlolohi;
   double acc1imhihilo,acc1imlohilo,acc1imhilolo,acc1imlololo;
   double acc2rehihihi,acc2relohihi,acc2rehilohi,acc2relolohi;
   double acc2rehihilo,acc2relohilo,acc2rehilolo,acc2relololo;
   double acc2imhihihi,acc2imlohihi,acc2imhilohi,acc2imlolohi;
   double acc2imhihilo,acc2imlohilo,acc2imhilolo,acc2imlololo;

   xrehihihi[0] = brehihihi[0]; xrelohihi[0] = brelohihi[0];
   xrehilohi[0] = brehilohi[0]; xrelolohi[0] = brelolohi[0];
   xrehihilo[0] = brehihilo[0]; xrelohilo[0] = brelohilo[0];
   xrehilolo[0] = brehilolo[0]; xrelololo[0] = brelololo[0];
   ximhihihi[0] = bimhihihi[0]; ximlohihi[0] = bimlohihi[0];
   ximhilohi[0] = bimhilohi[0]; ximlolohi[0] = bimlolohi[0];
   ximhihilo[0] = bimhihilo[0]; ximlohilo[0] = bimlohilo[0];
   ximhilolo[0] = bimhilolo[0]; ximlololo[0] = bimlololo[0];

   for(int i=1; i<dim; i++)
   {
      xrehihihi[i] = brehihihi[i]; xrelohihi[i] = brelohihi[i];
      xrehilohi[i] = brehilohi[i]; xrelolohi[i] = brelolohi[i];
      xrehihilo[i] = brehihilo[i]; xrelohilo[i] = brelohilo[i];
      xrehilolo[i] = brehilolo[i]; xrelololo[i] = brelololo[i];
      ximhihihi[i] = bimhihihi[i]; ximlohihi[i] = bimlohihi[i];
      ximhilohi[i] = bimhilohi[i]; ximlolohi[i] = bimlolohi[i];
      ximhihilo[i] = bimhihilo[i]; ximlohilo[i] = bimlohilo[i];
      ximhilolo[i] = bimhilolo[i]; ximlololo[i] = bimlololo[i];

      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         odf_mul(Lrehihihi[i][j],Lrelohihi[i][j],
                 Lrehilohi[i][j],Lrelolohi[i][j],
                 Lrehihilo[i][j],Lrelohilo[i][j],
                 Lrehilolo[i][j],Lrelololo[i][j],
                 xrehihihi[j],   xrelohihi[j],   xrehilohi[j],   xrelolohi[j],
                 xrehihilo[j],   xrelohilo[j],   xrehilolo[j],   xrelololo[j],
             &acc1rehihihi,  &acc1relohihi,  &acc1rehilohi,  &acc1relolohi,
             &acc1rehihilo,  &acc1relohilo,  &acc1rehilolo,  &acc1relololo);
         odf_mul(Limhihihi[i][j],Limlohihi[i][j],
                 Limhilohi[i][j],Limlolohi[i][j],
                 Limhihilo[i][j],Limlohilo[i][j],
                 Limhilolo[i][j],Limlololo[i][j],
                 ximhihihi[j],   ximlohihi[j],   ximhilohi[j],   ximlolohi[j],
                 ximhihilo[j],   ximlohilo[j],   ximhilolo[j],   ximlololo[j],
             &acc1imhihihi,  &acc1imlohihi,  &acc1imhilohi,  &acc1imlolohi,
             &acc1imhihilo,  &acc1imlohilo,  &acc1imhilolo,  &acc1imlololo);
         odf_mul(Limhihihi[i][j],Limlohihi[i][j],
                 Limhilohi[i][j],Limlolohi[i][j],
                 Limhihilo[i][j],Limlohilo[i][j],
                 Limhilolo[i][j],Limlololo[i][j],
                 xrehihihi[j],   xrelohihi[j],   xrehilohi[j],   xrelolohi[j],
                 xrehihilo[j],   xrelohilo[j],   xrehilolo[j],   xrelololo[j],
             &acc2rehihihi,  &acc2relohihi,  &acc2rehilohi,  &acc2relolohi,
             &acc2rehihilo,  &acc2relohilo,  &acc2rehilolo,  &acc2relololo);
         odf_mul(Lrehihihi[i][j],Lrelohihi[i][j],
                 Lrehilohi[i][j],Lrelolohi[i][j],
                 Lrehihilo[i][j],Lrelohilo[i][j],
                 Lrehilolo[i][j],Lrelololo[i][j],
                 ximhihihi[j],   ximlohihi[j],   ximhilohi[j],   ximlolohi[j],
                 ximhihilo[j],   ximlohilo[j],   ximhilolo[j],   ximlololo[j],
             &acc2imhihihi,  &acc2imlohihi,  &acc2imhilohi,  &acc2imlolohi,
             &acc2imhihilo,  &acc2imlohilo,  &acc2imhilolo,  &acc2imlololo);

         odf_dec(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
                 &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
               acc1rehihihi, acc1relohihi, acc1rehilohi, acc1relolohi,
               acc1rehihilo, acc1relohilo, acc1rehilolo, acc1relololo);
         odf_inc(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
                 &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
               acc1imhihihi, acc1imlohihi, acc1imhilohi, acc1imlolohi,
               acc1imhihilo, acc1imlohilo, acc1imhilolo, acc1imlololo);
         odf_dec(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
                 &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
               acc2rehihihi, acc2relohihi, acc2rehilohi, acc2relolohi,
               acc2rehihilo, acc2relohilo, acc2rehilolo, acc2relololo);
         odf_dec(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
                 &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
               acc2imhihihi, acc2imlohihi, acc2imhilohi, acc2imlolohi,
               acc2imhihilo, acc2imlohilo, acc2imhilolo, acc2imlololo);
      }
   }
}

void CPU_dbl8_factors_backward
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=dim-1; i>=0; i--)
   {
      xhihihi[i] = bhihihi[i]; xlohihi[i] = blohihi[i];
      xhilohi[i] = bhilohi[i]; xlolohi[i] = blolohi[i];
      xhihilo[i] = bhihilo[i]; xlohilo[i] = blohilo[i];
      xhilolo[i] = bhilolo[i]; xlololo[i] = blololo[i];

      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         odf_mul(Uhihihi[i][j],Ulohihi[i][j],Uhilohi[i][j],Ulolohi[i][j],
                 Uhihilo[i][j],Ulohilo[i][j],Uhilolo[i][j],Ulololo[i][j],
                 xhihihi[j],   xlohihi[j],   xhilohi[j],   xlolohi[j],
                 xhihilo[j],   xlohilo[j],   xhilolo[j],   xlololo[j],
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odf_dec(&xhihihi[i],&xlohihi[i],&xhilohi[i],&xlolohi[i],
                 &xhihilo[i],&xlohilo[i],&xhilolo[i],&xlololo[i],
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
      }
      // x[i] = x[i]/U[i][i];
      odf_div(xhihihi[i],   xlohihi[i],   xhilohi[i],   xlolohi[i],
              xhihilo[i],   xlohilo[i],   xhilolo[i],   xlololo[i],
              Uhihihi[i][i],Ulohihi[i][i],Uhilohi[i][i],Ulolohi[i][i],
              Uhihilo[i][i],Ulohilo[i][i],Uhilolo[i][i],Ulololo[i][i],
           &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
           &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
      xhihihi[i] = acchihihi; xlohihi[i] = acclohihi;
      xhilohi[i] = acchilohi; xlolohi[i] = acclolohi;
      xhihilo[i] = acchihilo; xlohilo[i] = acclohilo;
      xhilolo[i] = acchilolo; xlololo[i] = acclololo;
   }
}

void CPU_cmplx8_factors_backward
 ( int dim,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi,
   double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo,
   double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi,
   double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo,
   double *ximhilolo, double *ximlololo )
{
   double acc1rehihihi,acc1relohihi,acc1rehilohi,acc1relolohi;
   double acc1rehihilo,acc1relohilo,acc1rehilolo,acc1relololo;
   double acc1imhihihi,acc1imlohihi,acc1imhilohi,acc1imlolohi;
   double acc1imhihilo,acc1imlohilo,acc1imhilolo,acc1imlololo;
   double acc2rehihihi,acc2relohihi,acc2rehilohi,acc2relolohi;
   double acc2rehihilo,acc2relohilo,acc2rehilolo,acc2relololo;
   double acc2imhihihi,acc2imlohihi,acc2imhilohi,acc2imlolohi;
   double acc2imhihilo,acc2imlohilo,acc2imhilolo,acc2imlololo;
   double acc3rehihihi,acc3relohihi,acc3rehilohi,acc3relolohi;
   double acc3rehihilo,acc3relohilo,acc3rehilolo,acc3relololo;
   double acc3imhihihi,acc3imlohihi,acc3imhilohi,acc3imlolohi;
   double acc3imhihilo,acc3imlohilo,acc3imhilolo,acc3imlololo;
   double denhihihi,denlohihi,denhilohi,denlolohi;
   double denhihilo,denlohilo,denhilolo,denlololo;

   for(int i=dim-1; i>=0; i--)
   {
      xrehihihi[i] = brehihihi[i]; xrelohihi[i] = brelohihi[i];
      xrehilohi[i] = brehilohi[i]; xrelolohi[i] = brelolohi[i];
      ximhihilo[i] = bimhihilo[i]; ximlohilo[i] = bimlohilo[i];
      ximhilolo[i] = bimhilolo[i]; ximlololo[i] = bimlololo[i];
      xrehihihi[i] = brehihihi[i]; xrelohihi[i] = brelohihi[i];
      xrehilohi[i] = brehilohi[i]; xrelolohi[i] = brelolohi[i];
      ximhihilo[i] = bimhihilo[i]; ximlohilo[i] = bimlohilo[i];
      ximhilolo[i] = bimhilolo[i]; ximlololo[i] = bimlololo[i];

      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         odf_mul(Urehihihi[i][j],Urelohihi[i][j],
                 Urehilohi[i][j],Urelolohi[i][j],
                 Urehihilo[i][j],Urelohilo[i][j],
                 Urehilolo[i][j],Urelololo[i][j],
                 xrehihihi[j],   xrelohihi[j],   xrehilohi[j],   xrelolohi[j],
                 xrehihilo[j],   xrelohilo[j],   xrehilolo[j],   xrelololo[j],
             &acc1rehihihi,  &acc1relohihi,  &acc1rehilohi,  &acc1relolohi,
             &acc1rehihilo,  &acc1relohilo,  &acc1rehilolo,  &acc1relololo);
         odf_mul(Uimhihihi[i][j],Uimlohihi[i][j],
                 Uimhilohi[i][j],Uimlolohi[i][j],
                 Uimhihilo[i][j],Uimlohilo[i][j],
                 Uimhilolo[i][j],Uimlololo[i][j],
                 ximhihihi[j],   ximlohihi[j],   ximhilohi[j],   ximlolohi[j],
                 ximhihilo[j],   ximlohilo[j],   ximhilolo[j],   ximlololo[j],
             &acc1imhihihi,  &acc1imlohihi,  &acc1imhilohi,  &acc1imlolohi,
             &acc1imhihilo,  &acc1imlohilo,  &acc1imhilolo,  &acc1imlololo);
         odf_mul(Uimhihihi[i][j],Uimlohihi[i][j],
                 Uimhilohi[i][j],Uimlolohi[i][j],
                 Uimhihilo[i][j],Uimlohilo[i][j],
                 Uimhilolo[i][j],Uimlololo[i][j],
                 xrehihihi[j],   xrelohihi[j],   xrehilohi[j],   xrelolohi[j],
                 xrehihilo[j],   xrelohilo[j],   xrehilolo[j],   xrelololo[j],
             &acc2rehihihi,  &acc2relohihi,  &acc2rehilohi,  &acc2relolohi,
             &acc2rehihilo,  &acc2relohilo,  &acc2rehilolo,  &acc2relololo);
         odf_mul(Urehihihi[i][j],Urelohihi[i][j],
                 Urehilohi[i][j],Urelolohi[i][j],
                 Urehihilo[i][j],Urelohilo[i][j],
                 Urehilolo[i][j],Urelololo[i][j],
                 ximhihihi[j],   ximlohihi[j],   ximhilohi[j],   ximlolohi[j],
                 ximhihilo[j],   ximlohilo[j],   ximhilolo[j],   ximlololo[j],
             &acc2imhihihi,  &acc2imlohihi,  &acc2imhilohi,  &acc2imlolohi,
             &acc2imhihilo,  &acc2imlohilo,  &acc2imhilolo,  &acc2imlololo);

         odf_dec(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
                 &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
               acc1rehihihi, acc1relohihi, acc1rehilohi, acc1relolohi,
               acc1rehihilo, acc1relohilo, acc1rehilolo, acc1relololo);
         odf_inc(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
                 &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
               acc1imhihihi, acc1imlohihi, acc1imhilohi, acc1imlolohi,
               acc1imhihilo, acc1imlohilo, acc1imhilolo, acc1imlololo);
         odf_dec(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
                 &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
               acc2rehihihi, acc2relohihi, acc2rehilohi, acc2relolohi,
               acc2rehihilo, acc2relohilo, acc2rehilolo, acc2relololo);
         odf_dec(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
                 &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
               acc2imhihihi, acc2imlohihi, acc2imhilohi, acc2imlolohi,
               acc2imhihilo, acc2imlohilo, acc2imhilolo, acc2imlololo);
      }
      // x[i] = x[i]/U[i][i];
      odf_mul(Urehihihi[i][i],Urelohihi[i][i],Urehilohi[i][i],Urelolohi[i][i],
              Urehihilo[i][i],Urelohilo[i][i],Urehilolo[i][i],Urelololo[i][i],
              Urehihihi[i][i],Urelohihi[i][i],Urehilohi[i][i],Urelolohi[i][i],
              Urehihilo[i][i],Urelohilo[i][i],Urehilolo[i][i],Urelololo[i][i],
             &denhihihi,     &denlohihi,     &denhilohi,     &denlolohi,
             &denhihilo,     &denlohilo,     &denhilolo,     &denlololo);
      odf_mul(Uimhihihi[i][i],Uimlohihi[i][i],Uimhilohi[i][i],Uimlolohi[i][i],
              Uimhihilo[i][i],Uimlohilo[i][i],Uimhilolo[i][i],Uimlololo[i][i],
              Uimhihihi[i][i],Uimlohihi[i][i],Uimhilohi[i][i],Uimlolohi[i][i],
              Uimhihilo[i][i],Uimlohilo[i][i],Uimhilolo[i][i],Uimlololo[i][i],
          &acc1rehihihi,  &acc1relohihi,  &acc1rehilohi,  &acc1relolohi,
          &acc1rehihilo,  &acc1relohilo,  &acc1rehilolo,  &acc1relololo);
      odf_inc(&denhihihi,  &denlohihi,  &denhilohi,  &denlolohi,
              &denhihilo,  &denlohilo,  &denhilolo,  &denlololo,
            acc1rehihihi,acc1relohihi,acc1rehilohi,acc1relolohi,
            acc1rehihilo,acc1relohilo,acc1rehilolo,acc1relololo); // denominator
      odf_div(Urehihihi[i][i],Urelohihi[i][i],Urehilohi[i][i],Urelolohi[i][i],
              Urehihilo[i][i],Urelohilo[i][i],Urehilolo[i][i],Urelololo[i][i],
              denhihihi,      denlohihi,      denhilohi,      denlolohi,
              denhihilo,      denlohilo,      denhilolo,      denlololo,
          &acc1rehihihi,  &acc1relohihi,  &acc1rehilohi,  &acc1relolohi,
          &acc1rehihilo,  &acc1relohilo,  &acc1rehilolo,  &acc1relololo);
      // (acc1rehi,acc1relo) is real part of 1/U[i][i]
      odf_div(Uimhihihi[i][i],Uimlohihi[i][i],Uimhilohi[i][i],Uimlolohi[i][i],
              Uimhihilo[i][i],Uimlohilo[i][i],Uimhilolo[i][i],Uimlololo[i][i],
              denhihihi,      denlohihi,      denhilohi,      denlolohi,
              denhihilo,      denlohilo,      denhilolo,      denlololo,
          &acc1imhihihi,  &acc1imlohihi,  &acc1imhilohi,  &acc1imlolohi,
          &acc1imhihilo,  &acc1imlohilo,  &acc1imhilolo,  &acc1imlololo);
      odf_minus(&acc1imhihihi,&acc1imlohihi,&acc1imhilohi,&acc1imlolohi,
                &acc1imhihilo,&acc1imlohilo,&acc1imhilolo,&acc1imlololo);
      // (acc1imhi,acc1imlo) is imaginary part of 1/U[i][i]
      odf_mul(xrehihihi[i], xrelohihi[i], xrehilohi[i], xrelolohi[i],
              xrehihilo[i], xrelohilo[i], xrehilolo[i], xrelololo[i],
           acc1rehihihi, acc1relohihi, acc1rehilohi, acc1relolohi,
           acc1rehihilo, acc1relohilo, acc1rehilolo, acc1relololo,
          &acc2rehihihi,&acc2relohihi,&acc2rehilohi,&acc2relolohi,
          &acc2rehihilo,&acc2relohilo,&acc2rehilolo,&acc2relololo);
      odf_mul(ximhihihi[i], ximlohihi[i], ximhilohi[i], ximlolohi[i],
              ximhihilo[i], ximlohilo[i], ximhilolo[i], ximlololo[i],
           acc1imhihihi, acc1imlohihi, acc1imhilohi, acc1imlolohi,
           acc1imhihilo, acc1imlohilo, acc1imhilolo, acc1imlololo,
          &acc2imhihihi,&acc2imlohihi,&acc2imhilohi,&acc2imlolohi,
          &acc2imhihilo,&acc2imlohilo,&acc2imhilolo,&acc2imlololo);
      // acc2 stores the doubles for xre
      odf_mul(ximhihihi[i], ximlohihi[i], ximhilohi[i], ximlolohi[i],
              ximhihilo[i], ximlohilo[i], ximhilolo[i], ximlololo[i],
           acc1rehihihi, acc1relohihi, acc1rehilohi, acc1relolohi,
           acc1rehihilo, acc1relohilo, acc1rehilolo, acc1relololo,
          &acc3rehihihi,&acc3relohihi,&acc3rehilohi,&acc3relolohi,
          &acc3rehihilo,&acc3relohilo,&acc3rehilolo,&acc3relololo);
      odf_mul(xrehihihi[i], xrelohihi[i], xrehilohi[i], xrelolohi[i],
              xrehihilo[i], xrelohilo[i], xrehilolo[i], xrelololo[i],
           acc1imhihihi, acc1imlohihi, acc1imhilohi, acc1imlolohi,
           acc1imhihilo, acc1imlohilo, acc1imhilolo, acc1imlololo,
          &acc3imhihihi,&acc3imlohihi,&acc3imhilohi,&acc3imlolohi,
          &acc3imhihilo,&acc3imlohilo,&acc3imhilolo,&acc3imlololo);
      // acc3 stores the doubles for xim
      xrehihihi[i] = acc2rehihihi; xrelohihi[i] = acc2relohihi;
      xrehilolo[i] = acc2rehilohi; xrelolohi[i] = acc2relolohi;
      xrehihihi[i] = acc2rehihilo; xrelohilo[i] = acc2relohilo;
      xrehilolo[i] = acc2rehilolo; xrelololo[i] = acc2relololo;
      odf_dec(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
              &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
            acc2imhihihi, acc2imlohihi, acc2imhilohi, acc2imlolohi,
            acc2imhihilo, acc2imlohilo, acc2imhilolo, acc2imlololo);
      ximhihihi[i] = acc3rehihihi; ximlohihi[i] = acc3relohihi;
      ximhilohi[i] = acc3rehilohi; ximlolohi[i] = acc3relolohi;
      ximhihilo[i] = acc3rehihilo; ximlohilo[i] = acc3relohilo;
      ximhilolo[i] = acc3rehilolo; ximlololo[i] = acc3relololo;
      odf_inc(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
              &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
            acc3imhihihi, acc3imlohihi, acc3imhilohi, acc3imlolohi,
            acc3imhihilo, acc3imlohilo, acc3imhilolo, acc3imlololo);
   }
}

void CPU_dbl8_factors_lufac
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   int *pivots )
{
   double valmax,valtmp;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Ahihihi[j][j]); idxmax = j;     // find the pivot
      for(int i=j+1; i<dim; i++)
      {
         valtmp = fabs(Ahihihi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Ahihihi[idxmax][k];
            Ahihihi[idxmax][k] = Ahihihi[j][k];
            Ahihihi[j][k] = valtmp;
            valtmp = Alohihi[idxmax][k];
            Alohihi[idxmax][k] = Alohihi[j][k];
            Alohihi[j][k] = valtmp;
            valtmp = Ahilohi[idxmax][k];
            Ahilohi[idxmax][k] = Ahilohi[j][k];
            Ahilohi[j][k] = valtmp;
            valtmp = Alolohi[idxmax][k];
            Alolohi[idxmax][k] = Alolohi[j][k];
            Alolohi[j][k] = valtmp;
            valtmp = Ahihilo[idxmax][k];
            Ahihilo[idxmax][k] = Ahihilo[j][k];
            Ahihilo[j][k] = valtmp;
            valtmp = Alohilo[idxmax][k];
            Alohilo[idxmax][k] = Alohilo[j][k];
            Alohilo[j][k] = valtmp;
            valtmp = Ahilolo[idxmax][k];
            Ahilolo[idxmax][k] = Ahilolo[j][k];
            Ahilolo[j][k] = valtmp;
            valtmp = Alololo[idxmax][k];
            Alololo[idxmax][k] = Alololo[j][k];
            Alololo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         odf_div(Ahihihi[i][j],Alohihi[i][j],Ahilohi[i][j],Alolohi[i][j],
                 Ahihilo[i][j],Alohilo[i][j],Ahilolo[i][j],Alololo[i][j],
                 Ahihihi[j][j],Alohihi[j][j],Ahilohi[j][j],Alolohi[j][j],
                 Ahihilo[j][j],Alohilo[j][j],Ahilolo[j][j],Alololo[j][j],
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);

         Ahihihi[i][j] = acchihihi; Alohihi[i][j] = acclohihi;
         Ahilohi[i][j] = acchilohi; Alolohi[i][j] = acclolohi;
         Ahihilo[i][j] = acchihilo; Alohilo[i][j] = acclohilo;
         Ahilolo[i][j] = acchilolo; Alololo[i][j] = acclololo;

         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            odf_mul(Ahihihi[i][j],Alohihi[i][j],Ahilohi[i][j],Alolohi[i][j],
                    Ahihilo[i][j],Alohilo[i][j],Ahilolo[i][j],Alololo[i][j],
                    Ahihihi[j][k],Alohihi[j][k],Ahilohi[j][k],Alolohi[j][k],
                    Ahihilo[j][k],Alohilo[j][k],Ahilolo[j][k],Alololo[j][k],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
            odf_dec(&Ahihihi[i][k],&Alohihi[i][k],
                    &Ahilohi[i][k],&Alolohi[i][k],
                    &Ahihilo[i][k],&Alohilo[i][k],
                    &Ahilolo[i][k],&Alololo[i][k],
                   acchihihi,     acclohihi,     acchilohi,     acclolohi,
                   acchihilo,     acclohilo,     acchilolo,     acclololo);
         }
      }
   }
}

void CPU_cmplx8_factors_lufac
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   int *pivots )
{
   double valmax,valtmp;
   int idxmax,idxtmp;
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;
   double denhihihi,denlohihi,denhilohi,denlolohi;
   double denhihilo,denlohilo,denhilolo,denlololo;
   double acc3rehihihi,acc3relohihi,acc3rehilohi,acc3relolohi;
   double acc3rehihilo,acc3relohilo,acc3rehilolo,acc3relololo;
   double acc3imhihihi,acc3imlohihi,acc3imhilohi,acc3imlolohi;
   double acc3imhihilo,acc3imlohilo,acc3imhilolo,acc3imlololo;
   double acc4rehihihi,acc4relohihi,acc4rehilohi,acc4relolohi;
   double acc4rehihilo,acc4relohilo,acc4rehilolo,acc4relololo;
   double acc4imhihihi,acc4imlohihi,acc4imhilohi,acc4imlolohi;
   double acc4imhihilo,acc4imlohilo,acc4imhilolo,acc4imlololo;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Arehihihi[j][j]) + fabs(Aimhihihi[j][j]);
      idxmax = j;     // find the pivot
      for(int i=j+1; i<dim; i++)
      {
         valtmp = fabs(Arehihihi[i][j]) + fabs(Aimhihihi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Arehihihi[idxmax][k];
            Arehihihi[idxmax][k] = Arehihihi[j][k];
            Arehihihi[j][k] = valtmp;
            valtmp = Arelohihi[idxmax][k];
            Arelohihi[idxmax][k] = Arelohihi[j][k];
            Arelohihi[j][k] = valtmp;
            valtmp = Arehilohi[idxmax][k];
            Arehilohi[idxmax][k] = Arehilohi[j][k];
            Arehilohi[j][k] = valtmp;
            valtmp = Arelolohi[idxmax][k];
            Arelolohi[idxmax][k] = Arelolohi[j][k];
            Arelolohi[j][k] = valtmp;
            valtmp = Arehihilo[idxmax][k];
            Arehihilo[idxmax][k] = Arehihilo[j][k];
            Arehihilo[j][k] = valtmp;
            valtmp = Arelohilo[idxmax][k];
            Arelohilo[idxmax][k] = Arelohilo[j][k];
            Arelohilo[j][k] = valtmp;
            valtmp = Arehilolo[idxmax][k];
            Arehilolo[idxmax][k] = Arehilolo[j][k];
            Arehilolo[j][k] = valtmp;
            valtmp = Arelololo[idxmax][k];
            Arelololo[idxmax][k] = Arelololo[j][k];
            Arelololo[j][k] = valtmp;
            valtmp = Aimhihihi[idxmax][k];
            Aimhihihi[idxmax][k] = Aimhihihi[j][k];
            Aimhihihi[j][k] = valtmp;
            valtmp = Aimlohihi[idxmax][k];
            Aimlohihi[idxmax][k] = Aimlohihi[j][k];
            Aimlohihi[j][k] = valtmp;
            valtmp = Aimhilohi[idxmax][k];
            Aimhilohi[idxmax][k] = Aimhilohi[j][k];
            Aimhilohi[j][k] = valtmp;
            valtmp = Aimlolohi[idxmax][k];
            Aimlolohi[idxmax][k] = Aimlolohi[j][k];
            Aimlolohi[j][k] = valtmp;
            Aimhihilo[idxmax][k] = Aimhihilo[j][k];
            Aimhihilo[j][k] = valtmp;
            valtmp = Aimlohilo[idxmax][k];
            Aimlohilo[idxmax][k] = Aimlohilo[j][k];
            Aimlohilo[j][k] = valtmp;
            valtmp = Aimhilolo[idxmax][k];
            Aimhilolo[idxmax][k] = Aimhilolo[j][k];
            Aimhilolo[j][k] = valtmp;
            valtmp = Aimlololo[idxmax][k];
            Aimlololo[idxmax][k] = Aimlololo[j][k];
            Aimlololo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         odf_mul(Arehihihi[j][j],Arelohihi[j][j],
                 Arehilohi[j][j],Arelolohi[j][j],
                 Arehihilo[j][j],Arelohilo[j][j],
                 Arehilolo[j][j],Arelololo[j][j],
                 Arehihihi[j][j],Arelohihi[j][j],
                 Arehilohi[j][j],Arelolohi[j][j],
                 Arehihilo[j][j],Arelohilo[j][j],
                 Arehilolo[j][j],Arelololo[j][j],
                &denhihihi,     &denlohihi,     &denhilohi,     &denlolohi,
                &denhihilo,     &denlohilo,     &denhilolo,     &denlololo);
         odf_mul(Aimhihihi[j][j],Aimlohihi[j][j],
                 Aimhilohi[j][j],Aimlolohi[j][j],
                 Aimhihilo[j][j],Aimlohilo[j][j],
                 Aimhilolo[j][j],Aimlololo[j][j],
                 Aimhihihi[j][j],Aimlohihi[j][j],
                 Aimhilohi[j][j],Aimlolohi[j][j],
                 Aimhihilo[j][j],Aimlohilo[j][j],
                 Aimhilolo[j][j],Aimlololo[j][j],
               &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
               &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
         odf_inc(&denhihihi,&denlohihi,&denhilohi,&denlolohi,
                 &denhihilo,&denlohilo,&denhilolo,&denlololo,
                 acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
                 acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo); // denominator
         odf_div(Arehihihi[j][j],Arelohihi[j][j],
                 Arehilohi[j][j],Arelolohi[j][j],
                 Arehihilo[j][j],Arelohilo[j][j],
                 Arehilolo[j][j],Arelololo[j][j],
                 denhihihi,      denlohihi,      denhilohi,      denlolohi,
                 denhihilo,      denlohilo,      denhilolo,      denlololo,
               &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
               &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
         // (acc1hi,acc1lo) is real part of 1/A[j][j]
         odf_div(Aimhihihi[j][j],Aimlohihi[j][j],
                 Aimhilohi[j][j],Aimlolohi[j][j],
                 Aimhihilo[j][j],Aimlohilo[j][j],
                 Aimhilolo[j][j],Aimlololo[j][j],
                 denhihihi,      denlohihi,      denhilohi,      denlolohi,
                 denhihilo,      denlohilo,      denhilolo,      denlololo,
               &acc2hihihi,    &acc2lohihi,    &acc2hilohi,   &acc2lolohi,
               &acc2hihilo,    &acc2lohilo,    &acc2hilolo,   &acc2lololo);
         odf_minus(&acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
                   &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
         // (acc2hi,acc2lo) is imaginary part of 1/A[j][j]
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
                acc1hihihi,     acc1lohihi,     acc1hilohi,     acc1lolohi,
                acc1hihilo,     acc1lohilo,     acc1hilolo,     acc1lololo,
             &acc3rehihihi,  &acc3relohihi,  &acc3rehilohi,  &acc3relolohi,
             &acc3rehihilo,  &acc3relohilo,  &acc3rehilolo,  &acc3relololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
                acc2hihihi,     acc2lohihi,     acc2hilohi,     acc2lolohi,
                acc2hihilo,     acc2lohilo,     acc2hilolo,     acc2lololo,
             &acc3imhihihi,  &acc3imlohihi,  &acc3imhilohi,  &acc3imlolohi,
             &acc3imhihilo,  &acc3imlohilo,  &acc3imhilolo,  &acc3imlololo);
         // acc3 stores doubles for Arehi
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
                acc1hihihi,     acc1lohihi,     acc1hilohi,     acc1lolohi,
                acc1hihilo,     acc1lohilo,     acc1hilolo,     acc1lololo,
             &acc4rehihihi,  &acc4relohihi,  &acc4rehilohi,  &acc4relolohi,
             &acc4rehihilo,  &acc4relohilo,  &acc4rehilolo,  &acc4relololo);
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
                acc2hihihi     ,acc2lohihi,     acc2hilohi,     acc2lolohi,
                acc2hihilo     ,acc2lohilo,     acc2hilolo,     acc2lololo,
             &acc4imhihihi,  &acc4imlohihi,  &acc4imhilohi,  &acc4imlolohi,
             &acc4imhihilo,  &acc4imlohilo,  &acc4imhilolo,  &acc4imlololo);
         // acc4 stores doubles for Aimhi
         Arehihihi[i][j] = acc3rehihihi; Arelohihi[i][j] = acc3relohihi;
         Arehilolo[i][j] = acc3rehilohi; Arelolohi[i][j] = acc3relolohi;
         Arehihihi[i][j] = acc3rehihilo; Arelohilo[i][j] = acc3relohilo;
         Arehilolo[i][j] = acc3rehilolo; Arelololo[i][j] = acc3relololo;
         odf_dec(&Arehihihi[i][j],&Arelohihi[i][j],
                 &Arehilohi[i][j],&Arelolohi[i][j],
                 &Arehihilo[i][j],&Arelohilo[i][j],
                 &Arehilolo[i][j],&Arelololo[i][j],
               acc3imhihihi,    acc3imlohihi,  acc3imhilohi,  acc3imlolohi,
               acc3imhihilo,    acc3imlohilo,  acc3imhilolo,  acc3imlololo);
         Aimhihihi[i][j] = acc4rehihihi; Aimlohihi[i][j] = acc4relohihi;
         Aimhilohi[i][j] = acc4rehilohi; Aimlolohi[i][j] = acc4relolohi;
         Aimhihilo[i][j] = acc4rehihilo; Aimlohilo[i][j] = acc4relohilo;
         Aimhilolo[i][j] = acc4rehilolo; Aimlololo[i][j] = acc4relololo;
         odf_inc(&Aimhihihi[i][j],&Aimlohihi[i][j],
                 &Aimhilohi[i][j],&Aimlolohi[i][j],
                 &Aimhihilo[i][j],&Aimlohilo[i][j],
                 &Aimhilolo[i][j],&Aimlololo[i][j],
               acc4imhihihi,    acc4imlohihi,  acc4imhilohi,  acc4imlolohi,
               acc4imhihilo,    acc4imlohilo,  acc4imhilolo,  acc4imlololo);

         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                    Arehilohi[i][j],Arelolohi[i][j],
                    Arehihilo[i][j],Arelohilo[i][j],
                    Arehilolo[i][j],Arelololo[i][j],
                    Arehihihi[j][k],Arelohihi[j][k],
                    Arehilohi[j][k],Arelolohi[j][k],
                    Arehihilo[j][k],Arelohilo[j][k],
                    Arehilolo[j][k],Arelololo[j][k],
                &acc3rehihihi,  &acc3relohihi, &acc3rehilohi, &acc3relolohi,
                &acc3rehihilo,  &acc3relohilo, &acc3rehilolo, &acc3relololo);
            odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                    Aimhilohi[i][j],Aimlolohi[i][j],
                    Aimhihilo[i][j],Aimlohilo[i][j],
                    Aimhilolo[i][j],Aimlololo[i][j],
                    Aimhihihi[j][k],Aimlohihi[j][k],
                    Aimhilohi[j][k],Aimlolohi[j][k],
                    Aimhihilo[j][k],Aimlohilo[j][k],
                    Aimhilolo[j][k],Aimlololo[j][k],
                &acc3imhihihi,  &acc3imlohihi, &acc3imhilohi, &acc3imlolohi,
                &acc3imhihilo,  &acc3imlohilo, &acc3imhilolo, &acc3imlololo);
            odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                    Aimhilohi[i][j],Aimlolohi[i][j],
                    Aimhihilo[i][j],Aimlohilo[i][j],
                    Aimhilolo[i][j],Aimlololo[i][j],
                    Arehihihi[j][k],Arelohihi[j][k],
                    Arehilohi[j][k],Arelolohi[j][k],
                    Arehihilo[j][k],Arelohilo[j][k],
                    Arehilolo[j][k],Arelololo[j][k],
                &acc4rehihihi,  &acc4relohihi, &acc4rehilohi, &acc4relolohi,
                &acc4rehihilo,  &acc4relohilo, &acc4rehilolo, &acc4relololo);
            odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                    Arehilohi[i][j],Arelolohi[i][j],
                    Arehihilo[i][j],Arelohilo[i][j],
                    Arehilolo[i][j],Arelololo[i][j],
                    Aimhihihi[j][k],Aimlohihi[j][k],
                    Aimhilohi[j][k],Aimlolohi[j][k],
                    Aimhihilo[j][k],Aimlohilo[j][k],
                    Aimhilolo[j][k],Aimlololo[j][k],
                &acc4imhihihi,  &acc4imlohihi, &acc4imhilohi, &acc4imlolohi,
                &acc4imhihilo,  &acc4imlohilo, &acc4imhilolo, &acc4imlololo);

            odf_dec(&Arehihihi[i][k],&Arelohihi[i][k],
                    &Arehilohi[i][k],&Arelolohi[i][k],
                    &Arehihilo[i][k],&Arelohilo[i][k],
                    &Arehilolo[i][k],&Arelololo[i][k],
                    acc3rehihihi,acc3relohihi,acc3rehilohi,acc3relolohi,
                    acc3rehihilo,acc3relohilo,acc3rehilolo,acc3relololo);
            odf_inc(&Arehihihi[i][k],&Arelohihi[i][k],
                    &Arehilohi[i][k],&Arelolohi[i][k],
                    &Arehihilo[i][k],&Arelohilo[i][k],
                    &Arehilolo[i][k],&Arelololo[i][k],
                    acc3imhihihi,acc3imlohihi,acc3imhilohi,acc3imlolohi,
                    acc3imhihilo,acc3imlohilo,acc3imhilolo,acc3imlololo);
            odf_dec(&Aimhihihi[i][k],&Aimlohihi[i][k],
                    &Aimhilohi[i][k],&Aimlolohi[i][k],
                    &Aimhihilo[i][k],&Aimlohilo[i][k],
                    &Aimhilolo[i][k],&Aimlololo[i][k],
                    acc4rehihihi,acc4relohihi,acc4rehilohi,acc4relolohi,
                    acc4rehihilo,acc4relohilo,acc4rehilolo,acc4relololo);
            odf_dec(&Aimhihihi[i][k],&Aimlohihi[i][k],
                    &Aimhilohi[i][k],&Aimlolohi[i][k],
                    &Aimhihilo[i][k],&Aimlohilo[i][k],
                    &Aimhilolo[i][k],&Aimlololo[i][k],
                    acc4imhihihi,acc4imlohihi,acc4imhilohi,acc4imlolohi,
                    acc4imhihilo,acc4imlohilo,acc4imhilolo,acc4imlololo);
         }
      }
   }
}

void CPU_dbl8_factors_lusolve
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   int *pivots,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo )
{
   CPU_dbl8_factors_lufac
      (dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
           Ahihilo,Alohilo,Ahilolo,Alololo,pivots);

   for(int i=0; i<dim; i++) 
   {
      xhihihi[i] = bhihihi[pivots[i]];
      xlohihi[i] = blohihi[pivots[i]];
      xhilohi[i] = bhilohi[pivots[i]];
      xlolohi[i] = blolohi[pivots[i]];
      xhihilo[i] = bhihilo[pivots[i]];
      xlohilo[i] = blohilo[pivots[i]];
      xhilolo[i] = bhilolo[pivots[i]];
      xlololo[i] = blololo[pivots[i]];
   }
   CPU_dbl8_factors_forward
      (dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
           Ahihilo,Alohilo,Ahilolo,Alololo,
           xhihihi,xlohihi,xhilohi,xlolohi,
           xhihilo,xlohilo,xhilolo,xlololo,
           bhihihi,blohihi,bhilohi,blolohi,
           bhihilo,blohilo,bhilolo,blololo);

   CPU_dbl8_factors_backward
      (dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
           Ahihilo,Alohilo,Ahilolo,Alololo,
           bhihihi,blohihi,bhilohi,blolohi,
           bhihilo,blohilo,bhilolo,blololo,
           xhihihi,xlohihi,xhilohi,xlolohi,
           xhihilo,xlohilo,xhilolo,xlololo);
}

void CPU_cmplx8_factors_lusolve
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   int *pivots,
   double *brehihihi, double *brelohihi,
   double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo,
   double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi,
   double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo,
   double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi,
   double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo,
   double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi,
   double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo,
   double *ximhilolo, double *ximlololo )
{
   CPU_cmplx8_factors_lufac
      (dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
           Arehihilo,Arelohilo,Arehilolo,Arelololo,
           Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
           Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,pivots);

   for(int i=0; i<dim; i++) 
   {
      xrehihihi[i] = brehihihi[pivots[i]];
      xrelohihi[i] = brelohihi[pivots[i]];
      xrehilohi[i] = brehilohi[pivots[i]];
      xrelolohi[i] = brelolohi[pivots[i]];
      xrehihilo[i] = brehihilo[pivots[i]];
      xrelohilo[i] = brelohilo[pivots[i]];
      xrehilolo[i] = brehilolo[pivots[i]];
      xrelololo[i] = brelololo[pivots[i]];
      ximhihihi[i] = bimhihihi[pivots[i]];
      ximlohihi[i] = bimlohihi[pivots[i]];
      ximhilohi[i] = bimhilohi[pivots[i]];
      ximlolohi[i] = bimlolohi[pivots[i]];
      ximhihilo[i] = bimhihilo[pivots[i]];
      ximlohilo[i] = bimlohilo[pivots[i]];
      ximhilolo[i] = bimhilolo[pivots[i]];
      ximlololo[i] = bimlololo[pivots[i]];
   }
   CPU_cmplx8_factors_forward
      (dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
           Arehihilo,Arelohilo,Arehilolo,Arelololo,
           Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
           Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
           xrehihihi,xrelohihi,xrehilohi,xrelolohi,
           xrehihilo,xrelohilo,xrehilolo,xrelololo,
           ximhihihi,ximlohihi,ximhilohi,ximlolohi,
           ximhihilo,ximlohilo,ximhilolo,ximlololo,
           brehihihi,brelohihi,brehilohi,brelolohi,
           brehihilo,brelohilo,brehilolo,brelololo,
           bimhihihi,bimlohihi,bimhilohi,bimlolohi,
           bimhihilo,bimlohilo,bimhilolo,bimlololo);

   CPU_cmplx8_factors_backward
      (dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
           Arehihilo,Arelohilo,Arehilolo,Arelololo,
           Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
           Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
           brehihihi,brelohihi,brehilohi,brelolohi,
           brehihilo,brelohilo,brehilolo,brelololo,
           bimhihihi,bimlohihi,bimhilohi,bimlolohi,
           bimhihilo,bimlohilo,bimhilolo,bimlololo,
           xrehihihi,xrelohihi,xrehilohi,xrelolohi,
           xrehihilo,xrelohilo,xrehilolo,xrelololo,
           ximhihihi,ximlohihi,ximhilohi,ximlolohi,
           ximhihilo,ximlohilo,ximhilolo,ximlololo);
}

/*

void CPU_dbl4_factors_house
 ( int n,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
{
   double sigmahihi = 0.0;
   double sigmalohi = 0.0;
   double sigmahilo = 0.0;
   double sigmalolo = 0.0;
   double muhihi,mulohi,muhilo,mulolo;
   double v0p2hihi,v0p2lohi,v0p2hilo,v0p2lolo;
   double acchihi,acclohi,acchilo,acclolo;
   
   vhihi[0] = 1.0; vlohi[0] = 0.0;
   vhilo[0] = 0.0; vlolo[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      // sigma = sigma + x[i]*x[i];
      qdf_sqr(xhihi[i],xlohi[i],xhilo[i],xlolo[i],
           &acchihi,&acclohi,&acchilo,&acclolo);
      qdf_inc(&sigmahihi,&sigmalohi,&sigmahilo,&sigmalolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      vhihi[i] = xhihi[i]; vlohi[i] = xlohi[i];
      vhilo[i] = xhilo[i]; vlolo[i] = xlolo[i];
   }
   if((sigmahihi == 0.0) && (sigmalohi == 0.0) &&
      (sigmahilo == 0.0) && (sigmalolo == 0.0))
   {
      *betahihi = 0.0; *betalohi = 0.0;
      *betahilo = 0.0; *betalolo = 0.0;
   }
   else
   {
      // mu = sqrt(x[0]*x[0] + sigma);
      qdf_sqr(xhihi[0],xlohi[0],xhilo[0],xlolo[0],
           &acchihi,&acclohi,&acchilo,&acclolo);
      qdf_inc(&acchihi, &acclohi, &acchilo, &acclolo,
             sigmahihi,sigmalohi,sigmahilo,sigmalolo);
      qdf_sqrt(acchihi,acclohi,acchilo,acclolo,
               &muhihi,&mulohi,&muhilo,&mulolo);

      if(xhihi[0] <= 0.0)
      {
         // v[0] = x[0] - mu;
         qdf_sub(xhihi[0], xlohi[0], xhilo[0], xlolo[0],
                muhihi,   mulohi,   muhilo,   mulolo,
                &vhihi[0],&vlohi[0],&vhilo[0],&vlolo[0]);
      }
      else
      {
         // v[0] = -sigma/(x[0] + mu);
         qdf_add(xhihi[0],xlohi[0],xhilo[0],xlolo[0],
                muhihi,  mulohi,  muhilo,  mulolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_div(sigmahihi,sigmalohi,sigmahilo,sigmalolo,
                   acchihi,  acclohi,  acchilo,  acclolo,
                    &vhihi[0],&vlohi[0],&vhilo[0],&vlolo[0]);
         qdf_minus(&vhihi[0],&vlohi[0],&vhilo[0],&vlolo[0]);
      }
      // v0p2 = v[0]*v[0];
      qdf_sqr(vhihi[0], vlohi[0], vhilo[0], vlolo[0],
          &v0p2hihi,&v0p2lohi,&v0p2hilo,&v0p2lolo);
      // *beta = 2.0*v0p2/(sigma + v0p2);
      qdf_add(sigmahihi,sigmalohi,sigmahilo,sigmalolo,
               v0p2hihi, v0p2lohi, v0p2hilo, v0p2lolo,
               &acchihi, &acclohi, &acchilo, &acclolo);
      qdf_div(v0p2hihi,v0p2lohi,v0p2hilo,v0p2lolo,
               acchihi, acclohi, acchilo, acclolo,
              betahihi,betalohi,betahilo,betalolo);
      qdf_mlt_d(betahihi,betalohi,betahilo,betalolo,2.0);
      
      for(int i=1; i<n; i++) // v[i] = v[i]/v[0];
      {
         qdf_div(vhihi[i],vlohi[i],vhilo[i],vlolo[i],
                 vhihi[0],vlohi[0],vhilo[0],vlolo[0],
              &acchihi,&acclohi,&acchilo,&acclolo);

         vhihi[i] = acchihi; vlohi[i] = acclohi;
         vhilo[i] = acchilo; vlolo[i] = acclolo;
      }
      vhihi[0] = 1.0; vlohi[0] = 0.0;
      vhilo[0] = 0.0; vlolo[0] = 0.0;
   }
}

void CPU_cmplx4_factors_house 
( int n,
  double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
  double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
  double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
  double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
  double *betahihi, double *betalohi, double *betahilo, double *betalolo )
{
   double sigmahihi = 0.0;
   double sigmalohi = 0.0;
   double sigmahilo = 0.0;
   double sigmalolo = 0.0;
   double acchihi,acclohi,acchilo,acclolo;
   double muhihi,mulohi,muhilo,mulolo;
   double sqrx0hihi,sqrx0lohi,sqrx0hilo,sqrx0lolo;
   double x0radhihi,x0radlohi,x0radhilo,x0radlolo;
   double sqrv0hihi,sqrv0lohi,sqrv0hilo,sqrv0lolo;
   double inv0rehihi,inv0relohi,inv0rehilo,inv0relolo;
   double inv0imhihi,inv0imlohi,inv0imhilo,inv0imlolo;
   double zrehihi,zrelohi,zrehilo,zrelolo;
   double zimhihi,zimlohi,zimhilo,zimlolo;
   
   vrehihi[0] = 1.0; vrelohi[0] = 0.0;
   vrehilo[0] = 0.0; vrelolo[0] = 0.0;
   vimhihi[0] = 0.0; vimlohi[0] = 0.0;
   vimhilo[0] = 0.0; vimlolo[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      // sigma = sigma + xre[i]*xre[i] + xim[i]*xim[i];
      qdf_sqr(xrehihi[i],xrelohi[i],xrehilo[i],xrelolo[i],
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdf_inc(&sigmahihi,&sigmalohi,&sigmahilo,&sigmalolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      qdf_sqr(ximhihi[i],ximlohi[i],ximhilo[i],ximlolo[i],
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdf_inc(&sigmahihi,&sigmalohi,&sigmahilo,&sigmalolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      vrehihi[i] = xrehihi[i]; vrelohi[i] = xrelohi[i];
      vrehilo[i] = xrehilo[i]; vrelolo[i] = xrelolo[i];
      vimhihi[i] = ximhihi[i]; vimlohi[i] = ximlohi[i];
      vimhilo[i] = ximhilo[i]; vimlolo[i] = ximlolo[i];
   }
   if((sigmahihi == 0.0) && (sigmalohi == 0.0) &&
      (sigmahilo == 0.0) && (sigmalolo == 0.0))
   {
      *betahihi = 0.0; *betalohi = 0.0;
      *betahilo = 0.0; *betalolo = 0.0;
   }
   else
   {
      // sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      qdf_sqr(xrehihi[0],xrelohi[0],xrehilo[0],xrelolo[0],
           &sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo);
      qdf_sqr(ximhihi[0],ximlohi[0],ximhilo[0],ximlolo[0],
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdf_inc(&sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      // x0rad = sqrt(sqrx0);
      qdf_sqrt(sqrx0hihi, sqrx0lohi, sqrx0hilo, sqrx0lolo,
              &x0radhihi,&x0radlohi,&x0radhilo,&x0radlolo);
      // cout << "x0rad : " << x0radhi << "  " << x0radlo << endl;
      // mu = sqrt(sqrx0 + sigma); // norm of the vector x
      qdf_inc(&sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo,
               sigmahihi, sigmalohi, sigmahilo, sigmalolo);
      qdf_sqrt(sqrx0hihi,sqrx0lohi,sqrx0hilo,sqrx0lolo,
                 &muhihi,  &mulohi,  &muhilo,  &mulolo);

      if((x0radhihi == 0.0) && (x0radlohi == 0.0) &&
         (x0radhilo == 0.0) && (x0radlolo == 0.0))
      {
         vrehihi[0] = -muhihi; vrelohi[0] = -mulohi;
         vrehilo[0] = -muhilo; vrelolo[0] = -mulolo;
         vimhihi[0] = 0.0; vimlohi[0] = 0.0;
         vimhilo[0] = 0.0; vimlolo[0] = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         // mu = mu/x0rad;
         qdf_div(muhihi,   mulohi,   muhilo,   mulolo,
              x0radhihi,x0radlohi,x0radhilo,x0radlolo,
               &acchihi, &acclohi, &acchilo, &acclolo);
         muhihi = acchihi; mulohi = acclohi;
         muhilo = acchilo; mulolo = acclolo;
         // vre[0] = xre[0] - mu*xre[0];
         qdf_mul(muhihi,    mulohi,    muhilo,    mulolo,
                xrehihi[0],xrelohi[0],xrehilo[0],xrelolo[0],
               &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_sub(xrehihi[0], xrelohi[0], xrehilo[0], xrelolo[0],
                 acchihi,    acclohi,    acchilo,    acclolo,
                &vrehihi[0],&vrelohi[0],&vrehilo[0],&vrelolo[0]);
         // vim[0] = xim[0] - mu*xim[0];
         qdf_mul(muhihi,    mulohi,    muhilo,    mulolo,
                ximhihi[0],ximlohi[0],ximhilo[0],ximlolo[0],
               &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_sub(ximhihi[0], ximlohi[0], ximhilo[0], ximlolo[0],
                 acchihi,    acclohi,    acchilo,    acclolo,
                &vimhihi[0],&vimlohi[0],&vimhilo[0],&vimlolo[0]);
      }
      // cout << "mu : " << muhi << "  " << mulo << endl;

      // sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      qdf_sqr(vrehihi[0],vrelohi[0],vrehilo[0],vrelolo[0],
           &sqrv0hihi,&sqrv0lohi,&sqrv0hilo,&sqrv0lolo);
      // cout << "vre[0] : " << vrehi[0] << "  " << vrelo[0] << endl;
      // cout << "vim[0] : " << vimhi[0] << "  " << vimlo[0] << endl;
      qdf_sqr(vimhihi[0],vimlohi[0],vimhilo[0],vimlolo[0],
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdf_inc(&sqrv0hihi,&sqrv0lohi,&sqrv0hilo,&sqrv0lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      // cout << "sqrv0 : " << sqrv0hi << "  " << sqrv0lo << endl;
      // *beta = 2.0*sqrv0/(sigma + sqrv0);
      qdf_inc(&sigmahihi,&sigmalohi,&sigmahilo,&sigmalolo,
               sqrv0hihi, sqrv0lohi, sqrv0hilo, sqrv0lolo);
      qdf_div(sqrv0hihi,sqrv0lohi,sqrv0hilo,sqrv0lolo,
              sigmahihi,sigmalohi,sigmahilo,sigmalolo,
               betahihi, betalohi, betahilo, betalolo);
      qdf_mlt_d(betahihi,betalohi,betahilo,betalolo,2.0);
      // inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      qdf_div(vrehihi[0], vrelohi[0], vrehilo[0], vrelolo[0],
            sqrv0hihi,  sqrv0lohi,  sqrv0hilo,  sqrv0lolo,
          &inv0rehihi,&inv0relohi,&inv0rehilo,&inv0relolo);
      // inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]
      qdf_div(vimhihi[0], vimlohi[0], vimhilo[0], vimlolo[0],
            sqrv0hihi,  sqrv0lohi,  sqrv0hilo,  sqrv0lolo,
          &inv0imhihi,&inv0imlohi,&inv0imhilo,&inv0imlolo);
      qdf_minus(&inv0imhihi,&inv0imlohi,&inv0imhilo,&inv0imlolo);

      for(int i=1; i<n; i++)  // v[i] = v[i]/v[0]
      {
         // zre = vre[i]*inv0re - vim[i]*inv0im;
         qdf_mul(vrehihi[i],vrelohi[i],vrehilo[i],vrelolo[i],
              inv0rehihi,inv0relohi,inv0rehilo,inv0relolo,
                &zrehihi,  &zrelohi,  &zrehilo,  &zrelolo);
         qdf_mul(vimhihi[i],vimlohi[i],vimhilo[i],vimlolo[i],
              inv0imhihi,inv0imlohi,inv0imhilo,inv0imlolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_dec(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
                  acchihi, acclohi, acchilo, acclolo);
         // zim = vim[i]*inv0re + vre[i]*inv0im;
         qdf_mul(vimhihi[i],vimlohi[i],vimhilo[i],vimlolo[i],
              inv0rehihi,inv0relohi,inv0rehilo,inv0relolo,
                &zimhihi,  &zimlohi,  &zimhilo,  &zimlolo);
         qdf_mul(vrehihi[i],vrelohi[i],vrehilo[i],vrelolo[i],
              inv0imhihi,inv0imlohi,inv0imhilo,inv0imlolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_inc(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
                  acchihi, acclohi, acchilo, acclolo);
         vrehihi[i] = zrehihi; vrelohi[i] = zrelohi;
         vrehilo[i] = zrehilo; vrelolo[i] = zrelolo;
         vimhihi[i] = zimhihi; vimlohi[i] = zimlohi;
         vimhilo[i] = zimhilo; vimlolo[i] = zimlolo;
      }
      vrehihi[0] = 1.0; vrelohi[0] = 0.0;
      vrehilo[0] = 0.0; vrelolo[0] = 0.0;
      vimhihi[0] = 0.0; vimlohi[0] = 0.0;
      vimhilo[0] = 0.0; vimlolo[0] = 0.0;
   }
}

void CPU_dbl4_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double betahihi, double betalohi, double betahilo, double betalolo )
{
   double *whihi = new double[ncols-k];
   double *wlohi = new double[ncols-k];
   double *whilo = new double[ncols-k];
   double *wlolo = new double[ncols-k];
   double acchihi,acclohi,acchilo,acclolo;

   for(int j=k; j<ncols; j++)
   {
      whihi[j-k] = 0.0; wlohi[j-k] = 0.0;
      whilo[j-k] = 0.0; wlolo[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         qdf_mul(Rhihi[i][j],Rlohi[i][j],Rhilo[i][j],Rlolo[i][j],
                 vhihi[i-k], vlohi[i-k], vhilo[i-k], vlolo[i-k],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_inc(&whihi[j-k],&wlohi[j-k],&whilo[j-k],&wlolo[j-k],
                acchihi,    acclohi,    acchilo,    acclolo);
      }
      // w[j-k] = beta*w[j-k];
      qdf_mul(betahihi,  betalohi,  betahilo,  betalolo,
                 whihi[j-k],wlohi[j-k],whilo[j-k],wlolo[j-k],
              &acchihi,  &acclohi,  &acchilo,  &acclolo);

      whihi[j-k] = acchihi; wlohi[j-k] = acclohi;
      whilo[j-k] = acchilo; wlolo[j-k] = acclolo;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         qdf_mul(vhihi[i-k],vlohi[i-k],vhilo[i-k],vlolo[i-k],
                 whihi[j-k],wlohi[j-k],whilo[j-k],wlolo[j-k],
              &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_dec(&Rhihi[i][j],&Rlohi[i][j],&Rhilo[i][j],&Rlolo[i][j],
                acchihi,     acclohi,     acchilo,     acclolo);
      }

   free(whihi); free(wlohi); free(whilo); free(wlolo);
}

void CPU_cmplx4_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double betahihi, double betalohi, double betahilo, double betalolo )
{
   double *wrehihi = new double[ncols-k];
   double *wrelohi = new double[ncols-k];
   double *wrehilo = new double[ncols-k];
   double *wrelolo = new double[ncols-k];
   double *wimhihi = new double[ncols-k];
   double *wimlohi = new double[ncols-k];
   double *wimhilo = new double[ncols-k];
   double *wimlolo = new double[ncols-k];
   double zrehihi,zrelohi,zrehilo,zrelolo;
   double zimhihi,zimlohi,zimhilo,zimlolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int j=k; j<ncols; j++)
   {
      wrehihi[j-k] = 0.0; wrelohi[j-k] = 0.0;
      wrehilo[j-k] = 0.0; wrelolo[j-k] = 0.0;
      wimhihi[j-k] = 0.0; wimlohi[j-k] = 0.0;
      wimhilo[j-k] = 0.0; wimlolo[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         // Hermitian of R => flip sign of Rim
         // zre =   Rre[i][j]*vre[i-k] + Rim[i][j]*vim[i-k];
         qdf_mul(Rrehihi[i][j],Rrelohi[i][j],Rrehilo[i][j],Rrelolo[i][j],
                 vrehihi[i-k], vrelohi[i-k], vrehilo[i-k], vrelolo[i-k],
                &zrehihi,     &zrelohi,     &zrehilo,     &zrelolo);
         qdf_mul(Rimhihi[i][j],Rimlohi[i][j],Rimhilo[i][j],Rimlolo[i][j],
                 vimhihi[i-k], vimlohi[i-k], vimhilo[i-k], vimlolo[i-k],
                &acchihi,     &acclohi,     &acchilo,     &acclolo);
         qdf_inc(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
                  acchihi, acclohi, acchilo, acclolo);
         // zim = - Rim[i][j]*vre[i-k] + Rre[i][j]*vim[i-k];
         qdf_mul(Rrehihi[i][j],Rrelohi[i][j],Rrehilo[i][j],Rrelolo[i][j],
                 vimhihi[i-k], vimlohi[i-k], vimhilo[i-k], vimlolo[i-k],
                &zimhihi,     &zimlohi,     &zimhilo,     &zimlolo);
         qdf_mul(Rimhihi[i][j],Rimlohi[i][j],Rimhilo[i][j],Rimlolo[i][j],
                 vrehihi[i-k], vrelohi[i-k], vrehilo[i-k], vrelolo[i-k],
                &acchihi,     &acclohi,     &acchilo,     &acclolo);
         qdf_dec(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
                  acchihi, acclohi, acchilo, acclolo);
         // wre[j-k] = wre[j-k] + zre;
         qdf_inc(&wrehihi[j-k],&wrelohi[j-k],&wrehilo[j-k],&wrelolo[j-k],
                  zrehihi,      zrelohi,      zrehilo,      zrelolo);
         // wim[j-k] = wim[j-k] + zim;
         qdf_inc(&wimhihi[j-k],&wimlohi[j-k],&wimhilo[j-k],&wimlolo[j-k],
                  zimhihi,      zimlohi,      zimhilo,      zimlolo);
      }
      // wre[j-k] = beta*wre[j-k];
      qdf_mul(betahihi,    betalohi,    betahilo,    betalolo,
               wrehihi[j-k],wrelohi[j-k],wrehilo[j-k],wrelolo[j-k],
              &acchihi,    &acclohi,    &acchilo,    &acclolo);
      wrehihi[j-k] = acchihi; wrelohi[j-k] = acclohi;
      wrehilo[j-k] = acchilo; wrelolo[j-k] = acclolo;
      // wim[j-k] = beta*wim[j-k];
      qdf_mul(betahihi,    betalohi,    betahilo,    betalolo,
               wimhihi[j-k],wimlohi[j-k],wimhilo[j-k],wimlolo[j-k],
              &acchihi,    &acclohi,    &acchilo,    &acclolo);
      wimhihi[j-k] = acchihi; wimlohi[j-k] = acclohi;
      wimhilo[j-k] = acchilo; wimlolo[j-k] = acclolo;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         // Hermitian of w => flip sign of wim
         // zre = vre[i-k]*wre[j-k] + vim[i-k]*wim[j-k];
         qdf_mul(vrehihi[i-k],vrelohi[i-k],vrehilo[i-k],vrelolo[i-k],
                 wrehihi[j-k],wrelohi[j-k],wrehilo[j-k],wrelolo[j-k],
                &zrehihi,    &zrelohi,    &zrehilo,    &zrelolo);
         qdf_mul(vimhihi[i-k],vimlohi[i-k],vimhilo[i-k],vimlolo[i-k],
                 wimhihi[j-k],wimlohi[j-k],wimhilo[j-k],wimlolo[j-k],
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
         qdf_inc(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
                  acchihi, acclohi, acchilo, acclolo);
         // zim = vim[i-k]*wre[j-k] - vre[i-k]*wim[j-k];
         qdf_mul(vimhihi[i-k],vimlohi[i-k],vimhilo[i-k],vimlolo[i-k],
                 wrehihi[j-k],wrelohi[j-k],wrehilo[j-k],wrelolo[j-k],
                &zimhihi,    &zimlohi,    &zimhilo,    &zimlolo);
         qdf_mul(vrehihi[i-k],vrelohi[i-k],vrehilo[i-k],vrelolo[i-k],
                 wimhihi[j-k],wimlohi[j-k],wimhilo[j-k],wimlolo[j-k],
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
         qdf_dec(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
                  acchihi, acclohi, acchilo, acclolo);
         // Rre[i][j] = Rre[i][j] - zre;
         qdf_dec(&Rrehihi[i][j],&Rrelohi[i][j],&Rrehilo[i][j],&Rrelolo[i][j],
                  zrehihi,       zrelohi,       zrehilo,       zrelolo);
         // Rim[i][j] = Rim[i][j] - zim;
         qdf_dec(&Rimhihi[i][j],&Rimlohi[i][j],&Rimhilo[i][j],&Rimlolo[i][j],
                  zimhihi,       zimlohi,       zimhilo,       zimlolo);
      }

   free(wrehihi); free(wrelohi); free(wrehilo); free(wrelolo);
   free(wimhihi); free(wimlohi); free(wimhilo); free(wimlolo);
}

void CPU_dbl4_factors_rightQupdate
 ( int n, int k,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double betahihi, double betalohi, double betahilo, double betalolo )
{
   double *whihi = new double[n];
   double *wlohi = new double[n];
   double *whilo = new double[n];
   double *wlolo = new double[n];
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<n; i++)
   {
      whihi[i] = 0.0; wlohi[i] = 0.0;
      whilo[i] = 0.0; wlolo[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         qdf_mul(Qhihi[i][j],Qlohi[i][j],Qhilo[i][j],Qlolo[i][j],
                 vhihi[j-k], vlohi[j-k], vhilo[j-k], vlolo[j-k],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_inc(&whihi[i],&wlohi[i],&whilo[i],&wlolo[i],
                acchihi,  acclohi,  acchilo,  acclolo);
      }
      // w[i] = beta*w[i];
      qdf_mul(betahihi,betalohi,betahilo,betalolo,
                 whihi[i],wlohi[i],whilo[i],wlolo[i],
              &acchihi,&acclohi,&acchilo,&acclolo);
      whihi[i] = acchihi; wlohi[i] = acclohi;
      whilo[i] = acchilo; wlolo[i] = acclolo;
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         qdf_mul(whihi[i],  wlohi[i],  whilo[i],  wlolo[i],
                 vhihi[j-k],vlohi[j-k],vhilo[j-k],vlolo[j-k],
              &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_dec(&Qhihi[i][j],&Qlohi[i][j],&Qhilo[i][j],&Qlolo[i][j],
                acchihi,     acclohi,     acchilo,     acclolo);
      }

   free(whihi); free(wlohi); free(whilo); free(wlolo);
}

void CPU_cmplx4_factors_rightQupdate
 ( int n, int k,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double betahihi, double betalohi, double betahilo, double betalolo )
{
   double *wrehihi = new double[n];
   double *wrelohi = new double[n];
   double *wrehilo = new double[n];
   double *wrelolo = new double[n];
   double *wimhihi = new double[n];
   double *wimlohi = new double[n];
   double *wimhilo = new double[n];
   double *wimlolo = new double[n];
   double zrehihi,zrelohi,zrehilo,zrelolo;
   double zimhihi,zimlohi,zimhilo,zimlolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<n; i++)
   {
      wrehihi[i] = 0.0; wrelohi[i] = 0.0;
      wrehilo[i] = 0.0; wrelolo[i] = 0.0;
      wimhihi[i] = 0.0; wimlohi[i] = 0.0;
      wimhilo[i] = 0.0; wimlolo[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         // zre = Qre[i][j]*vre[j-k] - Qim[i][j]*vim[j-k];
         qdf_mul(Qrehihi[i][j],Qrelohi[i][j],Qrehilo[i][j],Qrelolo[i][j],
                 vrehihi[j-k], vrelohi[j-k], vrehilo[j-k], vrelolo[j-k],
                &zrehihi,     &zrelohi,     &zrehilo,     &zrelolo);
         qdf_mul(Qimhihi[i][j],Qimlohi[i][j],Qimhilo[i][j],Qimlolo[i][j],
                 vimhihi[j-k], vimlohi[j-k], vimhilo[j-k], vimlolo[j-k],
                &acchihi,     &acclohi,     &acchilo,     &acclolo);
         qdf_dec(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
                  acchihi, acclohi, acchilo, acclolo);
         // zim = Qim[i][j]*vre[j-k] + Qre[i][j]*vim[j-k];
         qdf_mul(Qimhihi[i][j],Qimlohi[i][j],Qimhilo[i][j],Qimlolo[i][j],
                 vrehihi[j-k], vrelohi[j-k], vrehilo[j-k], vrelolo[j-k],
                &zimhihi,     &zimlohi,     &zimhilo,     &zimlolo);
         qdf_mul(Qrehihi[i][j],Qrelohi[i][j],Qrehilo[i][j],Qrelolo[i][j],
                 vimhihi[j-k], vimlohi[j-k], vimhilo[j-k], vimlolo[j-k],
                &acchihi,     &acclohi,     &acchilo,     &acclolo);
         qdf_inc(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
                  acchihi, acclohi, acchilo, acclolo);
         // wre[i] = wre[i] + zre;
         qdf_inc(&wrehihi[i],&wrelohi[i],&wrehilo[i],&wrelolo[i],
                  zrehihi,    zrelohi,    zrehilo,    zrelolo);
         // wim[i] = wim[i] + zim;
         qdf_inc(&wimhihi[i],&wimlohi[i],&wimhilo[i],&wimlolo[i],
                  zimhihi,    zimlohi,    zimhilo,    zimlolo);
      }
      // wre[i] = beta*wre[i];
      qdf_mul(betahihi,  betalohi,  betahilo,  betalolo,
               wrehihi[i],wrelohi[i],wrehilo[i],wrelolo[i],
              &acchihi,  &acclohi,  &acchilo,  &acclolo);
      wrehihi[i] = acchihi; wrelohi[i] = acclohi;
      wrehilo[i] = acchilo; wrelolo[i] = acclolo;
      // wim[i] = beta*wim[i];
      qdf_mul(betahihi,  betalohi,  betahilo,  betalolo,
               wimhihi[i],wimlohi[i],wimhilo[i],wimlolo[i],
              &acchihi,  &acclohi,  &acchilo  ,&acclolo);
      wimhihi[i] = acchihi; wimlohi[i] = acclohi;
      wimhilo[i] = acchilo; wimlolo[i] = acclolo;
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         // Hermitian transpose => flip sign of vim
         // zre = wre[i]*vre[j-k] + wim[i]*vim[j-k];
         qdf_mul(wrehihi[i],  wrelohi[i],  wrehilo[i],  wrelolo[i],
                 vrehihi[j-k],vrelohi[j-k],vrehilo[j-k],vrelolo[j-k],
                 &zrehihi,   &zrelohi,    &zrehilo,    &zrelolo);
         qdf_mul(wimhihi[i],  wimlohi[i],  wimhilo[i],  wimlolo[i],
                 vimhihi[j-k],vimlohi[j-k],vimhilo[j-k],vimlolo[j-k],
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
         qdf_inc(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
                  acchihi, acclohi, acchilo, acclolo);
         // zim = wim[i]*vre[j-k] - wre[i]*vim[j-k];
         qdf_mul(wimhihi[i],  wimlohi[i],  wimhilo[i],  wimlolo[i],
                 vrehihi[j-k],vrelohi[j-k],vrehilo[j-k],vrelolo[j-k],
                &zimhihi,    &zimlohi,    &zimhilo,    &zimlolo);
         qdf_mul(wrehihi[i],  wrelohi[i],  wrehilo[i],  wrelolo[i],
                 vimhihi[j-k],vimlohi[j-k],vimhilo[j-k],vimlolo[j-k],
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
         qdf_dec(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
                  acchihi, acclohi, acchilo, acclolo);
         // Qre[i][j] = Qre[i][j] - zre;
         qdf_dec(&Qrehihi[i][j],&Qrelohi[i][j],&Qrehilo[i][j],&Qrelolo[i][j],
                  zrehihi,       zrelohi,       zrehilo,       zrelolo);
         // Qim[i][j] = Qim[i][j] - zim;
         qdf_dec(&Qimhihi[i][j],&Qimlohi[i][j],&Qimhilo[i][j],&Qimlolo[i][j],
                  zimhihi,       zimlohi,       zimhilo,       zimlolo);
      }

   free(wrehihi); free(wrelohi); free(wrehilo); free(wrelolo);
   free(wimhihi); free(wimlohi); free(wimhilo); free(wimlolo);
}

void CPU_dbl4_factors_houseqr
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo )
{
   double *xhihi = new double[nrows];
   double *xlohi = new double[nrows];
   double *xhilo = new double[nrows];
   double *xlolo = new double[nrows];
   double *vhihi = new double[nrows];
   double *vlohi = new double[nrows];
   double *vhilo = new double[nrows];
   double *vlolo = new double[nrows];
   double betahihi,betalohi,betahilo,betalolo;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qhihi[i][j] = 0.0; Qlohi[i][j] = 0.0;
         Qhilo[i][j] = 0.0; Qlolo[i][j] = 0.0;
      }
      Qhihi[i][i] = 1.0; Qlohi[i][i] = 0.0;
      Qhilo[i][i] = 0.0; Qlolo[i][i] = 0.0;

      for(int j=0; j<ncols; j++)
      {
         Rhihi[i][j] = Ahihi[i][j];
         Rlohi[i][j] = Alohi[i][j];
         Rhilo[i][j] = Ahilo[i][j];
         Rlolo[i][j] = Alolo[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xhihi[i-k] = Rhihi[i][k];
            xlohi[i-k] = Rlohi[i][k];
            xhilo[i-k] = Rhilo[i][k];
            xlolo[i-k] = Rlolo[i][k];
         }
         CPU_dbl4_factors_house
            (nrows-k,xhihi,xlohi,xhilo,xlolo,vhihi,vlohi,vhilo,vlolo,
             &betahihi,&betalohi,&betahilo,&betalolo);
         CPU_dbl4_factors_leftRupdate
            (nrows,ncols,k,Rhihi,Rlohi,Rhilo,Rlolo,vhihi,vlohi,vhilo,vlolo,
             betahihi,betalohi,betahilo,betalolo);
         CPU_dbl4_factors_rightQupdate
            (nrows,k,Qhihi,Qlohi,Qhilo,Qlolo,vhihi,vlohi,vhilo,vlolo,
             betahihi,betalohi,betahilo,betalolo);
      }
   }
   free(xhihi); free(xlohi); free(xhilo); free(xlolo);
   free(vhihi); free(vlohi); free(vhilo); free(vlolo);
}

void CPU_cmplx4_factors_houseqr
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo )
{
   double *xrehihi = new double[nrows];
   double *xrelohi = new double[nrows];
   double *xrehilo = new double[nrows];
   double *xrelolo = new double[nrows];
   double *ximhihi = new double[nrows];
   double *ximlohi = new double[nrows];
   double *ximhilo = new double[nrows];
   double *ximlolo = new double[nrows];
   double *vrehihi = new double[nrows];
   double *vrelohi = new double[nrows];
   double *vrehilo = new double[nrows];
   double *vrelolo = new double[nrows];
   double *vimhihi = new double[nrows];
   double *vimlohi = new double[nrows];
   double *vimhilo = new double[nrows];
   double *vimlolo = new double[nrows];
   double betahihi,betalohi,betahilo,betalolo;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qrehihi[i][j] = 0.0; Qrelohi[i][j] = 0.0;
         Qrehilo[i][j] = 0.0; Qrelolo[i][j] = 0.0;
         Qimhihi[i][j] = 0.0; Qimlohi[i][j] = 0.0;
         Qimhilo[i][j] = 0.0; Qimlolo[i][j] = 0.0;
      }
      Qrehihi[i][i] = 1.0; Qrelohi[i][i] = 0.0;
      Qrehilo[i][i] = 0.0; Qrelolo[i][i] = 0.0;
      Qimhihi[i][i] = 0.0; Qimlohi[i][i] = 0.0;
      Qimhilo[i][i] = 0.0; Qimlolo[i][i] = 0.0;

      for(int j=0; j<ncols; j++)
      {
         Rrehihi[i][j] = Arehihi[i][j];
         Rrelohi[i][j] = Arelohi[i][j];
         Rrehilo[i][j] = Arehilo[i][j];
         Rrelolo[i][j] = Arelolo[i][j];
         Rimhihi[i][j] = Aimhihi[i][j];
         Rimlohi[i][j] = Aimlohi[i][j];
         Rimhilo[i][j] = Aimhilo[i][j];
         Rimlolo[i][j] = Aimlolo[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xrehihi[i-k] = Rrehihi[i][k];
            xrelohi[i-k] = Rrelohi[i][k];
            xrehilo[i-k] = Rrehilo[i][k];
            xrelolo[i-k] = Rrelolo[i][k];
            ximhihi[i-k] = Rimhihi[i][k];
            ximlohi[i-k] = Rimlohi[i][k];
            ximhilo[i-k] = Rimhilo[i][k];
            ximlolo[i-k] = Rimlolo[i][k];
         }
         CPU_cmplx4_factors_house
            (nrows-k,
             xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
             vrehihi,vrelohi,vrehilo,vrelolo,vimhihi,vimlohi,vimhilo,vimlolo,
             &betahihi,&betalohi,&betahilo,&betalolo);
         CPU_cmplx4_factors_leftRupdate
            (nrows,ncols,k,
             Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
             vrehihi,vrelohi,vrehilo,vrelolo,vimhihi,vimlohi,vimhilo,vimlolo,
             betahihi,betalohi,betahilo,betalolo);
         CPU_cmplx4_factors_rightQupdate
            (nrows,k,
             Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
             vrehihi,vrelohi,vrehilo,vrelolo,vimhihi,vimlohi,vimhilo,vimlolo,
             betahihi,betalohi,betahilo,betalolo);
      }
   }
   free(xrehihi); free(xrelohi); free(xrehilo); free(xrelolo);
   free(ximhihi); free(ximlohi); free(ximhilo); free(ximlolo);
   free(vrehihi); free(vrelohi); free(vrehilo); free(vrelolo);
   free(vimhihi); free(vimlohi); free(vimhilo); free(vimlolo);
}
*/
