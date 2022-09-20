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
      xrehihilo[i] = brehihilo[i]; xrelohilo[i] = brelohilo[i];
      xrehilolo[i] = brehilolo[i]; xrelololo[i] = brelololo[i];
      ximhihihi[i] = bimhihihi[i]; ximlohihi[i] = bimlohihi[i];
      ximhilohi[i] = bimhilohi[i]; ximlolohi[i] = bimlolohi[i];
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
      xrehilohi[i] = acc2rehilohi; xrelolohi[i] = acc2relolohi;
      xrehihilo[i] = acc2rehihilo; xrelohilo[i] = acc2relohilo;
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
            valtmp = Aimhihilo[idxmax][k];
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
         Arehilohi[i][j] = acc3rehilohi; Arelolohi[i][j] = acc3relolohi;
         Arehihilo[i][j] = acc3rehihilo; Arelohilo[i][j] = acc3relohilo;
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

void CPU_dbl8_factors_house
 ( int n,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
   double sigmahihihi = 0.0;
   double sigmalohihi = 0.0;
   double sigmahilohi = 0.0;
   double sigmalolohi = 0.0;
   double sigmahihilo = 0.0;
   double sigmalohilo = 0.0;
   double sigmahilolo = 0.0;
   double sigmalololo = 0.0;
   double muhihihi,mulohihi,muhilohi,mulolohi;
   double muhihilo,mulohilo,muhilolo,mulololo;
   double v0p2hihihi,v0p2lohihi,v0p2hilohi,v0p2lolohi;
   double v0p2hihilo,v0p2lohilo,v0p2hilolo,v0p2lololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   
   vhihihi[0] = 1.0; vlohihi[0] = 0.0;
   vhilohi[0] = 0.0; vlolohi[0] = 0.0;
   vhihilo[0] = 0.0; vlohilo[0] = 0.0;
   vhilolo[0] = 0.0; vlololo[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      // sigma = sigma + x[i]*x[i];
      odf_sqr(xhihihi[i],xlohihi[i],xhilohi[i],xlolohi[i],
              xhihilo[i],xlohilo[i],xhilolo[i],xlololo[i],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);

      odf_inc(&sigmahihihi,&sigmalohihi,&sigmahilohi,&sigmalolohi,
              &sigmahihilo,&sigmalohilo,&sigmahilolo,&sigmalololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);

      vhihihi[i] = xhihihi[i]; vlohihi[i] = xlohihi[i];
      vhilohi[i] = xhilohi[i]; vlolohi[i] = xlolohi[i];
      vhihilo[i] = xhihilo[i]; vlohilo[i] = xlohilo[i];
      vhilolo[i] = xhilolo[i]; vlololo[i] = xlololo[i];
   }
   if((sigmahihihi == 0.0) && (sigmalohihi == 0.0) &&
      (sigmahilohi == 0.0) && (sigmalolohi == 0.0) &&
      (sigmahihilo == 0.0) && (sigmalohilo == 0.0) &&
      (sigmahilolo == 0.0) && (sigmalololo == 0.0))
   {
      *betahihihi = 0.0; *betalohihi = 0.0;
      *betahilohi = 0.0; *betalolohi = 0.0;
      *betahihilo = 0.0; *betalohilo = 0.0;
      *betahilolo = 0.0; *betalololo = 0.0;
   }
   else
   {
      // mu = sqrt(x[0]*x[0] + sigma);
      odf_sqr(xhihihi[0],xlohihi[0],xhilohi[0],xlolohi[0],
              xhihilo[0],xlohilo[0],xhilolo[0],xlololo[0],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);

      odf_inc(&acchihihi, &acclohihi, &acchilohi, &acclolohi,
              &acchihilo, &acclohilo, &acchilolo, &acclololo,
             sigmahihihi,sigmalohihi,sigmahilohi,sigmalolohi,
             sigmahihilo,sigmalohilo,sigmahilolo,sigmalololo);

      odf_sqrt(acchihihi,acclohihi,acchilohi,acclolohi,
               acchihilo,acclohilo,acchilolo,acclololo,
               &muhihihi,&mulohihi,&muhilohi,&mulolohi,
               &muhihilo,&mulohilo,&muhilolo,&mulololo);

      if(xhihihi[0] <= 0.0)
      {
         // v[0] = x[0] - mu;
         odf_sub(xhihihi[0], xlohihi[0], xhilohi[0], xlolohi[0],
                 xhihilo[0], xlohilo[0], xhilolo[0], xlololo[0],
                muhihihi,   mulohihi,   muhilohi,   mulolohi,
                muhihilo,   mulohilo,   muhilolo,   mulololo,
                &vhihihi[0],&vlohihi[0],&vhilohi[0],&vlolohi[0],
                &vhihilo[0],&vlohilo[0],&vhilolo[0],&vlololo[0]);
      }
      else
      {
         // v[0] = -sigma/(x[0] + mu);
         odf_add(xhihihi[0],xlohihi[0],xhilohi[0],xlolohi[0],
                 xhihilo[0],xlohilo[0],xhilolo[0],xlololo[0],
                muhihihi,  mulohihi,  muhilohi,  mulolohi,
                muhihilo,  mulohilo,  muhilolo,  mulololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);

         odf_div(sigmahihihi,sigmalohihi,sigmahilohi,sigmalolohi,
                 sigmahihilo,sigmalohilo,sigmahilolo,sigmalololo,
                   acchihihi,  acclohihi,  acchilohi,  acclolohi,
                   acchihilo,  acclohilo,  acchilolo,  acclololo,
                    &vhihihi[0],&vlohihi[0],&vhilohi[0],&vlolohi[0],
                    &vhihilo[0],&vlohilo[0],&vhilolo[0],&vlololo[0]);

         odf_minus(&vhihihi[0],&vlohihi[0],&vhilohi[0],&vlolohi[0],
                   &vhihilo[0],&vlohilo[0],&vhilolo[0],&vlololo[0]);
      }
      // v0p2 = v[0]*v[0];
      odf_sqr(vhihihi[0], vlohihi[0], vhilohi[0], vlolohi[0],
              vhihilo[0], vlohilo[0], vhilolo[0], vlololo[0],
          &v0p2hihihi,&v0p2lohihi,&v0p2hilohi,&v0p2lolohi,
          &v0p2hihilo,&v0p2lohilo,&v0p2hilolo,&v0p2lololo);
      // *beta = 2.0*v0p2/(sigma + v0p2);
      odf_add(sigmahihihi,sigmalohihi,sigmahilohi,sigmalolohi,
              sigmahihilo,sigmalohilo,sigmahilolo,sigmalololo,
               v0p2hihihi, v0p2lohihi, v0p2hilohi, v0p2lolohi,
               v0p2hihilo, v0p2lohilo, v0p2hilolo, v0p2lololo,
               &acchihihi, &acclohihi, &acchilohi, &acclolohi,
               &acchihilo, &acclohilo, &acchilolo, &acclololo);

      odf_div(v0p2hihihi,v0p2lohihi,v0p2hilohi,v0p2lolohi,
              v0p2hihilo,v0p2lohilo,v0p2hilolo,v0p2lololo,
               acchihihi, acclohihi, acchilohi, acclolohi,
               acchihilo, acclohilo, acchilolo, acclololo,
              betahihihi,betalohihi,betahilohi,betalolohi,
              betahihilo,betalohilo,betahilolo,betalololo);

      odf_mlt_d(betahihihi,betalohihi,betahilohi,betalolohi,
                betahihilo,betalohilo,betahilolo,betalololo,2.0);
      
      for(int i=1; i<n; i++) // v[i] = v[i]/v[0];
      {
         odf_div(vhihihi[i],vlohihi[i],vhilohi[i],vlolohi[i],
                 vhihilo[i],vlohilo[i],vhilolo[i],vlololo[i],
                 vhihihi[0],vlohihi[0],vhilohi[0],vlolohi[0],
                 vhihilo[0],vlohilo[0],vhilolo[0],vlololo[0],
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);

         vhihihi[i] = acchihihi; vlohihi[i] = acclohihi;
         vhilohi[i] = acchilohi; vlolohi[i] = acclolohi;
         vhihilo[i] = acchihilo; vlohilo[i] = acclohilo;
         vhilolo[i] = acchilolo; vlololo[i] = acclololo;
      }
      vhihihi[0] = 1.0; vlohihi[0] = 0.0;
      vhilohi[0] = 0.0; vlolohi[0] = 0.0;
      vhihilo[0] = 0.0; vlohilo[0] = 0.0;
      vhilolo[0] = 0.0; vlololo[0] = 0.0;
   }
}

void CPU_cmplx8_factors_house 
 ( int n,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
   double sigmahihihi = 0.0;
   double sigmalohihi = 0.0;
   double sigmahilohi = 0.0;
   double sigmalolohi = 0.0;
   double sigmahihilo = 0.0;
   double sigmalohilo = 0.0;
   double sigmahilolo = 0.0;
   double sigmalololo = 0.0;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double muhihihi,mulohihi,muhilohi,mulolohi;
   double muhihilo,mulohilo,muhilolo,mulololo;
   double sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi;
   double sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo;
   double x0radhihihi,x0radlohihi,x0radhilohi,x0radlolohi;
   double x0radhihilo,x0radlohilo,x0radhilolo,x0radlololo;
   double sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi;
   double sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo;
   double inv0rehihihi,inv0relohihi,inv0rehilohi,inv0relolohi;
   double inv0rehihilo,inv0relohilo,inv0rehilolo,inv0relololo;
   double inv0imhihihi,inv0imlohihi,inv0imhilohi,inv0imlolohi;
   double inv0imhihilo,inv0imlohilo,inv0imhilolo,inv0imlololo;
   double zrehihihi,zrelohihi,zrehilohi,zrelolohi;
   double zrehihilo,zrelohilo,zrehilolo,zrelololo;
   double zimhihihi,zimlohihi,zimhilohi,zimlolohi;
   double zimhihilo,zimlohilo,zimhilolo,zimlololo;
   
   vrehihihi[0] = 1.0; vrelohihi[0] = 0.0;
   vrehilohi[0] = 0.0; vrelolohi[0] = 0.0;
   vrehihilo[0] = 0.0; vrelohilo[0] = 0.0;
   vrehilolo[0] = 0.0; vrelololo[0] = 0.0;
   vimhihihi[0] = 0.0; vimlohihi[0] = 0.0;
   vimhilohi[0] = 0.0; vimlolohi[0] = 0.0;
   vimhihilo[0] = 0.0; vimlohilo[0] = 0.0;
   vimhilolo[0] = 0.0; vimlololo[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      // sigma = sigma + xre[i]*xre[i] + xim[i]*xim[i];
      odf_sqr(xrehihihi[i],xrelohihi[i],xrehilohi[i],xrelolohi[i],
              xrehihilo[i],xrelohilo[i],xrehilolo[i],xrelololo[i],
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odf_inc(&sigmahihihi,&sigmalohihi,&sigmahilohi,&sigmalolohi,
              &sigmahihilo,&sigmalohilo,&sigmahilolo,&sigmalololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      odf_sqr(ximhihihi[i],ximlohihi[i],ximhilohi[i],ximlolohi[i],
              ximhihilo[i],ximlohilo[i],ximhilolo[i],ximlololo[i],
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odf_inc(&sigmahihihi,&sigmalohihi,&sigmahilohi,&sigmalolohi,
              &sigmahihilo,&sigmalohilo,&sigmahilolo,&sigmalololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      vrehihihi[i] = xrehihihi[i]; vrelohihi[i] = xrelohihi[i];
      vrehilohi[i] = xrehilohi[i]; vrelolohi[i] = xrelolohi[i];
      vrehihilo[i] = xrehihilo[i]; vrelohilo[i] = xrelohilo[i];
      vrehilolo[i] = xrehilolo[i]; vrelololo[i] = xrelololo[i];
      vimhihihi[i] = ximhihihi[i]; vimlohihi[i] = ximlohihi[i];
      vimhilohi[i] = ximhilohi[i]; vimlolohi[i] = ximlolohi[i];
      vimhihilo[i] = ximhihilo[i]; vimlohilo[i] = ximlohilo[i];
      vimhilolo[i] = ximhilolo[i]; vimlololo[i] = ximlololo[i];
   }
   if((sigmahihihi == 0.0) && (sigmalohihi == 0.0) &&
      (sigmahilohi == 0.0) && (sigmalolohi == 0.0) &&
      (sigmahihilo == 0.0) && (sigmalohilo == 0.0) &&
      (sigmahilolo == 0.0) && (sigmalololo == 0.0))
   {
      *betahihihi = 0.0; *betalohihi = 0.0;
      *betahilohi = 0.0; *betalolohi = 0.0;
      *betahihilo = 0.0; *betalohilo = 0.0;
      *betahilolo = 0.0; *betalololo = 0.0;
   }
   else
   {
      // sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      odf_sqr(xrehihihi[0],xrelohihi[0],xrehilohi[0],xrelolohi[0],
              xrehihilo[0],xrelohilo[0],xrehilolo[0],xrelololo[0],
           &sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
           &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo);
      odf_sqr(ximhihihi[0],ximlohihi[0],ximhilohi[0],ximlolohi[0],
              ximhihilo[0],ximlohilo[0],ximhilolo[0],ximlololo[0],
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odf_inc(&sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      // x0rad = sqrt(sqrx0);
      odf_sqrt(sqrx0hihihi, sqrx0lohihi, sqrx0hilohi, sqrx0lolohi,
               sqrx0hihilo, sqrx0lohilo, sqrx0hilolo, sqrx0lololo,
              &x0radhihihi,&x0radlohihi,&x0radhilohi,&x0radlolohi,
              &x0radhihilo,&x0radlohilo,&x0radhilolo,&x0radlololo);
      // cout << "x0rad : " << x0radhi << "  " << x0radlo << endl;
      // mu = sqrt(sqrx0 + sigma); // norm of the vector x
      odf_inc(&sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo,
               sigmahihihi, sigmalohihi, sigmahilohi, sigmalolohi,
               sigmahihilo, sigmalohilo, sigmahilolo, sigmalololo);
      odf_sqrt(sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi,
               sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo,
                 &muhihihi,  &mulohihi,  &muhilohi,  &mulolohi,
                 &muhihilo,  &mulohilo,  &muhilolo,  &mulololo);

      if((x0radhihihi == 0.0) && (x0radlohihi == 0.0) &&
         (x0radhilohi == 0.0) && (x0radlolohi == 0.0) &&
         (x0radhihilo == 0.0) && (x0radlohilo == 0.0) &&
         (x0radhilolo == 0.0) && (x0radlololo == 0.0))
      {
         vrehihihi[0] = -muhihihi; vrelohihi[0] = -mulohihi;
         vrehilohi[0] = -muhilohi; vrelolohi[0] = -mulolohi;
         vrehihilo[0] = -muhihilo; vrelohilo[0] = -mulohilo;
         vrehilolo[0] = -muhilolo; vrelololo[0] = -mulololo;
         vimhihihi[0] = 0.0; vimlohihi[0] = 0.0;
         vimhilohi[0] = 0.0; vimlolohi[0] = 0.0;
         vimhihilo[0] = 0.0; vimlohilo[0] = 0.0;
         vimhilolo[0] = 0.0; vimlololo[0] = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         // mu = mu/x0rad;
         odf_div(muhihihi,   mulohihi,   muhilohi,   mulolohi,
                 muhihilo,   mulohilo,   muhilolo,   mulololo,
              x0radhihihi,x0radlohihi,x0radhilohi,x0radlolohi,
              x0radhihilo,x0radlohilo,x0radhilolo,x0radlololo,
               &acchihihi, &acclohihi, &acchilohi, &acclolohi,
               &acchihilo, &acclohilo, &acchilolo, &acclololo);
         muhihihi = acchihihi; mulohihi = acclohihi;
         muhilohi = acchilohi; mulolohi = acclolohi;
         muhihilo = acchihilo; mulohilo = acclohilo;
         muhilolo = acchilolo; mulololo = acclololo;
         // vre[0] = xre[0] - mu*xre[0];
         odf_mul(muhihihi,    mulohihi,    muhilohi,    mulolohi,
                 muhihilo,    mulohilo,    muhilolo,    mulololo,
                xrehihihi[0],xrelohihi[0],xrehilohi[0],xrelolohi[0],
                xrehihilo[0],xrelohilo[0],xrehilolo[0],xrelololo[0],
               &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
               &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_sub(xrehihihi[0], xrelohihi[0], xrehilohi[0], xrelolohi[0],
                 xrehihilo[0], xrelohilo[0], xrehilolo[0], xrelololo[0],
                 acchihihi,    acclohihi,    acchilohi,    acclolohi,
                 acchihilo,    acclohilo,    acchilolo,    acclololo,
                &vrehihihi[0],&vrelohihi[0],&vrehilohi[0],&vrelolohi[0],
                &vrehihilo[0],&vrelohilo[0],&vrehilolo[0],&vrelololo[0]);
         // vim[0] = xim[0] - mu*xim[0];
         odf_mul(muhihihi,    mulohihi,    muhilohi,    mulolohi,
                 muhihilo,    mulohilo,    muhilolo,    mulololo,
                ximhihihi[0],ximlohihi[0],ximhilohi[0],ximlolohi[0],
                ximhihilo[0],ximlohilo[0],ximhilolo[0],ximlololo[0],
               &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
               &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_sub(ximhihihi[0], ximlohihi[0], ximhilohi[0], ximlolohi[0],
                 ximhihilo[0], ximlohilo[0], ximhilolo[0], ximlololo[0],
                 acchihihi,    acclohihi,    acchilohi,    acclolohi,
                 acchihilo,    acclohilo,    acchilolo,    acclololo,
                &vimhihihi[0],&vimlohihi[0],&vimhilohi[0],&vimlolohi[0],
                &vimhihilo[0],&vimlohilo[0],&vimhilolo[0],&vimlololo[0]);
      }
      // cout << "mu : " << muhi << "  " << mulo << endl;

      // sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      odf_sqr(vrehihihi[0],vrelohihi[0],vrehilohi[0],vrelolohi[0],
              vrehihilo[0],vrelohilo[0],vrehilolo[0],vrelololo[0],
           &sqrv0hihihi,&sqrv0lohihi,&sqrv0hilohi,&sqrv0lolohi,
           &sqrv0hihilo,&sqrv0lohilo,&sqrv0hilolo,&sqrv0lololo);
      // cout << "vre[0] : " << vrehi[0] << "  " << vrelo[0] << endl;
      // cout << "vim[0] : " << vimhi[0] << "  " << vimlo[0] << endl;
      odf_sqr(vimhihihi[0],vimlohihi[0],vimhilohi[0],vimlolohi[0],
              vimhihilo[0],vimlohilo[0],vimhilolo[0],vimlololo[0],
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odf_inc(&sqrv0hihihi,&sqrv0lohihi,&sqrv0hilohi,&sqrv0lolohi,
              &sqrv0hihilo,&sqrv0lohilo,&sqrv0hilolo,&sqrv0lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      // cout << "sqrv0 : " << sqrv0hi << "  " << sqrv0lo << endl;
      // *beta = 2.0*sqrv0/(sigma + sqrv0);
      odf_inc(&sigmahihihi,&sigmalohihi,&sigmahilohi,&sigmalolohi,
              &sigmahihilo,&sigmalohilo,&sigmahilolo,&sigmalololo,
               sqrv0hihihi, sqrv0lohihi, sqrv0hilohi, sqrv0lolohi,
               sqrv0hihilo, sqrv0lohilo, sqrv0hilolo, sqrv0lololo);
      odf_div(sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi,
              sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo,
              sigmahihihi,sigmalohihi,sigmahilohi,sigmalolohi,
              sigmahihilo,sigmalohilo,sigmahilolo,sigmalololo,
               betahihihi, betalohihi, betahilohi, betalolohi,
               betahihilo, betalohilo, betahilolo, betalololo);
      odf_mlt_d(betahihihi,betalohihi,betahilohi,betalolohi,
                betahihilo,betalohilo,betahilolo,betalololo,2.0);
      // inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      odf_div(vrehihihi[0], vrelohihi[0], vrehilohi[0], vrelolohi[0],
              vrehihilo[0], vrelohilo[0], vrehilolo[0], vrelololo[0],
            sqrv0hihihi,  sqrv0lohihi,  sqrv0hilohi,  sqrv0lolohi,
            sqrv0hihilo,  sqrv0lohilo,  sqrv0hilolo,  sqrv0lololo,
          &inv0rehihihi,&inv0relohihi,&inv0rehilohi,&inv0relolohi,
          &inv0rehihilo,&inv0relohilo,&inv0rehilolo,&inv0relololo);
      // inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]
      odf_div(vimhihihi[0], vimlohihi[0], vimhilohi[0], vimlolohi[0],
              vimhihilo[0], vimlohilo[0], vimhilolo[0], vimlololo[0],
            sqrv0hihihi,  sqrv0lohihi,  sqrv0hilohi,  sqrv0lolohi,
            sqrv0hihilo,  sqrv0lohilo,  sqrv0hilolo,  sqrv0lololo,
          &inv0imhihihi,&inv0imlohihi,&inv0imhilohi,&inv0imlolohi,
          &inv0imhihilo,&inv0imlohilo,&inv0imhilolo,&inv0imlololo);
      odf_minus(&inv0imhihihi,&inv0imlohihi,&inv0imhilohi,&inv0imlolohi,
                &inv0imhihilo,&inv0imlohilo,&inv0imhilolo,&inv0imlololo);

      for(int i=1; i<n; i++)  // v[i] = v[i]/v[0]
      {
         // zre = vre[i]*inv0re - vim[i]*inv0im;
         odf_mul(vrehihihi[i],vrelohihi[i],vrehilohi[i],vrelolohi[i],
                 vrehihilo[i],vrelohilo[i],vrehilolo[i],vrelololo[i],
              inv0rehihihi,inv0relohihi,inv0rehilohi,inv0relolohi,
              inv0rehihilo,inv0relohilo,inv0rehilolo,inv0relololo,
                &zrehihihi,  &zrelohihi,  &zrehilohi,  &zrelolohi,
                &zrehihilo,  &zrelohilo,  &zrehilolo,  &zrelololo);
         odf_mul(vimhihihi[i],vimlohihi[i],vimhilohi[i],vimlolohi[i],
                 vimhihilo[i],vimlohilo[i],vimhilolo[i],vimlololo[i],
              inv0imhihihi,inv0imlohihi,inv0imhilohi,inv0imlolohi,
              inv0imhihilo,inv0imlohilo,inv0imhilolo,inv0imlololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_dec(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                 &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // zim = vim[i]*inv0re + vre[i]*inv0im;
         odf_mul(vimhihihi[i],vimlohihi[i],vimhilohi[i],vimlolohi[i],
                 vimhihilo[i],vimlohilo[i],vimhilolo[i],vimlololo[i],
              inv0rehihihi,inv0relohihi,inv0rehilohi,inv0relolohi,
              inv0rehihilo,inv0relohilo,inv0rehilolo,inv0relololo,
                &zimhihihi,  &zimlohihi,  &zimhilohi,  &zimlolohi,
                &zimhihilo,  &zimlohilo,  &zimhilolo,  &zimlololo);
         odf_mul(vrehihihi[i],vrelohihi[i],vrehilohi[i],vrelolohi[i],
                 vrehihilo[i],vrelohilo[i],vrehilolo[i],vrelololo[i],
              inv0imhihihi,inv0imlohihi,inv0imhilohi,inv0imlolohi,
              inv0imhihilo,inv0imlohilo,inv0imhilolo,inv0imlololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_inc(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                 &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         vrehihihi[i] = zrehihihi; vrelohihi[i] = zrelohihi;
         vrehilohi[i] = zrehilohi; vrelolohi[i] = zrelolohi;
         vrehihilo[i] = zrehihilo; vrelohilo[i] = zrelohilo;
         vrehilolo[i] = zrehilolo; vrelololo[i] = zrelololo;
         vimhihihi[i] = zimhihihi; vimlohihi[i] = zimlohihi;
         vimhilohi[i] = zimhilohi; vimlolohi[i] = zimlolohi;
         vimhihilo[i] = zimhihilo; vimlohilo[i] = zimlohilo;
         vimhilolo[i] = zimhilolo; vimlololo[i] = zimlololo;
      }
      vrehihihi[0] = 1.0; vrelohihi[0] = 0.0;
      vrehilohi[0] = 0.0; vrelolohi[0] = 0.0;
      vrehihilo[0] = 0.0; vrelohilo[0] = 0.0;
      vrehilolo[0] = 0.0; vrelololo[0] = 0.0;
      vimhihihi[0] = 0.0; vimlohihi[0] = 0.0;
      vimhilohi[0] = 0.0; vimlolohi[0] = 0.0;
      vimhihilo[0] = 0.0; vimlohilo[0] = 0.0;
      vimhilolo[0] = 0.0; vimlololo[0] = 0.0;
   }
}

void CPU_dbl8_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo )
{
   double *whihihi = new double[ncols-k];
   double *wlohihi = new double[ncols-k];
   double *whilohi = new double[ncols-k];
   double *wlolohi = new double[ncols-k];
   double *whihilo = new double[ncols-k];
   double *wlohilo = new double[ncols-k];
   double *whilolo = new double[ncols-k];
   double *wlololo = new double[ncols-k];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int j=k; j<ncols; j++)
   {
      whihihi[j-k] = 0.0; wlohihi[j-k] = 0.0;
      whilohi[j-k] = 0.0; wlolohi[j-k] = 0.0;
      whihilo[j-k] = 0.0; wlohilo[j-k] = 0.0;
      whilolo[j-k] = 0.0; wlololo[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         odf_mul(Rhihihi[i][j],Rlohihi[i][j],Rhilohi[i][j],Rlolohi[i][j],
                 Rhihilo[i][j],Rlohilo[i][j],Rhilolo[i][j],Rlololo[i][j],
                 vhihihi[i-k], vlohihi[i-k], vhilohi[i-k], vlolohi[i-k],
                 vhihilo[i-k], vlohilo[i-k], vhilolo[i-k], vlololo[i-k],
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odf_inc(&whihihi[j-k],&wlohihi[j-k],&whilohi[j-k],&wlolohi[j-k],
                 &whihilo[j-k],&wlohilo[j-k],&whilolo[j-k],&wlololo[j-k],
                acchihihi,    acclohihi,    acchilohi,    acclolohi,
                acchihilo,    acclohilo,    acchilolo,    acclololo);
      }
      // w[j-k] = beta*w[j-k];
      odf_mul(betahihihi,  betalohihi,  betahilohi,  betalolohi,
              betahihilo,  betalohilo,  betahilolo,  betalololo,
                 whihihi[j-k],wlohihi[j-k],whilohi[j-k],wlolohi[j-k],
                 whihilo[j-k],wlohilo[j-k],whilolo[j-k],wlololo[j-k],
              &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);

      whihihi[j-k] = acchihihi; wlohihi[j-k] = acclohihi;
      whilohi[j-k] = acchilohi; wlolohi[j-k] = acclolohi;
      whihilo[j-k] = acchihilo; wlohilo[j-k] = acclohilo;
      whilolo[j-k] = acchilolo; wlololo[j-k] = acclololo;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         odf_mul(vhihihi[i-k],vlohihi[i-k],vhilohi[i-k],vlolohi[i-k],
                 vhihilo[i-k],vlohilo[i-k],vhilolo[i-k],vlololo[i-k],
                 whihihi[j-k],wlohihi[j-k],whilohi[j-k],wlolohi[j-k],
                 whihilo[j-k],wlohilo[j-k],whilolo[j-k],wlololo[j-k],
              &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_dec(&Rhihihi[i][j],&Rlohihi[i][j],&Rhilohi[i][j],&Rlolohi[i][j],
                 &Rhihilo[i][j],&Rlohilo[i][j],&Rhilolo[i][j],&Rlololo[i][j],
                acchihihi,     acclohihi,     acchilohi,     acclolohi,
                acchihilo,     acclohilo,     acchilolo,     acclololo);
      }

   free(whihihi); free(wlohihi); free(whilohi); free(wlolohi);
   free(whihilo); free(wlohilo); free(whilolo); free(wlololo);
}

void CPU_cmplx8_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo )
{
   double *wrehihihi = new double[ncols-k];
   double *wrelohihi = new double[ncols-k];
   double *wrehilohi = new double[ncols-k];
   double *wrelolohi = new double[ncols-k];
   double *wrehihilo = new double[ncols-k];
   double *wrelohilo = new double[ncols-k];
   double *wrehilolo = new double[ncols-k];
   double *wrelololo = new double[ncols-k];
   double *wimhihihi = new double[ncols-k];
   double *wimlohihi = new double[ncols-k];
   double *wimhilohi = new double[ncols-k];
   double *wimlolohi = new double[ncols-k];
   double *wimhihilo = new double[ncols-k];
   double *wimlohilo = new double[ncols-k];
   double *wimhilolo = new double[ncols-k];
   double *wimlololo = new double[ncols-k];
   double zrehihihi,zrelohihi,zrehilohi,zrelolohi;
   double zrehihilo,zrelohilo,zrehilolo,zrelololo;
   double zimhihihi,zimlohihi,zimhilohi,zimlolohi;
   double zimhihilo,zimlohilo,zimhilolo,zimlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int j=k; j<ncols; j++)
   {
      wrehihihi[j-k] = 0.0; wrelohihi[j-k] = 0.0;
      wrehilohi[j-k] = 0.0; wrelolohi[j-k] = 0.0;
      wrehihilo[j-k] = 0.0; wrelohilo[j-k] = 0.0;
      wrehilolo[j-k] = 0.0; wrelololo[j-k] = 0.0;
      wimhihihi[j-k] = 0.0; wimlohihi[j-k] = 0.0;
      wimhilohi[j-k] = 0.0; wimlolohi[j-k] = 0.0;
      wimhihilo[j-k] = 0.0; wimlohilo[j-k] = 0.0;
      wimhilolo[j-k] = 0.0; wimlololo[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         // Hermitian of R => flip sign of Rim
         // zre =   Rre[i][j]*vre[i-k] + Rim[i][j]*vim[i-k];
         odf_mul(Rrehihihi[i][j],Rrelohihi[i][j],
                 Rrehilohi[i][j],Rrelolohi[i][j],
                 Rrehihilo[i][j],Rrelohilo[i][j],
                 Rrehilolo[i][j],Rrelololo[i][j],
                 vrehihihi[i-k], vrelohihi[i-k],
                 vrehilohi[i-k], vrelolohi[i-k],
                 vrehihilo[i-k], vrelohilo[i-k],
                 vrehilolo[i-k], vrelololo[i-k],
                &zrehihihi,     &zrelohihi,     &zrehilohi,     &zrelolohi,
                &zrehihilo,     &zrelohilo,     &zrehilolo,     &zrelololo);
         odf_mul(Rimhihihi[i][j],Rimlohihi[i][j],
                 Rimhilohi[i][j],Rimlolohi[i][j],
                 Rimhihilo[i][j],Rimlohilo[i][j],
                 Rimhilolo[i][j],Rimlololo[i][j],
                 vimhihihi[i-k], vimlohihi[i-k],
                 vimhilohi[i-k], vimlolohi[i-k],
                 vimhihilo[i-k], vimlohilo[i-k],
                 vimhilolo[i-k], vimlololo[i-k],
                &acchihihi,     &acclohihi,     &acchilohi,     &acclolohi,
                &acchihilo,     &acclohilo,     &acchilolo,     &acclololo);
         odf_inc(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                 &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // zim = - Rim[i][j]*vre[i-k] + Rre[i][j]*vim[i-k];
         odf_mul(Rrehihihi[i][j],Rrelohihi[i][j],
                 Rrehilohi[i][j],Rrelolohi[i][j],
                 Rrehihilo[i][j],Rrelohilo[i][j],
                 Rrehilolo[i][j],Rrelololo[i][j],
                 vimhihihi[i-k], vimlohihi[i-k],
                 vimhilohi[i-k], vimlolohi[i-k],
                 vimhihilo[i-k], vimlohilo[i-k],
                 vimhilolo[i-k], vimlololo[i-k],
                &zimhihihi,     &zimlohihi,     &zimhilohi,     &zimlolohi,
                &zimhihilo,     &zimlohilo,     &zimhilolo,     &zimlololo);
         odf_mul(Rimhihihi[i][j],Rimlohihi[i][j],
                 Rimhilohi[i][j],Rimlolohi[i][j],
                 Rimhihilo[i][j],Rimlohilo[i][j],
                 Rimhilolo[i][j],Rimlololo[i][j],
                 vrehihihi[i-k], vrelohihi[i-k],
                 vrehilohi[i-k], vrelolohi[i-k],
                 vrehihilo[i-k], vrelohilo[i-k],
                 vrehilolo[i-k], vrelololo[i-k],
                &acchihihi,     &acclohihi,     &acchilohi,     &acclolohi,
                &acchihilo,     &acclohilo,     &acchilolo,     &acclololo);
         odf_dec(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                 &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // wre[j-k] = wre[j-k] + zre;
         odf_inc(&wrehihihi[j-k],&wrelohihi[j-k],
                 &wrehilohi[j-k],&wrelolohi[j-k],
                 &wrehihilo[j-k],&wrelohilo[j-k],
                 &wrehilolo[j-k],&wrelololo[j-k],
                  zrehihihi,      zrelohihi,      zrehilohi,      zrelolohi,
                  zrehihilo,      zrelohilo,      zrehilolo,      zrelololo);
         // wim[j-k] = wim[j-k] + zim;
         odf_inc(&wimhihihi[j-k],&wimlohihi[j-k],
                 &wimhilohi[j-k],&wimlolohi[j-k],
                 &wimhihilo[j-k],&wimlohilo[j-k],
                 &wimhilolo[j-k],&wimlololo[j-k],
                  zimhihihi,      zimlohihi,      zimhilohi,      zimlolohi,
                  zimhihilo,      zimlohilo,      zimhilolo,      zimlololo);
      }
      // wre[j-k] = beta*wre[j-k];
      odf_mul(betahihihi,    betalohihi,    betahilohi,    betalolohi,
              betahihilo,    betalohilo,    betahilolo,    betalololo,
               wrehihihi[j-k],wrelohihi[j-k],wrehilohi[j-k],wrelolohi[j-k],
               wrehihilo[j-k],wrelohilo[j-k],wrehilolo[j-k],wrelololo[j-k],
              &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
              &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      wrehihihi[j-k] = acchihihi; wrelohihi[j-k] = acclohihi;
      wrehilohi[j-k] = acchilohi; wrelolohi[j-k] = acclolohi;
      wrehihilo[j-k] = acchihilo; wrelohilo[j-k] = acclohilo;
      wrehilolo[j-k] = acchilolo; wrelololo[j-k] = acclololo;
      // wim[j-k] = beta*wim[j-k];
      odf_mul(betahihihi,    betalohihi,    betahilohi,    betalolohi,
              betahihilo,    betalohilo,    betahilolo,    betalololo,
               wimhihihi[j-k],wimlohihi[j-k],wimhilohi[j-k],wimlolohi[j-k],
               wimhihilo[j-k],wimlohilo[j-k],wimhilolo[j-k],wimlololo[j-k],
              &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
              &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      wimhihihi[j-k] = acchihihi; wimlohihi[j-k] = acclohihi;
      wimhilohi[j-k] = acchilohi; wimlolohi[j-k] = acclolohi;
      wimhihilo[j-k] = acchihilo; wimlohilo[j-k] = acclohilo;
      wimhilolo[j-k] = acchilolo; wimlololo[j-k] = acclololo;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         // Hermitian of w => flip sign of wim
         // zre = vre[i-k]*wre[j-k] + vim[i-k]*wim[j-k];
         odf_mul(vrehihihi[i-k],vrelohihi[i-k],vrehilohi[i-k],vrelolohi[i-k],
                 vrehihilo[i-k],vrelohilo[i-k],vrehilolo[i-k],vrelololo[i-k],
                 wrehihihi[j-k],wrelohihi[j-k],wrehilohi[j-k],wrelolohi[j-k],
                 wrehihilo[j-k],wrelohilo[j-k],wrehilolo[j-k],wrelololo[j-k],
                &zrehihihi,    &zrelohihi,    &zrehilohi,    &zrelolohi,
                &zrehihilo,    &zrelohilo,    &zrehilolo,    &zrelololo);
         odf_mul(vimhihihi[i-k],vimlohihi[i-k],vimhilohi[i-k],vimlolohi[i-k],
                 vimhihilo[i-k],vimlohilo[i-k],vimhilolo[i-k],vimlololo[i-k],
                 wimhihihi[j-k],wimlohihi[j-k],wimhilohi[j-k],wimlolohi[j-k],
                 wimhihilo[j-k],wimlohilo[j-k],wimhilolo[j-k],wimlololo[j-k],
                &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
         odf_inc(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                 &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // zim = vim[i-k]*wre[j-k] - vre[i-k]*wim[j-k];
         odf_mul(vimhihihi[i-k],vimlohihi[i-k],vimhilohi[i-k],vimlolohi[i-k],
                 vimhihilo[i-k],vimlohilo[i-k],vimhilolo[i-k],vimlololo[i-k],
                 wrehihihi[j-k],wrelohihi[j-k],wrehilohi[j-k],wrelolohi[j-k],
                 wrehihilo[j-k],wrelohilo[j-k],wrehilolo[j-k],wrelololo[j-k],
                &zimhihihi,    &zimlohihi,    &zimhilohi,    &zimlolohi,
                &zimhihilo,    &zimlohilo,    &zimhilolo,    &zimlololo);
         odf_mul(vrehihihi[i-k],vrelohihi[i-k],vrehilohi[i-k],vrelolohi[i-k],
                 vrehihilo[i-k],vrelohilo[i-k],vrehilolo[i-k],vrelololo[i-k],
                 wimhihihi[j-k],wimlohihi[j-k],wimhilohi[j-k],wimlolohi[j-k],
                 wimhihilo[j-k],wimlohilo[j-k],wimhilolo[j-k],wimlololo[j-k],
                &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
         odf_dec(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                 &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // Rre[i][j] = Rre[i][j] - zre;
         odf_dec(&Rrehihihi[i][j],&Rrelohihi[i][j],
                 &Rrehilohi[i][j],&Rrelolohi[i][j],
                 &Rrehihilo[i][j],&Rrelohilo[i][j],
                 &Rrehilolo[i][j],&Rrelololo[i][j],
                  zrehihihi,       zrelohihi,     zrehilohi,     zrelolohi,
                  zrehihilo,       zrelohilo,     zrehilolo,     zrelololo);
         // Rim[i][j] = Rim[i][j] - zim;
         odf_dec(&Rimhihihi[i][j],&Rimlohihi[i][j],
                 &Rimhilohi[i][j],&Rimlolohi[i][j],
                 &Rimhihilo[i][j],&Rimlohilo[i][j],
                 &Rimhilolo[i][j],&Rimlololo[i][j],
                  zimhihihi,       zimlohihi,     zimhilohi,     zimlolohi,
                  zimhihilo,       zimlohilo,     zimhilolo,     zimlololo);
      }

   free(wrehihihi); free(wrelohihi); free(wrehilohi); free(wrelolohi);
   free(wrehihilo); free(wrelohilo); free(wrehilolo); free(wrelololo);
   free(wimhihihi); free(wimlohihi); free(wimhilohi); free(wimlolohi);
   free(wimhihilo); free(wimlohilo); free(wimhilolo); free(wimlololo);
}

void CPU_dbl8_factors_rightQupdate
 ( int n, int k,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo )
{
   double *whihihi = new double[n];
   double *wlohihi = new double[n];
   double *whilohi = new double[n];
   double *wlolohi = new double[n];
   double *whihilo = new double[n];
   double *wlohilo = new double[n];
   double *whilolo = new double[n];
   double *wlololo = new double[n];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<n; i++)
   {
      whihihi[i] = 0.0; wlohihi[i] = 0.0;
      whilohi[i] = 0.0; wlolohi[i] = 0.0;
      whihilo[i] = 0.0; wlohilo[i] = 0.0;
      whilolo[i] = 0.0; wlololo[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         odf_mul(Qhihihi[i][j],Qlohihi[i][j],Qhilohi[i][j],Qlolohi[i][j],
                 Qhihilo[i][j],Qlohilo[i][j],Qhilolo[i][j],Qlololo[i][j],
                 vhihihi[j-k], vlohihi[j-k], vhilohi[j-k], vlolohi[j-k],
                 vhihilo[j-k], vlohilo[j-k], vhilolo[j-k], vlololo[j-k],
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odf_inc(&whihihi[i],&wlohihi[i],&whilohi[i],&wlolohi[i],
                 &whihilo[i],&wlohilo[i],&whilolo[i],&wlololo[i],
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
      }
      // w[i] = beta*w[i];
      odf_mul(betahihihi,betalohihi,betahilohi,betalolohi,
              betahihilo,betalohilo,betahilolo,betalololo,
                 whihihi[i],wlohihi[i],whilohi[i],wlolohi[i],
                 whihilo[i],wlohilo[i],whilolo[i],wlololo[i],
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      whihihi[i] = acchihihi; wlohihi[i] = acclohihi;
      whilohi[i] = acchilohi; wlolohi[i] = acclolohi;
      whihilo[i] = acchihilo; wlohilo[i] = acclohilo;
      whilolo[i] = acchilolo; wlololo[i] = acclololo;
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         odf_mul(whihihi[i],  wlohihi[i],  whilohi[i],  wlolohi[i],
                 whihilo[i],  wlohilo[i],  whilolo[i],  wlololo[i],
                 vhihihi[j-k],vlohihi[j-k],vhilohi[j-k],vlolohi[j-k],
                 vhihilo[j-k],vlohilo[j-k],vhilolo[j-k],vlololo[j-k],
              &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_dec(&Qhihihi[i][j],&Qlohihi[i][j],&Qhilohi[i][j],&Qlolohi[i][j],
                 &Qhihilo[i][j],&Qlohilo[i][j],&Qhilolo[i][j],&Qlololo[i][j],
                acchihihi,     acclohihi,     acchilohi,     acclolohi,
                acchihilo,     acclohilo,     acchilolo,     acclololo);
      }

   free(whihihi); free(wlohihi); free(whilohi); free(wlolohi);
   free(whihilo); free(wlohilo); free(whilolo); free(wlololo);
}

void CPU_cmplx8_factors_rightQupdate
 ( int n, int k,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double betahihihi, double betalohihi,
   double betahilohi, double betalolohi,
   double betahihilo, double betalohilo,
   double betahilolo, double betalololo )
{
   double *wrehihihi = new double[n];
   double *wrelohihi = new double[n];
   double *wrehilohi = new double[n];
   double *wrelolohi = new double[n];
   double *wrehihilo = new double[n];
   double *wrelohilo = new double[n];
   double *wrehilolo = new double[n];
   double *wrelololo = new double[n];
   double *wimhihihi = new double[n];
   double *wimlohihi = new double[n];
   double *wimhilohi = new double[n];
   double *wimlolohi = new double[n];
   double *wimhihilo = new double[n];
   double *wimlohilo = new double[n];
   double *wimhilolo = new double[n];
   double *wimlololo = new double[n];
   double zrehihihi,zrelohihi,zrehilohi,zrelolohi;
   double zrehihilo,zrelohilo,zrehilolo,zrelololo;
   double zimhihihi,zimlohihi,zimhilohi,zimlolohi;
   double zimhihilo,zimlohilo,zimhilolo,zimlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<n; i++)
   {
      wrehihihi[i] = 0.0; wrelohihi[i] = 0.0;
      wrehilohi[i] = 0.0; wrelolohi[i] = 0.0;
      wrehihilo[i] = 0.0; wrelohilo[i] = 0.0;
      wrehilolo[i] = 0.0; wrelololo[i] = 0.0;
      wimhihihi[i] = 0.0; wimlohihi[i] = 0.0;
      wimhilohi[i] = 0.0; wimlolohi[i] = 0.0;
      wimhihilo[i] = 0.0; wimlohilo[i] = 0.0;
      wimhilolo[i] = 0.0; wimlololo[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         // zre = Qre[i][j]*vre[j-k] - Qim[i][j]*vim[j-k];
         odf_mul(Qrehihihi[i][j],Qrelohihi[i][j],
                 Qrehilohi[i][j],Qrelolohi[i][j],
                 Qrehihilo[i][j],Qrelohilo[i][j],
                 Qrehilolo[i][j],Qrelololo[i][j],
                 vrehihihi[j-k], vrelohihi[j-k],
                 vrehilohi[j-k], vrelolohi[j-k],
                 vrehihilo[j-k], vrelohilo[j-k],
                 vrehilolo[j-k], vrelololo[j-k],
                &zrehihihi,     &zrelohihi,     &zrehilohi,     &zrelolohi,
                &zrehihilo,     &zrelohilo,     &zrehilolo,     &zrelololo);
         odf_mul(Qimhihihi[i][j],Qimlohihi[i][j],
                 Qimhilohi[i][j],Qimlolohi[i][j],
                 Qimhihilo[i][j],Qimlohilo[i][j],
                 Qimhilolo[i][j],Qimlololo[i][j],
                 vimhihihi[j-k], vimlohihi[j-k],
                 vimhilohi[j-k], vimlolohi[j-k],
                 vimhihilo[j-k], vimlohilo[j-k],
                 vimhilolo[j-k], vimlololo[j-k],
                &acchihihi,     &acclohihi,     &acchilohi,     &acclolohi,
                &acchihilo,     &acclohilo,     &acchilolo,     &acclololo);
         odf_dec(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                 &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // zim = Qim[i][j]*vre[j-k] + Qre[i][j]*vim[j-k];
         odf_mul(Qimhihihi[i][j],Qimlohihi[i][j],
                 Qimhilohi[i][j],Qimlolohi[i][j],
                 Qimhihilo[i][j],Qimlohilo[i][j],
                 Qimhilolo[i][j],Qimlololo[i][j],
                 vrehihihi[j-k], vrelohihi[j-k],
                 vrehilohi[j-k], vrelolohi[j-k],
                 vrehihilo[j-k], vrelohilo[j-k],
                 vrehilolo[j-k], vrelololo[j-k],
                &zimhihihi,     &zimlohihi,     &zimhilohi,     &zimlolohi,
                &zimhihilo,     &zimlohilo,     &zimhilolo,     &zimlololo);
         odf_mul(Qrehihihi[i][j],Qrelohihi[i][j],
                 Qrehilohi[i][j],Qrelolohi[i][j],
                 Qrehihilo[i][j],Qrelohilo[i][j],
                 Qrehilolo[i][j],Qrelololo[i][j],
                 vimhihihi[j-k], vimlohihi[j-k],
                 vimhilohi[j-k], vimlolohi[j-k],
                 vimhihilo[j-k], vimlohilo[j-k],
                 vimhilolo[j-k], vimlololo[j-k],
                &acchihihi,     &acclohihi,     &acchilohi,     &acclolohi,
                &acchihilo,     &acclohilo,     &acchilolo,     &acclololo);
         odf_inc(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                 &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // wre[i] = wre[i] + zre;
         odf_inc(&wrehihihi[i],&wrelohihi[i],&wrehilohi[i],&wrelolohi[i],
                 &wrehihilo[i],&wrelohilo[i],&wrehilolo[i],&wrelololo[i],
                  zrehihihi,    zrelohihi,    zrehilohi,    zrelolohi,
                  zrehihilo,    zrelohilo,    zrehilolo,    zrelololo);
         // wim[i] = wim[i] + zim;
         odf_inc(&wimhihihi[i],&wimlohihi[i],&wimhilohi[i],&wimlolohi[i],
                 &wimhihilo[i],&wimlohilo[i],&wimhilolo[i],&wimlololo[i],
                  zimhihihi,    zimlohihi,    zimhilohi,    zimlolohi,
                  zimhihilo,    zimlohilo,    zimhilolo,    zimlololo);
      }
      // wre[i] = beta*wre[i];
      odf_mul(betahihihi,  betalohihi,  betahilohi,  betalolohi,
              betahihilo,  betalohilo,  betahilolo,  betalololo,
               wrehihihi[i],wrelohihi[i],wrehilohi[i],wrelolohi[i],
               wrehihilo[i],wrelohilo[i],wrehilolo[i],wrelololo[i],
              &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      wrehihihi[i] = acchihihi; wrelohihi[i] = acclohihi;
      wrehilohi[i] = acchilohi; wrelolohi[i] = acclolohi;
      wrehihilo[i] = acchihilo; wrelohilo[i] = acclohilo;
      wrehilolo[i] = acchilolo; wrelololo[i] = acclololo;
      // wim[i] = beta*wim[i];
      odf_mul(betahihihi,  betalohihi,  betahilohi,  betalolohi,
              betahihilo,  betalohilo,  betahilolo,  betalololo,
               wimhihihi[i],wimlohihi[i],wimhilohi[i],wimlolohi[i],
               wimhihilo[i],wimlohilo[i],wimhilolo[i],wimlololo[i],
              &acchihihi,  &acclohihi,  &acchilohi  ,&acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo  ,&acclololo);
      wimhihihi[i] = acchihihi; wimlohihi[i] = acclohihi;
      wimhilohi[i] = acchilohi; wimlolohi[i] = acclolohi;
      wimhihilo[i] = acchihilo; wimlohilo[i] = acclohilo;
      wimhilolo[i] = acchilolo; wimlololo[i] = acclololo;
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         // Hermitian transpose => flip sign of vim
         // zre = wre[i]*vre[j-k] + wim[i]*vim[j-k];
         odf_mul(wrehihihi[i],  wrelohihi[i],  wrehilohi[i],  wrelolohi[i],
                 wrehihilo[i],  wrelohilo[i],  wrehilolo[i],  wrelololo[i],
                 vrehihihi[j-k],vrelohihi[j-k],vrehilohi[j-k],vrelolohi[j-k],
                 vrehihilo[j-k],vrelohilo[j-k],vrehilolo[j-k],vrelololo[j-k],
                 &zrehihihi,   &zrelohihi,    &zrehilohi,    &zrelolohi,
                 &zrehihilo,   &zrelohilo,    &zrehilolo,    &zrelololo);
         odf_mul(wimhihihi[i],  wimlohihi[i],  wimhilohi[i],  wimlolohi[i],
                 wimhihilo[i],  wimlohilo[i],  wimhilolo[i],  wimlololo[i],
                 vimhihihi[j-k],vimlohihi[j-k],vimhilohi[j-k],vimlolohi[j-k],
                 vimhihilo[j-k],vimlohilo[j-k],vimhilolo[j-k],vimlololo[j-k],
                &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
         odf_inc(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
                 &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // zim = wim[i]*vre[j-k] - wre[i]*vim[j-k];
         odf_mul(wimhihihi[i],  wimlohihi[i],  wimhilohi[i],  wimlolohi[i],
                 wimhihilo[i],  wimlohilo[i],  wimhilolo[i],  wimlololo[i],
                 vrehihihi[j-k],vrelohihi[j-k],vrehilohi[j-k],vrelolohi[j-k],
                 vrehihilo[j-k],vrelohilo[j-k],vrehilolo[j-k],vrelololo[j-k],
                &zimhihihi,    &zimlohihi,    &zimhilohi,    &zimlolohi,
                &zimhihilo,    &zimlohilo,    &zimhilolo,    &zimlololo);
         odf_mul(wrehihihi[i],  wrelohihi[i],  wrehilohi[i],  wrelolohi[i],
                 wrehihilo[i],  wrelohilo[i],  wrehilolo[i],  wrelololo[i],
                 vimhihihi[j-k],vimlohihi[j-k],vimhilohi[j-k],vimlolohi[j-k],
                 vimhihilo[j-k],vimlohilo[j-k],vimhilolo[j-k],vimlololo[j-k],
                &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
         odf_dec(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
                 &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
         // Qre[i][j] = Qre[i][j] - zre;
         odf_dec(&Qrehihihi[i][j],&Qrelohihi[i][j],
                 &Qrehilohi[i][j],&Qrelolohi[i][j],
                 &Qrehihilo[i][j],&Qrelohilo[i][j],
                 &Qrehilolo[i][j],&Qrelololo[i][j],
                  zrehihihi,       zrelohihi,     zrehilohi,     zrelolohi,
                  zrehihilo,       zrelohilo,     zrehilolo,     zrelololo);
         // Qim[i][j] = Qim[i][j] - zim;
         odf_dec(&Qimhihihi[i][j],&Qimlohihi[i][j],
                 &Qimhilohi[i][j],&Qimlolohi[i][j],
                 &Qimhihilo[i][j],&Qimlohilo[i][j],
                 &Qimhilolo[i][j],&Qimlololo[i][j],
                  zimhihihi,       zimlohihi,     zimhilohi,     zimlolohi,
                  zimhihilo,       zimlohilo,     zimhilolo,     zimlololo);
      }

   free(wrehihihi); free(wrelohihi); free(wrehilohi); free(wrelolohi);
   free(wrehihilo); free(wrelohilo); free(wrehilolo); free(wrelololo);
   free(wimhihihi); free(wimlohihi); free(wimhilohi); free(wimlolohi);
   free(wimhihilo); free(wimlohilo); free(wimhilolo); free(wimlololo);
}

void CPU_dbl8_factors_houseqr
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo )
{
   double *xhihihi = new double[nrows];
   double *xlohihi = new double[nrows];
   double *xhilohi = new double[nrows];
   double *xlolohi = new double[nrows];
   double *xhihilo = new double[nrows];
   double *xlohilo = new double[nrows];
   double *xhilolo = new double[nrows];
   double *xlololo = new double[nrows];
   double *vhihihi = new double[nrows];
   double *vlohihi = new double[nrows];
   double *vhilohi = new double[nrows];
   double *vlolohi = new double[nrows];
   double *vhihilo = new double[nrows];
   double *vlohilo = new double[nrows];
   double *vhilolo = new double[nrows];
   double *vlololo = new double[nrows];
   double betahihihi,betalohihi,betahilohi,betalolohi;
   double betahihilo,betalohilo,betahilolo,betalololo;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qhihihi[i][j] = 0.0; Qlohihi[i][j] = 0.0;
         Qhilohi[i][j] = 0.0; Qlolohi[i][j] = 0.0;
         Qhihilo[i][j] = 0.0; Qlohilo[i][j] = 0.0;
         Qhilolo[i][j] = 0.0; Qlololo[i][j] = 0.0;
      }
      Qhihihi[i][i] = 1.0;

      for(int j=0; j<ncols; j++)
      {
         Rhihihi[i][j] = Ahihihi[i][j];
         Rlohihi[i][j] = Alohihi[i][j];
         Rhilohi[i][j] = Ahilohi[i][j];
         Rlolohi[i][j] = Alolohi[i][j];
         Rhihilo[i][j] = Ahihilo[i][j];
         Rlohilo[i][j] = Alohilo[i][j];
         Rhilolo[i][j] = Ahilolo[i][j];
         Rlololo[i][j] = Alololo[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xhihihi[i-k] = Rhihihi[i][k];
            xlohihi[i-k] = Rlohihi[i][k];
            xhilohi[i-k] = Rhilohi[i][k];
            xlolohi[i-k] = Rlolohi[i][k];
            xhihilo[i-k] = Rhihilo[i][k];
            xlohilo[i-k] = Rlohilo[i][k];
            xhilolo[i-k] = Rhilolo[i][k];
            xlololo[i-k] = Rlololo[i][k];
         }
         CPU_dbl8_factors_house
            (nrows-k,
             xhihihi,xlohihi,xhilohi,xlolohi,
             xhihilo,xlohilo,xhilolo,xlololo,
             vhihihi,vlohihi,vhilohi,vlolohi,
             vhihilo,vlohilo,vhilolo,vlololo,
             &betahihihi,&betalohihi,&betahilohi,&betalolohi,
             &betahihilo,&betalohilo,&betahilolo,&betalololo);
         CPU_dbl8_factors_leftRupdate
            (nrows,ncols,k,
             Rhihihi,Rlohihi,Rhilohi,Rlolohi,
             Rhihilo,Rlohilo,Rhilolo,Rlololo,
             vhihihi,vlohihi,vhilohi,vlolohi,
             vhihilo,vlohilo,vhilolo,vlololo,
             betahihihi,betalohihi,betahilohi,betalolohi,
             betahihilo,betalohilo,betahilolo,betalololo);
         CPU_dbl8_factors_rightQupdate
            (nrows,k,
             Qhihihi,Qlohihi,Qhilohi,Qlolohi,
             Qhihilo,Qlohilo,Qhilolo,Qlololo,
             vhihihi,vlohihi,vhilohi,vlolohi,
             vhihilo,vlohilo,vhilolo,vlololo,
             betahihihi,betalohihi,betahilohi,betalolohi,
             betahihilo,betalohilo,betahilolo,betalololo);
      }
   }
   free(xhihihi); free(xlohihi); free(xhilohi); free(xlolohi);
   free(xhihilo); free(xlohilo); free(xhilolo); free(xlololo);
   free(vhihihi); free(vlohihi); free(vhilohi); free(vlolohi);
   free(vhihilo); free(vlohilo); free(vhilolo); free(vlololo);
}

void CPU_cmplx8_factors_houseqr
 ( int nrows, int ncols,
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
   double **Rimhilolo, double **Rimlololo )
{
   double *xrehihihi = new double[nrows];
   double *xrelohihi = new double[nrows];
   double *xrehilohi = new double[nrows];
   double *xrelolohi = new double[nrows];
   double *xrehihilo = new double[nrows];
   double *xrelohilo = new double[nrows];
   double *xrehilolo = new double[nrows];
   double *xrelololo = new double[nrows];
   double *ximhihihi = new double[nrows];
   double *ximlohihi = new double[nrows];
   double *ximhilohi = new double[nrows];
   double *ximlolohi = new double[nrows];
   double *ximhihilo = new double[nrows];
   double *ximlohilo = new double[nrows];
   double *ximhilolo = new double[nrows];
   double *ximlololo = new double[nrows];
   double *vrehihihi = new double[nrows];
   double *vrelohihi = new double[nrows];
   double *vrehilohi = new double[nrows];
   double *vrelolohi = new double[nrows];
   double *vrehihilo = new double[nrows];
   double *vrelohilo = new double[nrows];
   double *vrehilolo = new double[nrows];
   double *vrelololo = new double[nrows];
   double *vimhihihi = new double[nrows];
   double *vimlohihi = new double[nrows];
   double *vimhilohi = new double[nrows];
   double *vimlolohi = new double[nrows];
   double *vimhihilo = new double[nrows];
   double *vimlohilo = new double[nrows];
   double *vimhilolo = new double[nrows];
   double *vimlololo = new double[nrows];
   double betahihihi,betalohihi,betahilohi,betalolohi;
   double betahihilo,betalohilo,betahilolo,betalololo;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qrehihihi[i][j] = 0.0; Qrelohihi[i][j] = 0.0;
         Qrehilohi[i][j] = 0.0; Qrelolohi[i][j] = 0.0;
         Qrehihilo[i][j] = 0.0; Qrelohilo[i][j] = 0.0;
         Qrehilolo[i][j] = 0.0; Qrelololo[i][j] = 0.0;
         Qimhihihi[i][j] = 0.0; Qimlohihi[i][j] = 0.0;
         Qimhilohi[i][j] = 0.0; Qimlolohi[i][j] = 0.0;
         Qimhihilo[i][j] = 0.0; Qimlohilo[i][j] = 0.0;
         Qimhilolo[i][j] = 0.0; Qimlololo[i][j] = 0.0;
      }
      Qrehihihi[i][i] = 1.0;

      for(int j=0; j<ncols; j++)
      {
         Rrehihihi[i][j] = Arehihihi[i][j];
         Rrelohihi[i][j] = Arelohihi[i][j];
         Rrehilohi[i][j] = Arehilohi[i][j];
         Rrelolohi[i][j] = Arelolohi[i][j];
         Rrehihilo[i][j] = Arehihilo[i][j];
         Rrelohilo[i][j] = Arelohilo[i][j];
         Rrehilolo[i][j] = Arehilolo[i][j];
         Rrelololo[i][j] = Arelololo[i][j];
         Rimhihihi[i][j] = Aimhihihi[i][j];
         Rimlohihi[i][j] = Aimlohihi[i][j];
         Rimhilohi[i][j] = Aimhilohi[i][j];
         Rimlolohi[i][j] = Aimlolohi[i][j];
         Rimhihilo[i][j] = Aimhihilo[i][j];
         Rimlohilo[i][j] = Aimlohilo[i][j];
         Rimhilolo[i][j] = Aimhilolo[i][j];
         Rimlololo[i][j] = Aimlololo[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xrehihihi[i-k] = Rrehihihi[i][k];
            xrelohihi[i-k] = Rrelohihi[i][k];
            xrehilohi[i-k] = Rrehilohi[i][k];
            xrelolohi[i-k] = Rrelolohi[i][k];
            xrehihilo[i-k] = Rrehihilo[i][k];
            xrelohilo[i-k] = Rrelohilo[i][k];
            xrehilolo[i-k] = Rrehilolo[i][k];
            xrelololo[i-k] = Rrelololo[i][k];
            ximhihihi[i-k] = Rimhihihi[i][k];
            ximlohihi[i-k] = Rimlohihi[i][k];
            ximhilohi[i-k] = Rimhilohi[i][k];
            ximlolohi[i-k] = Rimlolohi[i][k];
            ximhihilo[i-k] = Rimhihilo[i][k];
            ximlohilo[i-k] = Rimlohilo[i][k];
            ximhilolo[i-k] = Rimhilolo[i][k];
            ximlololo[i-k] = Rimlololo[i][k];
         }
         CPU_cmplx8_factors_house
            (nrows-k,
             xrehihihi,xrelohihi,xrehilohi,xrelolohi,
             xrehihilo,xrelohilo,xrehilolo,xrelololo,
             ximhihihi,ximlohihi,ximhilohi,ximlolohi,
             ximhihilo,ximlohilo,ximhilolo,ximlololo,
             vrehihihi,vrelohihi,vrehilohi,vrelolohi,
             vrehihilo,vrelohilo,vrehilolo,vrelololo,
             vimhihihi,vimlohihi,vimhilohi,vimlolohi,
             vimhihilo,vimlohilo,vimhilolo,vimlololo,
             &betahihihi,&betalohihi,&betahilohi,&betalolohi,
             &betahihilo,&betalohilo,&betahilolo,&betalololo);
         CPU_cmplx8_factors_leftRupdate
            (nrows,ncols,k,
             Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
             Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
             Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
             Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
             vrehihihi,vrelohihi,vrehilohi,vrelolohi,
             vrehihilo,vrelohilo,vrehilolo,vrelololo,
             vimhihihi,vimlohihi,vimhilohi,vimlolohi,
             vimhihilo,vimlohilo,vimhilolo,vimlololo,
             betahihihi,betalohihi,betahilohi,betalolohi,
             betahihilo,betalohilo,betahilolo,betalololo);
         CPU_cmplx8_factors_rightQupdate
            (nrows,k,
             Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
             Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
             Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
             Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
             vrehihihi,vrelohihi,vrehilohi,vrelolohi,
             vrehihilo,vrelohilo,vrehilolo,vrelololo,
             vimhihihi,vimlohihi,vimhilohi,vimlolohi,
             vimhihilo,vimlohilo,vimhilolo,vimlololo,
             betahihihi,betalohihi,betahilohi,betalolohi,
             betahihilo,betalohilo,betahilolo,betalololo);
      }
   }
   free(xrehihihi); free(xrelohihi); free(xrehilohi); free(xrelolohi);
   free(xrehihilo); free(xrelohilo); free(xrehilolo); free(xrelololo);
   free(ximhihihi); free(ximlohihi); free(ximhilohi); free(ximlolohi);
   free(ximhihilo); free(ximlohilo); free(ximhilolo); free(ximlololo);
   free(vrehihihi); free(vrelohihi); free(vrehilohi); free(vrelolohi);
   free(vrehihilo); free(vrelohilo); free(vrehilolo); free(vrelololo);
   free(vimhihihi); free(vimlohihi); free(vimhilohi); free(vimlolohi);
   free(vimhihilo); free(vimlohilo); free(vimhilolo); free(vimlololo);
}

void CPU_dbl8_factors_qrbs
 ( int nrows, int ncols,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *rhshihihi, double *rhslohihi, double *rhshilohi, double *rhslolohi,
   double *rhshihilo, double *rhslohilo, double *rhshilolo, double *rhslololo,
   double *solhihihi, double *sollohihi, double *solhilohi, double *sollolohi,
   double *solhihilo, double *sollohilo, double *solhilolo, double *sollololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<nrows; i++)   // compute Q^T*b, b is rhs
   {
      wrkvechihihi[i] = 0.0;
      wrkveclohihi[i] = 0.0;
      wrkvechilohi[i] = 0.0;
      wrkveclolohi[i] = 0.0;
      wrkvechihilo[i] = 0.0;
      wrkveclohilo[i] = 0.0;
      wrkvechilolo[i] = 0.0;
      wrkveclololo[i] = 0.0;

      for(int j=0; j<nrows; j++) // wrkvec[i] = wrkvec[i] + Q[j][i]*rhs[j];
      {
         odf_mul(Qhihihi[j][i],Qlohihi[j][i],Qhilohi[j][i],Qlolohi[j][i],
                 Qhihilo[j][i],Qlohilo[j][i],Qhilolo[j][i],Qlololo[j][i],
                 rhshihihi[j],rhslohihi[j],rhshilohi[j],rhslolohi[j],
                 rhshihilo[j],rhslohilo[j],rhshilolo[j],rhslololo[j],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&wrkvechihihi[i],&wrkveclohihi[i],
                 &wrkvechilohi[i],&wrkveclolohi[i],
                 &wrkvechihilo[i],&wrkveclohilo[i],
                 &wrkvechilolo[i],&wrkveclololo[i],
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo);
      }
   }
   CPU_dbl8_factors_backward
      (ncols,Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
       wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
       wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo);
}

void CPU_cmplx8_factors_qrbs
 ( int nrows, int ncols,
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
   double *rhsrehihihi, double *rhsrelohihi,
   double *rhsrehilohi, double *rhsrelolohi,
   double *rhsrehihilo, double *rhsrelohilo,
   double *rhsrehilolo, double *rhsrelololo,
   double *rhsimhihihi, double *rhsimlohihi,
   double *rhsimhilohi, double *rhsimlolohi,
   double *rhsimhihilo, double *rhsimlohilo,
   double *rhsimhilolo, double *rhsimlololo,
   double *solrehihihi, double *solrelohihi,
   double *solrehilohi, double *solrelolohi,
   double *solrehihilo, double *solrelohilo,
   double *solrehilolo, double *solrelololo,
   double *solimhihihi, double *solimlohihi,
   double *solimhilohi, double *solimlolohi,
   double *solimhihilo, double *solimlohilo,
   double *solimhilolo, double *solimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi; // accumulates product
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<nrows; i++)    // compute Q^H*b, b is rhs
   {
      wrkvecrehihihi[i] = 0.0; wrkvecrelohihi[i] = 0.0;
      wrkvecrehilohi[i] = 0.0; wrkvecrelolohi[i] = 0.0;
      wrkvecrehihilo[i] = 0.0; wrkvecrelohilo[i] = 0.0;
      wrkvecrehilolo[i] = 0.0; wrkvecrelololo[i] = 0.0;
      wrkvecimhihihi[i] = 0.0; wrkvecimlohihi[i] = 0.0;
      wrkvecimhilohi[i] = 0.0; wrkvecimlolohi[i] = 0.0;
      wrkvecimhihilo[i] = 0.0; wrkvecimlohilo[i] = 0.0;
      wrkvecimhilolo[i] = 0.0; wrkvecimlololo[i] = 0.0;

      for(int j=0; j<nrows; j++) // work with Hermitian transpose of Q
      {
         // accre =  Qre[j][i]*rhsre[j] + Qim[j][i]*rhsim[j];
         // wrkvecre[i] = wrkvecre[i] + accre;
         odf_mul(Qrehihihi[j][i],Qrelohihi[j][i],
                 Qrehilohi[j][i],Qrelolohi[j][i],
                 Qrehihilo[j][i],Qrelohilo[j][i],
                 Qrehilolo[j][i],Qrelololo[j][i],
                 rhsrehihihi[j],rhsrelohihi[j],rhsrehilohi[j],rhsrelolohi[j],
                 rhsrehihilo[j],rhsrelohilo[j],rhsrehilolo[j],rhsrelololo[j],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&wrkvecrehihihi[i],&wrkvecrelohihi[i],
                 &wrkvecrehilohi[i],&wrkvecrelolohi[i],
                 &wrkvecrehihilo[i],&wrkvecrelohilo[i],
                 &wrkvecrehilolo[i],&wrkvecrelololo[i],
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo);
         odf_mul(Qimhihihi[j][i],Qimlohihi[j][i],
                 Qimhilohi[j][i],Qimlolohi[j][i],
                 Qimhihilo[j][i],Qimlohilo[j][i],
                 Qimhilolo[j][i],Qimlololo[j][i],
                 rhsimhihihi[j],rhsimlohihi[j],rhsimhilohi[j],rhsimlolohi[j],
                 rhsimhihilo[j],rhsimlohilo[j],rhsimhilolo[j],rhsimlololo[j],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&wrkvecrehihihi[i],&wrkvecrelohihi[i],
                 &wrkvecrehilohi[i],&wrkvecrelolohi[i],
                 &wrkvecrehihilo[i],&wrkvecrelohilo[i],
                 &wrkvecrehilolo[i],&wrkvecrelololo[i],
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo);
         // accim = -Qim[j][i]*rhsre[j] + Qre[j][i]*rhsim[j];
         // wrkvecim[i] = wrkvecim[i] + accim;
         odf_mul(Qrehihihi[j][i],Qrelohihi[j][i],
                 Qrehilohi[j][i],Qrelolohi[j][i],
                 Qrehihilo[j][i],Qrelohilo[j][i],
                 Qrehilolo[j][i],Qrelololo[j][i],
                 rhsimhihihi[j],rhsimlohihi[j],rhsimhilohi[j],rhsimlolohi[j],
                 rhsimhihilo[j],rhsimlohilo[j],rhsimhilolo[j],rhsimlololo[j],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&wrkvecimhihihi[i],&wrkvecimlohihi[i],
                 &wrkvecimhilohi[i],&wrkvecimlolohi[i],
                 &wrkvecimhihilo[i],&wrkvecimlohilo[i],
                 &wrkvecimhilolo[i],&wrkvecimlololo[i],
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo);
         odf_mul(Qimhihihi[j][i],Qimlohihi[j][i],
                 Qimhilohi[j][i],Qimlolohi[j][i],
                 Qimhihilo[j][i],Qimlohilo[j][i],
                 Qimhilolo[j][i],Qimlololo[j][i],
                 rhsrehihihi[j],rhsrelohihi[j],rhsrehilohi[j],rhsrelolohi[j],
                 rhsrehihilo[j],rhsrelohilo[j],rhsrehilolo[j],rhsrelololo[j],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_dec(&wrkvecimhihihi[i],&wrkvecimlohihi[i],
                 &wrkvecimhilohi[i],&wrkvecimlolohi[i],
                 &wrkvecimhihilo[i],&wrkvecimlohilo[i],
                 &wrkvecimhilolo[i],&wrkvecimlololo[i],
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo);
      }
   }
   CPU_cmplx8_factors_backward
      (ncols,Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
             Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
             Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
             Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
       wrkvecrehihihi,wrkvecrelohihi,wrkvecrehilohi,wrkvecrelolohi,
       wrkvecrehihilo,wrkvecrelohilo,wrkvecrehilolo,wrkvecrelololo,
       wrkvecimhihihi,wrkvecimlohihi,wrkvecimhilohi,wrkvecimlolohi,
       wrkvecimhihilo,wrkvecimlohilo,wrkvecimhilolo,wrkvecimlololo,
       solrehihihi,solrelohihi,solrehilohi,solrelolohi,
       solrehihilo,solrelohilo,solrehilolo,solrelololo,
       solimhihihi,solimlohihi,solimhilohi,solimlolohi,
       solimhihilo,solimlohilo,solimhilolo,solimlololo);
}
