/* The file dbl4_factorizations.cpp defines functions specified in
 * the file dbl4_factorizations.h. */

#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"

// #include <iostream>
// using namespace std;

void CPU_dbl4_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo,
   double **Chihi, double **Clohi, double **Chilo, double **Clolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Chihi[i][j] = 0.0; Clohi[i][j] = 0.0;
         Chilo[i][j] = 0.0; Clolo[i][j] = 0.0;

         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            qdf_mul(Ahihi[i][k],Alohi[i][k],Ahilo[i][k],Alolo[i][k],
                    Bhihi[k][j],Blohi[k][j],Bhilo[k][j],Blolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&Chihi[i][j],&Clohi[i][j],&Chilo[i][j],&Clolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }
}

void CPU_cmplx4_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Brehihi, double **Brelohi, double **Brehilo, double **Brelolo,
   double **Bimhihi, double **Bimlohi, double **Bimhilo, double **Bimlolo,
   double **Crehihi, double **Crelohi, double **Crehilo, double **Crelolo,
   double **Cimhihi, double **Cimlohi, double **Cimhilo, double **Cimlolo )
{
   double zrehihi,zrelohi,zrehilo,zrelolo;
   double zimhihi,zimlohi,zimhilo,zimlolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Crehihi[i][j] = 0.0; Crelohi[i][j] = 0.0;
         Crehilo[i][j] = 0.0; Crelolo[i][j] = 0.0;
         Cimhihi[i][j] = 0.0; Cimlohi[i][j] = 0.0;
         Cimhilo[i][j] = 0.0; Cimlolo[i][j] = 0.0;

         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            // zre = Are[i][k]*Bre[k][j] - Aim[i][k]*Bim[k][j];
            qdf_mul(Arehihi[i][k],Arelohi[i][k],Arehilo[i][k],Arelolo[i][k],
                    Brehihi[k][j],Brelohi[k][j],Brehilo[k][j],Brelolo[k][j],
                    &zrehihi,&zrelohi,&zrehilo,&zrelolo);
            qdf_mul(Aimhihi[i][k],Aimlohi[i][k],Aimhilo[i][k],Aimlolo[i][k],
                    Bimhihi[k][j],Bimlohi[k][j],Bimhilo[k][j],Bimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_dec(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
                    acchihi,acclohi,acchilo,acclolo);
            // zim = Aim[i][k]*Bre[k][j] + Are[i][k]*Bim[k][j];
            qdf_mul(Aimhihi[i][k],Aimlohi[i][k],Aimhilo[i][k],Aimlolo[i][k],
                    Brehihi[k][j],Brelohi[k][j],Brehilo[k][j],Brelolo[k][j],
                    &zimhihi,&zimlohi,&zimhilo,&zimlolo);
            qdf_mul(Arehihi[i][k],Arelohi[i][k],Arehilo[i][k],Arelolo[i][k],
                    Bimhihi[k][j],Bimlohi[k][j],Bimhilo[k][j],Bimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
                    acchihi,acclohi,acchilo,acclolo);
            // Cre[i][j] = Cre[i][j] + zre;
            qdf_inc(&Crehihi[i][j],&Crelohi[i][j],
                    &Crehilo[i][j],&Crelolo[i][j],
                    zrehihi,zrelohi,zrehilo,zrelolo);
            // Cim[i][j] = Cim[i][j] + zim;
            qdf_inc(&Cimhihi[i][j],&Cimlohi[i][j],
                    &Cimhilo[i][j],&Cimlolo[i][j],
                    zimhihi,zimlohi,zimhilo,zimlolo);
         }
      }
}

void CPU_dbl4_factors_forward
 ( int dim,
   double **Lhihi, double **Llohi, double **Lhilo, double **Llolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   xhihi[0] = bhihi[0]; xlohi[0] = blohi[0];
   xhilo[0] = bhilo[0]; xlolo[0] = blolo[0];

   for(int i=1; i<dim; i++)
   {
      xhihi[i] = bhihi[i]; xlohi[i] = blohi[i];
      xhilo[i] = bhilo[i]; xlolo[i] = blolo[i];

      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         qdf_mul(Lhihi[i][j],Llohi[i][j],Lhilo[i][j],Llolo[i][j],
                 xhihi[j],   xlohi[j],   xhilo[j],   xlolo[j],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_dec(&xhihi[i],&xlohi[i],&xhilo[i],&xlolo[i],
                acchihi,  acclohi,  acchilo,  acclolo);
      }
   }
}

void CPU_cmplx4_factors_forward
 ( int dim,
   double **Lrehihi, double **Lrelohi, double **Lrehilo, double **Lrelolo,
   double **Limhihi, double **Limlohi, double **Limhilo, double **Limlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo )
{
   double acc1rehihi,acc1relohi,acc1rehilo,acc1relolo;
   double acc1imhihi,acc1imlohi,acc1imhilo,acc1imlolo;
   double acc2rehihi,acc2relohi,acc2rehilo,acc2relolo;
   double acc2imhihi,acc2imlohi,acc2imhilo,acc2imlolo;

   xrehihi[0] = brehihi[0]; xrelohi[0] = brelohi[0];
   xrehilo[0] = brehilo[0]; xrelolo[0] = brelolo[0];
   ximhihi[0] = bimhihi[0]; ximlohi[0] = bimlohi[0];
   ximhilo[0] = bimhilo[0]; ximlolo[0] = bimlolo[0];

   for(int i=1; i<dim; i++)
   {
      xrehihi[i] = brehihi[i]; xrelohi[i] = brelohi[i];
      xrehilo[i] = brehilo[i]; xrelolo[i] = brelolo[i];
      ximhihi[i] = bimhihi[i]; ximlohi[i] = bimlohi[i];
      ximhilo[i] = bimhilo[i]; ximlolo[i] = bimlolo[i];

      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         qdf_mul(Lrehihi[i][j],Lrelohi[i][j],Lrehilo[i][j],Lrelolo[i][j],
                 xrehihi[j],   xrelohi[j],   xrehilo[j],   xrelolo[j],
             &acc1rehihi,  &acc1relohi,  &acc1rehilo,  &acc1relolo);
         qdf_mul(Limhihi[i][j],Limlohi[i][j],Limhilo[i][j],Limlolo[i][j],
                 ximhihi[j],   ximlohi[j],   ximhilo[j],   ximlolo[j],
             &acc1imhihi,  &acc1imlohi,  &acc1imhilo,  &acc1imlolo);
         qdf_mul(Limhihi[i][j],Limlohi[i][j],Limhilo[i][j],Limlolo[i][j],
                 xrehihi[j],   xrelohi[j],   xrehilo[j],   xrelolo[j],
             &acc2rehihi,  &acc2relohi,  &acc2rehilo,  &acc2relolo);
         qdf_mul(Lrehihi[i][j],Lrelohi[i][j],Lrehilo[i][j],Lrelolo[i][j],
                 ximhihi[j],   ximlohi[j],   ximhilo[j],   ximlolo[j],
             &acc2imhihi,  &acc2imlohi,  &acc2imhilo,  &acc2imlolo);

         qdf_dec(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
               acc1rehihi, acc1relohi, acc1rehilo, acc1relolo);
         qdf_inc(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
               acc1imhihi, acc1imlohi, acc1imhilo, acc1imlolo);
         qdf_dec(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
               acc2rehihi, acc2relohi, acc2rehilo, acc2relolo);
         qdf_dec(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
               acc2imhihi, acc2imlohi, acc2imhilo, acc2imlolo);
      }
   }
}

void CPU_dbl4_factors_backward
 ( int dim,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=dim-1; i>=0; i--)
   {
      xhihi[i] = bhihi[i]; xlohi[i] = blohi[i];
      xhilo[i] = bhilo[i]; xlolo[i] = blolo[i];

      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         qdf_mul(Uhihi[i][j],Ulohi[i][j],Uhilo[i][j],Ulolo[i][j],
                 xhihi[j],   xlohi[j],   xhilo[j],   xlolo[j],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_dec(&xhihi[i],&xlohi[i],&xhilo[i],&xlolo[i],
                acchihi,  acclohi,  acchilo,  acclolo);
      }
      // x[i] = x[i]/U[i][i];
      qdf_div(xhihi[i],   xlohi[i],   xhilo[i],   xlolo[i],
              Uhihi[i][i],Ulohi[i][i],Uhilo[i][i],Ulolo[i][i],
           &acchihi,   &acclohi,   &acchilo,   &acclolo);
      xhihi[i] = acchihi; xlohi[i] = acclohi;
      xhilo[i] = acchilo; xlolo[i] = acclolo;
   }
}

void CPU_cmplx4_factors_backward
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo )
{
   double acc1rehihi,acc1relohi,acc1rehilo,acc1relolo;
   double acc1imhihi,acc1imlohi,acc1imhilo,acc1imlolo;
   double acc2rehihi,acc2relohi,acc2rehilo,acc2relolo;
   double acc2imhihi,acc2imlohi,acc2imhilo,acc2imlolo;
   double acc3rehihi,acc3relohi,acc3rehilo,acc3relolo;
   double acc3imhihi,acc3imlohi,acc3imhilo,acc3imlolo;
   double denhihi,denlohi,denhilo,denlolo;

   for(int i=dim-1; i>=0; i--)
   {
      xrehihi[i] = brehihi[i]; xrelohi[i] = brelohi[i];
      xrehilo[i] = brehilo[i]; xrelolo[i] = brelolo[i];
      ximhihi[i] = bimhihi[i]; ximlohi[i] = bimlohi[i];
      ximhilo[i] = bimhilo[i]; ximlolo[i] = bimlolo[i];

      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         qdf_mul(Urehihi[i][j],Urelohi[i][j],Urehilo[i][j],Urelolo[i][j],
                 xrehihi[j],   xrelohi[j],   xrehilo[j],   xrelolo[j],
             &acc1rehihi,  &acc1relohi,  &acc1rehilo,  &acc1relolo);
         qdf_mul(Uimhihi[i][j],Uimlohi[i][j],Uimhilo[i][j],Uimlolo[i][j],
                 ximhihi[j],   ximlohi[j],   ximhilo[j],   ximlolo[j],
             &acc1imhihi,  &acc1imlohi,  &acc1imhilo,  &acc1imlolo);
         qdf_mul(Uimhihi[i][j],Uimlohi[i][j],Uimhilo[i][j],Uimlolo[i][j],
                 xrehihi[j],   xrelohi[j],   xrehilo[j],   xrelolo[j],
             &acc2rehihi,  &acc2relohi,  &acc2rehilo,  &acc2relolo);
         qdf_mul(Urehihi[i][j],Urelohi[i][j],Urehilo[i][j],Urelolo[i][j],
                 ximhihi[j],   ximlohi[j],   ximhilo[j],   ximlolo[j],
             &acc2imhihi,  &acc2imlohi,  &acc2imhilo,  &acc2imlolo);
         qdf_dec(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
               acc1rehihi, acc1relohi, acc1rehilo, acc1relolo);
         qdf_inc(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
               acc1imhihi, acc1imlohi, acc1imhilo, acc1imlolo);
         qdf_dec(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
               acc2rehihi, acc2relohi, acc2rehilo, acc2relolo);
         qdf_dec(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
               acc2imhihi, acc2imlohi, acc2imhilo, acc2imlolo);
      }
      // x[i] = x[i]/U[i][i];
      qdf_mul(Urehihi[i][i],Urelohi[i][i],Urehilo[i][i],Urelolo[i][i],
              Urehihi[i][i],Urelohi[i][i],Urehilo[i][i],Urelolo[i][i],
             &denhihi,     &denlohi,     &denhilo,     &denlolo);
      qdf_mul(Uimhihi[i][i],Uimlohi[i][i],Uimhilo[i][i],Uimlolo[i][i],
              Uimhihi[i][i],Uimlohi[i][i],Uimhilo[i][i],Uimlolo[i][i],
          &acc1rehihi,  &acc1relohi,  &acc1rehilo,  &acc1relolo);
      qdf_inc(&denhihi,  &denlohi,  &denhilo,  &denlolo,
            acc1rehihi,acc1relohi,acc1rehilo,acc1relolo); // denominator
      qdf_div(Urehihi[i][i],Urelohi[i][i],Urehilo[i][i],Urelolo[i][i],
              denhihi,      denlohi,      denhilo,      denlolo,
          &acc1rehihi,  &acc1relohi,  &acc1rehilo,  &acc1relolo);
      // (acc1rehi,acc1relo) is real part of 1/U[i][i]
      qdf_div(Uimhihi[i][i],Uimlohi[i][i],Uimhilo[i][i],Uimlolo[i][i],
              denhihi,      denlohi,      denhilo,      denlolo,
          &acc1imhihi,  &acc1imlohi,  &acc1imhilo,  &acc1imlolo);
      qdf_minus(&acc1imhihi,&acc1imlohi,&acc1imhilo,&acc1imlolo);
      // (acc1imhi,acc1imlo) is imaginary part of 1/U[i][i]
      qdf_mul(xrehihi[i], xrelohi[i], xrehilo[i], xrelolo[i],
           acc1rehihi, acc1relohi, acc1rehilo, acc1relolo,
          &acc2rehihi,&acc2relohi,&acc2rehilo,&acc2relolo);
      qdf_mul(ximhihi[i], ximlohi[i], ximhilo[i], ximlolo[i],
           acc1imhihi, acc1imlohi, acc1imhilo, acc1imlolo,
          &acc2imhihi,&acc2imlohi,&acc2imhilo,&acc2imlolo);
      // acc2 stores the doubles for xre
      qdf_mul(ximhihi[i], ximlohi[i], ximhilo[i], ximlolo[i],
           acc1rehihi, acc1relohi, acc1rehilo, acc1relolo,
          &acc3rehihi,&acc3relohi,&acc3rehilo,&acc3relolo);
      qdf_mul(xrehihi[i],xrelohi[i],
              xrehilo[i],xrelolo[i],
              acc1imhihi,acc1imlohi,
              acc1imhilo,acc1imlolo,
              &acc3imhihi,&acc3imlohi,
              &acc3imhilo,&acc3imlolo);
      // acc3 stores the doubles for xim
      xrehihi[i] = acc2rehihi; xrelohi[i] = acc2relohi;
      xrehilo[i] = acc2rehilo; xrelolo[i] = acc2relolo;
      qdf_dec(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
            acc2imhihi, acc2imlohi, acc2imhilo, acc2imlolo);
      ximhihi[i] = acc3rehihi; ximlohi[i] = acc3relohi;
      ximhilo[i] = acc3rehilo; ximlolo[i] = acc3relolo;
      qdf_inc(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
            acc3imhihi, acc3imlohi, acc3imhilo, acc3imlolo);
   }
}

void CPU_dbl4_factors_lufac
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   int *pivots )
{
   double valmax,valtmp,acchihi,acclohi,acchilo,acclolo;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Ahihi[j][j]); idxmax = j;     // find the pivot
      for(int i=j+1; i<dim; i++)
      {
         valtmp = fabs(Ahihi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Ahihi[idxmax][k]; Ahihi[idxmax][k] = Ahihi[j][k];
            Ahihi[j][k] = valtmp;
            valtmp = Alohi[idxmax][k]; Alohi[idxmax][k] = Alohi[j][k];
            Alohi[j][k] = valtmp;
            valtmp = Ahilo[idxmax][k]; Ahilo[idxmax][k] = Ahilo[j][k];
            Ahilo[j][k] = valtmp;
            valtmp = Alolo[idxmax][k]; Alolo[idxmax][k] = Alolo[j][k];
            Alolo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         qdf_div(Ahihi[i][j],Alohi[i][j],Ahilo[i][j],Alolo[i][j],
                 Ahihi[j][j],Alohi[j][j],Ahilo[j][j],Alolo[j][j],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);

         Ahihi[i][j] = acchihi; Alohi[i][j] = acclohi;
         Ahilo[i][j] = acchilo; Alolo[i][j] = acclolo;

         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            qdf_mul(Ahihi[i][j],Alohi[i][j],Ahilo[i][j],Alolo[i][j],
                    Ahihi[j][k],Alohi[j][k],Ahilo[j][k],Alolo[j][k],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
            qdf_dec(&Ahihi[i][k],&Alohi[i][k],&Ahilo[i][k],&Alolo[i][k],
                   acchihi,     acclohi,     acchilo,     acclolo);
         }
      }
   }
}

void CPU_cmplx4_factors_lufac
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   int *pivots )
{
   double valmax,valtmp;
   int idxmax,idxtmp;
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double denhihi,denlohi,denhilo,denlolo;
   double acc3rehihi,acc3relohi,acc3rehilo,acc3relolo;
   double acc3imhihi,acc3imlohi,acc3imhilo,acc3imlolo;
   double acc4rehihi,acc4relohi,acc4rehilo,acc4relolo;
   double acc4imhihi,acc4imlohi,acc4imhilo,acc4imlolo;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Arehihi[j][j]) + fabs(Aimhihi[j][j]);
      idxmax = j;     // find the pivot
      for(int i=j+1; i<dim; i++)
      {
         valtmp = fabs(Arehihi[i][j]) + fabs(Aimhihi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Arehihi[idxmax][k]; Arehihi[idxmax][k] = Arehihi[j][k];
            Arehihi[j][k] = valtmp;
            valtmp = Arelohi[idxmax][k]; Arelohi[idxmax][k] = Arelohi[j][k];
            Arelohi[j][k] = valtmp;
            valtmp = Arehilo[idxmax][k]; Arehilo[idxmax][k] = Arehilo[j][k];
            Arehilo[j][k] = valtmp;
            valtmp = Arelolo[idxmax][k]; Arelolo[idxmax][k] = Arelolo[j][k];
            Arelolo[j][k] = valtmp;
            valtmp = Aimhihi[idxmax][k]; Aimhihi[idxmax][k] = Aimhihi[j][k];
            Aimhihi[j][k] = valtmp;
            valtmp = Aimlohi[idxmax][k]; Aimlohi[idxmax][k] = Aimlohi[j][k];
            Aimlohi[j][k] = valtmp;
            valtmp = Aimhilo[idxmax][k]; Aimhilo[idxmax][k] = Aimhilo[j][k];
            Aimhilo[j][k] = valtmp;
            valtmp = Aimlolo[idxmax][k]; Aimlolo[idxmax][k] = Aimlolo[j][k];
            Aimlolo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         qdf_mul(Arehihi[j][j],Arelohi[j][j],Arehilo[j][j],Arelolo[j][j],
                 Arehihi[j][j],Arelohi[j][j],Arehilo[j][j],Arelolo[j][j],
                &denhihi,     &denlohi,     &denhilo,     &denlolo);
         qdf_mul(Aimhihi[j][j],Aimlohi[j][j],Aimhilo[j][j],Aimlolo[j][j],
                 Aimhihi[j][j],Aimlohi[j][j],Aimhilo[j][j],Aimlolo[j][j],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdf_inc(&denhihi,&denlohi,&denhilo,&denlolo,
                 acc1hihi,acc1lohi,acc1hilo,acc1lolo); // denominator
         qdf_div(Arehihi[j][j],Arelohi[j][j],Arehilo[j][j],Arelolo[j][j],
                 denhihi,      denlohi,      denhilo,      denlolo,
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         // (acc1hi,acc1lo) is real part of 1/A[j][j]
         qdf_div(Aimhihi[j][j],Aimlohi[j][j],Aimhilo[j][j],Aimlolo[j][j],
                 denhihi,      denlohi,      denhilo,      denlolo,
               &acc2hihi,    &acc2lohi,    &acc2hilo,   &acc2lolo);
         qdf_minus(&acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
         // (acc2hi,acc2lo) is imaginary part of 1/A[j][j]
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
                acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo,
             &acc3rehihi,  &acc3relohi,  &acc3rehilo,  &acc3relolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
                acc2hihi,     acc2lohi,     acc2hilo,     acc2lolo,
             &acc3imhihi,  &acc3imlohi,  &acc3imhilo,  &acc3imlolo);
         // acc3 stores doubles for Arehi
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
                acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo,
             &acc4rehihi,  &acc4relohi,  &acc4rehilo,  &acc4relolo);
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
                acc2hihi     ,acc2lohi,     acc2hilo,     acc2lolo,
             &acc4imhihi,  &acc4imlohi,  &acc4imhilo,  &acc4imlolo);
         // acc4 stores doubles for Aimhi
         Arehihi[i][j] = acc3rehihi; Arelohi[i][j] = acc3relohi;
         Arehilo[i][j] = acc3rehilo; Arelolo[i][j] = acc3relolo;
         qdf_dec(&Arehihi[i][j],&Arelohi[i][j],&Arehilo[i][j],&Arelolo[i][j],
               acc3imhihi,    acc3imlohi,    acc3imhilo,    acc3imlolo);
         Aimhihi[i][j] = acc4rehihi; Aimlohi[i][j] = acc4relohi;
         Aimhilo[i][j] = acc4rehilo; Aimlolo[i][j] = acc4relolo;
         qdf_inc(&Aimhihi[i][j],&Aimlohi[i][j],&Aimhilo[i][j],&Aimlolo[i][j],
               acc4imhihi,    acc4imlohi,    acc4imhilo,    acc4imlolo);

         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
                    Arehihi[j][k],Arelohi[j][k],Arehilo[j][k],Arelolo[j][k],
                &acc3rehihi,  &acc3relohi,  &acc3rehilo,  &acc3relolo);
            qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
                    Aimhihi[j][k],Aimlohi[j][k],Aimhilo[j][k],Aimlolo[j][k],
                &acc3imhihi,  &acc3imlohi,  &acc3imhilo,  &acc3imlolo);
            qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
                    Arehihi[j][k],Arelohi[j][k],Arehilo[j][k],Arelolo[j][k],
                &acc4rehihi,  &acc4relohi,  &acc4rehilo,  &acc4relolo);
            qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
                    Aimhihi[j][k],Aimlohi[j][k],Aimhilo[j][k],Aimlolo[j][k],
                &acc4imhihi,  &acc4imlohi,  &acc4imhilo,  &acc4imlolo);

            qdf_dec(&Arehihi[i][k],&Arelohi[i][k],
                    &Arehilo[i][k],&Arelolo[i][k],
                    acc3rehihi,acc3relohi,acc3rehilo,acc3relolo);
            qdf_inc(&Arehihi[i][k],&Arelohi[i][k],
                    &Arehilo[i][k],&Arelolo[i][k],
                    acc3imhihi,acc3imlohi,acc3imhilo,acc3imlolo);
            qdf_dec(&Aimhihi[i][k],&Aimlohi[i][k],
                    &Aimhilo[i][k],&Aimlolo[i][k],
                    acc4rehihi,acc4relohi,acc4rehilo,acc4relolo);
            qdf_dec(&Aimhihi[i][k],&Aimlohi[i][k],
                    &Aimhilo[i][k],&Aimlolo[i][k],
                    acc4imhihi,acc4imlohi,acc4imhilo,acc4imlolo);
         }
      }
   }
}

void CPU_dbl4_factors_lusolve
 ( int dim,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   int *pivots,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   CPU_dbl4_factors_lufac(dim,Ahihi,Alohi,Ahilo,Alolo,pivots);
   for(int i=0; i<dim; i++) 
   {
      xhihi[i] = bhihi[pivots[i]];
      xlohi[i] = blohi[pivots[i]];
      xhilo[i] = bhilo[pivots[i]];
      xlolo[i] = blolo[pivots[i]];
   }
   CPU_dbl4_factors_forward
      (dim,Ahihi,Alohi,Ahilo,Alolo,
           xhihi,xlohi,xhilo,xlolo,
           bhihi,blohi,bhilo,blolo);
   CPU_dbl4_factors_backward
      (dim,Ahihi,Alohi,Ahilo,Alolo,
           bhihi,blohi,bhilo,blolo,
           xhihi,xlohi,xhilo,xlolo);
}

void CPU_cmplx4_factors_lusolve
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   int *pivots,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo )
{
   CPU_cmplx4_factors_lufac
      (dim,Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
       pivots);

   for(int i=0; i<dim; i++) 
   {
      xrehihi[i] = brehihi[pivots[i]];
      xrelohi[i] = brelohi[pivots[i]];
      xrehilo[i] = brehilo[pivots[i]];
      xrelolo[i] = brelolo[pivots[i]];
      ximhihi[i] = bimhihi[pivots[i]];
      ximlohi[i] = bimlohi[pivots[i]];
      ximhilo[i] = bimhilo[pivots[i]];
      ximlolo[i] = bimlolo[pivots[i]];
   }
   CPU_cmplx4_factors_forward
      (dim,Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
           xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
           brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo);

   CPU_cmplx4_factors_backward
      (dim,Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
           brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
           xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo);
}

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

void CPU_dbl4_factors_qrbs
 ( int nrows, int ncols,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *rhshihi, double *rhslohi, double *rhshilo, double *rhslolo,
   double *solhihi, double *sollohi, double *solhilo, double *sollolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<nrows; i++)   // compute Q^T*b, b is rhs
   {
      wrkvechihi[i] = 0.0;
      wrkveclohi[i] = 0.0;
      wrkvechilo[i] = 0.0;
      wrkveclolo[i] = 0.0;

      for(int j=0; j<nrows; j++) // wrkvec[i] = wrkvec[i] + Q[j][i]*rhs[j];
      {
         qdf_mul(Qhihi[j][i],Qlohi[j][i],Qhilo[j][i],Qlolo[j][i],
                 rhshihi[j],rhslohi[j],rhshilo[j],rhslolo[j],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&wrkvechihi[i],&wrkveclohi[i],&wrkvechilo[i],&wrkveclolo[i],
                 acchihi,acclohi,acchilo,acclolo);
      }
   }
   CPU_dbl4_factors_backward
      (ncols,Rhihi,Rlohi,Rhilo,Rlolo,
       wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo,
       solhihi,sollohi,solhilo,sollolo);
}

void CPU_cmplx4_factors_qrbs
 ( int nrows, int ncols,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo, 
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *rhsrehihi, double *rhsrelohi, double *rhsrehilo, double *rhsrelolo,
   double *rhsimhihi, double *rhsimlohi, double *rhsimhilo, double *rhsimlolo,
   double *solrehihi, double *solrelohi, double *solrehilo, double *solrelolo,
   double *solimhihi, double *solimlohi, double *solimhilo, double *solimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo )
{
   double acchihi,acclohi,acchilo,acclolo; // accumulates product

   for(int i=0; i<nrows; i++)    // compute Q^H*b, b is rhs
   {
      wrkvecrehihi[i] = 0.0; wrkvecrelohi[i] = 0.0;
      wrkvecrehilo[i] = 0.0; wrkvecrelolo[i] = 0.0;
      wrkvecimhihi[i] = 0.0; wrkvecimlohi[i] = 0.0;
      wrkvecimhilo[i] = 0.0; wrkvecimlolo[i] = 0.0;

      for(int j=0; j<nrows; j++) // work with Hermitian transpose of Q
      {
         // accre =  Qre[j][i]*rhsre[j] + Qim[j][i]*rhsim[j];
         // wrkvecre[i] = wrkvecre[i] + accre;
         qdf_mul(Qrehihi[j][i],Qrelohi[j][i],Qrehilo[j][i],Qrelolo[j][i],
                 rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                 &acchihi,&acclohi, &acchilo,&acclolo);
         qdf_inc(&wrkvecrehihi[i],&wrkvecrelohi[i],
                 &wrkvecrehilo[i],&wrkvecrelolo[i],
                 acchihi,acclohi,acchilo,acclolo);
         qdf_mul(Qimhihi[j][i],Qimlohi[j][i],Qimhilo[j][i],Qimlolo[j][i],
                 rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&wrkvecrehihi[i],&wrkvecrelohi[i],
                 &wrkvecrehilo[i],&wrkvecrelolo[i],
                 acchihi,acclohi,acchilo,acclolo);
         // accim = -Qim[j][i]*rhsre[j] + Qre[j][i]*rhsim[j];
         // wrkvecim[i] = wrkvecim[i] + accim;
         qdf_mul(Qrehihi[j][i],Qrelohi[j][i],Qrehilo[j][i],Qrelolo[j][i],
                 rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&wrkvecimhihi[i],&wrkvecimlohi[i],
                 &wrkvecimhilo[i],&wrkvecimlolo[i],
                 acchihi,acclohi,acchilo,acclolo);
         qdf_mul(Qimhihi[j][i],Qimlohi[j][i],Qimhilo[j][i],Qimlolo[j][i],
                 rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_dec(&wrkvecimhihi[i],&wrkvecimlohi[i],
                 &wrkvecimhilo[i],&wrkvecimlolo[i],
                 acchihi,acclohi,acchilo,acclolo);
      }
   }
   CPU_cmplx4_factors_backward
      (ncols,Rrehihi,Rrelohi,Rrehilo,Rrelolo,
             Rimhihi,Rimlohi,Rimhilo,Rimlolo,
       wrkvecrehihi,wrkvecrelohi,wrkvecrehilo,wrkvecrelolo,
       wrkvecimhihi,wrkvecimlohi,wrkvecimhilo,wrkvecimlolo,
       solrehihi,solrelohi,solrehilo,solrelolo,
       solimhihi,solimlohi,solimhilo,solimlolo);
}
