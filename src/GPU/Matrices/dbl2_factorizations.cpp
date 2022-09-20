/* The file dbl2_factorizations.cpp defines functions specified in
 * the file dbl2_factorizations.h. */

#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "dbl2_factorizations.h"

// #include <iostream>
// using namespace std;

void CPU_dbl2_factors_matmatmul
 ( int rows, int dim, int cols, double **Ahi, double **Alo,
   double **Bhi, double **Blo, double **Chi, double **Clo )
{
   double acchi,acclo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Chi[i][j] = 0.0;
         Clo[i][j] = 0.0;
         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            ddf_mul(Ahi[i][k],Alo[i][k],Bhi[k][j],Blo[k][j],&acchi,&acclo);
            ddf_inc(&Chi[i][j],&Clo[i][j],acchi,acclo);
         }
      }
}

void CPU_cmplx2_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Brehi, double **Brelo, double **Bimhi, double **Bimlo,
   double **Crehi, double **Crelo, double **Cimhi, double **Cimlo )
{
   double zrehi,zrelo,zimhi,zimlo,acchi,acclo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Crehi[i][j] = 0.0; Cimhi[i][j] = 0.0;
         Crelo[i][j] = 0.0; Cimlo[i][j] = 0.0;

         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            // zre = Are[i][k]*Bre[k][j] - Aim[i][k]*Bim[k][j];
            ddf_mul(Arehi[i][k],Arelo[i][k],Brehi[k][j],Brelo[k][j],
                    &zrehi,&zrelo);
            ddf_mul(Aimhi[i][k],Aimlo[i][k],Bimhi[k][j],Bimlo[k][j],
                    &acchi,&acclo);
            ddf_dec(&zrehi,&zrelo,acchi,acclo);
            // zim = Aim[i][k]*Bre[k][j] + Are[i][k]*Bim[k][j];
            ddf_mul(Aimhi[i][k],Aimlo[i][k],Brehi[k][j],Brelo[k][j],
                    &zimhi,&zimlo);
            ddf_mul(Arehi[i][k],Arelo[i][k],Bimhi[k][j],Bimlo[k][j],
                    &acchi,&acclo);
            ddf_inc(&zimhi,&zimlo,acchi,acclo);
            // Cre[i][j] = Cre[i][j] + zre;
            ddf_inc(&Crehi[i][j],&Crelo[i][j],zrehi,zrelo);
            // Cim[i][j] = Cim[i][j] + zim;
            ddf_inc(&Cimhi[i][j],&Cimlo[i][j],zimhi,zimlo);
         }
      }
}

void CPU_dbl2_factors_forward
 ( int dim, double **Lhi, double **Llo, double *bhi, double *blo,
   double *xhi, double *xlo )
{
   double acchi,acclo;

   xhi[0] = bhi[0];
   xlo[0] = blo[0];
   for(int i=1; i<dim; i++)
   {
      xhi[i] = bhi[i];
      xlo[i] = blo[i];
      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         ddf_mul(Lhi[i][j],Llo[i][j],xhi[j],xlo[j],&acchi,&acclo);
         ddf_dec(&xhi[i],&xlo[i],acchi,acclo);
      }
   }
}

void CPU_cmplx2_factors_forward
 ( int dim, double **Lrehi, double **Lrelo, double **Limhi, double **Limlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo )
{
   double acc1rehi,acc1relo,acc1imhi,acc1imlo;
   double acc2rehi,acc2relo,acc2imhi,acc2imlo;

   xrehi[0] = brehi[0]; xrelo[0] = brelo[0];
   ximhi[0] = bimhi[0]; ximlo[0] = bimlo[0];

   for(int i=1; i<dim; i++)
   {
      xrehi[i] = brehi[i]; xrelo[i] = brelo[i];
      ximhi[i] = bimhi[i]; ximlo[i] = bimlo[i];

      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         ddf_mul(Lrehi[i][j],Lrelo[i][j],xrehi[j],xrelo[j],
                 &acc1rehi,&acc1relo);
         ddf_mul(Limhi[i][j],Limlo[i][j],ximhi[j],ximlo[j],
                 &acc1imhi,&acc1imlo);
         ddf_mul(Limhi[i][j],Limlo[i][j],xrehi[j],xrelo[j],
                 &acc2rehi,&acc2relo);
         ddf_mul(Lrehi[i][j],Lrelo[i][j],ximhi[j],ximlo[j],
                 &acc2imhi,&acc2imlo);
         ddf_dec(&xrehi[i],&xrelo[i],acc1rehi,acc1relo);
         ddf_inc(&xrehi[i],&xrelo[i],acc1imhi,acc1imlo);
         ddf_dec(&ximhi[i],&ximlo[i],acc2rehi,acc2relo);
         ddf_dec(&ximhi[i],&ximlo[i],acc2imhi,acc2imlo);
      }
   }
}

void CPU_dbl2_factors_backward
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo )
{
   double acchi,acclo;

   for(int i=dim-1; i>=0; i--)
   {
      xhi[i] = bhi[i];
      xlo[i] = blo[i];
      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         ddf_mul(Uhi[i][j],Ulo[i][j],xhi[j],xlo[j],&acchi,&acclo);
         ddf_dec(&xhi[i],&xlo[i],acchi,acclo);
      }
      // x[i] = x[i]/U[i][i];
      ddf_div(xhi[i],xlo[i],Uhi[i][i],Ulo[i][i],&acchi,&acclo);
      xhi[i] = acchi; xlo[i] = acclo;
   }
}

void CPU_cmplx2_factors_backward
 ( int dim, double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo )
{
   double acc1rehi,acc1relo,acc1imhi,acc1imlo;
   double acc2rehi,acc2relo,acc2imhi,acc2imlo;
   double acc3rehi,acc3relo,acc3imhi,acc3imlo;
   double denhi,denlo;

   for(int i=dim-1; i>=0; i--)
   {
      xrehi[i] = brehi[i]; xrelo[i] = brelo[i];
      ximhi[i] = bimhi[i]; ximlo[i] = bimlo[i];

      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         ddf_mul(Urehi[i][j],Urelo[i][j],xrehi[j],xrelo[j],
                 &acc1rehi,&acc1relo);
         ddf_mul(Uimhi[i][j],Uimlo[i][j],ximhi[j],ximlo[j],
                 &acc1imhi,&acc1imlo);
         ddf_mul(Uimhi[i][j],Uimlo[i][j],xrehi[j],xrelo[j],
                 &acc2rehi,&acc2relo);
         ddf_mul(Urehi[i][j],Urelo[i][j],ximhi[j],ximlo[j],
                 &acc2imhi,&acc2imlo);
         ddf_dec(&xrehi[i],&xrelo[i],acc1rehi,acc1relo);
         ddf_inc(&xrehi[i],&xrelo[i],acc1imhi,acc1imlo);
         ddf_dec(&ximhi[i],&ximlo[i],acc2rehi,acc2relo);
         ddf_dec(&ximhi[i],&ximlo[i],acc2imhi,acc2imlo);
      }
      // x[i] = x[i]/U[i][i];
      ddf_mul(Urehi[i][i],Urelo[i][i],
              Urehi[i][i],Urelo[i][i],&denhi,&denlo);
      ddf_mul(Uimhi[i][i],Uimlo[i][i],
              Uimhi[i][i],Uimlo[i][i],&acc1rehi,&acc1relo);
      ddf_inc(&denhi,&denlo,acc1rehi,acc1relo); // denominator
      ddf_div(Urehi[i][i],Urelo[i][i],denhi,denlo,&acc1rehi,&acc1relo);
      // (acc1rehi,acc1relo) is real part of 1/U[i][i]
      ddf_div(Uimhi[i][i],Uimlo[i][i],denhi,denlo,&acc1imhi,&acc1imlo);
      ddf_minus(&acc1imhi,&acc1imlo);
      // (acc1imhi,acc1imlo) is imaginary part of 1/U[i][i]
      ddf_mul(xrehi[i],xrelo[i],acc1rehi,acc1relo,&acc2rehi,&acc2relo);
      ddf_mul(ximhi[i],ximlo[i],acc1imhi,acc1imlo,&acc2imhi,&acc2imlo);
      // acc2 stores the doubles for xre
      ddf_mul(ximhi[i],ximlo[i],acc1rehi,acc1relo,&acc3rehi,&acc3relo);
      ddf_mul(xrehi[i],xrelo[i],acc1imhi,acc1imlo,&acc3imhi,&acc3imlo);
      // acc3 stores the doubles for xim
      xrehi[i] = acc2rehi; xrelo[i] = acc2relo;
      ddf_dec(&xrehi[i],&xrelo[i],acc2imhi,acc2imlo);
      ximhi[i] = acc3rehi; ximlo[i] = acc3relo;
      ddf_inc(&ximhi[i],&ximlo[i],acc3imhi,acc3imlo);
   }
}

void CPU_dbl2_factors_lufac
 ( int dim, double **Ahi, double **Alo, int *pivots )
{
   double valmax,valtmp,acchi,acclo;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Ahi[j][j]); idxmax = j;     // find the pivot
      for(int i=j+1; i<dim; i++)
      {
         valtmp = fabs(Ahi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Ahi[idxmax][k];
            Ahi[idxmax][k] = Ahi[j][k];
            Ahi[j][k] = valtmp;
            valtmp = Alo[idxmax][k];
            Alo[idxmax][k] = Alo[j][k];
            Alo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         ddf_div(Ahi[i][j],Alo[i][j],Ahi[j][j],Alo[j][j],&acchi,&acclo);
         Ahi[i][j] = acchi;
         Alo[i][j] = acclo;
         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            ddf_mul(Ahi[i][j],Alo[i][j],Ahi[j][k],Alo[j][k],&acchi,&acclo);
            ddf_dec(&Ahi[i][k],&Alo[i][k],acchi,acclo);
         }
      }
   }
}

void CPU_cmplx2_factors_lufac
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   int *pivots )
{
   double valmax,valtmp;
   int idxmax,idxtmp;
   double acc1hi,acc1lo,acc2hi,acc2lo,denhi,denlo;
   double acc3rehi,acc3relo,acc3imhi,acc3imlo;
   double acc4rehi,acc4relo,acc4imhi,acc4imlo;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Arehi[j][j]) + fabs(Aimhi[j][j]);
      idxmax = j;     // find the pivot
      for(int i=j+1; i<dim; i++)
      {
         valtmp = fabs(Arehi[i][j]) + fabs(Aimhi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Arehi[idxmax][k];
            Arehi[idxmax][k] = Arehi[j][k];
            Arehi[j][k] = valtmp;
            valtmp = Arelo[idxmax][k];
            Arelo[idxmax][k] = Arelo[j][k];
            Arelo[j][k] = valtmp;
            valtmp = Aimhi[idxmax][k];
            Aimhi[idxmax][k] = Aimhi[j][k];
            Aimhi[j][k] = valtmp;
            valtmp = Aimlo[idxmax][k];
            Aimlo[idxmax][k] = Aimlo[j][k];
            Aimlo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         ddf_mul(Arehi[j][j],Arelo[j][j],
                 Arehi[j][j],Arelo[j][j],&denhi,&denlo);
         ddf_mul(Aimhi[j][j],Aimlo[j][j],
                 Aimhi[j][j],Aimlo[j][j],&acc1hi,&acc1lo);
         ddf_inc(&denhi,&denlo,acc1hi,acc1lo); // denominator
         ddf_div(Arehi[j][j],Arelo[j][j],denhi,denlo,&acc1hi,&acc1lo);
         // (acc1hi,acc1lo) is real part of 1/A[j][j]
         ddf_div(Aimhi[j][j],Aimlo[j][j],denhi,denlo,&acc2hi,&acc2lo);
         ddf_minus(&acc2hi,&acc2lo);
         // (acc2hi,acc2lo) is imaginary part of 1/A[j][j]
         ddf_mul(Arehi[i][j],Arelo[i][j],acc1hi,acc1lo,
                 &acc3rehi,&acc3relo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],acc2hi,acc2lo,
                 &acc3imhi,&acc3imlo);
         // acc3 stores doubles for Arehi
         ddf_mul(Aimhi[i][j],Aimlo[i][j],acc1hi,acc1lo,
                 &acc4rehi,&acc4relo);
         ddf_mul(Arehi[i][j],Arelo[i][j],acc2hi,acc2lo,
                 &acc4imhi,&acc4imlo);
         // acc4 stores doubles for Aimhi
         Arehi[i][j] = acc3rehi; Arelo[i][j] = acc3relo;
         ddf_dec(&Arehi[i][j],&Arelo[i][j],acc3imhi,acc3imlo);
         Aimhi[i][j] = acc4rehi; Aimlo[i][j] = acc4relo;
         ddf_inc(&Aimhi[i][j],&Aimlo[i][j],acc4imhi,acc4imlo);

         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            ddf_mul(Arehi[i][j],Arelo[i][j],Arehi[j][k],Arelo[j][k],
                    &acc3rehi,&acc3relo);
            ddf_mul(Aimhi[i][j],Aimlo[i][j],Aimhi[j][k],Aimlo[j][k],
                    &acc3imhi,&acc3imlo);
            ddf_mul(Aimhi[i][j],Aimlo[i][j],Arehi[j][k],Arelo[j][k],
                    &acc4rehi,&acc4relo);
            ddf_mul(Arehi[i][j],Arelo[i][j],Aimhi[j][k],Aimlo[j][k],
                    &acc4imhi,&acc4imlo);
            ddf_dec(&Arehi[i][k],&Arelo[i][k],acc3rehi,acc3relo);
            ddf_inc(&Arehi[i][k],&Arelo[i][k],acc3imhi,acc3imlo);
            ddf_dec(&Aimhi[i][k],&Aimlo[i][k],acc4rehi,acc4relo);
            ddf_dec(&Aimhi[i][k],&Aimlo[i][k],acc4imhi,acc4imlo);
         }
      }
   }
}

void CPU_dbl2_factors_lusolve
 ( int dim, double **Ahi, double **Alo, int *pivots,
   double *bhi, double *blo, double *xhi, double *xlo )
{
   CPU_dbl2_factors_lufac(dim,Ahi,Alo,pivots);
   for(int i=0; i<dim; i++) 
   {
      xhi[i] = bhi[pivots[i]];
      xlo[i] = blo[pivots[i]];
   }
   CPU_dbl2_factors_forward(dim,Ahi,Alo,xhi,xlo,bhi,blo);
   CPU_dbl2_factors_backward(dim,Ahi,Alo,bhi,blo,xhi,xlo);
}

void CPU_cmplx2_factors_lusolve
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   int *pivots, double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo )
{
   CPU_cmplx2_factors_lufac(dim,Arehi,Arelo,Aimhi,Aimlo,pivots);

   for(int i=0; i<dim; i++) 
   {
      xrehi[i] = brehi[pivots[i]];
      xrelo[i] = brelo[pivots[i]];
      ximhi[i] = bimhi[pivots[i]];
      ximlo[i] = bimlo[pivots[i]];
   }
   CPU_cmplx2_factors_forward
      (dim,Arehi,Arelo,Aimhi,Aimlo,
           xrehi,xrelo,ximhi,ximlo,
           brehi,brelo,bimhi,bimlo);

   CPU_cmplx2_factors_backward
      (dim,Arehi,Arelo,Aimhi,Aimlo,
           brehi,brelo,bimhi,bimlo,
           xrehi,xrelo,ximhi,ximlo);
}

void CPU_dbl2_factors_house
 ( int n, double *xhi, double *xlo, double *vhi, double *vlo,
   double *betahi, double *betalo )
{
   double sigmahi = 0.0;
   double sigmalo = 0.0;
   double muhi,v0p2hi,acchi;
   double mulo,v0p2lo,acclo;
   
   vhi[0] = 1.0;
   vlo[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      // sigma = sigma + x[i]*x[i];
      ddf_sqr(xhi[i],xlo[i],&acchi,&acclo);
      ddf_inc(&sigmahi,&sigmalo,acchi,acclo);
      vhi[i] = xhi[i];
      vlo[i] = xlo[i];
   }
   if((sigmahi == 0.0) && (sigmalo == 0.0))
   {
      *betahi = 0.0;
      *betalo = 0.0;
   }
   else
   {
      // mu = sqrt(x[0]*x[0] + sigma);
      ddf_sqr(xhi[0],xlo[0],&acchi,&acclo);
      ddf_inc(&acchi,&acclo,sigmahi,sigmalo);
      ddf_sqrt(acchi,acclo,&muhi,&mulo);

      if(xhi[0] <= 0.0)
      {
         // v[0] = x[0] - mu;
         ddf_sub(xhi[0],xlo[0],muhi,mulo,&vhi[0],&vlo[0]);
      }
      else
      {
         // v[0] = -sigma/(x[0] + mu);
         ddf_add(xhi[0],xlo[0],muhi,mulo,&acchi,&acclo);
         ddf_div(sigmahi,sigmalo,acchi,acclo,&vhi[0],&vlo[0]);
         ddf_minus(&vhi[0],&vlo[0]);
      }
      // v0p2 = v[0]*v[0];
      ddf_sqr(vhi[0],vlo[0],&v0p2hi,&v0p2lo);
      // *beta = 2.0*v0p2/(sigma + v0p2);
      ddf_add(sigmahi,sigmalo,v0p2hi,v0p2lo,&acchi,&acclo);
      ddf_div(v0p2hi,v0p2lo,acchi,acclo,betahi,betalo);
      ddf_mlt_d(betahi,betalo,2.0);
      
      for(int i=1; i<n; i++) // v[i] = v[i]/v[0];
      {
         ddf_div(vhi[i],vlo[i],vhi[0],vlo[0],&acchi,&acclo);
         vhi[i] = acchi;
         vlo[i] = acclo;
      }
      vhi[0] = 1.0;
      vlo[0] = 0.0;
   }
}

void CPU_cmplx2_factors_house 
( int n, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
         double *vrehi, double *vrelo, double *vimhi, double *vimlo,
  double *betahi, double *betalo )
{
   double sigmahi = 0.0;
   double sigmalo = 0.0;
   double acchi,acclo;
   double muhi,sqrx0hi,x0radhi,sqrv0hi;
   double mulo,sqrx0lo,x0radlo,sqrv0lo;
   double inv0rehi,inv0imhi,zrehi,zimhi;
   double inv0relo,inv0imlo,zrelo,zimlo;
   
   vrehi[0] = 1.0; vimhi[0] = 0.0;
   vrelo[0] = 0.0; vimlo[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      // sigma = sigma + xre[i]*xre[i] + xim[i]*xim[i];
      ddf_sqr(xrehi[i],xrelo[i],&acchi,&acclo);
      ddf_inc(&sigmahi,&sigmalo,acchi,acclo);
      ddf_sqr(ximhi[i],ximlo[i],&acchi,&acclo);
      ddf_inc(&sigmahi,&sigmalo,acchi,acclo);
      vrehi[i] = xrehi[i]; vimhi[i] = ximhi[i];
      vrelo[i] = xrelo[i]; vimlo[i] = ximlo[i];
   }
   if((sigmahi == 0.0) && (sigmalo == 0.0))
   {
      *betahi = 0.0;
      *betalo = 0.0;
   }
   else
   {
      // sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      ddf_sqr(xrehi[0],xrelo[0],&sqrx0hi,&sqrx0lo);
      ddf_sqr(ximhi[0],ximlo[0],&acchi,&acclo);
      ddf_inc(&sqrx0hi,&sqrx0lo,acchi,acclo);
      // x0rad = sqrt(sqrx0);
      ddf_sqrt(sqrx0hi,sqrx0lo,&x0radhi,&x0radlo);
      // cout << "x0rad : " << x0radhi << "  " << x0radlo << endl;
      // mu = sqrt(sqrx0 + sigma); // norm of the vector x
      ddf_inc(&sqrx0hi,&sqrx0lo,sigmahi,sigmalo);
      ddf_sqrt(sqrx0hi,sqrx0lo,&muhi,&mulo);

      if((x0radhi == 0.0) && (x0radlo == 0.0))
      {
         vrehi[0] = -muhi; vimhi[0] = 0.0;
         vrelo[0] = -mulo; vimlo[0] = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         // mu = mu/x0rad;
         ddf_div(muhi,mulo,x0radhi,x0radlo,&acchi,&acclo);
         muhi = acchi; mulo = acclo;
         // vre[0] = xre[0] - mu*xre[0];
         ddf_mul(muhi,mulo,xrehi[0],xrelo[0],&acchi,&acclo);
         ddf_sub(xrehi[0],xrelo[0],acchi,acclo,&vrehi[0],&vrelo[0]);
         // vim[0] = xim[0] - mu*xim[0];
         ddf_mul(muhi,mulo,ximhi[0],ximlo[0],&acchi,&acclo);
         ddf_sub(ximhi[0],ximlo[0],acchi,acclo,&vimhi[0],&vimlo[0]);
      }
      // cout << "mu : " << muhi << "  " << mulo << endl;

      // sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      ddf_sqr(vrehi[0],vrelo[0],&sqrv0hi,&sqrv0lo);
      // cout << "vre[0] : " << vrehi[0] << "  " << vrelo[0] << endl;
      // cout << "vim[0] : " << vimhi[0] << "  " << vimlo[0] << endl;
      ddf_sqr(vimhi[0],vimlo[0],&acchi,&acclo);
      ddf_inc(&sqrv0hi,&sqrv0lo,acchi,acclo);
      // cout << "sqrv0 : " << sqrv0hi << "  " << sqrv0lo << endl;
      // *beta = 2.0*sqrv0/(sigma + sqrv0);
      ddf_inc(&sigmahi,&sigmalo,sqrv0hi,sqrv0lo);
      ddf_div(sqrv0hi,sqrv0lo,sigmahi,sigmalo,betahi,betalo);
      ddf_mlt_d(betahi,betalo,2.0);
      // inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      ddf_div(vrehi[0],vrelo[0],sqrv0hi,sqrv0lo,&inv0rehi,&inv0relo);
      // inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]
      ddf_div(vimhi[0],vimlo[0],sqrv0hi,sqrv0lo,&inv0imhi,&inv0imlo);
      ddf_minus(&inv0imhi,&inv0imlo);

      for(int i=1; i<n; i++)  // v[i] = v[i]/v[0]
      {
         // zre = vre[i]*inv0re - vim[i]*inv0im;
         ddf_mul(vrehi[i],vrelo[i],inv0rehi,inv0relo,&zrehi,&zrelo);
         ddf_mul(vimhi[i],vimlo[i],inv0imhi,inv0imlo,&acchi,&acclo);
         ddf_dec(&zrehi,&zrelo,acchi,acclo);
         // zim = vim[i]*inv0re + vre[i]*inv0im;
         ddf_mul(vimhi[i],vimlo[i],inv0rehi,inv0relo,&zimhi,&zimlo);
         ddf_mul(vrehi[i],vrelo[i],inv0imhi,inv0imlo,&acchi,&acclo);
         ddf_inc(&zimhi,&zimlo,acchi,acclo);
         vrehi[i] = zrehi; vimhi[i] = zimhi;
         vrelo[i] = zrelo; vimlo[i] = zimlo;
      }
      vrehi[0] = 1.0; vimhi[0] = 0.0;
      vrelo[0] = 0.0; vimlo[0] = 0.0;
   }
}

void CPU_dbl2_factors_leftRupdate
 ( int nrows, int ncols, int k, double **Rhi, double **Rlo,
   double *vhi, double *vlo, double betahi, double betalo )
{
   double *whi = new double[ncols-k];
   double *wlo = new double[ncols-k];
   double acchi,acclo;

   for(int j=k; j<ncols; j++)
   {
      whi[j-k] = 0.0;
      wlo[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         ddf_mul(Rhi[i][j],Rlo[i][j],vhi[i-k],vlo[i-k],&acchi,&acclo);
         ddf_inc(&whi[j-k],&wlo[j-k],acchi,acclo);
      }
      // w[j-k] = beta*w[j-k];
      ddf_mul(betahi,betalo,whi[j-k],wlo[j-k],&acchi,&acclo);
      whi[j-k] = acchi;
      wlo[j-k] = acclo;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         ddf_mul(vhi[i-k],vlo[i-k],whi[j-k],wlo[j-k],&acchi,&acclo);
         ddf_dec(&Rhi[i][j],&Rlo[i][j],acchi,acclo);
      }

   free(whi); free(wlo);
}

void CPU_cmplx2_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double betahi, double betalo )
{
   double *wrehi = new double[ncols-k];
   double *wrelo = new double[ncols-k];
   double *wimhi = new double[ncols-k];
   double *wimlo = new double[ncols-k];
   double zrehi,zrelo,zimhi,zimlo,acchi,acclo;

   for(int j=k; j<ncols; j++)
   {
      wrehi[j-k] = 0.0; wimhi[j-k] = 0.0;
      wrelo[j-k] = 0.0; wimlo[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         // Hermitian of R => flip sign of Rim
         // zre =   Rre[i][j]*vre[i-k] + Rim[i][j]*vim[i-k];
         ddf_mul(Rrehi[i][j],Rrelo[i][j],vrehi[i-k],vrelo[i-k],&zrehi,&zrelo);
         ddf_mul(Rimhi[i][j],Rimlo[i][j],vimhi[i-k],vimlo[i-k],&acchi,&acclo);
         ddf_inc(&zrehi,&zrelo,acchi,acclo);
         // zim = - Rim[i][j]*vre[i-k] + Rre[i][j]*vim[i-k];
         ddf_mul(Rrehi[i][j],Rrelo[i][j],vimhi[i-k],vimlo[i-k],&zimhi,&zimlo);
         ddf_mul(Rimhi[i][j],Rimlo[i][j],vrehi[i-k],vrelo[i-k],&acchi,&acclo);
         ddf_dec(&zimhi,&zimlo,acchi,acclo);
         // wre[j-k] = wre[j-k] + zre;
         ddf_inc(&wrehi[j-k],&wrelo[j-k],zrehi,zrelo);
         // wim[j-k] = wim[j-k] + zim;
         ddf_inc(&wimhi[j-k],&wimlo[j-k],zimhi,zimlo);
      }
      // wre[j-k] = beta*wre[j-k];
      ddf_mul(betahi,betalo,wrehi[j-k],wrelo[j-k],&acchi,&acclo);
      wrehi[j-k] = acchi; wrelo[j-k] = acclo;
      // wim[j-k] = beta*wim[j-k];
      ddf_mul(betahi,betalo,wimhi[j-k],wimlo[j-k],&acchi,&acclo);
      wimhi[j-k] = acchi; wimlo[j-k] = acclo;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         // Hermitian of w => flip sign of wim
         // zre = vre[i-k]*wre[j-k] + vim[i-k]*wim[j-k];
         ddf_mul(vrehi[i-k],vrelo[i-k],wrehi[j-k],wrelo[j-k],&zrehi,&zrelo);
         ddf_mul(vimhi[i-k],vimlo[i-k],wimhi[j-k],wimlo[j-k],&acchi,&acclo);
         ddf_inc(&zrehi,&zrelo,acchi,acclo);
         // zim = vim[i-k]*wre[j-k] - vre[i-k]*wim[j-k];
         ddf_mul(vimhi[i-k],vimlo[i-k],wrehi[j-k],wrelo[j-k],&zimhi,&zimlo);
         ddf_mul(vrehi[i-k],vrelo[i-k],wimhi[j-k],wimlo[j-k],&acchi,&acclo);
         ddf_dec(&zimhi,&zimlo,acchi,acclo);
         // Rre[i][j] = Rre[i][j] - zre;
         ddf_dec(&Rrehi[i][j],&Rrelo[i][j],zrehi,zrelo);
         // Rim[i][j] = Rim[i][j] - zim;
         ddf_dec(&Rimhi[i][j],&Rimlo[i][j],zimhi,zimlo);
      }

   free(wrehi); free(wimhi);
   free(wrelo); free(wimlo);
}

void CPU_dbl2_factors_rightQupdate
 ( int n, int k, double **Qhi, double **Qlo,
   double *vhi, double *vlo, double betahi, double betalo )
{
   double *whi = new double[n];
   double *wlo = new double[n];
   double acchi,acclo;

   for(int i=0; i<n; i++)
   {
      whi[i] = 0.0;
      wlo[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         ddf_mul(Qhi[i][j],Qlo[i][j],vhi[j-k],vlo[j-k],&acchi,&acclo);
         ddf_inc(&whi[i],&wlo[i],acchi,acclo);
      }
      // w[i] = beta*w[i];
      ddf_mul(betahi,betalo,whi[i],wlo[i],&acchi,&acclo);
      whi[i] = acchi;
      wlo[i] = acclo;
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         ddf_mul(whi[i],wlo[i],vhi[j-k],vlo[j-k],&acchi,&acclo);
         ddf_dec(&Qhi[i][j],&Qlo[i][j],acchi,acclo);
      }

   free(whi); free(wlo);
}

void CPU_cmplx2_factors_rightQupdate
 ( int n, int k,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double betahi, double betalo )
{
   double *wrehi = new double[n];
   double *wrelo = new double[n];
   double *wimhi = new double[n];
   double *wimlo = new double[n];
   double zrehi,zrelo,zimhi,zimlo,acchi,acclo;

   for(int i=0; i<n; i++)
   {
      wrehi[i] = 0.0; wimhi[i] = 0.0;
      wrelo[i] = 0.0; wimlo[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         // zre = Qre[i][j]*vre[j-k] - Qim[i][j]*vim[j-k];
         ddf_mul(Qrehi[i][j],Qrelo[i][j],vrehi[j-k],vrelo[j-k],&zrehi,&zrelo);
         ddf_mul(Qimhi[i][j],Qimlo[i][j],vimhi[j-k],vimlo[j-k],&acchi,&acclo);
         ddf_dec(&zrehi,&zrelo,acchi,acclo);
         // zim = Qim[i][j]*vre[j-k] + Qre[i][j]*vim[j-k];
         ddf_mul(Qimhi[i][j],Qimlo[i][j],vrehi[j-k],vrelo[j-k],&zimhi,&zimlo);
         ddf_mul(Qrehi[i][j],Qrelo[i][j],vimhi[j-k],vimlo[j-k],&acchi,&acclo);
         ddf_inc(&zimhi,&zimlo,acchi,acclo);
         // wre[i] = wre[i] + zre;
         ddf_inc(&wrehi[i],&wrelo[i],zrehi,zrelo);
         // wim[i] = wim[i] + zim;
         ddf_inc(&wimhi[i],&wimlo[i],zimhi,zimlo);
      }
      // wre[i] = beta*wre[i];
      ddf_mul(betahi,betalo,wrehi[i],wrelo[i],&acchi,&acclo);
      wrehi[i] = acchi; wrelo[i] = acclo;
      // wim[i] = beta*wim[i];
      ddf_mul(betahi,betalo,wimhi[i],wimlo[i],&acchi,&acclo);
      wimhi[i] = acchi; wimlo[i] = acclo;
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         // Hermitian transpose => flip sign of vim
         // zre = wre[i]*vre[j-k] + wim[i]*vim[j-k];
         ddf_mul(wrehi[i],wrelo[i],vrehi[j-k],vrelo[j-k],&zrehi,&zrelo);
         ddf_mul(wimhi[i],wimlo[i],vimhi[j-k],vimlo[j-k],&acchi,&acclo);
         ddf_inc(&zrehi,&zrelo,acchi,acclo);
         // zim = wim[i]*vre[j-k] - wre[i]*vim[j-k];
         ddf_mul(wimhi[i],wimlo[i],vrehi[j-k],vrelo[j-k],&zimhi,&zimlo);
         ddf_mul(wrehi[i],wrelo[i],vimhi[j-k],vimlo[j-k],&acchi,&acclo);
         ddf_dec(&zimhi,&zimlo,acchi,acclo);
         // Qre[i][j] = Qre[i][j] - zre;
         ddf_dec(&Qrehi[i][j],&Qrelo[i][j],zrehi,zrelo);
         // Qim[i][j] = Qim[i][j] - zim;
         ddf_dec(&Qimhi[i][j],&Qimlo[i][j],zimhi,zimlo);
      }

   free(wrehi); free(wimhi);
   free(wrelo); free(wimlo);
}

void CPU_dbl2_factors_houseqr
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo )
{
   double *xhi = new double[nrows];
   double *xlo = new double[nrows];
   double *vhi = new double[nrows];
   double *vlo = new double[nrows];
   double betahi,betalo;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qhi[i][j] = 0.0;
         Qlo[i][j] = 0.0;
      }
      Qhi[i][i] = 1.0;
      Qlo[i][i] = 0.0;
      for(int j=0; j<ncols; j++)
      {
         Rhi[i][j] = Ahi[i][j];
         Rlo[i][j] = Alo[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xhi[i-k] = Rhi[i][k];
            xlo[i-k] = Rlo[i][k];
         }
         CPU_dbl2_factors_house(nrows-k,xhi,xlo,vhi,vlo,&betahi,&betalo);
         CPU_dbl2_factors_leftRupdate
            (nrows,ncols,k,Rhi,Rlo,vhi,vlo,betahi,betalo);
         CPU_dbl2_factors_rightQupdate(nrows,k,Qhi,Qlo,vhi,vlo,betahi,betalo);
      }
   }
   free(xhi); free(vhi);
   free(xlo); free(vlo);
}

void CPU_cmplx2_factors_houseqr
 ( int nrows, int ncols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo )
{
   double *xrehi = new double[nrows];
   double *xrelo = new double[nrows];
   double *ximhi = new double[nrows];
   double *ximlo = new double[nrows];
   double *vrehi = new double[nrows];
   double *vrelo = new double[nrows];
   double *vimhi = new double[nrows];
   double *vimlo = new double[nrows];
   double betahi,betalo;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qrehi[i][j] = 0.0; Qimhi[i][j] = 0.0;
         Qrelo[i][j] = 0.0; Qimlo[i][j] = 0.0;
      }
      Qrehi[i][i] = 1.0; Qimhi[i][i] = 0.0;
      Qrelo[i][i] = 0.0; Qimlo[i][i] = 0.0;

      for(int j=0; j<ncols; j++)
      {
         Rrehi[i][j] = Arehi[i][j]; Rimhi[i][j] = Aimhi[i][j];
         Rrelo[i][j] = Arelo[i][j]; Rimlo[i][j] = Aimlo[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xrehi[i-k] = Rrehi[i][k]; ximhi[i-k] = Rimhi[i][k];
            xrelo[i-k] = Rrelo[i][k]; ximlo[i-k] = Rimlo[i][k];
         }
         CPU_cmplx2_factors_house
            (nrows-k,xrehi,xrelo,ximhi,ximlo,
                     vrehi,vrelo,vimhi,vimlo,&betahi,&betalo);
         CPU_cmplx2_factors_leftRupdate
            (nrows,ncols,k,Rrehi,Rrelo,Rimhi,Rimlo,
                           vrehi,vrelo,vimhi,vimlo,betahi,betalo);
         CPU_cmplx2_factors_rightQupdate
            (nrows,k,Qrehi,Qrelo,Qimhi,Qimlo,
                     vrehi,vrelo,vimhi,vimlo,betahi,betalo);
      }
   }
   free(xrehi); free(ximhi); free(vrehi); free(vimhi);
   free(xrelo); free(ximlo); free(vrelo); free(vimlo);
}

void CPU_dbl2_factors_qrbs
 ( int nrows, int ncols,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *rhshi, double *rhslo, double *solhi, double *sollo,
   double *wrkvechi, double *wrkveclo )
{
   double acchi,acclo;

   for(int i=0; i<nrows; i++)   // compute Q^T*b, b is rhs
   {
      wrkvechi[i] = 0.0;
      wrkveclo[i] = 0.0;
      for(int j=0; j<nrows; j++) // wrkvec[i] = wrkvec[i] + Q[j][i]*rhs[j];
      {
         ddf_mul(Qhi[j][i],Qlo[j][i],rhshi[j],rhslo[j],&acchi,&acclo);
         ddf_inc(&wrkvechi[i],&wrkveclo[i],acchi,acclo);
      }
   }
   CPU_dbl2_factors_backward(ncols,Rhi,Rlo,wrkvechi,wrkveclo,solhi,sollo);
}

void CPU_cmplx2_factors_qrbs
 ( int nrows, int ncols,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo, 
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *rhsrehi, double *rhsrelo, double *rhsimhi, double *rhsimlo,
   double *solrehi, double *solrelo, double *solimhi, double *solimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo )
{
   double acchi,acclo; // accumulates product of two complex numbers

   for(int i=0; i<nrows; i++)    // compute Q^H*b, b is rhs
   {
      wrkvecrehi[i] = 0.0; wrkvecrelo[i] = 0.0;
      wrkvecimhi[i] = 0.0; wrkvecimlo[i] = 0.0;

      for(int j=0; j<nrows; j++) // work with Hermitian transpose of Q
      {
         // accre =  Qre[j][i]*rhsre[j] + Qim[j][i]*rhsim[j];
         // wrkvecre[i] = wrkvecre[i] + accre;
         ddf_mul(Qrehi[j][i],Qrelo[j][i],rhsrehi[j],rhsrelo[j],&acchi,&acclo);
         ddf_inc(&wrkvecrehi[i],&wrkvecrelo[i],acchi,acclo);
         ddf_mul(Qimhi[j][i],Qimlo[j][i],rhsimhi[j],rhsimlo[j],&acchi,&acclo);
         ddf_inc(&wrkvecrehi[i],&wrkvecrelo[i],acchi,acclo);
         // accim = -Qim[j][i]*rhsre[j] + Qre[j][i]*rhsim[j];
         // wrkvecim[i] = wrkvecim[i] + accim;
         ddf_mul(Qrehi[j][i],Qrelo[j][i],rhsimhi[j],rhsimlo[j],&acchi,&acclo);
         ddf_inc(&wrkvecimhi[i],&wrkvecimlo[i],acchi,acclo);
         ddf_mul(Qimhi[j][i],Qimlo[j][i],rhsrehi[j],rhsrelo[j],&acchi,&acclo);
         ddf_dec(&wrkvecimhi[i],&wrkvecimlo[i],acchi,acclo);
      }
   }
   CPU_cmplx2_factors_backward
      (ncols,Rrehi,Rrelo,Rimhi,Rimlo,wrkvecrehi,wrkvecrelo,
       wrkvecimhi,wrkvecimlo,solrehi,solrelo,solimhi,solimlo);
}
