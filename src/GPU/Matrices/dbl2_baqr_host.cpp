/* The file dbl2_baqr_host.cpp defines the functions specified in
 * the file dbl2_baqr_host.h. */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "double_double_functions.h"
#include "dbl2_factorizations.h"
#include "dbl2_baqr_host.h"

using namespace std;

void CPU_dbl2_blocked_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double **Vhi, double **Vlo, double **Whi, double **Wlo )
{
   double *vhi = new double[nrows];
   double *vlo = new double[nrows];
   double *phi = new double[ncols];
   double *plo = new double[ncols];
   double zihi,zilo,acchi,acclo;

   for(int i=0; i<nrows; i++)           //  W[0][i] = -B[0]*V[0][i];
   {
      ddf_mul(Bhi[0],Blo[0],Vhi[0][i],Vlo[0][i],&Whi[0][i],&Wlo[0][i]);
      ddf_minus(&Whi[0][i],&Wlo[0][i]);
   }

   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vhi[i] = Vhi[j][i];
         vlo[i] = Vlo[j][i];
      }
      for(int k=0; k<j; k++)
      {
         phi[k] = 0.0;                  // compute k-th component of Y^T*v
         plo[k] = 0.0;                  // over all rows of k-th column of Y
         for(int i=0; i<nrows; i++)     // p[k] = p[k] + V[k][i]*v[i]
         {
            ddf_mul(Vhi[k][i],Vlo[k][i],vhi[i],vlo[i],&acchi,&acclo);
            ddf_inc(&phi[k],&plo[k],acchi,acclo);
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zihi = 0.0;                    // compute i-th component of W*p
         zilo = 0.0;                    // over all rows of k-th column of W
         for(int k=0; k<j; k++)         // zi = zi + W[k][i]*p[k]
         {
            ddf_mul(Whi[k][i],Wlo[k][i],phi[k],plo[k],&acchi,&acclo);
            ddf_inc(&zihi,&zilo,acchi,acclo);
         }
         ddf_inc(&zihi,&zilo,vhi[i],vlo[i]);
         // W[j][i] = -B[j]*zi;
         ddf_mul(Bhi[j],Blo[j],zihi,zilo,&Whi[j][i],&Wlo[j][i]);
         ddf_minus(&Whi[j][i],&Wlo[j][i]);
      }
   }

   free(vhi); free(phi);
   free(vlo); free(plo);
}

void CPU_cmplx2_blocked_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double **Vrehi, double **Vrelo, double **Vimhi, double **Vimlo,
   double **Wrehi, double **Wrelo, double **Wimhi, double **Wimlo )
{
   double *vrehi = new double[nrows];
   double *vrelo = new double[nrows];
   double *vimhi = new double[nrows];
   double *vimlo = new double[nrows];
   double *prehi = new double[ncols];
   double *prelo = new double[ncols];
   double *pimhi = new double[ncols];
   double *pimlo = new double[ncols];
   double acchi,acclo;
   double zi_rehi,zi_relo,zi_imhi,zi_imlo;

   for(int i=0; i<nrows; i++) // W[0][i] = -B[0]*V[0][i]
   {
      ddf_mul(Bhi[0],Blo[0],Vrehi[0][i],Vrelo[0][i],
              &Wrehi[0][i],&Wrelo[0][i]);
      ddf_minus(&Wrehi[0][i],&Wrelo[0][i]);
      ddf_mul(Bhi[0],Blo[0],Vimhi[0][i],Vimlo[0][i],
              &Wimhi[0][i],&Wimlo[0][i]);
      ddf_minus(&Wimhi[0][i],&Wimlo[0][i]);
   }
   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vrehi[i] = Vrehi[j][i]; vimhi[i] = Vimhi[j][i];
         vrelo[i] = Vrelo[j][i]; vimlo[i] = Vimlo[j][i];
      }
      for(int k=0; k<j; k++)
      {
         prehi[k] = 0.0;                // compute k-th component of Y^H*v
         prelo[k] = 0.0;
         pimhi[k] = 0.0;
         pimlo[k] = 0.0;

         for(int i=0; i<nrows; i++)    // p[k] = p[k] + V[k][i]*v[i];
         {
            // accre =   Vre[k][i]*vre[i] + Vim[k][i]*vim[i];
            // pre[k] = pre[k] + accre;
            ddf_mul(Vrehi[k][i],Vrelo[k][i],vrehi[i],vrelo[i],
                    &acchi,&acclo);
            ddf_inc(&prehi[k],&prelo[k],acchi,acclo);
            ddf_mul(Vimhi[k][i],Vimlo[k][i],vimhi[i],vimlo[i],
                    &acchi,&acclo);
            ddf_inc(&prehi[k],&prelo[k],acchi,acclo);
            // accim = - Vim[k][i]*vre[i] + Vre[k][i]*vim[i];
            // pim[k] = pim[k] + accim;
            ddf_mul(Vimhi[k][i],Vimlo[k][i],vrehi[i],vrelo[i],
                    &acchi,&acclo);
            ddf_dec(&pimhi[k],&pimlo[k],acchi,acclo);
            ddf_mul(Vrehi[k][i],Vrelo[k][i],vimhi[i],vimlo[i],
                    &acchi,&acclo);
            ddf_inc(&pimhi[k],&pimlo[k],acchi,acclo);
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zi_rehi = 0.0;                 // compute i-th component of W*p
         zi_relo = 0.0;
         zi_imhi = 0.0;
         zi_imlo = 0.0;

         for(int k=0; k<j; k++)         // zi = zi + W[k][i]*p[k];
         {
            // accre = Wre[k][i]*pre[k] - Wim[k][i]*pim[k];
            // zi_re = zi_re + accre;
            ddf_mul(Wrehi[k][i],Wrelo[k][i],prehi[k],prelo[k],
                    &acchi,&acclo);
            ddf_inc(&zi_rehi,&zi_relo,acchi,acclo);
            ddf_mul(Wimhi[k][i],Wimlo[k][i],pimhi[k],pimlo[k],
                    &acchi,&acclo);
            ddf_dec(&zi_rehi,&zi_relo,acchi,acclo);
            // accim = Wim[k][i]*pre[k] + Wre[k][i]*pim[k];
            // zi_im = zi_im + accim;
            ddf_mul(Wimhi[k][i],Wimlo[k][i],prehi[k],prelo[k],
                    &acchi,&acclo);
            ddf_inc(&zi_imhi,&zi_imlo,acchi,acclo);
            ddf_mul(Wrehi[k][i],Wrelo[k][i],pimhi[k],pimlo[k],
                    &acchi,&acclo);
            ddf_inc(&zi_imhi,&zi_imlo,acchi,acclo);
         }
         // zi_re = zi_re + vre[i];
         ddf_inc(&zi_rehi,&zi_relo,vrehi[i],vrelo[i]);
         // zi_im = zi_im + vim[i];
         ddf_inc(&zi_imhi,&zi_imlo,vimhi[i],vimlo[i]);
         // Wre[j][i] = -B[j]*zi_re;
         ddf_mul(Bhi[j],Blo[j],zi_rehi,zi_relo,&Wrehi[j][i],&Wrelo[j][i]);
         ddf_minus(&Wrehi[j][i],&Wrelo[j][i]);
         // Wim[j][i] = -B[j]*zi_im;
         ddf_mul(Bhi[j],Blo[j],zi_imhi,zi_imlo,&Wimhi[j][i],&Wimlo[j][i]);
         ddf_minus(&Wimhi[j][i],&Wimlo[j][i]);
      }
   }
   free(vrehi); free(prehi);
   free(vrelo); free(prelo);
   free(vimhi); free(pimhi);
   free(vimlo); free(pimlo);
}

void CPU_dbl2_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx, double **Chi, double **Clo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo, bool verbose )
{
   const int rowoff = idx*szt;            // row offset for C
   const int rowdim = nrows - rowoff;     // number of rows in Y and W
   const int coloff = (idx+1)*szt;        // column offset for C
   const int coldim = ncols - coloff;     // number of columns in C

   if(verbose)
   {
      cout << "updating R ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff << endl;
   }

   double **YWThi = new double*[rowdim];
   double **YWTlo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      YWThi[i] = new double[rowdim];
      YWTlo[i] = new double[rowdim];
   }
   double **prdhi = new double*[rowdim];
   double **prdlo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      prdhi[i] = new double[coldim];
      prdlo[i] = new double[coldim];
   }
   double acchi,acclo;

   for(int i=0; i<rowdim; i++)        // compute Y*W^T
      for(int j=0; j<rowdim; j++)
      {
         YWThi[i][j] = 0.0;           // row i of Y with column j of W^T
         YWTlo[i][j] = 0.0;
         for(int k=0; k<szt; k++)     // YWT[i][j] += Y[k][i]*W[k][j]
         {
            ddf_mul(Yhi[k][i],Ylo[k][i],Whi[k][j],Wlo[k][j],&acchi,&acclo);
            ddf_inc(&YWThi[i][j],&YWTlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThi[i][j] << "  " << YWTlo[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         prdhi[i][j] = 0.0;
         prdlo[i][j] = 0.0;           // prd[i][j]
         for(int k=0; k<rowdim; k++)  // += YWT[i][k]*C[rowoff+k][coloff+j]
         {
            ddf_mul(YWThi[i][k],            YWTlo[i][k],
                      Chi[rowoff+k][coloff+j],Clo[rowoff+k][coloff+j],
                    &acchi,&acclo);
            ddf_inc(&prdhi[i][j],&prdlo[i][j],acchi,acclo);
         }
      }

   for(int i=0; i<rowdim; i++)      // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)   // C[rowoff+i][coloff+j] += prd[i][j]
      {
         ddf_inc( &Chi[rowoff+i][coloff+j],&Clo[rowoff+i][coloff+j],
                 prdhi[i][j],             prdlo[i][j]);
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWThi[i]); free(prdhi[i]);
      free(YWTlo[i]); free(prdlo[i]);
   }
   free(YWThi); free(prdhi);
   free(YWTlo); free(prdlo);
}

void CPU_cmplx2_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Crehi, double **Crelo, double **Cimhi, double **Cimlo,
   double **Yrehi, double **Yrelo, double **Yimhi, double **Yimlo,
   double **Wrehi, double **Wrelo, double **Wimhi, double **Wimlo,
   bool verbose )
{
   const int rowoff = idx*szt;             // row offset for C
   const int rowdim = nrows - rowoff;      // number of rows in Y and W
   const int coloff = (idx+1)*szt;         // column offset for C
   const int coldim = ncols - coloff;      // number of columns in C
   double acchi,acclo;

   if(verbose)
   {
      cout << "updating R ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff << endl;
   }

   double **YWTrehi = new double*[rowdim];
   double **YWTrelo = new double*[rowdim];
   double **YWTimhi = new double*[rowdim];
   double **YWTimlo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      YWTrehi[i] = new double[rowdim];
      YWTrelo[i] = new double[rowdim];
      YWTimhi[i] = new double[rowdim];
      YWTimlo[i] = new double[rowdim];
   }
   double **prdrehi = new double*[rowdim];
   double **prdrelo = new double*[rowdim];
   double **prdimhi = new double*[rowdim];
   double **prdimlo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      prdrehi[i] = new double[coldim];
      prdrelo[i] = new double[coldim];
      prdimhi[i] = new double[coldim];
      prdimlo[i] = new double[coldim];
   }
   for(int i=0; i<rowdim; i++)      // compute Y*W^H
      for(int j=0; j<rowdim; j++)
      {
         YWTrehi[i][j] = 0.0;       // row i of Y with column j of W^H
         YWTrelo[i][j] = 0.0;
         YWTimhi[i][j] = 0.0;
         YWTimlo[i][j] = 0.0;

         for(int k=0; k<szt; k++)   // YWT[i][j] = YWT[i][j] + Y[k][i]*W[k][j]
         {  
            // accre = Yre[k][i]*Wre[k][j] + Yim[k][i]*Wim[k][j];
            // YWTre[i][j] = YWTre[i][j] + accre;
            ddf_mul(Yrehi[k][i],Yrelo[k][i],Wrehi[k][j],Wrelo[k][j],
                    &acchi,&acclo);
            ddf_inc(&YWTrehi[i][j],&YWTrelo[i][j],acchi,acclo);
            ddf_mul(Yimhi[k][i],Yimlo[k][i],Wimhi[k][j],Wimlo[k][j],
                    &acchi,&acclo);
            ddf_inc(&YWTrehi[i][j],&YWTrelo[i][j],acchi,acclo);
            // accim = Yim[k][i]*Wre[k][j] - Yre[k][i]*Wim[k][j];
            // YWTim[i][j] = YWTim[i][j] + accim;
            ddf_mul(Yimhi[k][i],Yimlo[k][i],Wrehi[k][j],Wrelo[k][j],
                    &acchi,&acclo);
            ddf_inc(&YWTimhi[i][j],&YWTimlo[i][j],acchi,acclo);
            ddf_mul(Yrehi[k][i],Yrelo[k][i],Wimhi[k][j],Wimlo[k][j],
                    &acchi,&acclo);
            ddf_dec(&YWTimhi[i][j],&YWTimlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "YWT[" << i << "][" << j << "]re : "
                 << YWTrehi[i][j] << "  " << YWTrelo[i][j] << endl;
            cout << "YWT[" << i << "][" << j << "]im : "
                 << YWTimhi[i][j] << "  " << YWTimlo[i][j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^H)*C
      for(int j=0; j<coldim; j++)
      {
         prdrehi[i][j] = 0.0;
         prdrelo[i][j] = 0.0;
         prdimhi[i][j] = 0.0;
         prdimlo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++)
         {  // prd[i][j] = prd[i][j] + YWT[i][k]*C[rowoff+k][coloff+j];
            // accre = YWTre[i][k]*Cre[rowoff+k][coloff+j]
            //       - YWTim[i][k]*Cim[rowoff+k][coloff+j];
            // prdre[i][j] = prdre[i][j] + accre;
            ddf_mul(YWTrehi[i][k],            YWTrelo[i][k],
                      Crehi[rowoff+k][coloff+j],Crelo[rowoff+k][coloff+j],
                    &acchi,&acclo);
            ddf_inc(&prdrehi[i][j],&prdrelo[i][j],acchi,acclo);
            ddf_mul(YWTimhi[i][k],            YWTimlo[i][k],
                      Cimhi[rowoff+k][coloff+j],Cimlo[rowoff+k][coloff+j],
                    &acchi,&acclo);
            ddf_dec(&prdrehi[i][j],&prdrelo[i][j],acchi,acclo);
            // accim = YWTim[i][k]*Cre[rowoff+k][coloff+j]
            //       + YWTre[i][k]*Cim[rowoff+k][coloff+j];
            // prdim[i][j] = prdim[i][j] + accim;
            ddf_mul(YWTimhi[i][k],            YWTimlo[i][k],
                      Crehi[rowoff+k][coloff+j],Crelo[rowoff+k][coloff+j],
                    &acchi,&acclo);
            ddf_inc(&prdimhi[i][j],&prdimlo[i][j],acchi,acclo);
            ddf_mul(YWTrehi[i][k],            YWTrelo[i][k],
                      Cimhi[rowoff+k][coloff+j],Cimlo[rowoff+k][coloff+j],
                    &acchi,&acclo);
            ddf_inc(&prdimhi[i][j],&prdimlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      cout << endl;
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
         {
            cout << "prd[" << i << "][" << j << "]re : "
                 << prdrehi[i][j] << "  " << prdrelo[i][j] << endl;
            cout << "prd[" << i << "][" << j << "]im : "
                 << prdimlo[i][j] << "  " << prdimlo[i][j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)        // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         ddf_inc(&Crehi[rowoff+i][coloff+j],
                 &Crelo[rowoff+i][coloff+j],prdrehi[i][j],prdrelo[i][j]);
         ddf_inc(&Cimhi[rowoff+i][coloff+j],
                 &Cimlo[rowoff+i][coloff+j],prdimhi[i][j],prdimlo[i][j]);
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWTrehi[i]); free(prdrehi[i]);
      free(YWTrelo[i]); free(prdrelo[i]);
      free(YWTimhi[i]); free(prdimhi[i]);
      free(YWTimlo[i]); free(prdimlo[i]);
   }
   free(YWTrehi); free(prdrehi);
   free(YWTrelo); free(prdrelo);
   free(YWTimhi); free(prdimhi);
   free(YWTimlo); free(prdimlo);
}

void CPU_dbl2_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Qhi, double **Qlo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo, bool verbose )
{
   const int coloff = idx*szt;        // column offset for Q
   const int rowdim = dim - coloff;
   // the number of rows in Y and W is the number of columns to update

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;
   }
   double **WYThi = new double*[rowdim];
   double **WYTlo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      WYThi[i] = new double[rowdim];
      WYTlo[i] = new double[rowdim];
   }
   double **prdhi = new double*[dim];
   double **prdlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdhi[i] = new double[rowdim];
      prdlo[i] = new double[rowdim];
   }
   double acchi,acclo;

   for(int i=0; i<rowdim; i++)     // compute W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         WYThi[i][j] = 0.0;        // row i of W with column j of Y^T
         WYTlo[i][j] = 0.0;        // take k-th column of W
         for(int k=0; k<szt; k++)  // WYT[i][j] = WYT[i][j] + W[k][i]*Y[k][j]
         {
            ddf_mul(Whi[k][i],Wlo[k][i],Yhi[k][j],Ylo[k][j],&acchi,&acclo);
            ddf_inc(&WYThi[i][j],&WYTlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "Y[" << i << "][" << j << "] : "
                 << Yhi[j][i] << "  " << Ylo[j][i] << endl;

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "W[" << i << "][" << j << "] : "
                 << Whi[j][i] << "  " << Wlo[j][i] << endl;

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThi[i][j] << "  " << WYTlo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)          // prd = Q*W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         prdhi[i][j] = 0.0;
         prdlo[i][j] = 0.0;
         for(int k=0; k<rowdim; k++) // prd[i][j] += Q[i][coloff+k]*WYT[k][j]
         {
            ddf_mul(  Qhi[i][coloff+k],Qlo[i][coloff+k],
                    WYThi[k][j],     WYTlo[k][j],&acchi,&acclo);
            ddf_inc(&prdhi[i][j],&prdlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "QWYT[" << i << "][" << j << "] : "
                 << prdhi[i][j] << "  " << prdlo[i][j] << endl;

      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Qhi[i][coloff+j] << "  "
                 << Qlo[i][coloff+j] << endl;
   }
   for(int i=0; i<dim; i++)        // Q = Q + Q*W*Y^T
      for(int j=0; j<rowdim; j++)  // Q[i][coloff+j] += prd[i][j];
      {
         ddf_inc(&Qhi[i][coloff+j],&Qlo[i][coloff+j],prdhi[i][j],prdlo[i][j]);
      }

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Qhi[i][coloff+j] << "  "
                 << Qlo[i][coloff+j] << endl;
   }
   for(int i=0; i<rowdim; i++)
   {
      free(WYThi[i]);
      free(WYTlo[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdhi[i]);
      free(prdlo[i]);
   }
   free(WYThi); free(prdhi);
   free(WYTlo); free(prdlo);
}

void CPU_cmplx2_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Yrehi, double **Yrelo, double **Yimhi, double **Yimlo,
   double **Wrehi, double **Wrelo, double **Wimhi, double **Wimlo,
   bool verbose )
{
   const int coloff = idx*szt;        // column offset for Q
   const int rowdim = dim - coloff;
   // the number of rows in Y and W is the number of columns to update
   double acchi,acclo;

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;
   }
   double **WYTrehi = new double*[rowdim];
   double **WYTrelo = new double*[rowdim];
   double **WYTimhi = new double*[rowdim];
   double **WYTimlo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      WYTrehi[i] = new double[rowdim];
      WYTrelo[i] = new double[rowdim];
      WYTimhi[i] = new double[rowdim];
      WYTimlo[i] = new double[rowdim];
   }
   double **prdrehi = new double*[dim];
   double **prdrelo = new double*[dim];
   double **prdimhi = new double*[dim];
   double **prdimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrehi[i] = new double[rowdim];
      prdrelo[i] = new double[rowdim];
      prdimhi[i] = new double[rowdim];
      prdimlo[i] = new double[rowdim];
   }
   for(int i=0; i<rowdim; i++)        // compute W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         WYTrehi[i][j] = 0.0;         // row i of W with column j of Y^H
         WYTrelo[i][j] = 0.0;
         WYTimhi[i][j] = 0.0;
         WYTimlo[i][j] = 0.0;

         for(int k=0; k<szt; k++)     // WYT[i][j] += W[k][i]*Y[k][j]
         {
            // accre = Wre[k][i]*Yre[k][j] + Wim[k][i]*Yim[k][j];
            // WYTre[i][j] = WYTre[i][j] + accre;
            ddf_mul(Wrehi[k][i],Wrelo[k][i],Yrehi[k][j],Yrelo[k][j],
                    &acchi,&acclo);
            ddf_inc(&WYTrehi[i][j],&WYTrelo[i][j],acchi,acclo);
            ddf_mul(Wimhi[k][i],Wimlo[k][i],Yimhi[k][j],Yimlo[k][j],
                    &acchi,&acclo);
            ddf_inc(&WYTrehi[i][j],&WYTrelo[i][j],acchi,acclo);
            // accim = Wim[k][i]*Yre[k][j] - Wre[k][i]*Yim[k][j];
            // WYTim[i][j] = WYTim[i][j] + accim;
            ddf_mul(Wimhi[k][i],Wimlo[k][i],Yrehi[k][j],Yrelo[k][j],
                    &acchi,&acclo);
            ddf_inc(&WYTimhi[i][j],&WYTimlo[i][j],acchi,acclo);
            ddf_mul(Wrehi[k][i],Wrelo[k][i],Yimhi[k][j],Yimlo[k][j],
                    &acchi,&acclo);
            ddf_dec(&WYTimhi[i][j],&WYTimlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
         {
            cout << "Y[" << i << "][" << j << "]re : "
                 << Yrehi[j][i] << "  " << Yrelo[j][i] << endl;
            cout << "Y[" << i << "][" << j << "]im : "
                 << Yimhi[j][i] << "  " << Yimlo[j][i] << endl;
         }

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehi[j][i] << "  " << Wrelo[j][i] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhi[j][i] << "  " << Wimlo[j][i] << endl;
         }

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "WYT[" << i << "][" << j << "]re : "
                 << WYTrehi[i][j] << "  " << WYTrelo[i][j] << endl;
            cout << "WYT[" << i << "][" << j << "]im : "
                 << WYTimhi[i][j] << "  " << WYTimlo[i][j] << endl;
         }
   }
   for(int i=0; i<dim; i++)           // prd = Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         prdrehi[i][j] = 0.0;
         prdrelo[i][j] = 0.0;
         prdimhi[i][j] = 0.0;
         prdimlo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++) // prd[i][j] += Q[i][coloff+k]*WYT[k][j]
         {
            // accre = Qre[i][coloff+k]*WYTre[k][j]
            //       - Qim[i][coloff+k]*WYTim[k][j];
            // prdre[i][j] = prdre[i][j] + accre;
            ddf_mul(   Qrehi[i][coloff+k],Qrelo[i][coloff+k],
                     WYTrehi[k][j],     WYTrelo[k][j],&acchi,&acclo);
            ddf_inc(&prdrehi[i][j],&prdrelo[i][j],acchi,acclo);
            ddf_mul(   Qimhi[i][coloff+k],Qimlo[i][coloff+k],
                     WYTimhi[k][j],     WYTimlo[k][j],&acchi,&acclo);
            ddf_dec(&prdrehi[i][j],&prdrelo[i][j],acchi,acclo);
            // accim = Qim[i][coloff+k]*WYTre[k][j]
            //       + Qre[i][coloff+k]*WYTim[k][j];
            // prdim[i][j] = prdim[i][j] + accim;
            ddf_mul(   Qimhi[i][coloff+k],Qimlo[i][coloff+k],
                     WYTrehi[k][j],     WYTrelo[k][j],&acchi,&acclo);
            ddf_inc(&prdimhi[i][j],&prdimlo[i][j],acchi,acclo);
            ddf_mul(   Qrehi[i][coloff+k],Qrelo[i][coloff+k],
                     WYTimhi[k][j],     WYTimlo[k][j],&acchi,&acclo);
            ddf_inc(&prdimhi[i][j],&prdimlo[i][j],acchi,acclo);
         }
      }

   if(verbose)
   {
      cout << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "QWYT[" << i << "][" << j << "]re : "
                 << prdrehi[i][j] << "  " << prdrelo[i][j] << endl;
            cout << "QWYT[" << i << "][" << j << "]im : "
                 << prdimlo[i][j] << "  " << prdimlo[i][j] << endl;
         }

      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "Q[" << i << "][" << coloff+j << "]re : "
                 << Qrehi[i][coloff+j] << "  "
                 << Qrelo[i][coloff+j] << endl;
            cout << "Q[" << i << "][" << coloff+j << "]im : "
                 << Qimhi[i][coloff+j] << "  "
                 << Qimlo[i][coloff+j] << endl;
         }
   }
   for(int i=0; i<dim; i++)           // Q = Q + Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         // Qre[i][coloff+j] = Qre[i][coloff+j] + prdre[i][j];
         ddf_inc(&Qrehi[i][coloff+j],
                 &Qrelo[i][coloff+j],prdrehi[i][j],prdrelo[i][j]);
         // Qim[i][coloff+j] = Qim[i][coloff+j] + prdim[i][j];
         ddf_inc(&Qimhi[i][coloff+j],
                 &Qimlo[i][coloff+j],prdimhi[i][j],prdimlo[i][j]);
      }

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "Q[" << i << "][" << coloff+j << "]re : "
                 << Qrehi[i][coloff+j] << "  "
                 << Qrelo[i][coloff+j] << endl;
            cout << "Q[" << i << "][" << coloff+j << "]im : "
                 << Qimhi[i][coloff+j] << "  "
                 << Qimlo[i][coloff+j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)
   {
      free(WYTrehi[i]); free(WYTimhi[i]);
      free(WYTrelo[i]); free(WYTimlo[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrehi[i]); free(prdimhi[i]);
      free(prdrelo[i]); free(prdimlo[i]);
   }
   free(WYTrehi); free(prdrehi);
   free(WYTrelo); free(prdrelo);
   free(WYTimhi); free(prdimhi);
   free(WYTimlo); free(prdimlo);
}

void CPU_dbl2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *lapsec, bool verbose )
{
   double betahi,betalo;
   double *xhi = new double[nrows]; // input vector for house
   double *xlo = new double[nrows]; // low doubles of the input vector
   double *vhi = new double[nrows]; // Householder vector
   double *vlo = new double[nrows]; // low doubles of the Householder vector
   double *Bhi = new double[szt];   // high doubles of the betas
   double *Blo = new double[szt];   // low doubles of the betas
   double **Yhi = new double*[szt]; // Householder vectors in one block
   double **Ylo = new double*[szt]; // low doubles of Householder vectors
   double **Whi = new double*[szt]; // columns of W
   double **Wlo = new double*[szt]; // low doubles of the columns of W

   for(int j=0; j<szt; j++)
   {
      Yhi[j] = new double[nrows];
      Ylo[j] = new double[nrows];
      Whi[j] = new double[nrows];
      Wlo[j] = new double[nrows];
   }
   clock_t start = clock();

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
   int colidx,endcol,nrowscol;

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      for(int L=0; L<szt; L++)    // L runs over the columns in one block
      {
         colidx = k*szt + L;      // index of the current column
         endcol = (k+1)*szt;      // 1 + last column index in current block
         nrowscol = nrows-colidx; // number of rows in current column

         for(int i=colidx; i<nrows; i++)
         {
            xhi[i-colidx] = Rhi[i][colidx];
            xlo[i-colidx] = Rlo[i][colidx];
         }
         CPU_dbl2_factors_house(nrowscol,xhi,xlo,vhi,vlo,&betahi,&betalo);
         if(verbose)
         {
            cout << "beta[" << colidx << "] : "
                 << betahi << "  " << betalo << endl;
            for(int i=colidx; i<nrows; i++)
               cout << "v[" << i-colidx << "] : "
                    << vhi[i-colidx] << "  " << vlo[i-colidx] << endl;
            cout << "the R matrix :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rhi[i][j] << "  " << Rlo[i][j] << endl;
         }
         CPU_dbl2_factors_leftRupdate
            (nrows,endcol,colidx,Rhi,Rlo,vhi,vlo,betahi,betalo);
         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rhi[i][j] << "  " << Rlo[i][j] << endl;
         }
         Bhi[L] = betahi;
         Blo[L] = betalo;
         for(int i=0; i<L; i++)
         {
            Yhi[L][i] = 0.0;
            Ylo[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yhi[L][i+L] = vhi[i];
            Ylo[L][i+L] = vlo[i];
         }
      }
      CPU_dbl2_blocked_VB_to_W(nrows-k*szt,szt,Bhi,Blo,Yhi,Ylo,Whi,Wlo);
      if(k<nbt-1)
      {
         CPU_dbl2_blocked_leftRupdate
            (nrows,ncols,szt,k,Rhi,Rlo,Yhi,Ylo,Whi,Wlo,verbose);
      }
      CPU_dbl2_blocked_rightQupdate
         (nrows,szt,k,Qhi,Qlo,Yhi,Ylo,Whi,Wlo,verbose);
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(xhi); free(vhi); free(Bhi);
   free(xlo); free(vlo); free(Blo);

   for(int j=0; j<szt; j++)
   {
      free(Yhi[j]); free(Whi[j]);
      free(Ylo[j]); free(Wlo[j]);
   }
   free(Yhi); free(Whi);
   free(Ylo); free(Wlo);
}

void CPU_cmplx2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *lapsec, bool verbose )
{
   double betahi,betalo;
   double *xrehi = new double[nrows]; // real input vector for house
   double *xrelo = new double[nrows]; // low doubles of real input
   double *ximhi = new double[nrows]; // imaginary input vector for house
   double *ximlo = new double[nrows]; // low doubles of imaginary input
   double *vrehi = new double[nrows]; // real parts of a Householder vector
   double *vrelo = new double[nrows]; // low doubles of vre
   double *vimhi = new double[nrows]; // imag parts of a Householder vector
   double *vimlo = new double[nrows]; // low doubles of vim
   double *Bhi = new double[szt];     // the high doubles of the betas
   double *Blo = new double[szt];     // the low doubles of the betas
   double **Yrehi = new double*[szt]; // Householder vectors in one block
   double **Yrelo = new double*[szt]; // low doubles of Yre
   double **Yimhi = new double*[szt]; // imag parts of Householder vectors
   double **Yimlo = new double*[szt]; // low doubles of Yim
   double **Wrehi = new double*[szt]; // real parts of the columns of W
   double **Wrelo = new double*[szt]; // low doubles of Wre
   double **Wimhi = new double*[szt]; // imaginary parts of the columns of W
   double **Wimlo = new double*[szt]; // low doubles of Wim

   for(int j=0; j<szt; j++)
   {
      Yrehi[j] = new double[nrows];
      Yrelo[j] = new double[nrows];
      Yimhi[j] = new double[nrows];
      Yimlo[j] = new double[nrows];
      Wrehi[j] = new double[nrows];
      Wrelo[j] = new double[nrows];
      Wimhi[j] = new double[nrows];
      Wimlo[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qrehi[i][j] = 0.0; Qrelo[i][j] = 0.0;
         Qimhi[i][j] = 0.0; Qimlo[i][j] = 0.0;
      }
      Qrehi[i][i] = 1.0; Qrelo[i][i] = 0.0;
      Qimhi[i][i] = 0.0; Qimlo[i][i] = 0.0;
      for(int j=0; j<ncols; j++)
      {
         Rrehi[i][j] = Arehi[i][j]; Rrelo[i][j] = Arelo[i][j];
         Rimhi[i][j] = Aimhi[i][j]; Rimlo[i][j] = Aimlo[i][j];
      }
   }
   int colidx,endcol,nrowscol;

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      for(int L=0; L<szt; L++)    // L runs over the columns in one block
      {
         colidx = k*szt + L;      // index of the current column
         endcol = (k+1)*szt;      // 1 + last column index in current block
         nrowscol = nrows-colidx; // number of rows in current column

         for(int i=colidx; i<nrows; i++)
         {
            xrehi[i-colidx] = Rrehi[i][colidx];
            xrelo[i-colidx] = Rrelo[i][colidx];
            ximhi[i-colidx] = Rimhi[i][colidx];
            ximlo[i-colidx] = Rimlo[i][colidx];
         }
         CPU_cmplx2_factors_house
            (nrowscol,xrehi,xrelo,ximhi,ximlo,
                      vrehi,vrelo,vimhi,vimlo,&betahi,&betalo);
         if(verbose)
         {
            cout << "beta[" << colidx << "] : "
                 << betahi << "  " << betalo << endl;
            for(int i=colidx; i<nrows; i++)
            {
               cout << "vre[" << i-colidx << "] : "
                    << vrehi[i-colidx] << "  " << vrelo[i-colidx] << endl;
               cout << "vim[" << i-colidx << "] : "
                    << vimhi[i-colidx] << "  " << vimlo[i-colidx] << endl;
            }
         }
         CPU_cmplx2_factors_leftRupdate
            (nrows,endcol,colidx,Rrehi,Rrelo,Rimhi,Rimlo,
                                 vrehi,vrelo,vimhi,vimlo,betahi,betalo);
         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
               {
                  cout << "Rre[" << i << "][" << j << "] : "
                       << Rrehi[i][j] << "  " << Rrelo[i][j] << endl;
                  cout << "Rim[" << i << "][" << j << "] : "
                       << Rimhi[i][j] << "  " << Rimlo[i][j] << endl;
               }
         }
         Bhi[L] = betahi;
         Blo[L] = betalo;
         for(int i=0; i<L; i++)
         {
            Yrehi[L][i] = 0.0;
            Yrelo[L][i] = 0.0;
            Yimhi[L][i] = 0.0;
            Yimlo[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yrehi[L][i+L] = vrehi[i];
            Yrelo[L][i+L] = vrelo[i];
            Yimhi[L][i+L] = vimhi[i];
            Yimlo[L][i+L] = vimlo[i];
         }
      }
      CPU_cmplx2_blocked_VB_to_W
         (nrows-k*szt,szt,Bhi,Blo,Yrehi,Yrelo,Yimhi,Yimlo,
                                  Wrehi,Wrelo,Wimhi,Wimlo);
      if(k<nbt-1)
      {
         CPU_cmplx2_blocked_leftRupdate
            (nrows,ncols,szt,k,Rrehi,Rrelo,Rimhi,Rimlo,
                               Yrehi,Yrelo,Yimhi,Yimlo,
                               Wrehi,Wrelo,Wimhi,Wimlo,verbose);
      }
      CPU_cmplx2_blocked_rightQupdate
         (nrows,szt,k,Qrehi,Qrelo,Qimhi,Qimlo,
                      Yrehi,Yrelo,Yimhi,Yimlo,
                      Wrehi,Wrelo,Wimhi,Wimlo,verbose);
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(xrehi); free(vrehi); free(Bhi);
   free(xrelo); free(vrelo); free(Blo);
   free(ximhi); free(vimhi);
   free(ximlo); free(vimlo);

   for(int j=0; j<szt; j++)
   {
      free(Yrehi[j]); free(Wrehi[j]);
      free(Yrelo[j]); free(Wrelo[j]);
      free(Yimhi[j]); free(Wimhi[j]);
      free(Yimlo[j]); free(Wimlo[j]);
   }
   free(Yrehi); free(Wrehi);
   free(Yrelo); free(Wrelo);
   free(Yimhi); free(Wimhi);
   free(Yimlo); free(Wimlo);
}
