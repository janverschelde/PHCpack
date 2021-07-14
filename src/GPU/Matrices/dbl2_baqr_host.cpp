/* The file dbl2_baqr_host.cpp defines the functions specified in
 * the file dbl2_baqr_host.h. */

#include <iostream>
#include <cstdlib>
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

void CPU_dbl2_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx, double **Chi, double **Clo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo )
{
   const int rowdim = nrows - idx*szt;      // number of rows in Y and W
   const int coldim = ncols - (idx+1)*szt;  // number of columns in C
   const int rowoff = idx*szt;              // row offset for C
   const int coloff = (idx+1)*szt;          // column offset for C

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

void CPU_dbl2_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Qhi, double **Qlo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo, bool verbose )
{
   const int rowdim = dim - idx*szt;  // number of rows in Y and W
                                      // is number of columns to update
   const int coloff = idx*szt;        // column offset for Q

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
      cout << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "prd[" << i << "][" << j << "] : "
                 << prdhi[i][j] << "  " << prdlo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)        // Q = Q + Q*W*Y^T
      for(int j=0; j<rowdim; j++)  // Q[i][coloff+j] += prd[i][j];
      {
         ddf_inc(&Qhi[i][coloff+j],&Qlo[i][coloff+j],prdhi[i][j],prdlo[i][j]);
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

void CPU_dbl2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo, bool verbose )
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
         CPU_dbl2_factors_leftRupdate
            (nrows,endcol,colidx,Rhi,Rlo,vhi,vlo,betahi,betalo);
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
            (nrows,ncols,szt,k,Rhi,Rlo,Yhi,Ylo,Whi,Wlo);
      }
      CPU_dbl2_blocked_rightQupdate
         (nrows,szt,k,Qhi,Qlo,Yhi,Ylo,Whi,Wlo,verbose);
   }
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
