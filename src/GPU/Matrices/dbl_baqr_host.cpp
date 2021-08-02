/* The file dbl_baqr_host.cpp defines the functions specified in
 * the file dbl_baqr_host.h. */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "dbl_factorizations.h"
#include "dbl_baqr_host.h"

using namespace std;

void CPU_dbl_blocked_VB_to_W
 ( int nrows, int ncols, double *B, double **V, double **W )
{
   double *v = new double[nrows];
   double *p = new double[ncols];
   double zi;

   for(int i=0; i<nrows; i++) W[0][i] = -B[0]*V[0][i];

   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++) v[i] = V[j][i];
      for(int k=0; k<j; k++)
      {
         p[k] = 0.0;                    // compute k-th component of Y^T*v
         for(int i=0; i<nrows; i++)     // over all rows of k-th column of Y
            p[k] = p[k] + V[k][i]*v[i];
      }
      for(int i=0; i<nrows; i++)
      {
         zi = 0.0;                      // compute i-th component of W*p
         for(int k=0; k<j; k++)         // over all rows of k-th column of W
            zi = zi + W[k][i]*p[k];

         zi = zi + v[i];
         W[j][i] = -B[j]*zi;
      }
   }
   free(v); free(p);
}

void CPU_cmplx_blocked_VB_to_W
 ( int nrows, int ncols, double *B,
   double **Vre, double **Vim, double **Wre, double **Wim )
{
   double *vre = new double[nrows];
   double *vim = new double[nrows];
   double *pre = new double[ncols];
   double *pim = new double[ncols];
   double accre,accim,zi_re,zi_im;

   for(int i=0; i<nrows; i++)
   {
      Wre[0][i] = -B[0]*Vre[0][i];
      Wim[0][i] = -B[0]*Vim[0][i];
   }
   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vre[i] = Vre[j][i];
         vim[i] = Vim[j][i];
      }
      for(int k=0; k<j; k++)
      {
         pre[k] = 0.0;                  // compute k-th component of Y^H*v
         pim[k] = 0.0;
         for(int i=0; i<nrows; i++)     // over all rows of k-th column of Y
         {                              // p[k] = p[k] + V[k][i]*v[i];
            accre =   Vre[k][i]*vre[i] + Vim[k][i]*vim[i];
            accim = - Vim[k][i]*vre[i] + Vre[k][i]*vim[i];
            pre[k] = pre[k] + accre;
            pim[k] = pim[k] + accim;
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zi_re = 0.0;                   // compute i-th component of W*p
         zi_im = 0.0;
         for(int k=0; k<j; k++)         // over all rows of k-th column of W
         {                              // zi = zi + W[k][i]*p[k];
            accre = Wre[k][i]*pre[k] - Wim[k][i]*pim[k];
            accim = Wim[k][i]*pre[k] + Wre[k][i]*pim[k];
            zi_re = zi_re + accre;
            zi_im = zi_im + accim;
         }
         zi_re = zi_re + vre[i];
         zi_im = zi_im + vim[i];
         Wre[j][i] = -B[j]*zi_re;
         Wim[j][i] = -B[j]*zi_im;
      }
   }
   free(vre); free(pre);
   free(vim); free(pim);
}

void CPU_dbl_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **C, double **Y, double **W, bool verbose )
{
   const int rowoff = idx*szt;             // row offset for C
   const int rowdim = nrows - rowoff;      // number of rows in Y and W
   const int coloff = (idx+1)*szt;         // column offset for C
   const int coldim = ncols - coloff;      // number of columns in C

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

   double **YWT = new double*[rowdim];
   for(int i=0; i<rowdim; i++) YWT[i] = new double[rowdim];
   double **prd = new double*[rowdim];
   for(int i=0; i<rowdim; i++) prd[i] = new double[coldim];

   for(int i=0; i<rowdim; i++)        // compute Y*W^T
      for(int j=0; j<rowdim; j++)
      {
         YWT[i][j] = 0.0;             // row i of Y with column j of W^T
         for(int k=0; k<szt; k++)
            YWT[i][j] = YWT[i][j] + Y[k][i]*W[k][j];
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "YWT[" << i << "][" << j << "] : " << YWT[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         prd[i][j] = 0.0;
         for(int k=0; k<rowdim; k++)
            prd[i][j] = prd[i][j] + YWT[i][k]*C[rowoff+k][coloff+j];
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
            cout << "YWTC[" << i << "][" << j << "] : " << prd[i][j] << endl;
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
            cout << "C[" << rowoff+i << "][" << coloff+j << "] : "
                 << C[rowoff+i][coloff+j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)
        C[rowoff+i][coloff+j] = C[rowoff+i][coloff+j] + prd[i][j];
   if(verbose)
   {
      cout << "C after the update with YWTC :" << endl;
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
            cout << "C[" << rowoff+i << "][" << coloff+j << "] : "
                 << C[rowoff+i][coloff+j] << endl;
   }
   for(int i=0; i<rowdim; i++)
   {
      free(YWT[i]);
      free(prd[i]);
   }
   free(YWT); free(prd);
}

void CPU_cmplx_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Cre, double **Cim, double **Yre, double **Yim,
   double **Wre, double **Wim, bool verbose )
{
   const int rowoff = idx*szt;            // row offset for C
   const int rowdim = nrows - rowoff;     // number of rows in Y and W
   const int coloff = (idx+1)*szt;        // column offset for C
   const int coldim = ncols - coloff;     // number of columns in C
   double accre,accim;

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

   double **YWTre = new double*[rowdim];
   double **YWTim = new double*[rowdim];
   for(int i=0; i<rowdim; i++)
   {
      YWTre[i] = new double[rowdim];
      YWTim[i] = new double[rowdim];
   }
   double **prdre = new double*[rowdim];
   double **prdim = new double*[rowdim];
   for(int i=0; i<rowdim; i++)
   {
      prdre[i] = new double[coldim];
      prdim[i] = new double[coldim];
   }
   for(int i=0; i<rowdim; i++)      // compute Y*W^H
      for(int j=0; j<rowdim; j++)
      {
         YWTre[i][j] = 0.0;         // row i of Y with column j of W^H
         YWTim[i][j] = 0.0;
         for(int k=0; k<szt; k++)   // YWT[i][j] = YWT[i][j] + Y[k][i]*W[k][j]
         {  
            accre = Yre[k][i]*Wre[k][j] + Yim[k][i]*Wim[k][j];
            accim = Yim[k][i]*Wre[k][j] - Yre[k][i]*Wim[k][j];
            YWTre[i][j] = YWTre[i][j] + accre;
            YWTim[i][j] = YWTim[i][j] + accim;
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWTre[i][j] << "  "
                 << YWTim[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^H)*C
      for(int j=0; j<coldim; j++)
      {
         prdre[i][j] = 0.0;
         prdim[i][j] = 0.0;
         for(int k=0; k<rowdim; k++)
         {  // prd[i][j] = prd[i][j] + YWT[i][k]*C[rowoff+k][coloff+j];
            accre = YWTre[i][k]*Cre[rowoff+k][coloff+j]
                  - YWTim[i][k]*Cim[rowoff+k][coloff+j];
            accim = YWTim[i][k]*Cre[rowoff+k][coloff+j]
                  + YWTre[i][k]*Cim[rowoff+k][coloff+j];
            prdre[i][j] = prdre[i][j] + accre;
            prdim[i][j] = prdim[i][j] + accim;
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
            cout << "YWHC[" << i << "][" << j << "] : "
                 << prdre[i][j] << "  " << prdim[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         Cre[rowoff+i][coloff+j] = Cre[rowoff+i][coloff+j] + prdre[i][j];
         Cim[rowoff+i][coloff+j] = Cim[rowoff+i][coloff+j] + prdim[i][j];
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWTre[i]); free(prdre[i]);
      free(YWTim[i]); free(prdim[i]);
   }
   free(YWTre); free(prdre);
   free(YWTim); free(prdim);
}

void CPU_dbl_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Q, double **Y, double **W,
   bool verbose )
{
   const int coloff = idx*szt;       // column offset for Q
   const int rowdim = dim - coloff; 
   // the number of rows in Y and W is the number of columns to update

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;
   }

   double **WYT = new double*[rowdim];
   for(int i=0; i<rowdim; i++) WYT[i] = new double[rowdim];
   double **prd = new double*[dim];
   for(int i=0; i<dim; i++) prd[i] = new double[rowdim];

   for(int i=0; i<rowdim; i++)        // compute W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         WYT[i][j] = 0.0;             // row i of W with column j of Y^T
         for(int k=0; k<szt; k++)     // take k-th column of W
            WYT[i][j] = WYT[i][j] + W[k][i]*Y[k][j];
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "Y[" << i << "][" << j << "] : " << Y[j][i] << endl;

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "W[" << i << "][" << j << "] : " << W[j][i] << endl;

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "WYT[" << i << "][" << j << "] : " << WYT[i][j] << endl;
   }
   for(int i=0; i<dim; i++)           // prd = Q*W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         prd[i][j] = 0.0;
         for(int k=0; k<rowdim; k++)
            prd[i][j] = prd[i][j] + Q[i][coloff+k]*WYT[k][j];
      }

   if(verbose)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "QWYT[" << i << "][" << j << "] : " << prd[i][j] << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Q[i][coloff+j] << endl;
   }
   for(int i=0; i<dim; i++)           // Q = Q + Q*W*Y^T
      for(int j=0; j<rowdim; j++)
        Q[i][coloff+j] = Q[i][coloff+j] + prd[i][j];

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Q[i][coloff+j] << endl;
   }
   for(int i=0; i<rowdim; i++) free(WYT[i]);
   for(int i=0; i<dim; i++) free(prd[i]);
   free(WYT); free(prd);
}

void CPU_cmplx_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Qre, double **Qim,
   double **Yre, double **Yim, double **Wre, double **Wim,
   bool verbose )
{
   const int coloff = idx*szt;        // column offset for Q
   const int rowdim = dim - coloff;
   // the number of rows in Y and W is the number of columns to update
   double accre,accim;

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "Y[" << i << "][" << j << "] : "
                 << Yre[j][i] << "  " << Yim[j][i] << endl;

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "W[" << i << "][" << j << "] : "
                 << Wre[j][i] << "  " << Wim[j][i] << endl;
   }

   double **WYTre = new double*[rowdim];
   double **WYTim = new double*[rowdim];
   for(int i=0; i<rowdim; i++)
   {
      WYTre[i] = new double[rowdim];
      WYTim[i] = new double[rowdim];
   }
   double **prdre = new double*[dim];
   double **prdim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      prdre[i] = new double[rowdim];
      prdim[i] = new double[rowdim];
   }
   for(int i=0; i<rowdim; i++)        // compute W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         WYTre[i][j] = 0.0;           // row i of W with column j of Y^H
         WYTim[i][j] = 0.0;
         for(int k=0; k<szt; k++)     // take k-th column of W
         {  // WYT[i][j] = WYT[i][j] + W[k][i]*Y[k][j];
            accre = Wre[k][i]*Yre[k][j] + Wim[k][i]*Yim[k][j];
            accim = Wim[k][i]*Yre[k][j] - Wre[k][i]*Yim[k][j];
            WYTre[i][j] = WYTre[i][j] + accre;
            WYTim[i][j] = WYTim[i][j] + accim;
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYTre[i][j] << "  " << WYTim[i][j] << endl;
   }
   for(int i=0; i<dim; i++)           // prd = Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         prdre[i][j] = 0.0;
         prdim[i][j] = 0.0;
         for(int k=0; k<rowdim; k++)
         {  // prd[i][j] = prd[i][j] + Q[i][coloff+k]*WYT[k][j];
            accre = Qre[i][coloff+k]*WYTre[k][j]
                  - Qim[i][coloff+k]*WYTim[k][j];
            accim = Qim[i][coloff+k]*WYTre[k][j]
                  + Qre[i][coloff+k]*WYTim[k][j];
            prdre[i][j] = prdre[i][j] + accre;
            prdim[i][j] = prdim[i][j] + accim;
         }
      }

   if(verbose)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "QWYH[" << i << "][" << j << "] : "
                 << prdre[i][j] << "  " << prdim[i][j] << endl;
   }
   for(int i=0; i<dim; i++)           // Q = Q + Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         Qre[i][coloff+j] = Qre[i][coloff+j] + prdre[i][j];
         Qim[i][coloff+j] = Qim[i][coloff+j] + prdim[i][j];
      }

   for(int i=0; i<rowdim; i++)
   {
      free(WYTre[i]); free(WYTim[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdre[i]); free(prdim[i]);
   }
   free(WYTre); free(prdre);
   free(WYTim); free(prdim);
}

void CPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *lapsec,
   bool verbose )
{
   double beta;
   double *x = new double[nrows]; // input vector for house
   double *v = new double[nrows]; // Householder vector
   double *B = new double[szt];   // stores the betas
   double **Y = new double*[szt]; // Householder vectors in one block
   double **W = new double*[szt]; // columns of W

   for(int j=0; j<szt; j++)
   {
      Y[j] = new double[nrows];
      W[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++) Q[i][j] = 0.0;
      Q[i][i] = 1.0;
      for(int j=0; j<ncols; j++) R[i][j] = A[i][j];
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

         for(int i=colidx; i<nrows; i++) x[i-colidx] = R[i][colidx];

         CPU_dbl_factors_house(nrowscol,x,v,&beta);
         if(verbose)
         {
            cout << "beta[" << colidx << "] : " << beta << endl;
            for(int i=colidx; i<nrows; i++)
               cout << "v[" << i-colidx << "] : " << v[i-colidx] << endl;
         }
         CPU_dbl_factors_leftRupdate(nrows,endcol,colidx,R,v,beta,verbose);
         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;
         }
         B[L] = beta;
         for(int i=0; i<L; i++) Y[L][i] = 0.0;
         for(int i=0; i<nrowscol; i++) Y[L][i+L] = v[i];
      }
      CPU_dbl_blocked_VB_to_W(nrows-k*szt,szt,B,Y,W);
      if(k<nbt-1)
         CPU_dbl_blocked_leftRupdate(nrows,ncols,szt,k,R,Y,W,verbose);
      CPU_dbl_blocked_rightQupdate(nrows,szt,k,Q,Y,W,verbose);
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(x); free(v); free(B);
   for(int j=0; j<szt; j++)
   {
      free(Y[j]);
      free(W[j]);
   }
   free(Y); free(W);
}

void CPU_cmplx_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim, double *lapsec,
   bool verbose )
{
   double beta;
   double *xre = new double[nrows]; // real input vector for house
   double *xim = new double[nrows]; // imaginary input vector for house
   double *vre = new double[nrows]; // real parts of a Householder vector
   double *vim = new double[nrows]; // imaginary parts of a Householder vector
   double *B = new double[szt];     // stores the betas
   double **Yre = new double*[szt]; // Householder vectors in one block
   double **Yim = new double*[szt]; // imaginary parts of Householder vectors
   double **Wre = new double*[szt]; // real parts of the columns of W
   double **Wim = new double*[szt]; // imaginary parts of the columns of W

   for(int j=0; j<szt; j++)
   {
      Yre[j] = new double[nrows];
      Yim[j] = new double[nrows];
      Wre[j] = new double[nrows];
      Wim[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qre[i][j] = 0.0;
         Qim[i][j] = 0.0;
      }
      Qre[i][i] = 1.0;
      Qim[i][i] = 0.0;
      for(int j=0; j<ncols; j++)
      {
         Rre[i][j] = Are[i][j];
         Rim[i][j] = Aim[i][j];
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
            xre[i-colidx] = Rre[i][colidx];
            xim[i-colidx] = Rim[i][colidx];
         }
         CPU_cmplx_factors_house(nrowscol,xre,xim,vre,vim,&beta);
         if(verbose)
         {
            cout << "beta[" << colidx << "] : " << beta << endl;
            for(int i=colidx; i<nrows; i++)
               cout << "v[" << i-colidx << "] : "
                     << vre[i-colidx] << "  " << vim[i-colidx] << endl;
         }
         CPU_cmplx_factors_leftRupdate
            (nrows,endcol,colidx,Rre,Rim,vre,vim,beta);
         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rre[i][j] << "  " << Rim[i][j] << endl;
         }
         B[L] = beta;
         for(int i=0; i<L; i++)
         {
            Yre[L][i] = 0.0;
            Yim[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yre[L][i+L] = vre[i];
            Yim[L][i+L] = vim[i];
         }
      }
      CPU_cmplx_blocked_VB_to_W(nrows-k*szt,szt,B,Yre,Yim,Wre,Wim);
      if(k<nbt-1)
         CPU_cmplx_blocked_leftRupdate
            (nrows,ncols,szt,k,Rre,Rim,Yre,Yim,Wre,Wim,verbose);
      CPU_cmplx_blocked_rightQupdate
         (nrows,szt,k,Qre,Qim,Yre,Yim,Wre,Wim,verbose);
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(xre); free(vre); free(B);
   free(xim); free(vim);
   for(int j=0; j<szt; j++)
   {
      free(Yre[j]); free(Wre[j]);
      free(Yim[j]); free(Wim[j]);
   }
   free(Yre); free(Wre);
   free(Yim); free(Wim);
}
