/* The file dbl_baqr_host.cpp defines the functions specified in
 * the file dbl_baqr_host.h. */

#include <iostream>
#include <cstdlib>
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

void CPU_dbl_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **C, double **Y, double **W, bool verbose )
{
   const int rowdim = nrows - idx*szt;      // number of rows in Y and W
   const int coldim = ncols - (idx+1)*szt;  // number of columns in C
   const int rowoff = idx*szt;              // row offset for C
   const int coloff = (idx+1)*szt;          // column offset for C

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
      cout << endl;
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
            cout << "prd[" << i << "][" << j << "] : " << prd[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)
        C[rowoff+i][coloff+j] = C[rowoff+i][coloff+j] + prd[i][j];

   for(int i=0; i<rowdim; i++)
   {
      free(YWT[i]);
      free(prd[i]);
   }
   free(YWT); free(prd);
}

void CPU_dbl_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Q, double **Y, double **W,
   bool verbose )
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
      cout << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "prd[" << i << "][" << j << "] : " << prd[i][j] << endl;
   }
   for(int i=0; i<dim; i++)           // Q = Q + Q*W*Y^T
      for(int j=0; j<rowdim; j++)
        Q[i][coloff+j] = Q[i][coloff+j] + prd[i][j];

   for(int i=0; i<rowdim; i++) free(WYT[i]);
   for(int i=0; i<dim; i++) free(prd[i]);
   free(WYT); free(prd);
}

void CPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, bool verbose )
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
   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++) Q[i][j] = 0.0;
      Q[i][i] = 1.0;
      for(int j=0; j<ncols; j++) R[i][j] = A[i][j];
   }
   int colidx,endcol,nrowscol;

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      for(int L=0; L<szt; L++)    // L runs over the columns in one block
      {
         colidx = k*szt + L;      // index of the current column
         endcol = (k+1)*szt;      // 1 + last column index in current block
         nrowscol = nrows-colidx; // number of rows in current column

         for(int i=colidx; i<nrows; i++) x[i-colidx] = R[i][colidx];

         CPU_dbl_factors_house(nrowscol,x,v,&beta);
         CPU_dbl_factors_leftRupdate(nrows,endcol,colidx,R,v,beta);
         B[L] = beta;
         for(int i=0; i<L; i++) Y[L][i] = 0.0;
         for(int i=0; i<nrowscol; i++) Y[L][i+L] = v[i];
      }
      CPU_dbl_blocked_VB_to_W(nrows-k*szt,szt,B,Y,W);
      if(k<nbt-1)
         CPU_dbl_blocked_leftRupdate(nrows,ncols,szt,k,R,Y,W,verbose);
      CPU_dbl_blocked_rightQupdate(nrows,szt,k,Q,Y,W,verbose);
   }
   free(x); free(v); free(B);
   for(int j=0; j<szt; j++)
   {
      free(Y[j]);
      free(W[j]);
   }
   free(Y); free(W);
}
