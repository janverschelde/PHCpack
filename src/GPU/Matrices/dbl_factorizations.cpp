/* The file dbl_factorizations.cpp defines the functions specified in
 * the file dbl_factorizations.h. */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "dbl_factorizations.h"

using namespace std;

void CPU_dbl_factors_matmatmul
 ( int rows, int dim, int cols, double **A, double **B, double **C )
{
   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         C[i][j] = 0.0;
         for(int k=0; k<dim; k++)
            C[i][j] = C[i][j] + A[i][k]*B[k][j];
      }
}

void CPU_cmplx_factors_matmatmul
 ( int rows, int dim, int cols, double **Are, double **Aim,
   double **Bre, double **Bim, double **Cre, double **Cim )
{
   double zre,zim;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Cre[i][j] = 0.0;
         Cim[i][j] = 0.0;

         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            zre = Are[i][k]*Bre[k][j] - Aim[i][k]*Bim[k][j];
            zim = Aim[i][k]*Bre[k][j] + Are[i][k]*Bim[k][j];
            Cre[i][j] = Cre[i][j] + zre;
            Cim[i][j] = Cim[i][j] + zim;
         }
      }
}

void CPU_dbl_factors_forward ( int dim, double **L, double *b, double *x )
{
   x[0] = b[0];
   for(int i=1; i<dim; i++)
   {
      x[i] = b[i];
      for(int j=0; j<i; j++) x[i] = x[i] - L[i][j]*x[j];
   }
}

void CPU_cmplx_factors_forward
 ( int dim, double **Lre, double **Lim, double *bre, double *bim,
   double *xre, double *xim )
{
   double accre,accim;

   xre[0] = bre[0];
   xim[0] = bim[0];
   for(int i=1; i<dim; i++)
   {
      xre[i] = bre[i];
      xim[i] = bim[i];
      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         accre = Lre[i][j]*xre[j] - Lim[i][j]*xim[j];
         accim = Lim[i][j]*xre[j] + Lre[i][j]*xim[j];
         xre[i] = xre[i] - accre;
         xim[i] = xim[i] - accim;
      }
   }
}

void CPU_dbl_factors_backward ( int dim, double **U, double *b, double *x )
{
   for(int i=dim-1; i>=0; i--)
   {
      x[i] = b[i];
      for(int j=dim-1; j>i; j--) x[i] = x[i] - U[i][j]*x[j];
      x[i] = x[i]/U[i][i];
   }
}

void CPU_cmplx_factors_backward
 ( int dim, double **Ure, double **Uim, double *bre, double *bim,
   double *xre, double *xim )
{
   double accre,accim,det,zre,zim;

   for(int i=dim-1; i>=0; i--)
   {
      xre[i] = bre[i];
      xim[i] = bim[i];
      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         accre = Ure[i][j]*xre[j] - Uim[i][j]*xim[j];
         accim = Uim[i][j]*xre[j] + Ure[i][j]*xim[j];
         xre[i] = xre[i] - accre;
         xim[i] = xim[i] - accim;
      }
      // x[i] = x[i]/U[i][i];
      det = Ure[i][i]*Ure[i][i] + Uim[i][i]*Uim[i][i];
      accre = Ure[i][i]/det;   // accre is the real part of 1/U[i][i]
      accim = -Uim[i][i]/det;  // accim is imaginary part of 1/U[i][i]
      zre = xre[i]*accre - xim[i]*accim;
      zim = xim[i]*accre + xre[i]*accim;
      xre[i] = zre;
      xim[i] = zim;
   }
}

void CPU_dbl_factors_lufac ( int dim, double **A, int *pivots )
{
   double valmax,valtmp;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(A[j][j]); idxmax = j;     // find the pivot
      for(int i=j+1; i< dim; i++)
      {
         valtmp = fabs(A[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = A[idxmax][k];
            A[idxmax][k] = A[j][k];
            A[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         A[i][j] = A[i][j]/A[j][j];
         for(int k=j+1; k<dim; k++) 
            A[i][k] = A[i][k] - A[i][j]*A[j][k];
      }
   }
}

void CPU_cmplx_factors_lufac
 ( int dim, double **Are, double **Aim, int *pivots )
{
   double valmax,valtmp,accre,accim,det,zre,zim;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Are[j][j]) + fabs(Aim[j][j]); // find the pivot
      idxmax = j;
      for(int i=j+1; i< dim; i++)
      {
         valtmp = fabs(Are[i][j]) + fabs(Aim[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Are[idxmax][k];
            Are[idxmax][k] = Are[j][k];
            Are[j][k] = valtmp;
            valtmp = Aim[idxmax][k];
            Aim[idxmax][k] = Aim[j][k];
            Aim[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         det = Are[j][j]*Are[j][j] + Aim[j][j]*Aim[j][j];
         accre = Are[j][j]/det;
         accim = -Aim[j][j]/det;
         zre = Are[i][j]*accre - Aim[i][j]*accim;
         zim = Aim[i][j]*accre + Are[i][j]*accim;
         Are[i][j] = zre;
         Aim[i][j] = zim;
         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            accre = Are[i][j]*Are[j][k] - Aim[i][j]*Aim[j][k];
            accim = Aim[i][j]*Are[j][k] + Are[i][j]*Aim[j][k];
            Are[i][k] = Are[i][k] - accre;
            Aim[i][k] = Aim[i][k] - accim;
         }
      }
   }
}

void CPU_dbl_factors_lusolve
 ( int dim, double **A, int *pivots, double *b, double *x )
{
   CPU_dbl_factors_lufac(dim,A,pivots);
   for(int i=0; i<dim; i++) x[i] = b[pivots[i]];
   CPU_dbl_factors_forward(dim,A,x,b);
   CPU_dbl_factors_backward(dim,A,b,x);
}

void CPU_cmplx_factors_lusolve
 ( int dim, double **Are, double **Aim, int *pivots,
   double *bre, double *bim, double *xre, double *xim )
{
   CPU_cmplx_factors_lufac(dim,Are,Aim,pivots);
   for(int i=0; i<dim; i++)
   {
      xre[i] = bre[pivots[i]];
      xim[i] = bim[pivots[i]];
   }
   CPU_cmplx_factors_forward(dim,Are,Aim,xre,xim,bre,bim);
   CPU_cmplx_factors_backward(dim,Are,Aim,bre,bim,xre,xim);
}

void CPU_dbl_factors_house ( int n, double *x, double *v, double *beta )
{
   double sigma = 0.0;
   double mu,v0p2;
   
   v[0] = 1.0;
   for(int i=1; i<n; i++) 
   {
      sigma = sigma + x[i]*x[i];
      v[i] = x[i];
   }
   if(sigma == 0.0)
      *beta = 0.0;
   else
   {
      mu = sqrt(x[0]*x[0] + sigma);
      if(x[0] <= 0.0)
         v[0] = x[0] - mu;
      else
         v[0] = -sigma/(x[0] + mu);

      v0p2 = v[0]*v[0];
      *beta = 2.0*v0p2/(sigma + v0p2);
      
      for(int i=1; i<n; i++) v[i] = v[i]/v[0];
      v[0] = 1.0;
   }
}

void CPU_cmplx_factors_house 
 ( int n, double *xre, double *xim, double *vre, double *vim, double *beta )
{
   double sigma = 0.0;
   double mu,sqrx0,x0rad,sqrv0,inv0re,inv0im,zre,zim;
   
   vre[0] = 1.0;
   vim[0] = 0.0;

   for(int i=1; i<n; i++) 
   {
      sigma = sigma + xre[i]*xre[i] + xim[i]*xim[i];
      vre[i] = xre[i];
      vim[i] = xim[i];
   }
   if(sigma == 0.0)
      *beta = 0.0;
   else
   {
      sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      x0rad = sqrt(sqrx0);
      mu = sqrt(sqrx0 + sigma); // norm of the vector x

      if(x0rad == 0.0)
      {
         vre[0] = -mu;
         vim[0] = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         mu = mu/x0rad;
         vre[0] = xre[0] - mu*xre[0];
         vim[0] = xim[0] - mu*xim[0];
      }
      sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      *beta = 2.0*sqrv0/(sigma + sqrv0);

      inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]

      for(int i=1; i<n; i++)  // v[i] = v[i]/v[0]
      {
         zre = vre[i]*inv0re - vim[i]*inv0im;
         zim = vim[i]*inv0re + vre[i]*inv0im;
         vre[i] = zre;
         vim[i] = zim;
      }
      vre[0] = 1.0; vim[0] = 0.0;
   }
}

void CPU_dbl_factors_leftRupdate
 ( int nrows, int ncols, int k, double **R, double *v, double beta,
   bool verbose )
{
   double *w = new double[ncols-k];

   for(int j=k; j<ncols; j++)
   {
      w[j-k] = 0.0;
      for(int i=k; i<nrows; i++) w[j-k] = w[j-k] + R[i][j]*v[i-k];
      w[j-k] = beta*w[j-k];
   }
   if(verbose)
   {
      cout << "the vector w = beta*R^T*v : " << endl;
      for(int i=0; i<ncols-k; i++)
         cout << "w[" << i << "] : " << w[i] << endl;
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++)
         R[i][j] = R[i][j] - v[i-k]*w[j-k];

   free(w);
}

void CPU_cmplx_factors_leftRupdate
 ( int nrows, int ncols, int k, double **Rre, double **Rim,
   double *vre, double *vim, double beta )
{
   double *wre = new double[ncols-k];
   double *wim = new double[ncols-k];
   double zre,zim;

   for(int j=k; j<ncols; j++)
   {
      wre[j-k] = 0.0;
      wim[j-k] = 0.0;

      for(int i=k; i<nrows; i++) // w[j-k] = w[j-k] + R[i][j]*v[i-k];
      {
         zre =   Rre[i][j]*vre[i-k] + Rim[i][j]*vim[i-k]; // Hermitian of R
         zim = - Rim[i][j]*vre[i-k] + Rre[i][j]*vim[i-k]; // flip Rim sign
         wre[j-k] = wre[j-k] + zre;
         wim[j-k] = wim[j-k] + zim;
      }
      wre[j-k] = beta*wre[j-k];
      wim[j-k] = beta*wim[j-k];
   }
   for(int i=k; i<nrows; i++)
      for(int j=k; j<ncols; j++) // R[i][j] = R[i][j] - v[i-k]*w[j-k];
      {
         zre = vre[i-k]*wre[j-k] + vim[i-k]*wim[j-k]; // Hermitian of w
         zim = vim[i-k]*wre[j-k] - vre[i-k]*wim[j-k]; // flip wim sign
         Rre[i][j] = Rre[i][j] - zre;
         Rim[i][j] = Rim[i][j] - zim;
      }

   free(wre); free(wim);
}

void CPU_dbl_factors_rightQupdate
 ( int n, int k, double **Q, double *v, double beta )
{
   double *w = new double[n];

   for(int i=0; i<n; i++)
   {
      w[i] = 0.0;
      for(int j=k; j<n; j++) w[i] = w[i] + Q[i][j]*v[j-k];
      w[i] = beta*w[i];
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++)
         Q[i][j] = Q[i][j] - w[i]*v[j-k];

   free(w);
}

void CPU_cmplx_factors_rightQupdate
 ( int n, int k, double **Qre, double **Qim,
   double *vre, double *vim, double beta )
{
   double *wre = new double[n];
   double *wim = new double[n];
   double zre,zim;

   for(int i=0; i<n; i++)
   {
      wre[i] = 0.0;
      wim[i] = 0.0;

      for(int j=k; j<n; j++) // w[i] = w[i] + Q[i][j]*v[j-k];
      {
         zre = Qre[i][j]*vre[j-k] - Qim[i][j]*vim[j-k];
         zim = Qim[i][j]*vre[j-k] + Qre[i][j]*vim[j-k];
         wre[i] = wre[i] + zre;
         wim[i] = wim[i] + zim;
      }
      wre[i] = beta*wre[i];
      wim[i] = beta*wim[i];
   }
   for(int i=0; i<n; i++)
      for(int j=k; j<n; j++) // Q[i][j] = Q[i][j] - w[i]*v[j-k];
      {
         zre = wre[i]*vre[j-k] + wim[i]*vim[j-k]; // Hermitian transpose
         zim = wim[i]*vre[j-k] - wre[i]*vim[j-k]; // flip sign of vim
         Qre[i][j] = Qre[i][j] - zre;
         Qim[i][j] = Qim[i][j] - zim;
      }

   free(wre); free(wim);
}

void CPU_dbl_factors_houseqr
 ( int nrows, int ncols, double **A, double **Q, double **R )
{
   double *x = new double[nrows];
   double *v = new double[nrows];
   double beta;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++) Q[i][j] = 0.0;
      Q[i][i] = 1.0;
      for(int j=0; j<ncols; j++) R[i][j] = A[i][j];
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++) x[i-k] = R[i][k];
         CPU_dbl_factors_house(nrows-k,x,v,&beta);
         CPU_dbl_factors_leftRupdate(nrows,ncols,k,R,v,beta);
         CPU_dbl_factors_rightQupdate(nrows,k,Q,v,beta);
      }
   }
   free(x); free(v);
}

void CPU_cmplx_factors_houseqr
 ( int nrows, int ncols, double **Are, double **Aim,
   double **Qre, double **Qim, double **Rre, double **Rim )
{
   double *xre = new double[nrows];
   double *xim = new double[nrows];
   double *vre = new double[nrows];
   double *vim = new double[nrows];
   double beta;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qre[i][j] = 0.0;
         Qim[i][j] = 0.0;
      }
      Qre[i][i] = 1.0;

      for(int j=0; j<ncols; j++)
      {
         Rre[i][j] = Are[i][j];
         Rim[i][j] = Aim[i][j];
      }
   }
   for(int k=0; k<ncols; k++)
   {
      if(nrows - k > 0)
      {
         for(int i=k; i<nrows; i++)
         {
            xre[i-k] = Rre[i][k];
            xim[i-k] = Rim[i][k];
         }
         CPU_cmplx_factors_house(nrows-k,xre,xim,vre,vim,&beta);
         CPU_cmplx_factors_leftRupdate(nrows,ncols,k,Rre,Rim,vre,vim,beta);
         CPU_cmplx_factors_rightQupdate(nrows,k,Qre,Qim,vre,vim,beta);
      }
   }
   free(xre); free(xim); free(vre); free(vim);
}

void CPU_dbl_factors_qrbs
 ( int nrows, int ncols, double **Q, double **R,
   double *rhs, double *sol, double *wrkvec )
{
   for(int i=0; i<nrows; i++)   // compute Q^T*b, b is rhs
   {
      wrkvec[i] = 0.0;
      for(int j=0; j<nrows; j++)
         wrkvec[i] = wrkvec[i] + Q[j][i]*rhs[j];
   }
   CPU_dbl_factors_backward(ncols,R,wrkvec,sol);
}

void CPU_cmplx_factors_qrbs
 ( int nrows, int ncols,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *rhsre, double *rhsim, double *solre, double *solim,
   double *wrkvecre, double *wrkvecim )
{
   double accre,accim; // accumulates product of two complex numbers

   for(int i=0; i<nrows; i++)    // compute Q^H*b, b is rhs
   {
      wrkvecre[i] = 0.0;
      wrkvecim[i] = 0.0;

      for(int j=0; j<nrows; j++) // work with Hermitian transpose of Q
      {
         accre =  Qre[j][i]*rhsre[j] + Qim[j][i]*rhsim[j];
         accim = -Qim[j][i]*rhsre[j] + Qre[j][i]*rhsim[j];
         wrkvecre[i] = wrkvecre[i] + accre;
         wrkvecim[i] = wrkvecim[i] + accim;
      }
   }
   CPU_cmplx_factors_backward(ncols,Rre,Rim,wrkvecre,wrkvecim,solre,solim);
}
