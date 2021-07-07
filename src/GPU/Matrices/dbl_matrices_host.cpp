/* The file dbl_matrices_host.cpp defines the function specified
 * in dbl_matrices_host.h. */

#include <cstdlib>
#include <cmath>
#include "dbl_convolutions_host.h"
#include "dbl_matrices_host.h"

// for debugging and verbose mode
#include <iostream>

using namespace std;

void CPU_dbl_inner_product
 ( int dim, int deg, double **x, double **y, double *z )
{
   double *prod = new double[deg+1];

   for(int i=0; i<=deg; i++) z[i] = 0.0;

   for(int k=0; k<dim; k++)
   {
      CPU_dbl_product(deg,x[k],y[k],prod);
      for(int i=0; i<=deg; i++)
         z[i] = z[i] + prod[i];
   }
   free(prod);
}

void CPU_cmplx_inner_product
 ( int dim, int deg,
   double **xre, double **xim, double **yre, double **yim,
   double *zre, double *zim )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      zre[i] = 0.0; zim[i] = 0.0;
   }
   for(int k=0; k<dim; k++)
   {
      CPU_cmplx_product(deg,xre[k],xim[k],yre[k],yim[k],prodre,prodim);
      for(int i=0; i<=deg; i++)
      {
         zre[i] = zre[i] + prodre[i];
         zim[i] = zim[i] + prodim[i];
      }
   }
   free(prodre); free(prodim);
}

void CPU_dbl_matrix_vector_product
 ( int rows, int cols, int deg, double ***A, double **x, double **y )
{
   for(int k=0; k<rows; k++)
      CPU_dbl_inner_product(cols,deg,A[k],x,y[k]);
}

void CPU_cmplx_matrix_vector_product
 ( int rows, int cols, int deg, double ***Are, double ***Aim,
   double **xre, double **xim, double **yre, double **yim )
{
   for(int k=0; k<rows; k++)
      CPU_cmplx_inner_product(cols,deg,Are[k],Aim[k],xre,xim,yre[k],yim[k]);
}

void CPU_dbl_matrix_matrix_product
 ( int rows, int dim, int cols, int deg, 
   double ***A, double ***B, double ***C )
{
   double *prod = new double[deg+1];

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
      {
         for(int d=0; d<=deg; d++) C[i][j][d] = 0.0;

         for(int k=0; k<dim; k++)
         {
            CPU_dbl_product(deg,A[i][k],B[k][j],prod);
            for(int d=0; d<=deg; d++)
               C[i][j][d] = C[i][j][d] + prod[d];
         }
      }
   }
   free(prod);
}

void CPU_cmplx_matrix_matrix_product
 ( int rows, int dim, int cols, int deg, 
   double ***Are, double ***Aim, double ***Bre, double ***Bim,
   double ***Cre, double ***Cim )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
      {
         for(int d=0; d<=deg; d++)
         {
            Cre[i][j][d] = 0.0; Cim[i][j][d] = 0.0;
         }
         for(int k=0; k<dim; k++)
         {
            CPU_cmplx_product
               (deg,Are[i][k],Aim[i][k],Bre[k][j],Bim[k][j],prodre,prodim);
            for(int d=0; d<=deg; d++)
            {
               Cre[i][j][d] = Cre[i][j][d] + prodre[d];
               Cim[i][j][d] = Cim[i][j][d] + prodim[d];
            }
         }
      }
   }
   free(prodre); free(prodim);
}

void CPU_dbl_lower_solver
 ( int dim, int deg, double ***L, double **b, double **x, bool verbose )
{
   double *prod = new double[deg+1]; // for products
   // double *accu = new double[deg+1]; // accumulates products

   // CPU_dbl_inverse(deg,L[0][0],prod);      // prod = 1/L[0][0]
   // CPU_dbl_product(deg,b[0],prod,x[0]);    // x[0] = b[0]/L[0][0]

   if(verbose) cout << "Level of operations in forward substitution :"
                    << endl;

   if(verbose) cout << "x[0] = b[0], level 0" << endl;

   for(int k=0; k<=deg; k++) x[0][k] = b[0][k];

   for(int i=1; i<dim; i++)
   {
      if(verbose) cout << "x[" << i << "] = b[" << i << "], level 0" << endl;

      for(int k=0; k<=deg; k++) x[i][k] = b[i][k];

      for(int j=0; j<i; j++)
      {
         if(verbose) cout << "prod = L[" << i << "][" << j << "]*x[" << j
                          << "], level " << 2*(j+1)-1 << endl;

         CPU_dbl_product(deg,L[i][j],x[j],prod);

         if(verbose) cout << "x[" << i << "] = x[" << i
                          << "] - prod, level " << 2*(j+1) << endl;

         for(int k=0; k<=deg; k++) x[i][k] = x[i][k] - prod[k];
      }
      // CPU_dbl_inverse(deg,L[i][i],prod);    // prod = 1/L[i][i]
      // CPU_dbl_product(deg,accu,prod,x[i]);  // x[i] = acc/L[i][i]
   }
   free(prod); // free(accu);
}

void CPU_cmplx_lower_solver
 ( int dim, int deg, double ***Lre, double ***Lim,
   double **bre, double **bim, double **xre, double **xim, bool verbose )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];

   if(verbose) cout << "Level of operations in forward substitution :"
                    << endl;

   if(verbose) cout << "x[0] = b[0], level 0" << endl;

   for(int k=0; k<=deg; k++)
   {
      xre[0][k] = bre[0][k];
      xim[0][k] = bim[0][k];
   }
   for(int i=1; i<dim; i++)
   {
      if(verbose) cout << "x[" << i << "] = b[" << i << "], level 0" << endl;

      for(int k=0; k<=deg; k++)
      {
         xre[i][k] = bre[i][k];
         xim[i][k] = bim[i][k];
      }
      for(int j=0; j<i; j++)
      {
         if(verbose) cout << "prod = L[" << i << "][" << j << "]*x[" << j
                          << "], level " << 2*(j+1)-1 << endl;

         CPU_cmplx_product
            (deg,Lre[i][j],Lim[i][j],xre[j],xim[j],prodre,prodim);

         if(verbose) cout << "x[" << i << "] = x[" << i
                          << "] - prod, level " << 2*(j+1) << endl;

         for(int k=0; k<=deg; k++)
         {
            xre[i][k] = xre[i][k] - prodre[k];
            xim[i][k] = xim[i][k] - prodim[k];
         }
      }
   }
   free(prodre); free(prodim);
}

void CPU_dbl_upper_solver
 ( int dim, int deg, double ***U, double **b, double **x, bool verbose )
{
   double *prod = new double[deg+1];
   double *work = new double[deg+1];
   int lvl;

   if(verbose) cout << "Level of operations in backward substitution :"
                    << endl;

   for(int i=dim-1; i>=0; i--)
   {
      if(verbose) cout << "prod = b[" << i << "], level 0" << endl;

      for(int k=0; k<=deg; k++) prod[k] = b[i][k];

      lvl = 1;

      for(int j=dim-1; j>i; j--)
      {
         if(verbose) cout << "work = U[" << i << "][" << j << "]*x[" << j
                          << "], level " << ++lvl << endl;

         CPU_dbl_product(deg,U[i][j],x[j],work);

         if(verbose) cout << "prod = prod - work, level = "
                          << ++lvl << endl;

         for(int k=0; k<=deg; k++) prod[k] = prod[k] - work[k];

         lvl = lvl + 1;
      }
      if(verbose) cout << "work = 1/U[" << i << "][" << i
                       << "], level 0" << endl;

      CPU_dbl_inverse(deg,U[i][i],work);

      if(verbose) cout << "x[" << i << "] = prod/U["
                       << i << "][" << i << "], level "
                       << lvl << endl;

      CPU_dbl_product(deg,work,prod,x[i]);
   }
   free(prod); free(work);
}

void CPU_cmplx_upper_solver
 ( int dim, int deg, double ***Ure, double ***Uim,
   double **bre, double **bim, double **xre, double **xim, bool verbose )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];
   double *workre = new double[deg+1];
   double *workim = new double[deg+1];
   int lvl;

   if(verbose) cout << "Level of operations in backward substitution :"
                    << endl;

   for(int i=dim-1; i>=0; i--)
   {
      if(verbose) cout << "prod = b[" << i << "], level 0" << endl;

      for(int k=0; k<=deg; k++)
      {
         prodre[k] = bre[i][k];
         prodim[k] = bim[i][k];
      }
      lvl = 1;

      for(int j=dim-1; j>i; j--)
      {
         if(verbose) cout << "work = U[" << i << "][" << j << "]*x[" << j
                          << "], level " << ++lvl << endl;

         CPU_cmplx_product
            (deg,Ure[i][j],Uim[i][j],xre[j],xim[j],workre,workim);

         if(verbose) cout << "prod = prod - work, level = "
                          << ++lvl << endl;

         for(int k=0; k<=deg; k++)
         {
            prodre[k] = prodre[k] - workre[k];
            prodim[k] = prodim[k] - workim[k];
         }
         lvl = lvl + 1;
      }
      if(verbose) cout << "work = 1/U[" << i << "][" << i
                       << "], level 0" << endl;

      CPU_cmplx_inverse(deg,Ure[i][i],Uim[i][i],workre,workim);

      if(verbose) cout << "x[" << i << "] = prod/U["
                       << i << "][" << i << "], level "
                       << lvl << endl;

      CPU_cmplx_product(deg,workre,workim,prodre,prodim,xre[i],xim[i]);
   }
   free(prodre); free(workre);
   free(prodim); free(workim);
}

void CPU_dbl_upper_lead_solver
 ( int dim, double **U, double *b, double *x )
{
   int idx = dim-1;
   x[idx] = b[idx]/U[idx][idx];

   for(int i=idx-1; i>=0; i--)
   {
      x[i] = b[i];
      for(int j=i+1; j<dim; j++) x[i] = x[i] - U[i][j]*x[j];
      x[i] = x[i]/U[i][i];
   }
}

void CPU_cmplx_upper_lead_solver
 ( int dim, double **Ure, double **Uim, double *bre, double *bim,
   double *xre, double *xim )
{
   double zre,zim,det,vre,vim;
   int idx = dim-1;

   det = Ure[idx][idx]*Ure[idx][idx] + Uim[idx][idx]*Uim[idx][idx];
   vre = Ure[idx][idx]/det;   // vim is the real part of 1/U[idx][idx]
   vim = -Uim[idx][idx]/det;  // vim is imaginary part of 1/U[idx][idx]
   xre[idx] = bre[idx]*vre - bim[idx]*vim;
   xim[idx] = bim[idx]*vre + bre[idx]*vim; // x[idx] = b[idx]/U[idx][idx];

   for(int i=idx-1; i>=0; i--)
   {
      xre[i] = bre[i]; xim[i] = bim[i];
      for(int j=i+1; j<dim; j++) // x[i] = x[i] - U[i][j]*x[j];
      {
         zre = Ure[i][j]*xre[j] - Uim[i][j]*xim[j];
         zim = Uim[i][j]*xre[j] + Ure[i][j]*xim[j];
         xre[i] = xre[i] - zre;
         xim[i] = xim[i] - zim;
      }
      det = Ure[i][i]*Ure[i][i] + Uim[i][i]*Uim[i][i];
      vre = Ure[i][i]/det;   // vim is the real part of 1/U[i][i]
      vim = -Uim[i][i]/det;  // vim is imaginary part of 1/U[i][i]
      zre = xre[i]*vre - xim[i]*vim;
      zim = xim[i]*vre + xre[i]*vim;
      xre[i] = zre; xim[i] = zim;     // x[i] = x[i]/U[i][i];
   }
}

void CPU_dbl_upper_linearized_solver
 ( int dim, int deg, double ***U, double **b, double **x )
{
   CPU_dbl_upper_lead_solver(dim,U[0],b[0],x[0]);

   double *work = new double[dim];
   double prod;

   for(int i=1; i<deg; i++)
   {
      for(int j=0; j<dim; j++) work[j] = b[i][j];
      for(int j=1; j<=i; j++)
      {                            // multiply U[j] with x[i-j]
         for(int k=0; k<dim; k++)  // multiply k-th row of U[j] with x[i-j]
         {
            prod = 0.0;                   // accumulate over dim columns
            for(int L=0; L<dim; L++)
               prod = prod + U[j][k][L]*x[i-j][L];
            work[k] = work[k] - prod;
         }
      }
      CPU_dbl_upper_lead_solver(dim,U[0],work,x[i]);
   }
   free(work);
}

void CPU_cmplx_upper_linearized_solver
 ( int dim, int deg, double ***Ure, double ***Uim,
   double **bre, double **bim, double **xre, double **xim )
{
   CPU_cmplx_upper_lead_solver
      (dim,Ure[0],Uim[0],bre[0],bim[0],xre[0],xim[0]);

   double *workre = new double[dim];
   double *workim = new double[dim];
   double prodre,prodim,zre,zim;

   for(int i=1; i<deg; i++)
   {
      for(int j=0; j<dim; j++)
      {
         workre[j] = bre[i][j];
         workim[j] = bim[i][j];
      }
      for(int j=1; j<=i; j++)
      {                            // multiply U[j] with x[i-j]
         for(int k=0; k<dim; k++)  // multiply k-th row of U[j] with x[i-j]
         {
            prodre = 0.0;                   // accumulate over dim columns
            prodim = 0.0;
            for(int L=0; L<dim; L++)  // prod = prod + U[j][k][L]*x[i-j][L]
            {
               zre = Ure[j][k][L]*xre[i-j][L] - Uim[j][k][L]*xim[i-j][L];
               zim = Ure[j][k][L]*xim[i-j][L] + Uim[j][k][L]*xre[i-j][L];
               prodre = prodre + zre;
               prodim = prodim + zim;
            }
            workre[k] = workre[k] - prodre;
            workim[k] = workim[k] - prodim;
         }
      }
      CPU_cmplx_upper_lead_solver
         (dim,Ure[0],Uim[0],workre,workim,xre[i],xim[i]);
   }
   free(workre); free(workim);
}

void CPU_dbl_lufac
 ( int dim, int deg, double ***A, int *pivots, bool verbose )
{
   double valmax,valtmp;
   int idxmax,idxtmp;
   double *work = new double[deg+1];
   double *prod = new double[deg+1];

   for(int j=0; j<dim; j++) pivots[j] = j;
   
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(A[j][j][0]); idxmax = j;
      for(int i=j+1; i<dim; i++)
      {                                // find the pivot
         valtmp = fabs(A[i][j][0]);
         if(valtmp > valmax) 
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                  // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            for(int i=0; i<=deg; i++)
            {
               valtmp = A[idxmax][k][i];
               A[idxmax][k][i] = A[j][k][i];
               A[j][k][i] = valtmp;
            }
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)
      {                                        // A[i][j] = A[i][j]/A[j][j]
         CPU_dbl_inverse(deg,A[j][j],work);  
         CPU_dbl_product(deg,A[i][j],work,prod);
         for(int k=0; k<=deg; k++) A[i][j][k] = prod[k];

         for(int k=j+1; k<dim; k++)  // A[i][k] = A[i][k] - A[i][j]*A[j][k]
         {
            CPU_dbl_product(deg,A[i][j],A[j][k],prod);
            for(int d=0; d<=deg; d++) A[i][k][d] = A[i][k][d] - prod[d];
         }
      }
   }
   free(work); free(prod);
}

void CPU_cmplx_lufac
 ( int dim, int deg, double ***Are, double ***Aim, int *pivots, bool verbose )
{
   double valmax,valtmp;
   int idxmax,idxtmp;
   double *workre = new double[deg+1];
   double *workim = new double[deg+1];
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];

   for(int j=0; j<dim; j++) pivots[j] = j;
   
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Are[j][j][0]) + fabs(Aim[j][j][0]);
      idxmax = j;
      for(int i=j+1; i<dim; i++)
      {                                // find the pivot
         valtmp = fabs(Are[i][j][0]) + fabs(Aim[i][j][0]);
         if(valtmp > valmax) 
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                  // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            for(int i=0; i<=deg; i++)
            {
               valtmp = Are[idxmax][k][i];
               Are[idxmax][k][i] = Are[j][k][i];
               Are[j][k][i] = valtmp;
               valtmp = Aim[idxmax][k][i];
               Aim[idxmax][k][i] = Aim[j][k][i];
               Aim[j][k][i] = valtmp;
            }
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)
      {                                        // A[i][j] = A[i][j]/A[j][j]
         CPU_cmplx_inverse
            (deg,Are[j][j],Aim[j][j],workre,workim);  
         CPU_cmplx_product
            (deg,Are[i][j],Aim[i][j],workre,workim,prodre,prodim);

         for(int k=0; k<=deg; k++)
         {
            Are[i][j][k] = prodre[k];
            Aim[i][j][k] = prodim[k];
         }
         for(int k=j+1; k<dim; k++)  // A[i][k] = A[i][k] - A[i][j]*A[j][k]
         {
            CPU_cmplx_product
               (deg,Are[i][j],Aim[i][j],Are[j][k],Aim[j][k],prodre,prodim);

            for(int d=0; d<=deg; d++)
            {
               Are[i][k][d] = Are[i][k][d] - prodre[d];
               Aim[i][k][d] = Aim[i][k][d] - prodim[d];
            }
         }
      }
   }
   free(workre); free(workim);
   free(prodre); free(prodim);
}

void CPU_dbl_lu_solver
 ( int dim, int deg, double ***A, int *pivots, double **b, double **x,
   bool verbose )
{
   CPU_dbl_lufac(dim,deg,A,pivots,verbose);
  
   for(int i=0; i<dim; i++)        // permute b according to the pivots
   {
      for(int d=0; d<=deg; d++) x[i][d] = b[pivots[i]][d];
   }
   CPU_dbl_lower_solver(dim,deg,A,x,b,verbose); // forward substitution
   CPU_dbl_upper_solver(dim,deg,A,b,x,verbose); // backward substitution
}

void CPU_cmplx_lu_solver
 ( int dim, int deg, double ***Are, double ***Aim, int *pivots,
   double **bre, double **bim, double **xre, double **xim, bool verbose )
{
   CPU_cmplx_lufac(dim,deg,Are,Aim,pivots,verbose);
  
   for(int i=0; i<dim; i++)        // permute b according to the pivots
   {
      for(int d=0; d<=deg; d++)
      {
         xre[i][d] = bre[pivots[i]][d];
         xim[i][d] = bim[pivots[i]][d];
      }
   }
   CPU_cmplx_lower_solver(dim,deg,Are,Aim,xre,xim,bre,bim,verbose);
   CPU_cmplx_upper_solver(dim,deg,Are,Aim,bre,bim,xre,xim,verbose);
}
