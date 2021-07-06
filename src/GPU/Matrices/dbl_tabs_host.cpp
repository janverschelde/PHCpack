/* The file dbl_tabs_host.cpp defines the functions specified in
 * the file dbl_tabs_host.h. */

#include <cstdlib>
#include "dbl_matrices_host.h"
#include "dbl_tabs_host.h"

void CPU_dbl_backsubs ( int dim, double **U, double *b, double *x )
{
   CPU_dbl_upper_lead_solver(dim,U,b,x);
}

void CPU_dbl_upper_inverse ( int dim, double **U, double **invU )
{
   double *col = new double[dim];
   double *rhs = new double[dim];

   for(int i=0; i<dim; i++) rhs[i] = 0.0;

   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhs[j] = 1.0;
      CPU_dbl_backsubs(dim,U,rhs,col);
      for(int i=0; i<dim; i++) invU[i][j] = col[i];
      rhs[j] = 0.0;
   }
   free(rhs); free(col);
}

void CPU_dbl_matmatmul ( int dim, double **A, double **F )
{
   double **result = new double*[dim];
   for(int i=0; i<dim; i++) result[i] = new double[dim];

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         result[i][j] = 0.0;
         for(int k=0; k<dim; k++)
            result[i][j] = result[i][j] + F[i][k]*A[k][j];
      }

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) A[i][j] = result[i][j];

   for(int i=0; i<dim; i++) free(result[i]);
   free(result);
}

void CPU_dbl_upper_tiled_solver
 ( int dim, int szt, int nbt, double **U, double *b, double *x )
{
   double **T = new double*[szt];
   double **invT = new double*[szt];
   for(int i=0; i<szt; i++)
   {
      T[i] = new double[szt];
      invT[i] = new double[szt];
   }
   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++) T[i][j] = U[idx+i][idx+j];

   CPU_dbl_upper_inverse(szt,T,invT);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++) U[idx+i][idx+j] = invT[i][j];

   for(int i=0; i<szt; i++)
   {
      x[idx+i] = 0.0;
      for(int j=0; j<szt; j++) x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
   }
   double *wb = new double[szt];    // work space for b
   double **wT = new double*[szt];  // work space for a tile
   for(int i=0; i<szt; i++)
      wT[i] = new double[szt];

   double prod;

   for(int k=1; k<nbt; k++)
   {
      idx = idx - szt;

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++) T[i][j] = U[idx+i][idx+j];

      CPU_dbl_upper_inverse(szt,T,invT);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++) U[idx+i][idx+j] = invT[i][j];

      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prod = 0.0;
         for(int j=0; j<szt; j++) prod = prod + invT[i][j]*b[idx+j];
         wb[i] = prod;
      }
      for(int L=1; L<=k; L++)   // update wb as many times as k
      {
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++) wT[i][j] = U[idx+i][idx+L*szt+j];

         CPU_dbl_matmatmul(szt,wT,invT); // multiply tile with invT

         for(int i=0; i<szt; i++) // update wb
         {
            prod = 0.0;
            for(int j=0; j<szt; j++) prod = prod + wT[i][j]*x[idx+L*szt+j];
            wb[i] = wb[i] - prod;
         }
      }
      for(int i=0; i<szt; i++) x[idx+i] = wb[i];
   }

   for(int i=0; i<szt; i++)
   {
      free(T[i]); free(invT[i]); free(wT[i]);
   }
   free(T); free(invT); free(wb); free(wT);
}
