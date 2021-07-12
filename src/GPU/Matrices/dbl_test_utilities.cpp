// The file dbl_test_utilities.cpp defines the function with prototypes in
// the file dbl_test_utilities.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random_matrices.h"
#include "dbl_factorizations.h"
#include "dbl_test_utilities.h"

using namespace std;

double dbl_Difference_Sum ( int n, double *x, double *y )
{
   double result = 0.0;

   for(int i=0; i<n; i++) result = result + abs(x[i] - y[i]);

   return result;
}

double cmplx_Difference_Sum
 ( int n, double *xre, double *xim, double *yre, double *yim )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xre[i] - yre[i]) + abs(xim[i] - yim[i]);

   return result;
}

double dbl_Column_Sum ( int dim, int col, double **A )
{
   double result = 0.0;

   for(int i=0; i<dim; i++) result = result + abs(A[i][col]);

   return result;
}

double cmplx_Column_Sum ( int dim, int col, double **Are, double **Aim )
{
   double result = 0.0;

   for(int i=0; i<dim; i++)
      result = result + abs(Are[i][col]) + abs(Aim[i][col]);

   return result;
}

double dbl_Max_Column_Sum ( int dim, double **A )
{
   double result = dbl_Column_Sum(dim,0,A);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = dbl_Column_Sum(dim,j,A);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double cmplx_Max_Column_Sum ( int dim, double **Are, double **Aim )
{
   double result = cmplx_Column_Sum(dim,0,Are,Aim);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = cmplx_Column_Sum(dim,j,Are,Aim);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double dbl_condition ( int dim, double **A, double **invA )
{
   double Amaxcolsum = dbl_Max_Column_Sum(dim,A);
   double invAmaxcolsum = dbl_Max_Column_Sum(dim,invA);

   return Amaxcolsum*invAmaxcolsum;
}

double cmplx_condition
 ( int dim, double **Are, double **Aim, double **invAre, double **invAim )
{
   double Amaxcolsum = cmplx_Max_Column_Sum(dim,Are,Aim);
   double invAmaxcolsum = cmplx_Max_Column_Sum(dim,invAre,invAim);

   return Amaxcolsum*invAmaxcolsum;
}

double dbl_Matrix_Difference_Sum ( int n, double **A, double **B )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(A[i][j] - B[i][j]);

   return result;
}

double cmplx_Matrix_Difference_Sum
 ( int n, double **Are, double **Aim, double **Bre, double **Bim )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Are[i][j] - Bre[i][j])
                         + abs(Aim[i][j] - Bim[i][j]);

   return result;
}

double dbl_Diagonal_Difference_Sum
 ( int nbt, int szt, double **A, double **B )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(A[offset+i][offset+j]
                                - B[offset+i][offset+j]);
   }
   return result;
}

double cmplx_Diagonal_Difference_Sum
 ( int nbt, int szt, double **Are, double **Aim,
   double **Bre, double **Bim )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Are[offset+i][offset+j]
                                - Bre[offset+i][offset+j])
                            + abs(Aim[offset+i][offset+j]
                                - Bim[offset+i][offset+j]);
   }
   return result;
}

void dbl_random_upper_factor ( int dim, double **A )
{
   random_dbl_matrix(dim,dim,A);

   int *pivots = new int[dim];

   CPU_dbl_factors_lufac(dim,A,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++) A[i][j] = 0.0;

   free(pivots);
}

void cmplx_random_upper_factor ( int dim, double **Are, double **Aim )
{
   random_cmplx_matrix(dim,dim,Are,Aim);

   int *pivots = new int[dim];

   CPU_cmplx_factors_lufac(dim,Are,Aim,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Are[i][j] = 0.0;
         Aim[i][j] = 0.0;
      }

   free(pivots);
}
