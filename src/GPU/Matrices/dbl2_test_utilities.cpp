// The file dbl2_test_utilities.cpp defines the functions specified in
// the file dbl2_test_utilities.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "double_double_functions.h"
#include "random2_matrices.h"
#include "dbl2_factorizations.h"
#include "dbl2_test_utilities.h"

using namespace std;

double dbl2_Difference_Sum
 ( int n, double *xhi, double *xlo, double *yhi, double *ylo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xhi[i] - yhi[i]) + abs(xlo[i] - ylo[i]);

   return result;
}

double cmplx2_Difference_Sum
 ( int n, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
          double *yrehi, double *yrelo, double *yimhi, double *yimlo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xrehi[i] - yrehi[i])
                      + abs(xrelo[i] - yrelo[i])
                      + abs(ximhi[i] - yimhi[i])
                      + abs(ximlo[i] - yimlo[i]);

   return result;
}

double dbl2_Column_Sum ( int dim, int col, double **Ahi, double **Alo )
{
   double resulthi = 0.0;
   double resultlo = 0.0;

   for(int i=0; i<dim; i++)
      ddf_inc(&resulthi,&resultlo,abs(Ahi[i][col]),abs(Alo[i][col]));

   return resulthi;
}

double cmplx2_Column_Sum
 ( int dim, int col,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo )
{
   double resultrehi = 0.0;
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;

   for(int i=0; i<dim; i++)
   {
      ddf_inc(&resultrehi,&resultrelo,abs(Arehi[i][col]),abs(Arelo[i][col]));
      ddf_inc(&resultimhi,&resultimlo,abs(Aimhi[i][col]),abs(Aimlo[i][col]));
   }
   return (resultrehi + resultimhi);
}

double dbl2_Max_Column_Sum ( int dim, double **Ahi, double **Alo )
{
   double result = dbl2_Column_Sum(dim,0,Ahi,Alo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = dbl2_Column_Sum(dim,j,Ahi,Alo);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double cmplx2_Max_Column_Sum
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo )
{
   double result = cmplx2_Column_Sum(dim,0,Arehi,Arelo,Aimhi,Aimlo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = cmplx2_Column_Sum(dim,j,Arehi,Arelo,Aimhi,Aimlo);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double dbl2_condition
 ( int dim, double **Ahi, double **Alo, double **invAhi, double **invAlo )
{
   double Amaxcolsum = dbl2_Max_Column_Sum(dim,Ahi,Alo);
   double invAmaxcolsum = dbl2_Max_Column_Sum(dim,invAhi,invAlo);

   return Amaxcolsum*invAmaxcolsum;
}

double cmplx2_condition
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **invArehi, double **invArelo,
   double **invAimhi, double **invAimlo )
{
   double Amaxcolsum = cmplx2_Max_Column_Sum(dim,Arehi,Arelo,Aimhi,Aimlo);
   double invAmaxcolsum = cmplx2_Max_Column_Sum(dim,invArehi,invArelo,
                                                    invAimhi,invAimlo);

   return Amaxcolsum*invAmaxcolsum;
}

double dbl2_Matrix_Difference_Sum
 ( int n, double **Ahi, double **Alo, double **Bhi, double **Blo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Ahi[i][j] - Bhi[i][j])
                         + abs(Alo[i][j] - Blo[i][j]);

   return result;
}

double cmplx2_Matrix_Difference_Sum
 ( int n, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Brehi, double **Brelo, double **Bimhi, double **Bimlo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Arehi[i][j] - Brehi[i][j])
                         + abs(Arelo[i][j] - Brelo[i][j])
                         + abs(Aimhi[i][j] - Bimhi[i][j])
                         + abs(Aimlo[i][j] - Bimlo[i][j]);

   return result;
}

double dbl2_Diagonal_Difference_Sum
 ( int nbt, int szt, double **Ahi, double **Alo,
   double **Bhi, double **Blo )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Ahi[offset+i][offset+j]
                                - Bhi[offset+i][offset+j])
                            + abs(Alo[offset+i][offset+j]
                                - Blo[offset+i][offset+j]);
   }
   return result;
}

double cmplx2_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Brehi, double **Brelo, double **Bimhi, double **Bimlo )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Arehi[offset+i][offset+j]
                                - Brehi[offset+i][offset+j])
                            + abs(Arelo[offset+i][offset+j]
                                - Brelo[offset+i][offset+j])
                            + abs(Aimhi[offset+i][offset+j]
                                - Bimhi[offset+i][offset+j])
                            + abs(Aimlo[offset+i][offset+j]
                                - Bimlo[offset+i][offset+j]);
   }
   return result;
}

void dbl2_random_upper_factor ( int dim, double **Ahi, double **Alo )
{
   random_dbl2_matrix(dim,dim,Ahi,Alo);

   int *pivots = new int[dim];

   CPU_dbl2_factors_lufac(dim,Ahi,Alo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Ahi[i][j] = 0.0;
         Alo[i][j] = 0.0;
      }

   free(pivots);
}

void cmplx2_random_upper_factor
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo )
{
   random_cmplx2_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);

   int *pivots = new int[dim];

   CPU_cmplx2_factors_lufac(dim,Arehi,Arelo,Aimhi,Aimlo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Arehi[i][j] = 0.0; Arelo[i][j] = 0.0;
         Aimhi[i][j] = 0.0; Aimlo[i][j] = 0.0;
      }

   free(pivots);
}
