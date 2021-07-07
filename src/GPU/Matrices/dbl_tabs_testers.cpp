// The file dbl_tabs_testers.cpp defines the function with prototypes
// in dbl_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random_matrices.h"
#include "dbl_factorizations.h"
#include "dbl_tabs_host.h"
#include "dbl_tabs_kernels.h"
#include "dbl_tabs_testers.h"

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

void test_real_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   // random_dbl_upper_matrix(dim,dim,A);
   dbl_random_upper_factor(dim,A);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhs[i] << endl;
   }
   double **invA_h = new double*[dim];
   for(int i=0; i<dim; i++) invA_h[i] = new double[dim];

   CPU_dbl_upper_inverse(dim,A,invA_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invA_h[i][j] << endl;
   }
   double **invA_d = new double*[dim];
   for(int i=0; i<dim; i++) invA_d[i] = new double[dim];

   GPU_dbl_upper_inverse(dim,A,invA_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invA_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : "
        << dbl_Matrix_Difference_Sum(dim,invA_h,invA_d) << endl;

   double *x = new double[dim];
   for(int i=0; i<dim; i++)
   {
      x[i] = 0.0;
      for(int j=0; j<dim; j++) x[i] = x[i] + invA_h[i][j]*rhs[j];
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << x[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : " << dbl_Difference_Sum(dim,sol,x) << endl;
   cout << "Condition number : "
        << dbl_condition(dim,A,invA_h) << endl;
}

void test_cmplx_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Are = new double*[dim];
   double **Aim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Are[i] = new double[dim];
      Aim[i] = new double[dim];
   }
   // random_cmplx_upper_matrix(dim,dim,Are,Aim);
   cmplx_random_upper_factor(dim,Are,Aim);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are[i][j] << "  " << Aim[i][j] << endl;
   }
   double *solre = new double[dim];
   double *solim = new double[dim];
   for(int i=0; i<dim; i++) 
   {
      solre[i] = 1.0;
      solim[i] = 0.0;
   }
   double *rhsre = new double[dim];
   double *rhsim = new double[dim];
   double accre,accim;

   for(int i=0; i<dim; i++)
   {
      rhsre[i] = 0.0;
      rhsim[i] = 0.0;
      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         accre = Are[i][j]*solre[j] - Aim[i][j]*solim[j];
         accim = Aim[i][j]*solre[j] + Are[i][j]*solim[j];
         rhsre[i] = rhsre[i] + accre;
         rhsim[i] = rhsim[i] + accim;
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhsre[i] << "  " << rhsim[i] << endl;
   }
   double **invAre_h = new double*[dim];
   double **invAim_h = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAre_h[i] = new double[dim];
      invAim_h[i] = new double[dim];
   }
   CPU_cmplx_upper_inverse(dim,Are,Aim,invAre_h,invAim_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invAre_h[i][j] << "  "
                 << invAim_h[i][j] << endl;
   }
   double **invAre_d = new double*[dim];
   double **invAim_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAre_d[i] = new double[dim];
      invAim_d[i] = new double[dim];
   }
   GPU_cmplx_upper_inverse(dim,Are,Aim,invAre_d,invAim_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAre_d[i][j] << "  " << invAim_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : "
        << cmplx_Matrix_Difference_Sum
              (dim,invAre_h,invAim_h,invAre_d,invAim_d) << endl;

   double *xre = new double[dim];
   double *xim = new double[dim];
   for(int i=0; i<dim; i++)
   {
      xre[i] = 0.0;
      xim[i] = 0.0;
      for(int j=0; j<dim; j++) // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         accre = invAre_h[i][j]*rhsre[j] - invAim_h[i][j]*rhsim[j];
         accim = invAim_h[i][j]*rhsre[j] + invAre_h[i][j]*rhsim[j];
         xre[i] = xre[i] + accre;
         xim[i] = xim[i] + accim;
      }
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xre[i] << "  " << xim[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : "
        << cmplx_Difference_Sum(dim,solre,solim,xre,xim) << endl;
   cout << "Condition number : "
        << cmplx_condition(dim,Are,Aim,invAre_h,invAim_h) << endl;
}

void test_real_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   // random_dbl_upper_matrix(dim,dim,A);
   dbl_random_upper_factor(dim,A);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhs[i] << endl;
   }
   double *x = new double[dim];
   double *x_d = new double[dim];
   double *rhs_d = new double[dim];
   double **A_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      rhs_d[i] = rhs[i];
      A_d[i] = new double[dim];
      for(int j=0; j<dim; j++) A_d[i][j] = A[i][j];
   }

   CPU_dbl_upper_tiled_solver(dim,sizetile,numtiles,A,rhs,x);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << A[i][j] << endl;
   }

   GPU_dbl_upper_tiled_solver(dim,sizetile,numtiles,A_d,rhs_d,x_d);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << A_d[i][j] << endl;
   }

   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << dbl_Diagonal_Difference_Sum(numtiles,sizetile,A,A_d) << endl;

   if(verbose > 0)
   {
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << x[i] << endl;
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << x_d[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors : "
        << dbl_Difference_Sum(dim,sol,x) << endl;
   cout << "   Sum of GPU errors : "
        << dbl_Difference_Sum(dim,sol,x_d) << endl;
}

void test_cmplx_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Are = new double*[dim];
   double **Aim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Are[i] = new double[dim];
      Aim[i] = new double[dim];
   }
   // random_cmplx_upper_matrix(dim,dim,Are,Aim);
   cmplx_random_upper_factor(dim,Are,Aim);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are[i][j] << "  " << Aim[i][j] << endl;
   }
   double *solre = new double[dim];
   double *solim = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solre[i] = 1.0;
      solim[i] = 0.0;
   }
   double *rhsre = new double[dim];
   double *rhsim = new double[dim];
   double accre,accim;

   for(int i=0; i<dim; i++)
   {
      rhsre[i] = 0.0;
      rhsim[i] = 0.0;
      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         accre = Are[i][j]*solre[j] - Aim[i][j]*solim[j];
         accim = Aim[i][j]*solre[j] + Are[i][j]*solim[j];
         rhsre[i] = rhsre[i] + accre;
         rhsim[i] = rhsim[i] + accim;
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhsre[i] << "  " << rhsim[i] << endl;
   }
   double *xre = new double[dim];
   double *xim = new double[dim];
   double *xre_d = new double[dim];
   double *xim_d = new double[dim];
   double *rhsre_d = new double[dim];
   double *rhsim_d = new double[dim];
   double **Are_d = new double*[dim];
   double **Aim_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      rhsre_d[i] = rhsre[i];
      rhsim_d[i] = rhsim[i];
      Are_d[i] = new double[dim];
      Aim_d[i] = new double[dim];
      for(int j=0; j<dim; j++)
      {
         Are_d[i][j] = Are[i][j];
         Aim_d[i][j] = Aim[i][j];
      }
   }
   CPU_cmplx_upper_tiled_solver
      (dim,sizetile,numtiles,Are,Aim,rhsre,rhsim,xre,xim);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are[i][j] << "  " << Aim[i][j] << endl;
   }
   GPU_cmplx_upper_tiled_solver
      (dim,sizetile,numtiles,Are_d,Aim_d,rhsre_d,rhsim_d,xre_d,xim_d);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are_d[i][j] << "  " << Aim_d[i][j] << endl;
   }

   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << cmplx_Diagonal_Difference_Sum
              (numtiles,sizetile,Are,Aim,Are_d,Aim_d) << endl;

   if(verbose > 0)
   {
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << xre[i] << "  " << xim[i] << endl;
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << xre_d[i] << "  " << xim_d[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors : "
        << cmplx_Difference_Sum(dim,solre,solim,xre,xim) << endl;
   cout << "   Sum of GPU errors : "
        << cmplx_Difference_Sum(dim,solre,solim,xre_d,xim_d) << endl;
}
