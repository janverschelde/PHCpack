// The file dbl_tabs_testers.cpp defines the function with prototypes in
// the file dbl_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random_matrices.h"
#include "dbl_factorizations.h"
#include "dbl_tabs_host.h"
#include "dbl_tabs_kernels.h"
#include "dbl_test_utilities.h"
#include "write_dbl_bstimeflops.h"
#include "dbl_tabs_testers.h"
#include "dbl_data_files.h"

using namespace std;

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
/*
 * The code below makes an identity matrix for A.
 *
  for(int i=0; i<dim; i++)
  {
     for(int j=0; j<dim; j++) A[i][j] = 0.0;
     A[i][i] = 1.0;
  }
*/
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

   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_dbl_upper_inverse(dim,A,invA_h,&timelapsed_h);

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

   cout << "-> GPU computes the inverse ..." << endl;

   GPU_dbl_upper_inverse(dim,A,invA_d,&elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invA_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
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
   cout << "   Sum of errors on solution : "
        << dbl_Difference_Sum(dim,sol,x) << endl;
   cout << "Condition number : "
        << dbl_condition(dim,A,invA_h) << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
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
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;
 
   CPU_cmplx_upper_inverse(dim,Are,Aim,invAre_h,invAim_h,&timelapsed_h);

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
   cout << "-> GPU computes the inverse ..." << endl;

   GPU_cmplx_upper_inverse
      (dim,Are,Aim,invAre_d,invAim_d,&elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAre_d[i][j] << "  " << invAim_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
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
   cout << "   Sum of errors on solution : "
        << cmplx_Difference_Sum(dim,solre,solim,xre,xim) << endl;
   cout << "Condition number : "
        << cmplx_condition(dim,Are,Aim,invAre_h,invAim_h) << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
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

   cout << "Generate a random matrix (1 = yes, 0 = read matrix) : ";
   int rndmat; cin >> rndmat;

   if(rndmat == 1)
      cout << "-> generating a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
   else
      cout << "-> reading a random upper triangular matrix of dimension "
           << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   // random_dbl_upper_matrix(dim,dim,A);
   if(rndmat == 1)
      dbl_random_upper_factor(dim,A);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      dbl_read_matrix(filename,dim,A);
   }
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
   double timelapsed_h,timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_dbl_upper_tiled_solver(dim,sizetile,numtiles,A,rhs,x,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << A[i][j] << endl;
   }
   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_dbl_upper_tiled_solver
      (dim,sizetile,numtiles,A_d,rhs_d,x_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

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
   cout << "   Sum of CPU errors on solution : "
        << dbl_Difference_Sum(dim,sol,x) << endl;
   cout << "   Sum of GPU errors on solution : "
        << dbl_Difference_Sum(dim,sol,x_d) << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl_bstimeflops
      (sizetile,numtiles,0,invlapsed,mullapsed,sublapsed,elapsedms,
       timelapsed_d,addcnt,mulcnt,divcnt);

   for(int i=0; i<dim; i++)
   {
      free(A[i]); free(A_d[i]);
   }
   free(A); free(A_d); free(sol);
   free(rhs); free(rhs_d); free(x); free(x_d);
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

   cout << "Generate a random matrix (1 = yes, 0 = read matrix) : ";
   int rndmat; cin >> rndmat;

   if(rndmat == 1)
      cout << "-> generating a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
   else
      cout << "-> reading a random upper triangular matrix of dimension "
           << dim << " ..." << endl;

   double **Are = new double*[dim];
   double **Aim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Are[i] = new double[dim];
      Aim[i] = new double[dim];
   }
   // random_cmplx_upper_matrix(dim,dim,Are,Aim);
   if(rndmat == 1)
      cmplx_random_upper_factor(dim,Are,Aim);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      dbl_read_matrix(filename,dim,Are);
   }
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
   double timelapsed_h,timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_cmplx_upper_tiled_solver
      (dim,sizetile,numtiles,Are,Aim,rhsre,rhsim,xre,xim,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are[i][j] << "  " << Aim[i][j] << endl;
   }
   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_cmplx_upper_tiled_solver
      (dim,sizetile,numtiles,Are_d,Aim_d,rhsre_d,rhsim_d,xre_d,xim_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

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
   cout << "   Sum of CPU errors on solution : "
        << cmplx_Difference_Sum(dim,solre,solim,xre,xim) << endl;
   cout << "   Sum of GPU errors on solution : "
        << cmplx_Difference_Sum(dim,solre,solim,xre_d,xim_d) << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl_bstimeflops
      (sizetile,numtiles,1,invlapsed,mullapsed,sublapsed,elapsedms,
       timelapsed_d,addcnt,mulcnt,divcnt);

   for(int i=0; i<dim; i++)
   {
      free(Are[i]); free(Are_d[i]);
      free(Aim[i]); free(Aim_d[i]);
   }
   free(Are); free(Aim); free(Are_d); free(Aim_d);
   free(solre); free(solim);
   free(rhsre); free(rhsim); free(rhsre_d); free(rhsim_d);
   free(xre); free(xim); free(xre_d); free(xim_d);
}
