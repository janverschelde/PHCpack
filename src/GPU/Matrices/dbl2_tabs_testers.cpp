// The file dbl2_tabs_testers.cpp defines the functions specified in
// the file dbl2_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "double_double_functions.h"
#include "random2_matrices.h"
#include "dbl2_factorizations.h"
#include "dbl2_tabs_host.h"
#include "dbl2_tabs_kernels.h"
#include "dbl_test_utilities.h"
#include "dbl2_test_utilities.h"
#include "write_dbl2_bstimeflops.h"
#include "dbl_tabs_testers.h"
#include "dbl_data_files.h"

using namespace std;

void test_real2_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }
   // random_dbl2_upper_matrix(dim,dim,Ahi,Alo);
   dbl2_random_upper_factor(dim,Ahi,Alo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   double *solhi = new double[dim];
   double *sollo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solhi[i] = 1.0;
      sollo[i] = 0.0;
   }
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];
   double acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
      for(int j=0; j<dim; j++)  // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhshi[i] << "  " << rhslo[i] << endl;
   }
   double **invAhi_h = new double*[dim];
   double **invAlo_h = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAhi_h[i] = new double[dim];
      invAlo_h[i] = new double[dim];
   }
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_dbl2_upper_inverse(dim,Ahi,Alo,invAhi_h,invAlo_h,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invAhi_h[i][j] << "  " << invAlo_h[i][j] << endl;
   }
   double **invAhi_d = new double*[dim];
   double **invAlo_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      invAhi_d[i] = new double[dim];
      invAlo_d[i] = new double[dim];
   }
   cout << "-> GPU computes the inverse ..." << endl;

   GPU_dbl2_upper_inverse
      (dim,Ahi,Alo,invAhi_d,invAlo_d,&elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAhi_d[i][j] << "  " << invAlo_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
        << dbl2_Matrix_Difference_Sum
              (dim,invAhi_h,invAlo_h,invAhi_d,invAlo_d)
        << endl;

   double *xhi = new double[dim];
   double *xlo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      xhi[i] = 0.0;
      xlo[i] = 0.0;
      for(int j=0; j<dim; j++)   // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         ddf_mul(invAhi_h[i][j],invAlo_h[i][j],rhshi[j],rhslo[j],
                 &acchi,&acclo);
         ddf_inc(&xhi[i],&xlo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhi[i] << "  " << xlo[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on solution : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi,xlo) << endl;
   cout << "Condition number : "
        << dbl2_condition(dim,Ahi,Alo,invAhi_h,invAlo_h) << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
}

void test_cmplx2_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Arehi = new double*[dim];
   double **Arelo = new double*[dim];
   double **Aimhi = new double*[dim];
   double **Aimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Arehi[i] = new double[dim];
      Arelo[i] = new double[dim];
      Aimhi[i] = new double[dim];
      Aimlo[i] = new double[dim];
   }
   // random_cmplx2_upper_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);
   cmplx2_random_upper_factor(dim,Arehi,Arelo,Aimhi,Aimlo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi[i][j] << "  " << Arelo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
         }
   }
   double *solrehi = new double[dim];
   double *solrelo = new double[dim];
   double *solimhi = new double[dim];
   double *solimlo = new double[dim];

   for(int i=0; i<dim; i++) 
   {
      solrehi[i] = 1.0; solrelo[i] = 0.0;
      solimhi[i] = 0.0; solimlo[i] = 0.0;
   }
   double *rhsrehi = new double[dim];
   double *rhsrelo = new double[dim];
   double *rhsimhi = new double[dim];
   double *rhsimlo = new double[dim];
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;

   for(int i=0; i<dim; i++)
   {
      rhsrehi[i] = 0.0; rhsrelo[i] = 0.0;
      rhsimhi[i] = 0.0; rhsimlo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Arehi[i][j],Arelo[i][j],solrehi[j],solrelo[j],
                 &acc1hi,&acc1lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solimhi[j],solimlo[j],
                 &acc2hi,&acc2lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solrehi[j],solrelo[j],
                 &acc3hi,&acc3lo);
         ddf_mul(Arehi[i][j],Arelo[i][j],solimhi[j],solimlo[j],
                 &acc4hi,&acc4lo);
         ddf_inc(&rhsrehi[i],&rhsrelo[i],acc1hi,acc1lo);
         ddf_dec(&rhsrehi[i],&rhsrelo[i],acc2hi,acc2lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc3hi,acc3lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc4hi,acc4lo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehi[i] << "  " << rhsrelo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhi[i] << "  " << rhsimlo[i] << endl;
      }
   }
   double **invArehi_h = new double*[dim];
   double **invArelo_h = new double*[dim];
   double **invAimhi_h = new double*[dim];
   double **invAimlo_h = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invArehi_h[i] = new double[dim];
      invArelo_h[i] = new double[dim];
      invAimhi_h[i] = new double[dim];
      invAimlo_h[i] = new double[dim];
   }
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_cmplx2_upper_inverse
      (dim,   Arehi,     Arelo,     Aimhi,     Aimlo,
           invArehi_h,invArelo_h,invAimhi_h,invAimlo_h,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "invA_h[" << i << "][" << j << "]re : "
                 << invArehi_h[i][j] << "  " << invArelo_h[i][j] << endl;
            cout << "invA_h[" << i << "][" << j << "]im : "
                 << invAimhi_h[i][j] << "  " << invAimlo_h[i][j] << endl;
         }
   }

   double **invArehi_d = new double*[dim];
   double **invArelo_d = new double*[dim];
   double **invAimhi_d = new double*[dim];
   double **invAimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invArehi_d[i] = new double[dim];
      invArelo_d[i] = new double[dim];
      invAimhi_d[i] = new double[dim];
      invAimlo_d[i] = new double[dim];
   }
   cout << "-> GPU computes the inverse ..." << endl;

   GPU_cmplx2_upper_inverse
      (dim,   Arehi,     Arelo,     Aimhi,     Aimlo,
           invArehi_d,invArelo_d,invAimhi_d,invAimlo_d,
       &elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "invA_d[" << i << "][" << j << "]re : "
                 << invArehi_d[i][j] << "  " << invArelo_d[i][j] << endl;
            cout << "invA_d[" << i << "][" << j << "]im : "
                 << invAimhi_d[i][j] << "  " << invAimlo_d[i][j] << endl;
         }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
        << cmplx2_Matrix_Difference_Sum
              (dim,invArehi_h,invArelo_h,invAimhi_h,invAimlo_h,
                   invArehi_d,invArelo_d,invAimhi_d,invAimlo_d) << endl;

   double *xrehi = new double[dim];
   double *xrelo = new double[dim];
   double *ximhi = new double[dim];
   double *ximlo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      xrehi[i] = 0.0; xrelo[i] = 0.0;
      ximhi[i] = 0.0; ximlo[i] = 0.0;

      for(int j=0; j<dim; j++) // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         ddf_mul(invArehi_h[i][j],invArelo_h[i][j],rhsrehi[j],rhsrelo[j],
                 &acc1hi,&acc1lo);
         ddf_mul(invAimhi_h[i][j],invAimlo_h[i][j],rhsimhi[j],rhsimlo[j],
                 &acc2hi,&acc2lo);
         ddf_mul(invAimhi_h[i][j],invAimlo_h[i][j],rhsrehi[j],rhsrelo[j],
                 &acc3hi,&acc3lo);
         ddf_mul(invArehi_h[i][j],invArelo_h[i][j],rhsimhi[j],rhsimlo[j],
                 &acc4hi,&acc4lo);
         ddf_inc(&xrehi[i],&xrelo[i],acc1hi,acc1lo);
         ddf_dec(&xrehi[i],&xrelo[i],acc2hi,acc2lo);
         ddf_inc(&ximhi[i],&ximlo[i],acc3hi,acc3lo);
         ddf_inc(&ximhi[i],&ximlo[i],acc4hi,acc4lo);
      }
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehi[i] << "  " << xrelo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhi[i] << "  " << ximlo[i] << endl;
      }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on solution : "
        << cmplx2_Difference_Sum
              (dim,solrehi,solrelo,solimhi,solimlo,
                     xrehi,  xrelo,  ximhi,  ximlo) << endl;
   cout << "Condition number : "
        << cmplx2_condition(dim,   Arehi,     Arelo,     Aimhi,     Aimlo,
                                invArehi_h,invArelo_h,invAimhi_h,invAimlo_h)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
}

void test_real2_upper_tiling ( void )
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

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++) 
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
      for(int j=0; j<dim; j++) Alo[i][j] = 0.0;
   }
   // random_dbl_upper_matrix(dim,dim,A);
   // dbl2_random_upper_factor(dim,Ahi,Alo);
   if(rndmat == 1)
      dbl_random_upper_factor(dim,Ahi);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      dbl_read_matrix(filename,dim,Ahi);
   }
   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   double *solhi = new double[dim];
   double *sollo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solhi[i] = 1.0;
      sollo[i] = 0.0;
   }
   double *xhi = new double[dim];
   double *xlo = new double[dim];
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];
   double acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   double *xhi_d = new double[dim];
   double *xlo_d = new double[dim];
   double *rhshi_d = new double[dim];
   double *rhslo_d = new double[dim];
   double **Ahi_d = new double*[dim];
   double **Alo_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      rhshi_d[i] = rhshi[i];
      rhslo_d[i] = rhslo[i];
      Ahi_d[i] = new double[dim];
      Alo_d[i] = new double[dim];
      for(int j=0; j<dim; j++)
      {
         Ahi_d[i][j] = Ahi[i][j];
         Alo_d[i][j] = Alo[i][j];
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshi[i] << "  " << rhslo[i] << endl;
   }
   double timelapsed_h,timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_dbl2_upper_tiled_solver
      (dim,sizetile,numtiles,Ahi,Alo,rhshi,rhslo,xhi,xlo,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_dbl2_upper_tiled_solver
      (dim,sizetile,numtiles,Ahi_d,Alo_d,rhshi_d,rhslo_d,xhi_d,xlo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi_d[i][j] << "  " << Alo_d[i][j] << endl;
   }

   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << dbl2_Diagonal_Difference_Sum(numtiles,sizetile,Ahi,Alo,Ahi_d,Alo_d)
        << endl;

   if(verbose > 0)
   {
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhi[i] << "  " << xlo[i] << endl;
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhi_d[i] << "  " << xlo_d[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi,xlo) << endl;
   cout << "   Sum of GPU errors on solution : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi_d,xlo_d) << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl2_bstimeflops
     (sizetile,numtiles,0,invlapsed,mullapsed,sublapsed,elapsedms,
      timelapsed_d,addcnt,mulcnt,divcnt);

   for(int i=0; i<dim; i++)
   {
      free(Ahi[i]); free(Ahi_d[i]);
      free(Alo[i]); free(Alo_d[i]);
   }
   free(Ahi); free(Ahi_d); free(solhi);
   free(Alo); free(Alo_d); free(sollo);
   free(rhshi); free(rhshi_d); free(xhi); free(xhi_d);
   free(rhslo); free(rhslo_d); free(xlo); free(xlo_d);
}

void test_cmplx2_upper_tiling ( void )
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

   double **Arehi = new double*[dim];
   double **Arelo = new double*[dim];
   double **Aimhi = new double*[dim];
   double **Aimlo = new double*[dim];

   for(int i=0; i<dim; i++) 
   {
      Arehi[i] = new double[dim];
      Arelo[i] = new double[dim];
      for(int j=0; j<dim; j++) Arelo[i][j] = 0.0;
      Aimhi[i] = new double[dim];
      Aimlo[i] = new double[dim];
      for(int j=0; j<dim; j++) Aimlo[i][j] = 0.0;
   }

   // random_dbl_upper_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);
   // cmplx2_random_upper_factor(dim,Arehi,Arelo,Aimhi,Aimlo);
   if(rndmat == 1)
      cmplx_random_upper_factor(dim,Arehi,Aimhi);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      cmplx_read_matrix(filename,dim,Arehi,Aimhi);
   }
   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi[i][j] << "  " << Arelo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
         }
   }
   double *solrehi = new double[dim];
   double *solrelo = new double[dim];
   double *solimhi = new double[dim];
   double *solimlo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solrehi[i] = 1.0; solrelo[i] = 0.0;
      solimhi[i] = 0.0; solimlo[i] = 0.0;
   }
   double *xrehi = new double[dim];
   double *xrelo = new double[dim];
   double *ximhi = new double[dim];
   double *ximlo = new double[dim];
   double *rhsrehi = new double[dim];
   double *rhsrelo = new double[dim];
   double *rhsimhi = new double[dim];
   double *rhsimlo = new double[dim];
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;

   for(int i=0; i<dim; i++)
   {
      rhsrehi[i] = 0.0; rhsrelo[i] = 0.0;
      rhsimhi[i] = 0.0; rhsimlo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Arehi[i][j],Arelo[i][j],solrehi[j],solrelo[j],
                 &acc1hi,&acc1lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solimhi[j],solimlo[j],
                 &acc2hi,&acc2lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solrehi[j],solrelo[j],
                 &acc3hi,&acc3lo);
         ddf_mul(Arehi[i][j],Arelo[i][j],solimhi[j],solimlo[j],
                 &acc4hi,&acc4lo);
         ddf_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
         ddf_inc(&rhsrehi[i],&rhsrelo[i],acc1hi,acc1lo);
         ddf_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc3hi,acc3lo);
      }
   }
   double *xrehi_d = new double[dim];
   double *xrelo_d = new double[dim];
   double *ximhi_d = new double[dim];
   double *ximlo_d = new double[dim];
   double *rhsrehi_d = new double[dim];
   double *rhsrelo_d = new double[dim];
   double *rhsimhi_d = new double[dim];
   double *rhsimlo_d = new double[dim];
   double **Arehi_d = new double*[dim];
   double **Arelo_d = new double*[dim];
   double **Aimhi_d = new double*[dim];
   double **Aimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhsrehi_d[i] = rhsrehi[i];
      rhsrelo_d[i] = rhsrelo[i];
      rhsimhi_d[i] = rhsimhi[i];
      rhsimlo_d[i] = rhsimlo[i];
      Arehi_d[i] = new double[dim];
      Arelo_d[i] = new double[dim];
      Aimhi_d[i] = new double[dim];
      Aimlo_d[i] = new double[dim];

      for(int j=0; j<dim; j++)
      {
         Arehi_d[i][j] = Arehi[i][j];
         Arelo_d[i][j] = Arelo[i][j];
         Aimhi_d[i][j] = Aimhi[i][j];
         Aimlo_d[i][j] = Aimlo[i][j];
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehi[i] << "  " << rhsrelo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhi[i] << "  " << rhsimlo[i] << endl;
      }
   }
   double timelapsed_h,timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_cmplx2_upper_tiled_solver
      (dim,sizetile,numtiles,Arehi,Arelo,Aimhi,Aimlo,
       rhsrehi,rhsrelo,rhsimhi,rhsimlo,xrehi,xrelo,ximhi,ximlo,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi[i][j] << "  " << Arelo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
         }
   }
   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_cmplx2_upper_tiled_solver
      (dim,sizetile,numtiles,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
       rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
         xrehi_d,  xrelo_d,  ximhi_d,  ximlo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi_d[i][j] << "  " << Arelo_d[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi_d[i][j] << "  " << Aimlo_d[i][j] << endl;
         }
   }

   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << cmplx2_Diagonal_Difference_Sum
             (numtiles,sizetile,Arehi,  Arelo,  Aimhi,  Aimlo,
                                Arehi_d,Arelo_d,Aimhi_d,Aimlo_d)
        << endl;

   if(verbose > 0)
   {
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehi[i] << "  " << xrelo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhi[i] << "  " << ximlo[i] << endl;
      }
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehi_d[i] << "  " << xrelo_d[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhi_d[i] << "  " << ximlo_d[i] << endl;
      }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << cmplx2_Difference_Sum(dim,solrehi,solrelo,solimhi,solimlo,
                                       xrehi,  xrelo,  ximhi,  ximlo)
        << endl;
   cout << "   Sum of GPU errors on solution : "
        << cmplx2_Difference_Sum(dim,solrehi,solrelo,solimhi,solimlo,
                                       xrehi_d,xrelo_d,ximhi_d,ximlo_d)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl2_bstimeflops
     (sizetile,numtiles,1,invlapsed,mullapsed,sublapsed,elapsedms,
      timelapsed_d,addcnt,mulcnt,divcnt);

   for(int i=0; i<dim; i++)
   {
      free(Arehi[i]); free(Arehi_d[i]);
      free(Arelo[i]); free(Arelo_d[i]);
      free(Aimhi[i]); free(Aimhi_d[i]);
      free(Aimlo[i]); free(Aimlo_d[i]);
   }
   free(Arehi); free(Aimhi); free(Arehi_d); free(Aimhi_d);
   free(Arelo); free(Aimlo); free(Arelo_d); free(Aimlo_d);
   free(solrehi); free(solimhi);
   free(solrelo); free(solimlo);
   free(rhsrehi); free(rhsimhi); free(rhsrehi_d); free(rhsimhi_d);
   free(rhsrelo); free(rhsimlo); free(rhsrelo_d); free(rhsimlo_d);
   free(xrehi); free(ximhi); free(xrehi_d); free(ximhi_d);
   free(xrelo); free(ximlo); free(xrelo_d); free(ximlo_d);
}
