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
   CPU_dbl2_upper_inverse(dim,Ahi,Alo,invAhi_h,invAlo_h);

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
   GPU_dbl2_upper_inverse(dim,Ahi,Alo,invAhi_d,invAlo_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAhi_d[i][j] << "  " << invAlo_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors : "
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
   cout << "   Sum of errors : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi,xlo) << endl;
   cout << "Condition number : "
        << dbl2_condition(dim,Ahi,Alo,invAhi_h,invAlo_h) << endl;
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
   CPU_cmplx2_upper_inverse
      (dim,   Arehi,     Arelo,     Aimhi,     Aimlo,
           invArehi_h,invArelo_h,invAimhi_h,invAimlo_h);

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
   GPU_cmplx2_upper_inverse
      (dim,   Arehi,     Arelo,     Aimhi,     Aimlo,
           invArehi_d,invArelo_d,invAimhi_d,invAimlo_d);

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
   cout << "   Sum of errors : "
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
   cout << "   Sum of errors : "
        << cmplx2_Difference_Sum
              (dim,solrehi,solrelo,solimhi,solimlo,
                     xrehi,  xrelo,  ximhi,  ximlo) << endl;
   cout << "Condition number : "
        << cmplx2_condition(dim,   Arehi,     Arelo,     Aimhi,     Aimlo,
                                invArehi_h,invArelo_h,invAimhi_h,invAimlo_h)
        << endl;
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

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++) 
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }

   // random_dbl_upper_matrix(dim,dim,A);
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

   CPU_dbl2_upper_tiled_solver
      (dim,sizetile,numtiles,Ahi,Alo,rhshi,rhslo,xhi,xlo);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }

   GPU_dbl2_upper_tiled_solver
      (dim,sizetile,numtiles,Ahi_d,Alo_d,rhshi_d,rhslo_d,xhi_d,xlo_d);

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
   cout << "   Sum of errors : "
        << dbl2_Difference_Sum(dim,solhi,sollo,xhi,xlo) << endl;
}
