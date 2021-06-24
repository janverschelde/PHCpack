// The file dbl_vectors_testers.cpp defines the function with prototypes
// in dbl_vectors_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector_types.h>
#include "random_numbers.h"
#include "random_series.h"
#include "random_matrices.h"
#include "dbl_matrices_host.h"
#include "dbl_matrices_kernels.h"
#include "dbl_vectors_testers.h"

using namespace std;

void test_real_inner_product ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   double *x = new double[dim];
   double **px = new double*[dim];
   double **mx = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      px[i] = new double[deg+1];
      mx[i] = new double[deg+1];
   }
   random_dbl_series_vectors(dim,deg,x,px,mx);

   cout << scientific << setprecision(16);

   for(int k=0; k<dim; k++)
   {
      cout << "a random x : " << x[k] << endl;
      cout << "series for exp(+x) :" << endl; 
      for(int i=0; i<=deg; i++) cout << px[k][i] << endl;
      cout << "series for exp(-x) :" << endl; 
      for(int i=0; i<=deg; i++) cout << mx[k][i] << endl;
   }

   double *ip_h = new double[deg+1];
   double *ip_d = new double[deg+1];

   CPU_dbl_inner_product(dim,deg,px,mx,ip_h);

   cout << "the inner product computed by the CPU :" << endl;
   for(int i=0; i<=deg; i++) cout << ip_h[i] << endl;

   GPU_dbl_inner_product(deg+1,dim,deg,px,mx,ip_d,1);

   cout << "the inner product computed by the GPU :" << endl;
   for(int i=0; i<=deg; i++) cout << ip_d[i] << endl;
}

void test_cmplx_inner_product ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   double *xre = new double[dim];
   double *xim = new double[dim];
   double **pxre = new double*[dim];
   double **pxim = new double*[dim];
   double **mxre = new double*[dim];
   double **mxim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      pxre[i] = new double[deg+1];
      pxim[i] = new double[deg+1];
      mxre[i] = new double[deg+1];
      mxim[i] = new double[deg+1];
   }
   random_cmplx_series_vectors
      (dim,deg,xre,xim,pxre,pxim,mxre,mxim);

   cout << scientific << setprecision(16);

   for(int k=0; k<dim; k++)
   {
      cout << "a random x : "
           << xre[k] << "  " << xim[k] << endl;
      cout << "series for exp(+x) :" << endl; 
      for(int i=0; i<=deg; i++)
         cout << pxre[k][i] << "  " << pxim[k][i] << endl;
      cout << "series for exp(-x) :" << endl; 
      for(int i=0; i<=deg; i++)
         cout << mxre[k][i] << "  " << mxim[k][i] << endl;
   }

   double *ipre = new double[deg+1];
   double *ipim = new double[deg+1];

   CPU_cmplx_inner_product(dim,deg,pxre,pxim,mxre,mxim,ipre,ipim);

   cout << "the inner product :" << endl;
   for(int i=0; i<=deg; i++)
      cout << ipre[i] << "  " << ipim[i] << endl;
}

void test_real_matrix_vector_product ( void )
{
   cout << "Give the number of rows : ";
   int nbrows; cin >> nbrows;

   cout << "Give the number of columns : ";
   int nbcols; cin >> nbcols;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   double **rnd = new double*[nbrows];
   double ***mat = new double**[nbrows];
   for(int i=0; i<nbrows; i++)
   {
      rnd[i] = new double[nbcols];
      mat[i] = new double*[nbcols];
      for(int j=0; j<nbcols; j++)
         mat[i][j] = new double[deg+1];
   }
   random_dbl_series_matrix(nbrows,nbcols,deg,rnd,mat);

   cout << scientific << setprecision(16);

   for(int i=0; i<nbrows; i++)
      for(int j=0; j<nbcols; j++)
      {
         cout << "A[" << i << "][" << j << "] is exp("
              << rnd[i][j] << ") :" << endl;
         for(int k=0; k<=deg; k++) cout << mat[i][j][k] << endl;
      }
   
   double **x = new double*[nbcols];
   for(int i=0; i<nbcols; i++) x[i] = new double[deg+1];
   double **y = new double*[nbrows];
   for(int i=0; i<nbrows; i++) y[i] = new double[deg+1];

   for(int i=0; i<nbcols; i++)
      dbl_exponential(deg,-rnd[nbrows-1][i],x[i]);

   CPU_dbl_matrix_vector_product(nbrows,nbcols,deg,mat,x,y);

   for(int i=0; i<nbrows; i++)
   {
      cout << "y[" << i << "] :" << endl;
      for(int k=0; k<=deg; k++) cout << y[i][k] << endl;
   }
}

void test_cmplx_matrix_vector_product ( void )
{
   cout << "Give the number of rows : ";
   int nbrows; cin >> nbrows;

   cout << "Give the number of columns : ";
   int nbcols; cin >> nbcols;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   double **rndre = new double*[nbrows];
   double **rndim = new double*[nbrows];
   double ***matre = new double**[nbrows];
   double ***matim = new double**[nbrows];
   for(int i=0; i<nbrows; i++)
   {
      rndre[i] = new double[nbcols];
      rndim[i] = new double[nbcols];
      matre[i] = new double*[nbcols];
      matim[i] = new double*[nbcols];
      for(int j=0; j<nbcols; j++)
      {
         matre[i][j] = new double[deg+1];
         matim[i][j] = new double[deg+1];
      }
   }
   random_cmplx_series_matrix(nbrows,nbcols,deg,rndre,rndim,matre,matim);

   cout << scientific << setprecision(16);

   for(int i=0; i<nbrows; i++)
      for(int j=0; j<nbcols; j++)
      {
         cout << "A[" << i << "][" << j << "] is exp("
              << rndre[i][j] << ", "
              << rndim[i][j] << ") :" << endl;
         for(int k=0; k<=deg; k++)
            cout << matre[i][j][k] << "  " << matim[i][j][k] << endl;
      }
   
   double **xre = new double*[nbcols];
   double **xim = new double*[nbcols];
   for(int i=0; i<nbcols; i++)
   {
      xre[i] = new double[deg+1];
      xim[i] = new double[deg+1];
   }
   double **yre = new double*[nbrows];
   double **yim = new double*[nbrows];
   for(int i=0; i<nbrows; i++)
   {
      yre[i] = new double[deg+1];
      yim[i] = new double[deg+1];
   }
   for(int i=0; i<nbcols; i++)
      cmplx_exponential
         (deg,-rndre[nbrows-1][i],-rndim[nbrows-1][i],xre[i],xim[i]);

   CPU_cmplx_matrix_vector_product
      (nbrows,nbcols,deg,matre,matim,xre,xim,yre,yim);

   for(int i=0; i<nbrows; i++)
   {
      cout << "y[" << i << "] :" << endl;
      for(int k=0; k<=deg; k++)
         cout << yre[i][k] << "  " << yim[i][k] << endl;
   }
}
