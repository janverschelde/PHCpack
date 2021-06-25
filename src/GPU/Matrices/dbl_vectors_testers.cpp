// The file dbl_vectors_testers.cpp defines the function with prototypes
// in dbl_vectors_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
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

   cout << "Give the verbose level : ";
   int vrblvl; cin >> vrblvl;

   double *x = new double[dim];
   double **px = new double*[dim];
   double **mx = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      px[i] = new double[deg+1];
      mx[i] = new double[deg+1];
   }
   random_dbl_series_vectors(dim,deg,x,px,mx);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);

      for(int k=0; k<dim; k++)
      {
         cout << "a random x : " << x[k] << endl;
         cout << "series for exp(+x) :" << endl; 
         for(int i=0; i<=deg; i++) cout << px[k][i] << endl;
         cout << "series for exp(-x) :" << endl; 
         for(int i=0; i<=deg; i++) cout << mx[k][i] << endl;
      }
   }
   double *ip_h = new double[deg+1];
   double *ip_d = new double[deg+1];

   CPU_dbl_inner_product(dim,deg,px,mx,ip_h);

   if(vrblvl > 1)
   {
      cout << "the inner product computed by the CPU :" << endl;
      for(int i=0; i<=deg; i++) cout << ip_h[i] << endl;
   }
   if(vrblvl > 0)
      GPU_dbl_inner_product(deg+1,dim,deg,px,mx,ip_d,0,true);
   else
      GPU_dbl_inner_product(deg+1,dim,deg,px,mx,ip_d,0,false);

   if(vrblvl > 1)
   {
      cout << "the inner product computed by the GPU :" << endl;
      for(int i=0; i<=deg; i++) cout << ip_d[i] << endl;
   }
   double sumerr = 0.0;

   for(int i=0; i<=deg; i++)
      sumerr = sumerr + abs(ip_h[i] - ip_d[i]); 
   cout << scientific << setprecision(2);
   cout << "Sum of all errors : " << sumerr << endl;
}

void test_cmplx_inner_product ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrblvl; cin >> vrblvl;

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

   if(vrblvl > 1)
   {
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
   }
   double *ipre_h = new double[deg+1];
   double *ipim_h = new double[deg+1];
   double *ipre_d = new double[deg+1];
   double *ipim_d = new double[deg+1];

   CPU_cmplx_inner_product(dim,deg,pxre,pxim,mxre,mxim,ipre_h,ipim_h);

   if(vrblvl > 1)
   {
      cout << "the inner product computed by the CPU :" << endl;
      for(int i=0; i<=deg; i++)
         cout << ipre_h[i] << "  " << ipim_h[i] << endl;
   }
   if(vrblvl > 0)
      GPU_cmplx_inner_product
         (deg+1,dim,deg,pxre,pxim,mxre,mxim,ipre_d,ipim_d,0,true);
   else
      GPU_cmplx_inner_product
         (deg+1,dim,deg,pxre,pxim,mxre,mxim,ipre_d,ipim_d,0,false);

   if(vrblvl > 1)
   {
      cout << "the inner product computed by the GPU :" << endl;
      for(int i=0; i<=deg; i++)
         cout << ipre_d[i] << "  " << ipim_d[i] << endl;
   }
   double sumerr = 0.0;
   for(int i=0; i<=deg; i++)
      sumerr = sumerr + abs(ipre_h[i] - ipre_d[i])
                      + abs(ipim_h[i] - ipim_d[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of all errors : " << sumerr << endl;
}

void test_real_matrix_vector_product ( void )
{
   cout << "Give the number of rows : ";
   int nbrows; cin >> nbrows;

   cout << "Give the number of columns : ";
   int nbcols; cin >> nbcols;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrblvl; cin >> vrblvl;

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

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);

      for(int i=0; i<nbrows; i++)
         for(int j=0; j<nbcols; j++)
         {
            cout << "A[" << i << "][" << j << "] is exp("
                 << rnd[i][j] << ") :" << endl;
            for(int k=0; k<=deg; k++) cout << mat[i][j][k] << endl;
         }
   }
   double **x = new double*[nbcols];
   for(int i=0; i<nbcols; i++) x[i] = new double[deg+1];
   double **y_h = new double*[nbrows];
   for(int i=0; i<nbrows; i++) y_h[i] = new double[deg+1];
   double **y_d = new double*[nbrows];
   for(int i=0; i<nbrows; i++) y_d[i] = new double[deg+1];

   for(int i=0; i<nbcols; i++)
      dbl_exponential(deg,-rnd[nbrows-1][i],x[i]);

   CPU_dbl_matrix_vector_product(nbrows,nbcols,deg,mat,x,y_h);

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrows; i++)
      {
         cout << "y[" << i << "] :" << endl;
         for(int k=0; k<=deg; k++) cout << y_h[i][k] << endl;
      }
   }
   if(vrblvl > 0)
      GPU_dbl_matrix_vector_product
         (deg+1,nbrows,nbcols,deg,mat,x,y_d,0,true);
   else
      GPU_dbl_matrix_vector_product
         (deg+1,nbrows,nbcols,deg,mat,x,y_d,0,false);

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrows; i++)
      {
         cout << "y[" << i << "] :" << endl;
         for(int k=0; k<=deg; k++) cout << y_d[i][k] << endl;
      }
   }
   double sumerr = 0.0;

   for(int i=0; i<nbrows; i++)
      for(int j=0; j<=deg; j++)
         sumerr = sumerr + abs(y_h[i][j] - y_d[i][j]); 

   cout << scientific << setprecision(2);
   cout << "Sum of all errors : " << sumerr << endl;
}

void test_cmplx_matrix_vector_product ( void )
{
   cout << "Give the number of rows : ";
   int nbrows; cin >> nbrows;

   cout << "Give the number of columns : ";
   int nbcols; cin >> nbcols;

   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrblvl; cin >> vrblvl;

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

   if(vrblvl > 1)
   {
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
   }
   double **xre = new double*[nbcols];
   double **xim = new double*[nbcols];
   for(int i=0; i<nbcols; i++)
   {
      xre[i] = new double[deg+1];
      xim[i] = new double[deg+1];
   }
   double **yre_h = new double*[nbrows];
   double **yim_h = new double*[nbrows];
   double **yre_d = new double*[nbrows];
   double **yim_d = new double*[nbrows];
   for(int i=0; i<nbrows; i++)
   {
      yre_h[i] = new double[deg+1];
      yim_h[i] = new double[deg+1];
      yre_d[i] = new double[deg+1];
      yim_d[i] = new double[deg+1];
   }
   for(int i=0; i<nbcols; i++)
      cmplx_exponential
         (deg,-rndre[nbrows-1][i],-rndim[nbrows-1][i],xre[i],xim[i]);

   CPU_cmplx_matrix_vector_product
      (nbrows,nbcols,deg,matre,matim,xre,xim,yre_h,yim_h);

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrows; i++)
      {
         cout << "y[" << i << "] :" << endl;
         for(int k=0; k<=deg; k++)
            cout << yre_h[i][k] << "  " << yim_h[i][k] << endl;
      }
   }
   if(vrblvl > 0)
      GPU_cmplx_matrix_vector_product
         (deg+1,nbrows,nbcols,deg,matre,matim,xre,xim,yre_d,yim_d,0,true);
   else
      GPU_cmplx_matrix_vector_product
         (deg+1,nbrows,nbcols,deg,matre,matim,xre,xim,yre_d,yim_d,0,false);

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrows; i++)
      {
         cout << "y[" << i << "] :" << endl;
         for(int k=0; k<=deg; k++)
            cout << yre_d[i][k] << "  " << yim_d[i][k] << endl;
      }
   }
   double sumerr = 0.0;

   for(int i=0; i<nbrows; i++)
      for(int j=0; j<=deg; j++)
         sumerr = sumerr + abs(yre_h[i][j] - yre_d[i][j])
                         + abs(yim_h[i][j] - yim_d[i][j]); 

   cout << scientific << setprecision(2);
   cout << "Sum of all errors : " << sumerr << endl;
}
