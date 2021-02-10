/* Tests operations on series vectors in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "random_series.h"
#include "random_matrices.h"
#include "dbl_matrices_host.h"

using namespace std;

void test_real_inner_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the inner product on random real data. */

void test_cmplx_inner_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the inner product on random complex data. */

void test_real_matrix_vector_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the matrix vector product on random real data. */

void test_cmplx_matrix_vector_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the matrix vector product on random complex data. */
  
int main ( void )
{
   cout << "testing a real inner product ..." << endl;
   test_real_inner_product();
   cout << "testing a complex inner product ..." << endl;
   test_cmplx_inner_product();
   cout << "testing a real matrix-vector product ..." << endl;
   test_real_matrix_vector_product();
   cout << "testing a complex matrix-vector product ..." << endl;
   test_cmplx_matrix_vector_product();

   return 0;
}

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

   double *ip = new double[deg+1];

   real_inner_product(dim,deg,px,mx,ip);

   cout << "the inner product :" << endl;
   for(int i=0; i<=deg; i++) cout << ip[i] << endl;
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

   cmplx_inner_product(dim,deg,pxre,pxim,mxre,mxim,ipre,ipim);

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

   real_matrix_vector_product(nbrows,nbcols,deg,mat,x,y);
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

   cmplx_matrix_vector_product
      (nbrows,nbcols,deg,matre,matim,xre,xim,yre,yim);
   for(int i=0; i<nbrows; i++)
   {
      cout << "y[" << i << "] :" << endl;
      for(int k=0; k<=deg; k++)
         cout << yre[i][k] << "  " << yim[i][k] << endl;
   }
}
