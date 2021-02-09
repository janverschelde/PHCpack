/* Tests operations on series vectors in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "random_series.h"
#include "random_matrices.h"
#include "dbl_convolutions_host.h"

using namespace std;

void real_inner_product
 ( int dim, int deg, double **x, double **y, double *z );
/*
 * DESCRIPTION :
 *   Computes the product of two real vectors x and y of power series
 *   and assigns the result to z.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors x and y;
 *   deg      truncation degree of the series;
 *   x        dim series truncated to degree deg;
 *   y        dim series truncated to degree deg;
 *   z        space for deg+1 doubles.
 *
 * ON RETURN :
 *   z        the sum of all x[k]*y[k] for k from 0 to dim-1,
 *            as a power series truncated to the degree deg. */

void cmplx_inner_product
 ( int dim, int deg,
   double **xre, double **xim, double **yre, double **yim,
   double *zre, double *zim );
/*
 * DESCRIPTION :
 *   Computes the product of two complex vectors x and y of power series
 *   and assigns the result to z.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors x and y;
 *   deg      truncation degree of the series;
 *   xre      real parts of dim series truncated to degree deg;
 *   xim      imaginary parts of dim series truncated to degree deg;
 *   yre      real parts of dim series truncated to degree deg;
 *   yim      imaginary parts of dim series truncated to degree deg;
 *   zre      space for deg+1 doubles;
 *   zim      space for deg+1 doubles.
 *
 * ON RETURN :
 *   zre      real parts of the sum of all x[k]*y[k] for k from 0 to dim-1,
 *            as a power series truncated to the degree deg;
 *   zim      imaginary parts of the sum of all x[k]*y[k] for k from 0
 *            to dim-1, as a power series truncated to the degree deg. */

void real_matrix_vector_product
 ( int rows, int cols, int deg, double ***A, double **x, double **y );
/*
 * DESCRIPTION :
 *   Computes the product y of the matrix A with x on real data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrix A
 *            and the dimension of y;
 *   cols     the number of columns in the matrix A
 *            and the dimension of x;
 *   deg      truncation degree of the series;
 *   A        matrix of dimensions rows and cols,
 *            of power series truncated at the degree deg;
 *   x        vector of dimension cols
 *            of power series truncated at the degree deg;
 *   y        space allocated for a vector of dimension rows 
 *            of power series truncated at the degree deg.
 *
 * ON RETURN :
 *   y        product of A with x. */

void cmplx_matrix_vector_product
 ( int rows, int cols, int deg, double ***Are, double ***Aim,
   double **xre, double **xim, double **yre, double **yim );
/*
 * DESCRIPTION :
 *   Computes the product y of the matrix A with x on complex data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrix A
 *            and the dimension of y;
 *   cols     the number of columns in the matrix A
 *            and the dimension of x;
 *   deg      truncation degree of the series;
 *   Are      real parts of a matrix of dimensions rows and cols,
 *            of power series truncated at the degree deg;
 *   Aim      imaginary parts of a matrix of dimensions rows and cols,
 *            of power series truncated at the degree deg;
 *   xre      real parts of a vector of dimension cols
 *            of power series truncated at the degree deg;
 *   xim      imaginary parts of a vector of dimension cols
 *            of power series truncated at the degree deg;
 *   yre      space allocated for a vector of dimension rows for
 *            the real parts of series truncated at the degree deg;
 *   yim      space allocated for a vector of dimension rows for
 *            the imaginary parts of series truncated at the degree deg.
 *
 * ON RETURN :
 *   yre      real parts of the product of A with x;
 *   yim      imaginary parts of the product of A with x. */

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

void real_inner_product
 ( int dim, int deg, double **x, double **y, double *z )
{
   double *prod = new double[deg+1];

   for(int i=0; i<=deg; i++) z[i] = 0.0;

   for(int k=0; k<dim; k++)
   {
      CPU_dbl_product(deg,x[k],y[k],prod);
      for(int i=0; i<=deg; i++)
         z[i] = z[i] + prod[i];
   }
   free(prod);
}

void cmplx_inner_product
 ( int dim, int deg,
   double **xre, double **xim, double **yre, double **yim,
   double *zre, double *zim )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      zre[i] = 0.0; zim[i] = 0.0;
   }
   for(int k=0; k<dim; k++)
   {
      CPU_cmplx_product(deg,xre[k],xim[k],yre[k],yim[k],prodre,prodim);
      for(int i=0; i<=deg; i++)
      {
         zre[i] = zre[i] + prodre[i];
         zim[i] = zim[i] + prodim[i];
      }
   }
   free(prodre); free(prodim);
}

void real_matrix_vector_product
 ( int rows, int cols, int deg, double ***A, double **x, double **y )
{
   for(int k=0; k<rows; k++)
      real_inner_product(cols,deg,A[k],x,y[k]);
}

void cmplx_matrix_vector_product
 ( int rows, int cols, int deg, double ***Are, double ***Aim,
   double **xre, double **xim, double **yre, double **yim )
{
   for(int k=0; k<rows; k++)
      cmplx_inner_product(cols,deg,Are[k],Aim[k],xre,xim,yre[k],yim[k]);
}
