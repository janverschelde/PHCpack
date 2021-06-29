/* Tests operations on series vectors in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "random_matrices.h"
#include "dbl_linearization.h"

using namespace std;

int test_dbl_vector_linearization ( int dim, int deg );
/*
 * Tests the linearization of a random vector of dimension dim,
 * with series truncated at degree deg, for real data. */

int test_cmplx_vector_linearization ( int dim, int deg );
/*
 * Tests the linearization of a random vector of dimension dim,
 * with series truncated at degree deg, for complex data. */

int test_dbl_matrix_linearization ( int nrows, int ncols, int deg );
/*
 * Tests the linearization of a random nrows-by-ncols matrix
 * with series truncated at degree deg, for real data. */

int test_cmplx_matrix_linearization ( int nrows, int ncols, int deg );
/*
 * Tests the linearization of a random nrows-by_ncols matrix,
 * with series truncated at degree deg, for complex data. */

int main ( void )
{
   srand(time(NULL));

   cout << "Give the dimension of the vector : ";
   int dim; cin >> dim;

   cout << "Give the truncation degree of the series : ";
   int deg; cin >> deg;

   cout << "Testing on real data ..." << endl;
   test_dbl_vector_linearization(dim,deg);

   cout << "Testing on complex data ..." << endl;
   test_cmplx_vector_linearization(dim,deg);

   cout << "Testing on matrices ..." << endl;
   cout << "-> give the number of rows : ";
   int nrows; cin >> nrows;
   cout << "-> give the number of columns : ";
   int ncols; cin >> ncols;
   cout << "-> give the truncation degree : "; cin >> deg;

   test_dbl_matrix_linearization(nrows,ncols,deg);

   cout << "Testing on complex data ..." << endl;
   test_cmplx_matrix_linearization(nrows,ncols,deg);

   return 0;
}

int test_dbl_vector_linearization ( int dim, int deg )
{
   double *x = new double[dim];
   double **v = new double*[dim];

   for(int i=0; i<dim; i++) v[i] = new double[deg+1];

   random_dbl_series_vector(dim,deg,x,v,false);

   cout << scientific << setprecision(16);

   cout << "A vector of random series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "v[" << i << "] :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << v[i][j];
         if(j == 0)
            cout << endl;
         else if(j == 1)
            cout << "*t" << endl;
         else
            cout << "*t^" << j << endl;
      }
   }
   double **w = new double*[deg+1];

   for(int i=0; i<=deg; i++) w[i] = new double[dim];

   dbl_linear_series_vector(dim,deg,v,w);

   cout << "The coefficient vectors of the series :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << "w[" << i << "] :" << endl;
      for(int j=0; j<dim; j++)
         cout << w[i][j] << endl;
   }
   return 0;
}

int test_cmplx_vector_linearization ( int dim, int deg )
{
   double *xre = new double[dim];
   double *xim = new double[dim];
   double **vre = new double*[dim];
   double **vim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      vre[i] = new double[deg+1];
      vim[i] = new double[deg+1];
   }
   random_cmplx_series_vector(dim,deg,xre,xim,vre,vim,false);

   cout << scientific << setprecision(16);

   cout << "A vector of random series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "v[" << i << "] :" << endl;
      for(int j=0; j<=deg; j++)
      {
         if(j > 0) cout << "+ ( ";
         cout << vre[i][j] << "  " << vim[i][j];
         if(j == 0)
            cout << endl;
         else if(j == 1)
            cout << " )*t" << endl;
         else
            cout << " )*t^" << j << endl;
      }
   }
   double **wre = new double*[deg+1];
   double **wim = new double*[deg+1];

   for(int i=0; i<=deg; i++)
   {
      wre[i] = new double[dim];
      wim[i] = new double[dim];
   }
   cmplx_linear_series_vector(dim,deg,vre,vim,wre,wim);

   cout << "The coefficient vectors of the series :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << "w[" << i << "] :" << endl;
      for(int j=0; j<dim; j++)
         cout << wre[i][j] << "  " << wim[i][j] << endl;
   }
   return 0;
}

int test_dbl_matrix_linearization ( int nrows, int ncols, int deg )
{
   double **rnd = new double*[nrows];
   double ***mat = new double**[nrows];

   for(int i=0; i<nrows; i++)
   {
      rnd[i] = new double[ncols];
      mat[i] = new double*[ncols];
      for(int j=0; j<ncols; j++)
         mat[i][j] = new double[deg+1];
   }
   random_dbl_series_matrix(nrows,ncols,deg,rnd,mat,false);

   cout << scientific << setprecision(16);

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         cout << "A[" << i << "][" << j << "] :" << endl;
         for(int k=0; k<=deg; k++)
         {
            cout << mat[i][j][k];
            if(k == 0)
               cout << endl;
            else if(k == 1)
               cout << "*t" << endl;
            else
               cout << "*t^" << k << endl;
         }
      }

   double ***linmat = new double**[deg+1];

   for(int k=0; k<=deg; k++)
   {
      linmat[k] = new double*[nrows];
      for(int i=0; i<nrows; i++)
         linmat[k][i] = new double[ncols];
   }
   dbl_linear_series_matrix(nrows,ncols,deg,mat,linmat);

   cout << "The coefficient matrices of the series :" << endl;
   for(int k=0; k<=deg; k++)
   {
      cout << "B[" << k << "] :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "[" << i << "][" << j << "] : "
                 << linmat[k][i][j] << endl;
   }

   return 0;
}

int test_cmplx_matrix_linearization ( int nrows, int ncols, int deg )
{
   double **rndre = new double*[nrows];
   double **rndim = new double*[nrows];
   double ***Are = new double**[nrows];
   double ***Aim = new double**[nrows];

   for(int i=0; i<nrows; i++)
   {
      rndre[i] = new double[ncols];
      rndim[i] = new double[ncols];
      Are[i] = new double*[ncols];
      Aim[i] = new double*[ncols];
      for(int j=0; j<ncols; j++)
      {
         Are[i][j] = new double[deg+1];
         Aim[i][j] = new double[deg+1];
      }
   }
   random_cmplx_series_matrix(nrows,ncols,deg,rndre,rndim,Are,Aim,false);

   cout << scientific << setprecision(16);

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         cout << "A[" << i << "][" << j << "] :" << endl;
         for(int k=0; k<=deg; k++)
         {
            if(k > 0) cout << "+ ( ";
            cout << Are[i][j][k] << "  " << Aim[i][j][k];
            if(k == 0)
               cout << endl;
            else if(k == 1)
               cout << " )*t" << endl;
            else
               cout << " )*t^" << k << endl;
         }
      }

   double ***Bre = new double**[deg+1];
   double ***Bim = new double**[deg+1];

   for(int k=0; k<=deg; k++)
   {
      Bre[k] = new double*[nrows];
      Bim[k] = new double*[nrows];
      for(int i=0; i<nrows; i++)
      {
         Bre[k][i] = new double[ncols];
         Bim[k][i] = new double[ncols];
      }
   }
   cmplx_linear_series_matrix(nrows,ncols,deg,Are,Aim,Bre,Bim);

   cout << "The coefficient matrices of the series :" << endl;
   for(int k=0; k<=deg; k++)
   {
      cout << "B[" << k << "] :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "[" << i << "][" << j << "] : "
                 << Bre[k][i][j] << "  " << Bim[k][i][j] << endl;
   }

   return 0;
}
