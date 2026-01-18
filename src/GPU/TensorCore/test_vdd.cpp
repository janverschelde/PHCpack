/* Tests the collection of function on vectored double double arithmetic. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "double_double.h"
#include "random2_vectors.h"
#include "double_double_functions.h"
#include "vectored_double_doubles.h"
#include "random2_matrices.h"

int test_quarter_double_double ( void );
/*
 * Generates a random positive double double, quarters the high and low
 * parts, and then checks if their sum equals the original double double. */

int test_vectored_dd_product ( int dim );
/*
 * Generates two random vectors of double doubles of size dim,
 * and compares their inner producted with the vectored inner product. */

int test_vectored_dd_matmatmul ( int nrows, int ncols, int nrc ); 
/*
 * Generates two random double double matrices of dimension
 * nrows-by-nrc for A, nrc-by-ncols for B, and then tests
 * the matrix matrix multiplication of A with B. */

using namespace std;

int main ( void )
{
   int fail;

   srand(time(NULL));

   fail = test_quarter_double_double();

   if(fail == 1)
      cout << "\nTest on quarter double double failed?!!!\n\n";
   else
      cout << "\nTest on quarter double double succeeded.\n\n";

   cout << "Give the dimension : ";
   int dim; cin >> dim;

   fail = test_vectored_dd_product(dim);

   if(fail == 1)
      cout << "\nTest on vectored double double product failed?!!!\n\n";
   else
      cout << "\nTest on vectored double double product succeeded.\n\n";

   cout << "Give #rows of the product A*B : ";
   int m; cin >> m;
   cout << "Give #columns of the product A*B : ";
   int n; cin >> n;
   cout << "Give #columns of A, #rows of B : ";
   int k; cin >> k;

   fail = test_vectored_dd_matmatmul(m, n, k);

   return 0;
}

int test_quarter_double_double ( void )
{
   int fail = 0;

   double x[2];
   double xhi0,xhi1,xhi2,xhi3;
   double xlo0,xlo1,xlo2,xlo3;
   double y[2];
   double e[2];

   random_double_double(&x[0], &x[1]);

   if(x[0] < 0.0) x[0] = -x[0]; // both high and low part
   if(x[1] < 0.0) x[1] = -x[1]; // must be positive

   cout << scientific << setprecision(16);

   cout << "x : "; dd_write(x, 32); cout << endl;

   quarter_double_double
      (x[0], x[1], &xhi0, &xhi1, &xhi2, &xhi3, &xlo0, &xlo1, &xlo2, &xlo3);

   cout << "xhi0 : " << xhi0 << endl;
   cout << "xhi1 : " << xhi1 << endl;
   cout << "xhi2 : " << xhi2 << endl;
   cout << "xhi3 : " << xhi3 << endl;
   cout << "xlo0 : " << xlo0 << endl;
   cout << "xlo1 : " << xlo1 << endl;
   cout << "xlo2 : " << xlo2 << endl;
   cout << "xlo3 : " << xlo3 << endl;

   to_double_double
      (xhi0, xhi1, xhi2, xhi3, xlo0, xlo1, xlo2, xlo3, &y[0], &y[1]);

   cout << "y : "; dd_write(y, 32); cout << endl;

   ddf_sub(x[0], x[1], y[0], y[1], &e[0], &e[1]);

   cout << "e : "; dd_write(e, 32); cout << endl;

   fail = not(e[0] == 0.0)
        + not(e[1] == 0.0);

   return fail;
}

int test_vectored_dd_product ( int dim )
{
   int fail = 0;

   double xhi[dim],xlo[dim],yhi[dim],ylo[dim];
   double x0[dim],x1[dim],x2[dim],x3[dim],x4[dim],x5[dim],x6[dim],x7[dim];
   double y0[dim],y1[dim],y2[dim],y3[dim],y4[dim],y5[dim],y6[dim],y7[dim];
   double prd[2],vpd[2],err[2];
   double s0,s1,s2,s3,s4,s5,s6,s7;

   for(int i=0; i<dim; i++)
   {
      random_double_double(&xhi[i], &xlo[i]);
      if(xhi[i] < 0.0) xhi[i] = -xhi[i];
      if(xlo[i] < 0.0) xlo[i] = -xlo[i];
      random_double_double(&yhi[i], &ylo[i]);
      if(yhi[i] < 0.0) yhi[i] = -yhi[i];
      if(ylo[i] < 0.0) ylo[i] = -ylo[i];
   }
   cout << scientific << setprecision(16);

   cout << "double double vector x :" << endl;
   dd_write_vector(dim, xhi, xlo);
   cout << "double double vector y :" << endl;
   dd_write_vector(dim, yhi, ylo);

   double_double_product(dim, xhi, xlo, yhi, ylo, &prd[0], &prd[1]);

   quarter_dd_vector(dim, xhi, xlo, x0, x1, x2, x3, x4, x5, x6, x7);
   quarter_dd_vector(dim, yhi, ylo, y0, y1, y2, y3, y4, y5, y6, y7);

   vectored_dd_product
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7,
       &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7);

   to_double_double(s0, s1, s2, s3, s4, s5, s6, s7, &vpd[0], &vpd[1]);

   cout << "dd x*y : "; dd_write(prd, 32); cout << endl;
   cout << "vd x*y : "; dd_write(vpd, 32); cout << endl;

   ddf_sub(prd[0], prd[1], vpd[0], vpd[1], &err[0], &err[1]);

   if(err[0] < 0.0) ddf_minus(&err[0], &err[1]);

   cout << " error : "; dd_write(err, 32); cout << endl;

   fail = (fabs(err[0]) > 1.0E-28);

   return fail;
}

int test_vectored_dd_matmatmul ( int nrows, int ncols, int nrc )
{
   int fail = 0;

   double **Chi = new double*[nrows];
   double **Clo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Chi[i] = new double[ncols];
      Clo[i] = new double[ncols];
   }
   double **Ahi = new double*[nrows];
   double **Alo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahi[i] = new double[nrc];
      Alo[i] = new double[nrc];
   }
   random_dbl2_matrix(nrows, nrc, Ahi, Alo);

   cout << scientific << setprecision(16);

   cout << "A random " << nrows << "-by-" << nrc << " matrix A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
         cout << "A[" << i << "][" << j << "] : "
              << Ahi[i][j] << "  " << Alo[i][j] << endl;

   double **Bhi = new double*[nrc];
   double **Blo = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhi[i] = new double[ncols];
      Blo[i] = new double[ncols];
   }
   random_dbl2_matrix(nrc, ncols, Bhi, Blo);

   cout << "A random " << nrc << "-by-" << ncols << " matrix B :" << endl;
   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
         cout << "B[" << i << "][" << j << "] : "
              << Bhi[i][j] << "  " << Blo[i][j] << endl;

   double_double_matmatmul(nrows, ncols, nrc, Ahi, Alo, Bhi, Blo, Chi, Clo);

   cout << "the product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "C[" << i << "][" << j << "] : "
              << Chi[i][j] << "  " << Clo[i][j] << endl;

   return fail;
}
