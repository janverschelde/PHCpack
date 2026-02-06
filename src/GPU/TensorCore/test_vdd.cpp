/* Tests the collection of function on vectored double double arithmetic. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "double_double.h"
#include "random2_vectors.h"
#include "double_double_functions.h"
#include "splitting_doubles.h"
#include "vectored_double_doubles.h"
#include "random2_matrices.h"
#include "double_matrix_multiplications.h"

int test_quarter_double_double ( void );
/*
 * Generates a random positive double double, quarters the high and low
 * parts, and then checks if their sum equals the original double double. */

int test_vectored_dd_product ( int dim );
/*
 * Generates two random vectors of double doubles of size dim,
 * and compares their inner product with the vectored inner product. */

int test_vectored_dd_matmatmul ( int nrows, int ncols, int nrc ); 
/*
 * Generates two random double double matrices of dimension
 * nrows-by-nrc for A, nrc-by-ncols for B, and then tests
 * the matrix matrix multiplication of A with B. */

using namespace std;

int main ( void )
{
   cout << "Give the seed for the random numbers (0 by default) : ";
   int seed; cin >> seed;

   if(seed == 0) seed = time(NULL);
   srand(seed);

   int fail = test_quarter_double_double();

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

   if(fail == 1)
      cout << "\nTest on vectored double double matmatmul failed?!!!\n\n";
   else
      cout << "\nTest on vectored double double matmatmul succeeded.\n\n";

   cout << "Seed used : " << seed << endl;

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

   make_dd_exponent_zero(&x[0], &x[1], 2);

   cout << scientific << setprecision(16);

   cout << "x : "; dd_write(x, 32); cout << endl;

   quarter_double_double
      (x[0], x[1], &xhi0, &xhi1, &xhi2, &xhi3, &xlo0, &xlo1, &xlo2, &xlo3, 2);

   cout << "xhi0 : " << xhi0 << endl;
   cout << "xhi1 : " << xhi1 << endl;
   cout << "xhi2 : " << xhi2 << endl;
   cout << "xhi3 : " << xhi3 << endl;
   cout << "xlo0 : " << xlo0 << endl;
   cout << "xlo1 : " << xlo1 << endl;
   cout << "xlo2 : " << xlo2 << endl;
   cout << "xlo3 : " << xlo3 << endl;

   if(is_dd_quarter_balanced
        (xhi0, xhi1, xhi2, xhi3, xlo0, xlo1, xlo2, xlo3, 2))
      cout << "The quarters are balanced." << endl;
   else
      cout << "The quarters are NOT balanced!?" << endl;

   to_double_double8sum
      (xhi0, xhi1, xhi2, xhi3, xlo0, xlo1, xlo2, xlo3, &y[0], &y[1]);

   cout << "y : "; dd_write(y, 32); cout << endl;

   ddf_sub(x[0], x[1], y[0], y[1], &e[0], &e[1]);

   cout << "e : "; dd_write(e, 3); cout << endl;

   cout << "xhi0 : "; write_52double(xhi0);
   cout << "xhi1 : "; write_52double(xhi1);
   cout << "xhi2 : "; write_52double(xhi2);
   cout << "xhi3 : "; write_52double(xhi3);
   cout << "xlo0 : "; write_52double(xlo0);
   cout << "xlo1 : "; write_52double(xlo1);
   cout << "xlo2 : "; write_52double(xlo2);
   cout << "xlo3 : "; write_52double(xlo3);

   fail = not(e[0] == 0.0)
        + not(e[1] == 0.0);

   return fail;
}

int test_vectored_dd_product ( int dim )
{
   int fail = 0;

   double xhi[dim],xlo[dim],yhi[dim],ylo[dim];
   double x0[dim],x1[dim],x2[dim],x3[dim],x4[dim],x5[dim],x6[dim],x7[dim];
   double x8[dim],x9[dim],x10[dim],x11[dim];
   double y0[dim],y1[dim],y2[dim],y3[dim],y4[dim],y5[dim],y6[dim],y7[dim];
   double y8[dim],y9[dim],y10[dim],y11[dim];
   double prd[2],rpd[2],vpd8[2],vpd12[2],vpd12split[2];
   double err0[2],err8sum[2],err12sum[2],err12split[2];
   double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11;
   double s4a,s5a,s6a,s7a,s4b,s5b,s6b,s7b;

   for(int i=0; i<dim; i++)
   {
      random_double_double(&xhi[i], &xlo[i]);
      if(xhi[i] < 0.0) xhi[i] = -xhi[i];
      if(xlo[i] < 0.0) xlo[i] = -xlo[i];
      make_dd_exponent_zero(&xhi[i], &xlo[i]);
      random_double_double(&yhi[i], &ylo[i]);
      if(yhi[i] < 0.0) yhi[i] = -yhi[i];
      if(ylo[i] < 0.0) ylo[i] = -ylo[i];
      make_dd_exponent_zero(&yhi[i], &ylo[i]);
   }
   cout << scientific << setprecision(16);

   cout << "double double vector x :" << endl;
   dd_write_vector(dim, xhi, xlo);
   cout << "double double vector y :" << endl;
   dd_write_vector(dim, yhi, ylo);

   double_double_product(dim, xhi, xlo, yhi, ylo, &prd[0], &prd[1]);
   recursive_dd_product(dim, xhi, xlo, yhi, ylo, &rpd[0], &rpd[1]);

   quarter_dd_vector(dim, xhi, xlo, x0, x1, x2, x3, x4, x5, x6, x7, 1);
   quarter_dd_vector(dim, yhi, ylo, y0, y1, y2, y3, y4, y5, y6, y7, 1);

   vectored_dd_product8sum
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7,
       &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7);

   cout << "the 8 sums :" << endl;
   cout << "s0 : "; write_52double(s0);
   cout << "s1 : "; write_52double(s1);
   cout << "s2 : "; write_52double(s2);
   cout << "s3 : "; write_52double(s3);
   cout << "s4 : "; write_52double(s4);
   cout << "s5 : "; write_52double(s5);
   cout << "s6 : "; write_52double(s6);
   cout << "s7 : "; write_52double(s7);

   to_double_double8sum(s0, s1, s2, s3, s4, s5, s6, s7, &vpd8[0], &vpd8[1]);

   vectored_dd_product12sum
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7,
       &s0, &s1, &s2, &s3, &s4a, &s5a, &s6a, &s7a, &s4b, &s5b, &s6b, &s7b);

   cout << "the 12 sums :" << endl;
   cout << " s0 : "; write_52double(s0);
   cout << " s1 : "; write_52double(s1);
   cout << " s2 : "; write_52double(s2);
   cout << " s3 : "; write_52double(s3);
   cout << "s4a : "; write_52double(s4a);
   cout << "s4b : "; write_52double(s4b);
   cout << "s5a : "; write_52double(s5a);
   cout << "s5b : "; write_52double(s5b);
   cout << "s6a : "; write_52double(s6a);
   cout << "s6b : "; write_52double(s6b);
   cout << "s7a : "; write_52double(s7a);
   cout << "s7b : "; write_52double(s7b);

   to_double_double12sum
      (s0, s1, s2, s3, s4a, s5a, s6a, s7a, s4b, s5b, s6b, s7b, 
       &vpd12[0], &vpd12[1]);

   split_dd_vector
      (dim, xhi, xlo, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, 1);
   split_dd_vector
      (dim, yhi, ylo, y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, 1);

   vectored_dd_product
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11,
            y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11,
       &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &s10, &s11);

   cout << "the 12 sums of the 12-split :" << endl;
   cout << " s0 : "; write_52double(s0);
   cout << " s1 : "; write_52double(s1);
   cout << " s2 : "; write_52double(s2);
   cout << " s3 : "; write_52double(s3);
   cout << " s4 : "; write_52double(s4);
   cout << " s5 : "; write_52double(s5);
   cout << " s6 : "; write_52double(s6);
   cout << " s7 : "; write_52double(s7);
   cout << " s8 : "; write_52double(s8);
   cout << " s9 : "; write_52double(s9);
   cout << "s10 : "; write_52double(s10);
   cout << "s11 : "; write_52double(s11);

   to_double_double12sum
      (s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, 
       &vpd12split[0], &vpd12split[1]);

   cout << " dd x*y : "; dd_write(prd, 32); cout << endl;
   cout << "rdd x*y : "; dd_write(rpd, 32); cout << endl;
   cout << " v8 x*y : "; dd_write(vpd8, 32); cout << endl;
   cout << "v12 x*y : "; dd_write(vpd12, 32); cout << endl;
   cout << "p12 x*y : "; dd_write(vpd12split, 32); cout << endl;

   ddf_sub(rpd[0], rpd[1], prd[0], prd[1], &err0[0], &err0[1]);
   ddf_sub(rpd[0], rpd[1], vpd8[0], vpd8[1], &err8sum[0], &err8sum[1]);
   ddf_sub(rpd[0], rpd[1], vpd12[0], vpd12[1], &err12sum[0], &err12sum[1]);
   ddf_sub(rpd[0], rpd[1], vpd12split[0], vpd12split[1],
           &err12split[0], &err12split[1]);

   if(err0[0] < 0.0) ddf_minus(&err0[0], &err0[1]);
   if(err8sum[0] < 0.0) ddf_minus(&err8sum[0], &err8sum[1]);
   if(err12sum[0] < 0.0) ddf_minus(&err12sum[0], &err12sum[1]);
   if(err12split[0] < 0.0) ddf_minus(&err12split[0], &err12split[1]);

   cout << "   plain error : "; dd_write(err0, 3); cout << endl;
   cout << "   error 8-sum : "; dd_write(err8sum, 3); cout << endl;
   cout << "  error 12-sum : "; dd_write(err12sum, 3); cout << endl;
   cout << "error 12-split : "; dd_write(err12split, 3); cout << endl;

   fail = (int(err8sum[0]) > err0[0]/1000) and
          (int(err12sum[0]) > err0[0]/1000) and
          (int(err12split[0]) > err0[0]/1000);

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

   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
      {
         if(Ahi[i][j] < 0.0) Ahi[i][j] = -Ahi[i][j];
         if(Alo[i][j] < 0.0) Alo[i][j] = -Alo[i][j];
      }

   cout << scientific << setprecision(16);

/*
   cout << "A random " << nrows << "-by-" << nrc << " matrix A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
         cout << "A[" << i << "][" << j << "] : "
              << Ahi[i][j] << "  " << Alo[i][j] << endl;
 */
   double **Bhi = new double*[nrc];
   double **Blo = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhi[i] = new double[ncols];
      Blo[i] = new double[ncols];
   }
   random_dbl2_matrix(nrc, ncols, Bhi, Blo);

   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
      {
         if(Bhi[i][j] < 0.0) Bhi[i][j] = -Bhi[i][j];
         if(Blo[i][j] < 0.0) Blo[i][j] = -Blo[i][j];
      }
/*
   cout << "A random " << nrc << "-by-" << ncols << " matrix B :" << endl;
   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
         cout << "B[" << i << "][" << j << "] : "
              << Bhi[i][j] << "  " << Blo[i][j] << endl;
 */
   double_double_matmatmul(nrows, ncols, nrc, Ahi, Alo, Bhi, Blo, Chi, Clo);
/*
   cout << "the product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "C[" << i << "][" << j << "] : "
              << Chi[i][j] << "  " << Clo[i][j] << endl;
 */
   double **Ahi0 = new double*[nrows];
   double **Ahi1 = new double*[nrows];
   double **Ahi2 = new double*[nrows];
   double **Ahi3 = new double*[nrows];
   double **Alo0 = new double*[nrows];
   double **Alo1 = new double*[nrows];
   double **Alo2 = new double*[nrows];
   double **Alo3 = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahi0[i] = new double[nrc];
      Ahi1[i] = new double[nrc];
      Ahi2[i] = new double[nrc];
      Ahi3[i] = new double[nrc];
      Alo0[i] = new double[nrc];
      Alo1[i] = new double[nrc];
      Alo2[i] = new double[nrc];
      Alo3[i] = new double[nrc];
   }
   quarter_dd_matrix
      (nrows, nrc, Ahi, Alo,
       Ahi0, Ahi1, Ahi2, Ahi3, Alo0, Alo1, Alo2, Alo3);

   double **Bhi0 = new double*[nrc];
   double **Bhi1 = new double*[nrc];
   double **Bhi2 = new double*[nrc];
   double **Bhi3 = new double*[nrc];
   double **Blo0 = new double*[nrc];
   double **Blo1 = new double*[nrc];
   double **Blo2 = new double*[nrc];
   double **Blo3 = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhi0[i] = new double[ncols];
      Bhi1[i] = new double[ncols];
      Bhi2[i] = new double[ncols];
      Bhi3[i] = new double[ncols];
      Blo0[i] = new double[ncols];
      Blo1[i] = new double[ncols];
      Blo2[i] = new double[ncols];
      Blo3[i] = new double[ncols];
   }
   quarter_dd_matrix
      (nrc, ncols, Bhi, Blo,
       Bhi0, Bhi1, Bhi2, Bhi3, Blo0, Blo1, Blo2, Blo3);

   double **Thi0 = new double*[ncols];
   double **Thi1 = new double*[ncols];
   double **Thi2 = new double*[ncols];
   double **Thi3 = new double*[ncols];
   double **Tlo0 = new double*[ncols];
   double **Tlo1 = new double*[ncols];
   double **Tlo2 = new double*[ncols];
   double **Tlo3 = new double*[ncols];

   for(int i=0; i<ncols; i++)
   {
      Thi0[i] = new double[nrc];
      Thi1[i] = new double[nrc];
      Thi2[i] = new double[nrc];
      Thi3[i] = new double[nrc];
      Tlo0[i] = new double[nrc];
      Tlo1[i] = new double[nrc];
      Tlo2[i] = new double[nrc];
      Tlo3[i] = new double[nrc];
   }
   transpose_dd_quarters
      (nrc, ncols, Bhi0, Bhi1, Bhi2, Bhi3, Blo0, Blo1, Blo2, Blo3,
                   Thi0, Thi1, Thi2, Thi3, Tlo0, Tlo1, Tlo2, Tlo3);

   double **Chi0 = new double*[nrows];
   double **Chi1 = new double*[nrows];
   double **Chi2 = new double*[nrows];
   double **Chi3 = new double*[nrows];
   double **Clo0 = new double*[nrows];
   double **Clo1 = new double*[nrows];
   double **Clo2 = new double*[nrows];
   double **Clo3 = new double*[nrows];
   double **Clo0a = new double*[nrows];
   double **Clo1a = new double*[nrows];
   double **Clo2a = new double*[nrows];
   double **Clo3a = new double*[nrows];
   double **Clo0b = new double*[nrows];
   double **Clo1b = new double*[nrows];
   double **Clo2b = new double*[nrows];
   double **Clo3b = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Chi0[i] = new double[ncols];
      Chi1[i] = new double[ncols];
      Chi2[i] = new double[ncols];
      Chi3[i] = new double[ncols];
      Clo0[i] = new double[ncols];
      Clo1[i] = new double[ncols];
      Clo2[i] = new double[ncols];
      Clo3[i] = new double[ncols];
      Clo0a[i] = new double[ncols];
      Clo1a[i] = new double[ncols];
      Clo2a[i] = new double[ncols];
      Clo3a[i] = new double[ncols];
      Clo0b[i] = new double[ncols];
      Clo1b[i] = new double[ncols];
      Clo2b[i] = new double[ncols];
      Clo3b[i] = new double[ncols];
   }
   vectored_dd_matmatmul8sum
      (nrows, ncols, nrc,
       Ahi0, Ahi1, Ahi2, Ahi3, Alo0, Alo1, Alo2, Alo3,
       Thi0, Thi1, Thi2, Thi3, Tlo0, Tlo1, Tlo2, Tlo3,
       Chi0, Chi1, Chi2, Chi3, Clo0, Clo1, Clo2, Clo3);

   double **Vhi = new double*[nrows];
   double **Vlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Vhi[i] = new double[ncols];
      Vlo[i] = new double[ncols];
   }
   to_double_double8sum_matrix
      (nrows, ncols,
       Chi0, Chi1, Chi2, Chi3, Clo0, Clo1, Clo2, Clo3, Vhi, Vlo);

   vectored_dd_matmatmul12sum
      (nrows, ncols, nrc,
       Ahi0, Ahi1, Ahi2, Ahi3, Alo0, Alo1, Alo2, Alo3,
       Thi0, Thi1, Thi2, Thi3, Tlo0, Tlo1, Tlo2, Tlo3,
       Chi0, Chi1, Chi2, Chi3,
       Clo0a, Clo1a, Clo2a, Clo3a, Clo0b, Clo1b, Clo2b, Clo3b);

   double **Whi = new double*[nrows];
   double **Wlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Whi[i] = new double[ncols];
      Wlo[i] = new double[ncols];
   }
   to_double_double12sum_matrix
      (nrows, ncols,
       Chi0, Chi1, Chi2, Chi3,
       Clo0a, Clo1a, Clo2a, Clo3a, Clo0b, Clo1b, Clo2b, Clo3b, Whi, Wlo);
/*
   cout << "the vectored product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "V[" << i << "][" << j << "] : "
              << Vhi[i][j] << "  " << Vlo[i][j] << endl;
 */
   double err8[2],err12[2],acc[2];
   double maxerr8 = 0.0;
   double maxerr12 = 0.0;
   err8[0] = 0.0; err8[1] = 0.0;
   err12[0] = 0.0; err12[1] = 0.0;

   cout << scientific << setprecision(3);

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         ddf_sub(Chi[i][j], Clo[i][j], Vhi[i][j], Vlo[i][j],
                 &acc[0], &acc[1]);
         if(acc[0] < 0.0) ddf_minus(&acc[0], &acc[1]);
         ddf_inc(&err8[0], &err8[1], acc[0], acc[1]);
         if(abs(acc[0]) > maxerr8) 
         {
            maxerr8 = abs(acc[0]);
            cout << "increase max error 8-sum to " << maxerr8
                 << ", at i = " << i << ", j = " << j << endl;
         }
         ddf_sub(Chi[i][j], Clo[i][j], Whi[i][j], Wlo[i][j],
                 &acc[0], &acc[1]);
         if(acc[0] < 0.0) ddf_minus(&acc[0], &acc[1]);
         ddf_inc(&err12[0], &err12[1], acc[0], acc[1]);
         if(abs(acc[0]) > maxerr12) 
         {
            maxerr12 = abs(acc[0]);
            cout << "increase max error 12-sum to " << maxerr12
                 << ", at i = " << i << ", j = " << j << endl;
         }
      }

   cout << "-> error  8-sum : "; dd_write(err8, 3); cout << endl;
   cout << "-> error 12-sum : "; dd_write(err12, 3); cout << endl;
   cout << scientific << setprecision(3);
   cout << "-> max error  8-sum : " << maxerr8 << endl;
   cout << "-> max error 12-sum : " << maxerr12 << endl;

   fail = (maxerr8 > 1.0E-28);

   if(fail == 1) return fail; // no point to continue ...

   return fail; // return anyway ...

   double **cA = new double*[8*nrows];
   for(int i=0; i<8*nrows; i++) cA[i] = new double[8*nrc];

   dd_convolute_quarters
      (nrows, nrc, Ahi0, Ahi1, Ahi2, Ahi3, Alo0, Alo1, Alo2, Alo3, cA);

   cout << "the convoluted quartered matrix A :" << endl;
   for(int i=0; i<8*nrows; i++)
      for(int j=0; j<8*nrc; j++)
         cout << "cA[" << i << "][" << j << "] : " << cA[i][j] << endl;

   double **sB = new double*[8*nrc];
   for(int i=0; i<8*nrc; i++) sB[i] = new double[ncols];

   dd_stack_quarters
      (nrc, ncols, Bhi0, Bhi1, Bhi2, Bhi3, Blo0, Blo1, Blo2, Blo3, sB);

   cout << "the stacked quartered matrix B :" << endl;
   for(int i=0; i<8*nrc; i++)
      for(int j=0; j<ncols; j++)
         cout << "sB[" << i << "][" << j << "] : " << sB[i][j] << endl;

   double **qC = new double*[8*nrows];
   for(int i=0; i<8*nrows; i++) qC[i] = new double[ncols];

   double_indexed_matrix_multiplication(8*nrows, ncols, 8*nrc, cA, sB, qC);

   cout << "the quartered product C :" << endl;
   for(int i=0; i<8*nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "qC[" << i << "][" << j << "] : " << qC[i][j] << endl;

   double **Dhi0 = new double*[nrows];
   double **Dhi1 = new double*[nrows];
   double **Dhi2 = new double*[nrows];
   double **Dhi3 = new double*[nrows];
   double **Dlo0 = new double*[nrows];
   double **Dlo1 = new double*[nrows];
   double **Dlo2 = new double*[nrows];
   double **Dlo3 = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Dhi0[i] = new double[ncols];
      Dhi1[i] = new double[ncols];
      Dhi2[i] = new double[ncols];
      Dhi3[i] = new double[ncols];
      Dlo0[i] = new double[ncols];
      Dlo1[i] = new double[ncols];
      Dlo2[i] = new double[ncols];
      Dlo3[i] = new double[ncols];
   }
   extract_dd_quarters
      (nrows, ncols, qC, Dhi0, Dhi1, Dhi2, Dhi3, Dlo0, Dlo1, Dlo2, Dlo3);

   double error = 0.0;

   cout << "comparing the extracted quarters ..." << endl;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         cout << "Chi0[" << i << "][" << j << "] : " << Chi0[i][j] << endl
              << "Dhi0[" << i << "][" << j << "] : " << Dhi0[i][j] << endl;
         error = error + abs(Chi0[i][j] - Dhi0[i][j]);
         cout << "Chi1[" << i << "][" << j << "] : " << Chi1[i][j] << endl
              << "Dhi1[" << i << "][" << j << "] : " << Dhi1[i][j] << endl;
         error = error + abs(Chi1[i][j] - Dhi1[i][j]);
         cout << "Chi2[" << i << "][" << j << "] : " << Chi2[i][j] << endl
              << "Dhi2[" << i << "][" << j << "] : " << Dhi2[i][j] << endl;
         error = error + abs(Chi2[i][j] - Dhi2[i][j]);
         cout << "Chi3[" << i << "][" << j << "] : " << Chi3[i][j] << endl
              << "Dhi3[" << i << "][" << j << "] : " << Dhi3[i][j] << endl;
         error = error + abs(Chi3[i][j] - Dhi3[i][j]);
         cout << "Clo0[" << i << "][" << j << "] : " << Clo0[i][j] << endl
              << "Dlo0[" << i << "][" << j << "] : " << Dlo0[i][j] << endl;
         error = error + abs(Clo0[i][j] - Dlo0[i][j]);
         cout << "Clo1[" << i << "][" << j << "] : " << Clo1[i][j] << endl
              << "Dlo1[" << i << "][" << j << "] : " << Dlo1[i][j] << endl;
         error = error + abs(Clo1[i][j] - Dlo1[i][j]);
         cout << "Clo2[" << i << "][" << j << "] : " << Clo2[i][j] << endl
              << "Dlo2[" << i << "][" << j << "] : " << Dlo2[i][j] << endl;
         error = error + abs(Clo2[i][j] - Dlo2[i][j]);
         cout << "Clo3[" << i << "][" << j << "] : " << Clo3[i][j] << endl
              << "Dlo3[" << i << "][" << j << "] : " << Dlo3[i][j] << endl;
         error = error + abs(Clo3[i][j] - Dlo3[i][j]);
      }

   cout << scientific << setprecision(3)
        << "sum of all errors : " << error << endl; 

   fail = (error > 1.0E-28);

   if(fail == 1) return fail; // no point to continue ...

   double *cAs = new double[8*nrows*8*nrc]; // single indexed cA
   double2single_row_major(8*nrows, 8*nrc, cA, cAs);

   cout << scientific << setprecision(16);

   error = 0.0;

   cout << "the single indexed convoluted quartered matrix A :" << endl;
   for(int i=0, idx=0; i<8*nrows; i++)
      for(int j=0; j<8*nrc; j++)
      {
         cout << " cA[" << i << "][" << j << "] : " << cA[i][j] << endl
              << "cAs[" << idx << "]    : " << cAs[idx] << endl;
         error = error + abs(cA[i][j] - cAs[idx++]);
      }

   cout << scientific << setprecision(3)
        << "sum of errors : " << error << endl;

   fail = (error > 1.0E-28);

   if(fail == 1) return fail; // no point to continue ...

   double **sBT = new double*[ncols];
   for(int i=0; i<ncols; i++) sBT[i] = new double[8*nrc];
   transpose_rows_columns(8*nrc, ncols, sB, sBT);

   cout << scientific << setprecision(16);

   cout << "Transpose of the stacked matrix sB :" << endl;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<8*nrc; j++)
         cout << "sBT[" << i << "][" << j << "] : " << sBT[i][j] << endl;

   cout << "converting into single indexed matrix ..." << endl;
  
   double *sBs = new double[8*ncols*nrc];
   double2single_column_major(8*nrc, ncols, sBT, sBs);

   double *qCs = new double[8*nrows*ncols];

   cout << "running a single indexed matrix matrix multiplication ..." << endl;

   single_indexed_matrix_multiplication(8*nrows, ncols, 8*nrc, cAs, sBs, qCs);

   double **qC2 = new double*[8*nrows];
   for(int i=0; i<8*nrows; i++) qC2[i] = new double[ncols];

   cout << "converting product to double indexed matrix ..." << endl;

   single2double_row_major(8*nrows, ncols, qCs, qC2);

   cout << scientific << setprecision(16);

   error = 0.0;

   for(int i=0; i<8*nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         cout << " qC[" << i << "][" << j << "] : " << qC[i][j] << endl
              << "qC2[" << i << "][" << j << "] : " << qC2[i][j] << endl;
         error = error + abs(qC[i][j] - qC2[i][j]);
      }

   cout << scientific << setprecision(3)
        << "sum of all errors : " << error << endl; 

   fail = (error > 1.0E-28);

   return fail;
}
