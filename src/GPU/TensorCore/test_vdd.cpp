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

   if(fail == 1)
      cout << "\nTest on vectored double double matmatmul failed?!!!\n\n";
   else
      cout << "\nTest on vectored double double matmatmul succeeded.\n\n";

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

   cout << "e : "; dd_write(e, 3); cout << endl;

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

   cout << " error : "; dd_write(err, 3); cout << endl;

   fail = (abs(err[0]) > 1.0E-28);

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

   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
      {
         if(Bhi[i][j] < 0.0) Bhi[i][j] = -Bhi[i][j];
         if(Blo[i][j] < 0.0) Blo[i][j] = -Blo[i][j];
      }

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
   }
   vectored_dd_matmatmul
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
   to_double_double_matrix
      (nrows, ncols,
       Chi0, Chi1, Chi2, Chi3, Clo0, Clo1, Clo2, Clo3, Vhi, Vlo);

   cout << "the vectored product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "V[" << i << "][" << j << "] : "
              << Vhi[i][j] << "  " << Vlo[i][j] << endl;

   double err[2],acc[2];
   err[0] = 0.0; err[1] = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         ddf_sub(Chi[i][j], Clo[i][j], Vhi[i][j], Vlo[i][j],
                 &acc[0], &acc[1]);
         ddf_inc(&err[0], &err[1], acc[0], acc[1]);
      }

   if(err[0] < 0.0) ddf_minus(&err[0], &err[1]);

   cout << "-> error : "; dd_write(err, 3); cout << endl;

   fail = (abs(err[0]) > 1.0E-28);

   if(fail == 1) return fail; // no point to continue ...

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

   return fail;
}
