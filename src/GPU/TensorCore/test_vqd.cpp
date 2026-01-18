/* Tests the collection of functions on vectored quad double arithmetic.  */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "quad_double.h"
#include "random4_vectors.h"
#include "random4_matrices.h"
#include "quad_double_functions.h"
#include "vectored_quad_doubles.h"

int test_quarter_quad_double ( void );
/*
 * Generates a random positive quad double, quarters the parts,
 * and then checks if their sum equals the original quad double.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_qd_product ( int dim );
/*
 * Generates two random vectors of quad doubles of size dim,
 * and compares their inner producted with the vectored inner product.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_qd_matmatmul ( int nrows, int ncols, int nrc ); 
/*
 * Generates two random quad double matrices of dimension
 * nrows-by-nrc for A, nrc-by-ncols for B, and then tests
 * the matrix matrix multiplication of A with B. */

using namespace std;

int main ( void )
{
   int fail;

   srand(time(NULL));

   fail = test_quarter_quad_double();

   if(fail == 1)
      cout << "\nTest on quarter quad double failed?!!!\n\n";
   else
      cout << "\nTest on quarter quad double succeeded.\n\n";

   cout << "Give the dimension : ";
   int dim; cin >> dim;

   fail = test_vectored_qd_product(dim);

   if(fail == 1)
      cout << "\nTest on vectored quad double product failed?!!!\n\n";
   else
      cout << "\nTest on vectored quad double product succeeded.\n\n";

   cout << "Give #rows of the product A*B : ";
   int m; cin >> m;
   cout << "Give #columns of the product A*B : ";
   int n; cin >> n;
   cout << "Give #columns of A, #rows of B : ";
   int k; cin >> k;

   fail = test_vectored_qd_matmatmul(m, n, k);

   return 0;
}

int test_quarter_quad_double ( void )
{
   int fail = 0;

   double x[4];
   double xhihi0,xhihi1,xhihi2,xhihi3;
   double xlohi0,xlohi1,xlohi2,xlohi3;
   double xhilo0,xhilo1,xhilo2,xhilo3;
   double xlolo0,xlolo1,xlolo2,xlolo3;
   double y[4];
   double e[4];

   random_quad_double(&x[0], &x[1], &x[2], &x[3]);

   if(x[0] < 0.0) x[0] = -x[0]; // all parts must be positive
   if(x[1] < 0.0) x[1] = -x[1];
   if(x[2] < 0.0) x[2] = -x[2];
   if(x[3] < 0.0) x[3] = -x[3];

   cout << scientific << setprecision(16);

   cout << "x : "; qd_write(x, 64); cout << endl;

   quarter_quad_double
      (x[0], x[1], x[2], x[3],
       &xhihi0, &xhihi1, &xhihi2, &xhihi3,
       &xlohi0, &xlohi1, &xlohi2, &xlohi3,
       &xhilo0, &xhilo1, &xhilo2, &xhilo3,
       &xlolo0, &xlolo1, &xlolo2, &xlolo3);

   cout << "xhihi0 : " << xhihi0 << endl;
   cout << "xhihi1 : " << xhihi1 << endl;
   cout << "xhihi2 : " << xhihi2 << endl;
   cout << "xhihi3 : " << xhihi3 << endl;
   cout << "xlohi0 : " << xlohi0 << endl;
   cout << "xlohi1 : " << xlohi1 << endl;
   cout << "xlohi2 : " << xlohi2 << endl;
   cout << "xlohi3 : " << xlohi3 << endl;
   cout << "xhilo0 : " << xhilo0 << endl;
   cout << "xhilo1 : " << xhilo1 << endl;
   cout << "xhilo2 : " << xhilo2 << endl;
   cout << "xhilo3 : " << xhilo3 << endl;
   cout << "xlolo0 : " << xlolo0 << endl;
   cout << "xlolo1 : " << xlolo1 << endl;
   cout << "xlolo2 : " << xlolo2 << endl;
   cout << "xlolo3 : " << xlolo3 << endl;

   to_quad_double
      (xhihi0, xhihi1, xhihi2, xhihi3, xlohi0, xlohi1, xlohi2, xlohi3,
       xhilo0, xhilo1, xhilo2, xhilo3, xlolo0, xlolo1, xlolo2, xlolo3,
       &y[0], &y[1], &y[2], &y[3]);

   cout << "y : "; qd_write(y, 64); cout << endl;

   qdf_sub(x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3],
           &e[0], &e[1], &e[2], &e[3]);

   cout << "e : "; qd_write(e, 64); cout << endl;

   fail = not(e[0] == 0.0)
        + not(e[1] == 0.0)
        + not(e[2] == 0.0)
        + not(e[3] == 0.0);

   return fail;
}

int test_vectored_qd_product ( int dim )
{
   int fail = 0;

   double xhihi[dim],xlohi[dim],xhilo[dim],xlolo[dim];
   double yhihi[dim],ylohi[dim],yhilo[dim],ylolo[dim];
   double x0[dim],x1[dim],x2[dim],x3[dim],x4[dim],x5[dim],x6[dim],x7[dim];
   double x8[dim],x9[dim],xA[dim],xB[dim],xC[dim],xD[dim],xE[dim],xF[dim];
   double y0[dim],y1[dim],y2[dim],y3[dim],y4[dim],y5[dim],y6[dim],y7[dim];
   double y8[dim],y9[dim],yA[dim],yB[dim],yC[dim],yD[dim],yE[dim],yF[dim];
   double prd[4],vpd[4],err[4];
   double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF;

   for(int i=0; i<dim; i++)
   {
      random_quad_double(&xhihi[i], &xlohi[i], &xhilo[i], &xlolo[i]);
      if(xhihi[i] < 0.0) xhihi[i] = -xhihi[i];
      if(xlohi[i] < 0.0) xlohi[i] = -xlohi[i];
      if(xhilo[i] < 0.0) xhilo[i] = -xhilo[i];
      if(xlolo[i] < 0.0) xlolo[i] = -xlolo[i];
      random_quad_double(&yhihi[i], &ylohi[i], &yhilo[i], &ylolo[i]);
      if(yhihi[i] < 0.0) yhihi[i] = -yhihi[i];
      if(ylohi[i] < 0.0) ylohi[i] = -ylohi[i];
      if(yhilo[i] < 0.0) yhilo[i] = -yhilo[i];
      if(ylolo[i] < 0.0) ylolo[i] = -ylolo[i];
   }
   cout << scientific << setprecision(16);

   cout << "quad double vector x :"; cout << endl;
   qd_write_vector(dim, xhihi, xlohi, xhilo, xlolo);
   cout << "quad double vector y :"; cout << endl;
   qd_write_vector(dim, yhihi, ylohi, yhilo, ylolo);

   quad_double_product
      (dim, xhihi, xlohi, xhilo, xlolo, yhihi, ylohi, yhilo, ylolo,
       &prd[0], &prd[1], &prd[2], &prd[3]);

   quarter_qd_vector
      (dim, xhihi, xlohi, xhilo, xlolo,
       x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xA, xB, xC, xD, xE, xF);
   quarter_qd_vector
      (dim, yhihi, ylohi, yhilo, ylolo,
       y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, yA, yB, yC, yD, yE, yF);

   vectored_qd_product
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xA, xB, xC, xD, xE, xF,
            y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, yA, yB, yC, yD, yE, yF,
       &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7,
       &s8, &s9, &sA, &sB, &sC, &sD, &sE, &sF);

   to_quad_double
      (s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF,
       &vpd[0], &vpd[1], &vpd[2], &vpd[3]);

   cout << "qd x*y : "; qd_write(prd, 64); cout << endl;
   cout << "vd x*y : "; qd_write(vpd, 64); cout << endl;

   qdf_sub(prd[0], prd[1], prd[2], prd[3], vpd[0], vpd[1], vpd[2], vpd[3],
           &err[0], &err[1], &err[2], &err[3]);

   if(err[0] < 0.0) qdf_minus(&err[0], &err[1], &err[2], &err[3]);

   cout << " error : "; qd_write(err, 64); cout << endl;

   fail = (fabs(err[0]) > 1.0E-58);

   return fail;
}

int test_vectored_qd_matmatmul ( int nrows, int ncols, int nrc )
{
   int fail = 0;

   double **Chihi = new double*[nrows];
   double **Clohi = new double*[nrows];
   double **Chilo = new double*[nrows];
   double **Clolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Chihi[i] = new double[ncols];
      Clohi[i] = new double[ncols];
      Chilo[i] = new double[ncols];
      Clolo[i] = new double[ncols];
   }
   double **Ahihi = new double*[nrows];
   double **Alohi = new double*[nrows];
   double **Ahilo = new double*[nrows];
   double **Alolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihi[i] = new double[nrc];
      Alohi[i] = new double[nrc];
      Ahilo[i] = new double[nrc];
      Alolo[i] = new double[nrc];
   }
   random_dbl4_matrix(nrows, nrc, Ahihi, Alohi, Ahilo, Alolo);

   cout << scientific << setprecision(16);

   cout << "A random " << nrows << "-by-" << nrc << " matrix A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
         cout << "A[" << i << "][" << j << "] : "
              << Ahihi[i][j] << "  " << Alohi[i][j] << endl
              << "          "
              << Ahilo[i][j] << "  " << Alolo[i][j] << endl;

   double **Bhihi = new double*[nrc];
   double **Blohi = new double*[nrc];
   double **Bhilo = new double*[nrc];
   double **Blolo = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhihi[i] = new double[ncols];
      Blohi[i] = new double[ncols];
      Bhilo[i] = new double[ncols];
      Blolo[i] = new double[ncols];
   }
   random_dbl4_matrix(nrc, ncols, Bhihi, Blohi, Bhilo, Blolo);

   cout << "A random " << nrc << "-by-" << ncols << " matrix B :" << endl;
   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
         cout << "B[" << i << "][" << j << "] : "
              << Bhihi[i][j] << "  " << Blohi[i][j] << endl
              << "          "
              << Bhilo[i][j] << "  " << Blolo[i][j] << endl;

   quad_double_matmatmul
      (nrows, ncols, nrc, Ahihi, Alohi, Ahilo, Alolo,
       Bhihi, Blohi, Bhilo, Blolo, Chihi, Clohi, Chilo, Clolo);

   cout << "the product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "C[" << i << "][" << j << "] : "
              << Chihi[i][j] << "  " << Clohi[i][j] << endl
              << "          "
              << Chilo[i][j] << "  " << Clolo[i][j] << endl;

   return fail;
}
