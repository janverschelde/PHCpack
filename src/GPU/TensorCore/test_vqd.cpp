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
 * and compares their inner product with the vectored inner product.
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

   if(fail == 1)
      cout << "\nTest on vectored quad double matmatmul failed?!!!\n\n";
   else
      cout << "\nTest on vectored quad double matmatmul succeeded.\n\n";

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

   cout << "e : "; qd_write(e, 3); cout << endl;

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

   cout << " error : "; qd_write(err, 3); cout << endl;

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

   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
      {
         if(Ahihi[i][j] < 0.0) Ahihi[i][j] = -Ahihi[i][j];
         if(Alohi[i][j] < 0.0) Alohi[i][j] = -Alohi[i][j];
         if(Ahilo[i][j] < 0.0) Ahilo[i][j] = -Ahilo[i][j];
         if(Alolo[i][j] < 0.0) Alolo[i][j] = -Alolo[i][j];
      }

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

   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
      {
         if(Bhihi[i][j] < 0.0) Bhihi[i][j] = -Bhihi[i][j];
         if(Blohi[i][j] < 0.0) Blohi[i][j] = -Blohi[i][j];
         if(Bhilo[i][j] < 0.0) Bhilo[i][j] = -Bhilo[i][j];
         if(Blolo[i][j] < 0.0) Blolo[i][j] = -Blolo[i][j];
      }

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

   double **Ahihi0 = new double*[nrows];
   double **Ahihi1 = new double*[nrows];
   double **Ahihi2 = new double*[nrows];
   double **Ahihi3 = new double*[nrows];
   double **Alohi0 = new double*[nrows];
   double **Alohi1 = new double*[nrows];
   double **Alohi2 = new double*[nrows];
   double **Alohi3 = new double*[nrows];
   double **Ahilo0 = new double*[nrows];
   double **Ahilo1 = new double*[nrows];
   double **Ahilo2 = new double*[nrows];
   double **Ahilo3 = new double*[nrows];
   double **Alolo0 = new double*[nrows];
   double **Alolo1 = new double*[nrows];
   double **Alolo2 = new double*[nrows];
   double **Alolo3 = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihi0[i] = new double[nrc];
      Ahihi1[i] = new double[nrc];
      Ahihi2[i] = new double[nrc];
      Ahihi3[i] = new double[nrc];
      Alohi0[i] = new double[nrc];
      Alohi1[i] = new double[nrc];
      Alohi2[i] = new double[nrc];
      Alohi3[i] = new double[nrc];
      Ahilo0[i] = new double[nrc];
      Ahilo1[i] = new double[nrc];
      Ahilo2[i] = new double[nrc];
      Ahilo3[i] = new double[nrc];
      Alolo0[i] = new double[nrc];
      Alolo1[i] = new double[nrc];
      Alolo2[i] = new double[nrc];
      Alolo3[i] = new double[nrc];
   }
   quarter_qd_matrix
      (nrows, nrc, Ahihi, Alohi, Ahilo, Alolo,
       Ahihi0, Ahihi1, Ahihi2, Ahihi3, Alohi0, Alohi1, Alohi2, Alohi3,
       Ahilo0, Ahilo1, Ahilo2, Ahilo3, Alolo0, Alolo1, Alolo2, Alolo3);

   double **Bhihi0 = new double*[nrc];
   double **Bhihi1 = new double*[nrc];
   double **Bhihi2 = new double*[nrc];
   double **Bhihi3 = new double*[nrc];
   double **Blohi0 = new double*[nrc];
   double **Blohi1 = new double*[nrc];
   double **Blohi2 = new double*[nrc];
   double **Blohi3 = new double*[nrc];
   double **Bhilo0 = new double*[nrc];
   double **Bhilo1 = new double*[nrc];
   double **Bhilo2 = new double*[nrc];
   double **Bhilo3 = new double*[nrc];
   double **Blolo0 = new double*[nrc];
   double **Blolo1 = new double*[nrc];
   double **Blolo2 = new double*[nrc];
   double **Blolo3 = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhihi0[i] = new double[ncols];
      Bhihi1[i] = new double[ncols];
      Bhihi2[i] = new double[ncols];
      Bhihi3[i] = new double[ncols];
      Blohi0[i] = new double[ncols];
      Blohi1[i] = new double[ncols];
      Blohi2[i] = new double[ncols];
      Blohi3[i] = new double[ncols];
      Bhilo0[i] = new double[ncols];
      Bhilo1[i] = new double[ncols];
      Bhilo2[i] = new double[ncols];
      Bhilo3[i] = new double[ncols];
      Blolo0[i] = new double[ncols];
      Blolo1[i] = new double[ncols];
      Blolo2[i] = new double[ncols];
      Blolo3[i] = new double[ncols];
   }
   quarter_qd_matrix
      (nrc, ncols, Bhihi, Blohi, Bhilo, Blolo,
       Bhihi0, Bhihi1, Bhihi2, Bhihi3, Blohi0, Blohi1, Blohi2, Blohi3,
       Bhilo0, Bhilo1, Bhilo2, Bhilo3, Blolo0, Blolo1, Blolo2, Blolo3);

   double **Thihi0 = new double*[ncols];
   double **Thihi1 = new double*[ncols];
   double **Thihi2 = new double*[ncols];
   double **Thihi3 = new double*[ncols];
   double **Tlohi0 = new double*[ncols];
   double **Tlohi1 = new double*[ncols];
   double **Tlohi2 = new double*[ncols];
   double **Tlohi3 = new double*[ncols];
   double **Thilo0 = new double*[ncols];
   double **Thilo1 = new double*[ncols];
   double **Thilo2 = new double*[ncols];
   double **Thilo3 = new double*[ncols];
   double **Tlolo0 = new double*[ncols];
   double **Tlolo1 = new double*[ncols];
   double **Tlolo2 = new double*[ncols];
   double **Tlolo3 = new double*[ncols];

   for(int i=0; i<ncols; i++)
   {
      Thihi0[i] = new double[nrc];
      Thihi1[i] = new double[nrc];
      Thihi2[i] = new double[nrc];
      Thihi3[i] = new double[nrc];
      Tlohi0[i] = new double[nrc];
      Tlohi1[i] = new double[nrc];
      Tlohi2[i] = new double[nrc];
      Tlohi3[i] = new double[nrc];
      Thilo0[i] = new double[nrc];
      Thilo1[i] = new double[nrc];
      Thilo2[i] = new double[nrc];
      Thilo3[i] = new double[nrc];
      Tlolo0[i] = new double[nrc];
      Tlolo1[i] = new double[nrc];
      Tlolo2[i] = new double[nrc];
      Tlolo3[i] = new double[nrc];
   }
   transpose_qd_quarters
      (nrc, ncols,
       Bhihi0, Bhihi1, Bhihi2, Bhihi3, Blohi0, Blohi1, Blohi2, Blohi3,
       Bhilo0, Bhilo1, Bhilo2, Bhilo3, Blolo0, Blolo1, Blolo2, Blolo3,
       Thihi0, Thihi1, Thihi2, Thihi3, Tlohi0, Tlohi1, Tlohi2, Tlohi3,
       Thilo0, Thilo1, Thilo2, Thilo3, Tlolo0, Tlolo1, Tlolo2, Tlolo3);

   double **Chihi0 = new double*[nrows];
   double **Chihi1 = new double*[nrows];
   double **Chihi2 = new double*[nrows];
   double **Chihi3 = new double*[nrows];
   double **Clohi0 = new double*[nrows];
   double **Clohi1 = new double*[nrows];
   double **Clohi2 = new double*[nrows];
   double **Clohi3 = new double*[nrows];
   double **Chilo0 = new double*[nrows];
   double **Chilo1 = new double*[nrows];
   double **Chilo2 = new double*[nrows];
   double **Chilo3 = new double*[nrows];
   double **Clolo0 = new double*[nrows];
   double **Clolo1 = new double*[nrows];
   double **Clolo2 = new double*[nrows];
   double **Clolo3 = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Chihi0[i] = new double[ncols];
      Chihi1[i] = new double[ncols];
      Chihi2[i] = new double[ncols];
      Chihi3[i] = new double[ncols];
      Clohi0[i] = new double[ncols];
      Clohi1[i] = new double[ncols];
      Clohi2[i] = new double[ncols];
      Clohi3[i] = new double[ncols];
      Chilo0[i] = new double[ncols];
      Chilo1[i] = new double[ncols];
      Chilo2[i] = new double[ncols];
      Chilo3[i] = new double[ncols];
      Clolo0[i] = new double[ncols];
      Clolo1[i] = new double[ncols];
      Clolo2[i] = new double[ncols];
      Clolo3[i] = new double[ncols];
   }
   vectored_qd_matmatmul
      (nrows, ncols, nrc,
       Ahihi0, Ahihi1, Ahihi2, Ahihi3, Alohi0, Alohi1, Alohi2, Alohi3,
       Ahilo0, Ahilo1, Ahilo2, Ahilo3, Alolo0, Alolo1, Alolo2, Alolo3,
       Thihi0, Thihi1, Thihi2, Thihi3, Tlohi0, Tlohi1, Tlohi2, Tlohi3,
       Thilo0, Thilo1, Thilo2, Thilo3, Tlolo0, Tlolo1, Tlolo2, Tlolo3,
       Chihi0, Chihi1, Chihi2, Chihi3, Clohi0, Clohi1, Clohi2, Clohi3,
       Chilo0, Chilo1, Chilo2, Chilo3, Clolo0, Clolo1, Clolo2, Clolo3);

   double **Vhihi = new double*[nrows];
   double **Vlohi = new double*[nrows];
   double **Vhilo = new double*[nrows];
   double **Vlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Vhihi[i] = new double[ncols];
      Vlohi[i] = new double[ncols];
      Vhilo[i] = new double[ncols];
      Vlolo[i] = new double[ncols];
   }
   to_quad_double_matrix
      (nrows, ncols,
       Chihi0, Chihi1, Chihi2, Chihi3, Clohi0, Clohi1, Clohi2, Clohi3,
       Chilo0, Chilo1, Chilo2, Chilo3, Clolo0, Clolo1, Clolo2, Clolo3,
       Vhihi, Vlohi, Vhilo, Vlolo);

   cout << "the vectored product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "V[" << i << "][" << j << "] : "
              << Vhihi[i][j] << "  " << Vlohi[i][j] << endl
              << "          "
              << Vhilo[i][j] << "  " << Vlolo[i][j] << endl;

   double err[4],acc[4];
   err[0] = 0.0; err[1] = 0.0;
   err[2] = 0.0; err[3] = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         qdf_sub(Chihi[i][j], Clohi[i][j], Chilo[i][j], Clolo[i][j],
                 Vhihi[i][j], Vlohi[i][j], Vhilo[i][j], Vlolo[i][j],
                 &acc[0], &acc[1], &acc[2], &acc[3]);
         qdf_inc(&err[0], &err[1], &err[2], &err[3],
                 acc[0], acc[1], acc[2], acc[3]);
      }

   if(err[0] < 0.0) qdf_minus(&err[0], &err[1], &err[2], &err[3]);

   cout << "-> error : "; qd_write(err, 3); cout << endl;

   fail = (abs(err[0]) > 1.0E-58);

   return fail;
}
