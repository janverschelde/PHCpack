/* Tests the collection of functions on vectored octo double arithmetic.  */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "octo_double.h"
#include "random8_vectors.h"
#include "random8_matrices.h"
#include "octo_double_functions.h"
#include "vectored_octo_doubles.h"

int test_quarter_octo_double ( void );
/*
 * Generates a random positive octo double, quarters the parts,
 * and then checks if their sum equals the original octo double.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_od_product ( int dim );
/*
 * Generates two random vectors of octo doubles of size dim,
 * and compares their inner product with the vectored inner product.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_od_matmatmul ( int nrows, int ncols, int nrc ); 
/*
 * Generates two random octo double matrices of dimension
 * nrows-by-nrc for A, nrc-by-ncols for B, and then tests
 * the matrix matrix multiplication of A with B. */

using namespace std;

int main ( void )
{
   int fail;

   srand(time(NULL));

   fail = test_quarter_octo_double();

   if(fail == 1)
      cout << "\nTest on quarter octo double failed?!!!\n\n";
   else
      cout << "\nTest on quarter octo double succeeded.\n\n";

   cout << "Give the dimension : ";
   int dim; cin >> dim;

   fail = test_vectored_od_product(dim);

   if(fail == 1)
      cout << "\nTest on vectored octo double product failed?!!!\n\n";
   else
      cout << "\nTest on vectored octo double product succeeded.\n\n";

   cout << "Give #rows of the product A*B : ";
   int m; cin >> m;
   cout << "Give #columns of the product A*B : ";
   int n; cin >> n;
   cout << "Give #columns of A, #rows of B : ";
   int k; cin >> k;

   fail = test_vectored_od_matmatmul(m, n, k);

   return 0;
}

int test_quarter_octo_double ( void )
{
   int fail = 0;

   double x[8];
   double xhihihi0,xhihihi1,xhihihi2,xhihihi3;
   double xlohihi0,xlohihi1,xlohihi2,xlohihi3;
   double xhilohi0,xhilohi1,xhilohi2,xhilohi3;
   double xlolohi0,xlolohi1,xlolohi2,xlolohi3;
   double xhihilo0,xhihilo1,xhihilo2,xhihilo3;
   double xlohilo0,xlohilo1,xlohilo2,xlohilo3;
   double xhilolo0,xhilolo1,xhilolo2,xhilolo3;
   double xlololo0,xlololo1,xlololo2,xlololo3;
   double y[8];
   double e[8];

   random_octo_double
      (&x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &x[6], &x[7]);

   if(x[0] < 0.0) x[0] = -x[0]; // all parts must be positive
   if(x[1] < 0.0) x[1] = -x[1];
   if(x[2] < 0.0) x[2] = -x[2];
   if(x[3] < 0.0) x[3] = -x[3];
   if(x[4] < 0.0) x[4] = -x[4];
   if(x[5] < 0.0) x[5] = -x[5];
   if(x[6] < 0.0) x[6] = -x[6];
   if(x[7] < 0.0) x[7] = -x[7];

   cout << scientific << setprecision(16);

   cout << "x :" << endl; od_write_doubles(x); cout << endl;

   quarter_octo_double
      (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
       &xhihihi0, &xhihihi1, &xhihihi2, &xhihihi3,
       &xlohihi0, &xlohihi1, &xlohihi2, &xlohihi3,
       &xhilohi0, &xhilohi1, &xhilohi2, &xhilohi3,
       &xlolohi0, &xlolohi1, &xlolohi2, &xlolohi3,
       &xhihilo0, &xhihilo1, &xhihilo2, &xhihilo3,
       &xlohilo0, &xlohilo1, &xlohilo2, &xlohilo3,
       &xhilolo0, &xhilolo1, &xhilolo2, &xhilolo3,
       &xlololo0, &xlololo1, &xlololo2, &xlololo3);

   cout << "xhihihi0 : " << xhihihi0 << endl;
   cout << "xhihihi1 : " << xhihihi1 << endl;
   cout << "xhihihi2 : " << xhihihi2 << endl;
   cout << "xhihihi3 : " << xhihihi3 << endl;
   cout << "xlohihi0 : " << xlohihi0 << endl;
   cout << "xlohihi1 : " << xlohihi1 << endl;
   cout << "xlohihi2 : " << xlohihi2 << endl;
   cout << "xlohihi3 : " << xlohihi3 << endl;
   cout << "xhilohi0 : " << xhilohi0 << endl;
   cout << "xhilohi1 : " << xhilohi1 << endl;
   cout << "xhilohi2 : " << xhilohi2 << endl;
   cout << "xhilohi3 : " << xhilohi3 << endl;
   cout << "xlolohi0 : " << xlolohi0 << endl;
   cout << "xlolohi1 : " << xlolohi1 << endl;
   cout << "xlolohi2 : " << xlolohi2 << endl;
   cout << "xlolohi3 : " << xlolohi3 << endl;
   cout << "xhihilo0 : " << xhihilo0 << endl;
   cout << "xhihilo1 : " << xhihilo1 << endl;
   cout << "xhihilo2 : " << xhihilo2 << endl;
   cout << "xhihilo3 : " << xhihilo3 << endl;
   cout << "xlohilo0 : " << xlohilo0 << endl;
   cout << "xlohilo1 : " << xlohilo1 << endl;
   cout << "xlohilo2 : " << xlohilo2 << endl;
   cout << "xlohilo3 : " << xlohilo3 << endl;
   cout << "xhilolo0 : " << xhilolo0 << endl;
   cout << "xhilolo1 : " << xhilolo1 << endl;
   cout << "xhilolo2 : " << xhilolo2 << endl;
   cout << "xhilolo3 : " << xhilolo3 << endl;
   cout << "xlololo0 : " << xlololo0 << endl;
   cout << "xlololo1 : " << xlololo1 << endl;
   cout << "xlololo2 : " << xlololo2 << endl;
   cout << "xlololo3 : " << xlololo3 << endl;

   to_octo_double
      (xhihihi0, xhihihi1, xhihihi2, xhihihi3,
       xlohihi0, xlohihi1, xlohihi2, xlohihi3,
       xhilohi0, xhilohi1, xhilohi2, xhilohi3,
       xlolohi0, xlolohi1, xlolohi2, xlolohi3,
       xhihilo0, xhihilo1, xhihilo2, xhihilo3,
       xlohilo0, xlohilo1, xlohilo2, xlohilo3,
       xhilolo0, xhilolo1, xhilolo2, xhilolo3,
       xlololo0, xlololo1, xlololo2, xlololo3,
       &y[0], &y[1], &y[2], &y[3], &y[4], &y[5], &y[6], &y[7]);

   cout << "x :"; cout << endl; od_write_doubles(x); cout << endl;
   cout << "y :"; cout << endl; od_write_doubles(y); cout << endl;

   odf_sub(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
           y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
           &e[0], &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7]);

   cout << "e :"; cout << endl; od_write_doubles(e); cout << endl;

   fail = not(e[0] == 0.0) + not(e[1] == 0.0)
        + not(e[2] == 0.0) + not(e[3] == 0.0)
        + not(e[4] == 0.0) + not(e[5] == 0.0)
        + not(e[6] == 0.0) + not(e[7] == 0.0);

   return fail;
}

int test_vectored_od_product ( int dim )
{
   int fail = 0;

   double xhihihi[dim],xlohihi[dim],xhilohi[dim],xlolohi[dim];
   double xhihilo[dim],xlohilo[dim],xhilolo[dim],xlololo[dim];
   double yhihihi[dim],ylohihi[dim],yhilohi[dim],ylolohi[dim];
   double yhihilo[dim],ylohilo[dim],yhilolo[dim],ylololo[dim];

   double x0[dim],x1[dim],x2[dim],x3[dim],x4[dim],x5[dim],x6[dim],x7[dim];
   double x8[dim],x9[dim],x10[dim],x11[dim],x12[dim],x13[dim],x14[dim];
   double x15[dim],x16[dim],x17[dim],x18[dim],x19[dim],x20[dim],x21[dim];
   double x22[dim],x23[dim],x24[dim],x25[dim],x26[dim],x27[dim],x28[dim];
   double x29[dim],x30[dim],x31[dim];
   double y0[dim],y1[dim],y2[dim],y3[dim],y4[dim],y5[dim],y6[dim],y7[dim];
   double y8[dim],y9[dim],y10[dim],y11[dim],y12[dim],y13[dim],y14[dim];
   double y15[dim],y16[dim],y17[dim],y18[dim],y19[dim],y20[dim],y21[dim];
   double y22[dim],y23[dim],y24[dim],y25[dim],y26[dim],y27[dim],y28[dim];
   double y29[dim],y30[dim],y31[dim];

   double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;
   double s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31;
   double prd[8],vpd[8],err[8];

   for(int i=0; i<dim; i++)
   {
      random_octo_double
         (&xhihihi[i], &xlohihi[i], &xhilohi[i], &xlolohi[i],
          &xhihilo[i], &xlohilo[i], &xhilolo[i], &xlololo[i]);

      if(xhihihi[i] < 0.0) xhihihi[i] = -xhihihi[i];
      if(xlohihi[i] < 0.0) xlohihi[i] = -xlohihi[i];
      if(xhilohi[i] < 0.0) xhilohi[i] = -xhilohi[i];
      if(xlolohi[i] < 0.0) xlolohi[i] = -xlolohi[i];
      if(xhihilo[i] < 0.0) xhihilo[i] = -xhihilo[i];
      if(xlohilo[i] < 0.0) xlohilo[i] = -xlohilo[i];
      if(xhilolo[i] < 0.0) xhilolo[i] = -xhilolo[i];
      if(xlololo[i] < 0.0) xlololo[i] = -xlololo[i];

      random_octo_double
         (&yhihihi[i], &ylohihi[i], &yhilohi[i], &ylolohi[i],
          &yhihilo[i], &ylohilo[i], &yhilolo[i], &ylololo[i]);

      if(yhihihi[i] < 0.0) yhihihi[i] = -yhihihi[i];
      if(ylohihi[i] < 0.0) ylohihi[i] = -ylohihi[i];
      if(yhilohi[i] < 0.0) yhilohi[i] = -yhilohi[i];
      if(ylolohi[i] < 0.0) ylolohi[i] = -ylolohi[i];
      if(yhihilo[i] < 0.0) yhihilo[i] = -yhihilo[i];
      if(ylohilo[i] < 0.0) ylohilo[i] = -ylohilo[i];
      if(yhilolo[i] < 0.0) yhilolo[i] = -yhilolo[i];
      if(ylololo[i] < 0.0) ylololo[i] = -ylololo[i];
   }
   cout << "octo double vector x :"; cout << endl;
   od_write_vector
      (dim, xhihihi, xlohihi, xhilohi, xlolohi,
            xhihilo, xlohilo, xhilolo, xlololo);
   cout << "octo double vector y :"; cout << endl;
   od_write_vector
      (dim, yhihihi, ylohihi, yhilohi, ylolohi,
            yhihilo, ylohilo, yhilolo, ylololo);

   octo_double_product
      (dim,
       xhihihi, xlohihi, xhilohi, xlolohi, xhihilo, xlohilo, xhilolo, xlololo,
       yhihihi, ylohihi, yhilohi, ylolohi, yhihilo, ylohilo, yhilolo, ylololo,
       &prd[0], &prd[1], &prd[2], &prd[3], &prd[4], &prd[5], &prd[6], &prd[7]);

   quarter_od_vector
      (dim,
       xhihihi, xlohihi, xhilohi, xlolohi, xhihilo, xlohilo, xhilolo, xlololo,
       x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15,
       x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29,
       x30, x31);
   quarter_od_vector
      (dim,
       yhihihi, ylohihi, yhilohi, ylolohi, yhihilo, ylohilo, yhilolo, ylololo,
       y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15,
       y16, y17, y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29,
       y30, y31);

   vectored_od_product
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,
       x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28,
       x29, x30, x31, y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12,
       y13, y14, y15, y16, y17, y18, y19, y20, y21, y22, y23, y24, y25, y26,
       y27, y28, y29, y30, y31, &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8,
       &s9, &s10, &s11, &s12, &s13, &s14, &s15, &s16, &s17, &s18, &s19, &s20,
       &s21, &s22, &s23, &s24, &s25, &s26, &s27, &s28, &s29, &s30, &s31);

   to_octo_double
      (s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15,
       s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29,
       s30, s31, &vpd[0], &vpd[1], &vpd[2], &vpd[3], &vpd[4], &vpd[5],
       &vpd[6], &vpd[7]);

   cout << "od x*y :"; cout << endl; od_write_doubles(prd); cout << endl;
   cout << "vd x*y :"; cout << endl; od_write_doubles(vpd); cout << endl;

   odf_sub(prd[0], prd[1], prd[2], prd[3], prd[4], prd[5], prd[6], prd[7],
           vpd[0], vpd[1], vpd[2], vpd[3], vpd[4], vpd[5], vpd[6], vpd[7],
           &err[0], &err[1], &err[2], &err[3],
           &err[4], &err[5], &err[6], &err[7]);

   if(err[0] < 0.0)
      odf_minus(&err[0], &err[1], &err[2], &err[3],
                &err[4], &err[5], &err[6], &err[7]);

   cout << " error :"; cout << endl; od_write_doubles(err); cout << endl;

   fail = (abs(err[0]) > 1.0E-120);

   return fail;
}

int test_vectored_od_matmatmul ( int nrows, int ncols, int nrc )
{
   int fail = 0;

   double **Chihihi = new double*[nrows];
   double **Clohihi = new double*[nrows];
   double **Chilohi = new double*[nrows];
   double **Clolohi = new double*[nrows];
   double **Chihilo = new double*[nrows];
   double **Clohilo = new double*[nrows];
   double **Chilolo = new double*[nrows];
   double **Clololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Chihihi[i] = new double[ncols];
      Clohihi[i] = new double[ncols];
      Chilohi[i] = new double[ncols];
      Clolohi[i] = new double[ncols];
      Chihilo[i] = new double[ncols];
      Clohilo[i] = new double[ncols];
      Chilolo[i] = new double[ncols];
      Clololo[i] = new double[ncols];
   }
   double **Ahihihi = new double*[nrows];
   double **Alohihi = new double*[nrows];
   double **Ahilohi = new double*[nrows];
   double **Alolohi = new double*[nrows];
   double **Ahihilo = new double*[nrows];
   double **Alohilo = new double*[nrows];
   double **Ahilolo = new double*[nrows];
   double **Alololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihihi[i] = new double[nrc];
      Alohihi[i] = new double[nrc];
      Ahilohi[i] = new double[nrc];
      Alolohi[i] = new double[nrc];
      Ahihilo[i] = new double[nrc];
      Alohilo[i] = new double[nrc];
      Ahilolo[i] = new double[nrc];
      Alololo[i] = new double[nrc];
   }
   random_dbl8_matrix
      (nrows, nrc,
       Ahihihi, Alohihi, Ahilohi, Alolohi,Ahihilo, Alohilo, Ahilolo, Alololo);

   cout << scientific << setprecision(16);

   cout << "A random " << nrows << "-by-" << nrc << " matrix A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
         cout << "A[" << i << "][" << j << "] : "
              << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
              << "          "
              << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
              << "          "
              << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
              << "          "
              << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;

   double **Bhihihi = new double*[nrc];
   double **Blohihi = new double*[nrc];
   double **Bhilohi = new double*[nrc];
   double **Blolohi = new double*[nrc];
   double **Bhihilo = new double*[nrc];
   double **Blohilo = new double*[nrc];
   double **Bhilolo = new double*[nrc];
   double **Blololo = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhihihi[i] = new double[ncols];
      Blohihi[i] = new double[ncols];
      Bhilohi[i] = new double[ncols];
      Blolohi[i] = new double[ncols];
      Bhihilo[i] = new double[ncols];
      Blohilo[i] = new double[ncols];
      Bhilolo[i] = new double[ncols];
      Blololo[i] = new double[ncols];
   }
   random_dbl8_matrix
      (nrc, ncols,
       Bhihihi, Blohihi, Bhilohi, Blolohi, Bhihilo, Blohilo, Bhilolo, Blololo);

   cout << "A random " << nrc << "-by-" << ncols << " matrix B :" << endl;
   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
         cout << "B[" << i << "][" << j << "] : "
              << Bhihihi[i][j] << "  " << Blohihi[i][j] << endl
              << "          "
              << Bhilohi[i][j] << "  " << Blolohi[i][j] << endl
              << "          "
              << Bhihilo[i][j] << "  " << Blohilo[i][j] << endl
              << "          "
              << Bhilolo[i][j] << "  " << Blololo[i][j] << endl;

   octo_double_matmatmul
      (nrows, ncols, nrc,
       Ahihihi, Alohihi, Ahilohi, Alolohi, Ahihilo, Alohilo, Ahilolo, Alololo,
       Bhihihi, Blohihi, Bhilohi, Blolohi, Bhihilo, Blohilo, Bhilolo, Blololo,
       Chihihi, Clohihi, Chilohi, Clolohi, Chihilo, Clohilo, Chilolo, Clololo);

   cout << "the product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "C[" << i << "][" << j << "] : "
              << Chihihi[i][j] << "  " << Clohihi[i][j] << endl
              << "          "
              << Chilohi[i][j] << "  " << Clolohi[i][j] << endl
              << "          "
              << Chihilo[i][j] << "  " << Clohilo[i][j] << endl
              << "          "
              << Chilolo[i][j] << "  " << Clololo[i][j] << endl;

   return fail;
}
