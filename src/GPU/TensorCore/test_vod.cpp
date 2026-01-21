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

   if(fail == 1)
      cout << "\nTest on vectored octo double matmatmul failed?!!!\n\n";
   else
      cout << "\nTest on vectored octo double matmatmul succeeded.\n\n";

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
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrc; j++)
      {
         if(Ahihihi[i][j] < 0.0) Ahihihi[i][j] = -Ahihihi[i][j];
         if(Alohihi[i][j] < 0.0) Alohihi[i][j] = -Alohihi[i][j];
         if(Ahilohi[i][j] < 0.0) Ahilohi[i][j] = -Ahilohi[i][j];
         if(Alolohi[i][j] < 0.0) Alolohi[i][j] = -Alolohi[i][j];
         if(Ahihilo[i][j] < 0.0) Ahihilo[i][j] = -Ahihilo[i][j];
         if(Alohilo[i][j] < 0.0) Alohilo[i][j] = -Alohilo[i][j];
         if(Ahilolo[i][j] < 0.0) Ahilolo[i][j] = -Ahilolo[i][j];
         if(Alololo[i][j] < 0.0) Alololo[i][j] = -Alololo[i][j];
      }

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

   for(int i=0; i<nrc; i++)
      for(int j=0; j<ncols; j++)
      {
         if(Bhihihi[i][j] < 0.0) Bhihihi[i][j] = -Bhihihi[i][j];
         if(Blohihi[i][j] < 0.0) Blohihi[i][j] = -Blohihi[i][j];
         if(Bhilohi[i][j] < 0.0) Bhilohi[i][j] = -Bhilohi[i][j];
         if(Blolohi[i][j] < 0.0) Blolohi[i][j] = -Blolohi[i][j];
         if(Bhihilo[i][j] < 0.0) Bhihilo[i][j] = -Bhihilo[i][j];
         if(Blohilo[i][j] < 0.0) Blohilo[i][j] = -Blohilo[i][j];
         if(Bhilolo[i][j] < 0.0) Bhilolo[i][j] = -Bhilolo[i][j];
         if(Blololo[i][j] < 0.0) Blololo[i][j] = -Blololo[i][j];
      }

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

   double **Ahihihi0 = new double*[nrows];
   double **Ahihihi1 = new double*[nrows];
   double **Ahihihi2 = new double*[nrows];
   double **Ahihihi3 = new double*[nrows];
   double **Alohihi0 = new double*[nrows];
   double **Alohihi1 = new double*[nrows];
   double **Alohihi2 = new double*[nrows];
   double **Alohihi3 = new double*[nrows];
   double **Ahilohi0 = new double*[nrows];
   double **Ahilohi1 = new double*[nrows];
   double **Ahilohi2 = new double*[nrows];
   double **Ahilohi3 = new double*[nrows];
   double **Alolohi0 = new double*[nrows];
   double **Alolohi1 = new double*[nrows];
   double **Alolohi2 = new double*[nrows];
   double **Alolohi3 = new double*[nrows];
   double **Ahihilo0 = new double*[nrows];
   double **Ahihilo1 = new double*[nrows];
   double **Ahihilo2 = new double*[nrows];
   double **Ahihilo3 = new double*[nrows];
   double **Alohilo0 = new double*[nrows];
   double **Alohilo1 = new double*[nrows];
   double **Alohilo2 = new double*[nrows];
   double **Alohilo3 = new double*[nrows];
   double **Ahilolo0 = new double*[nrows];
   double **Ahilolo1 = new double*[nrows];
   double **Ahilolo2 = new double*[nrows];
   double **Ahilolo3 = new double*[nrows];
   double **Alololo0 = new double*[nrows];
   double **Alololo1 = new double*[nrows];
   double **Alololo2 = new double*[nrows];
   double **Alololo3 = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihihi0[i] = new double[nrc];
      Ahihihi1[i] = new double[nrc];
      Ahihihi2[i] = new double[nrc];
      Ahihihi3[i] = new double[nrc];
      Alohihi0[i] = new double[nrc];
      Alohihi1[i] = new double[nrc];
      Alohihi2[i] = new double[nrc];
      Alohihi3[i] = new double[nrc];
      Ahilohi0[i] = new double[nrc];
      Ahilohi1[i] = new double[nrc];
      Ahilohi2[i] = new double[nrc];
      Ahilohi3[i] = new double[nrc];
      Alolohi0[i] = new double[nrc];
      Alolohi1[i] = new double[nrc];
      Alolohi2[i] = new double[nrc];
      Alolohi3[i] = new double[nrc];
      Ahihilo0[i] = new double[nrc];
      Ahihilo1[i] = new double[nrc];
      Ahihilo2[i] = new double[nrc];
      Ahihilo3[i] = new double[nrc];
      Alohilo0[i] = new double[nrc];
      Alohilo1[i] = new double[nrc];
      Alohilo2[i] = new double[nrc];
      Alohilo3[i] = new double[nrc];
      Ahilolo0[i] = new double[nrc];
      Ahilolo1[i] = new double[nrc];
      Ahilolo2[i] = new double[nrc];
      Ahilolo3[i] = new double[nrc];
      Alololo0[i] = new double[nrc];
      Alololo1[i] = new double[nrc];
      Alololo2[i] = new double[nrc];
      Alololo3[i] = new double[nrc];
   }
   quarter_od_matrix
      (nrows, nrc,
       Ahihihi, Alohihi, Ahilohi, Alolohi, Ahihilo, Alohilo, Ahilolo, Alololo,
       Ahihihi0, Ahihihi1, Ahihihi2, Ahihihi3,
       Alohihi0, Alohihi1, Alohihi2, Alohihi3,
       Ahilohi0, Ahilohi1, Ahilohi2, Ahilohi3,
       Alolohi0, Alolohi1, Alolohi2, Alolohi3,
       Ahihilo0, Ahihilo1, Ahihilo2, Ahihilo3,
       Alohilo0, Alohilo1, Alohilo2, Alohilo3,
       Ahilolo0, Ahilolo1, Ahilolo2, Ahilolo3,
       Alololo0, Alololo1, Alololo2, Alololo3);

   double **Bhihihi0 = new double*[nrc];
   double **Bhihihi1 = new double*[nrc];
   double **Bhihihi2 = new double*[nrc];
   double **Bhihihi3 = new double*[nrc];
   double **Blohihi0 = new double*[nrc];
   double **Blohihi1 = new double*[nrc];
   double **Blohihi2 = new double*[nrc];
   double **Blohihi3 = new double*[nrc];
   double **Bhilohi0 = new double*[nrc];
   double **Bhilohi1 = new double*[nrc];
   double **Bhilohi2 = new double*[nrc];
   double **Bhilohi3 = new double*[nrc];
   double **Blolohi0 = new double*[nrc];
   double **Blolohi1 = new double*[nrc];
   double **Blolohi2 = new double*[nrc];
   double **Blolohi3 = new double*[nrc];
   double **Bhihilo0 = new double*[nrc];
   double **Bhihilo1 = new double*[nrc];
   double **Bhihilo2 = new double*[nrc];
   double **Bhihilo3 = new double*[nrc];
   double **Blohilo0 = new double*[nrc];
   double **Blohilo1 = new double*[nrc];
   double **Blohilo2 = new double*[nrc];
   double **Blohilo3 = new double*[nrc];
   double **Bhilolo0 = new double*[nrc];
   double **Bhilolo1 = new double*[nrc];
   double **Bhilolo2 = new double*[nrc];
   double **Bhilolo3 = new double*[nrc];
   double **Blololo0 = new double*[nrc];
   double **Blololo1 = new double*[nrc];
   double **Blololo2 = new double*[nrc];
   double **Blololo3 = new double*[nrc];

   for(int i=0; i<nrc; i++)
   {
      Bhihihi0[i] = new double[ncols];
      Bhihihi1[i] = new double[ncols];
      Bhihihi2[i] = new double[ncols];
      Bhihihi3[i] = new double[ncols];
      Blohihi0[i] = new double[ncols];
      Blohihi1[i] = new double[ncols];
      Blohihi2[i] = new double[ncols];
      Blohihi3[i] = new double[ncols];
      Bhilohi0[i] = new double[ncols];
      Bhilohi1[i] = new double[ncols];
      Bhilohi2[i] = new double[ncols];
      Bhilohi3[i] = new double[ncols];
      Blolohi0[i] = new double[ncols];
      Blolohi1[i] = new double[ncols];
      Blolohi2[i] = new double[ncols];
      Blolohi3[i] = new double[ncols];
      Bhihilo0[i] = new double[ncols];
      Bhihilo1[i] = new double[ncols];
      Bhihilo2[i] = new double[ncols];
      Bhihilo3[i] = new double[ncols];
      Blohilo0[i] = new double[ncols];
      Blohilo1[i] = new double[ncols];
      Blohilo2[i] = new double[ncols];
      Blohilo3[i] = new double[ncols];
      Bhilolo0[i] = new double[ncols];
      Bhilolo1[i] = new double[ncols];
      Bhilolo2[i] = new double[ncols];
      Bhilolo3[i] = new double[ncols];
      Blololo0[i] = new double[ncols];
      Blololo1[i] = new double[ncols];
      Blololo2[i] = new double[ncols];
      Blololo3[i] = new double[ncols];
   }
   quarter_od_matrix
      (nrc, ncols,
       Bhihihi, Blohihi, Bhilohi, Blolohi,
       Bhihilo, Blohilo, Bhilolo, Blololo,
       Bhihihi0, Bhihihi1, Bhihihi2, Bhihihi3,
       Blohihi0, Blohihi1, Blohihi2, Blohihi3,
       Bhilohi0, Bhilohi1, Bhilohi2, Bhilohi3,
       Blolohi0, Blolohi1, Blolohi2, Blolohi3,
       Bhihilo0, Bhihilo1, Bhihilo2, Bhihilo3,
       Blohilo0, Blohilo1, Blohilo2, Blohilo3,
       Bhilolo0, Bhilolo1, Bhilolo2, Bhilolo3,
       Blololo0, Blololo1, Blololo2, Blololo3);

   double **Thihihi0 = new double*[ncols];
   double **Thihihi1 = new double*[ncols];
   double **Thihihi2 = new double*[ncols];
   double **Thihihi3 = new double*[ncols];
   double **Tlohihi0 = new double*[ncols];
   double **Tlohihi1 = new double*[ncols];
   double **Tlohihi2 = new double*[ncols];
   double **Tlohihi3 = new double*[ncols];
   double **Thilohi0 = new double*[ncols];
   double **Thilohi1 = new double*[ncols];
   double **Thilohi2 = new double*[ncols];
   double **Thilohi3 = new double*[ncols];
   double **Tlolohi0 = new double*[ncols];
   double **Tlolohi1 = new double*[ncols];
   double **Tlolohi2 = new double*[ncols];
   double **Tlolohi3 = new double*[ncols];
   double **Thihilo0 = new double*[ncols];
   double **Thihilo1 = new double*[ncols];
   double **Thihilo2 = new double*[ncols];
   double **Thihilo3 = new double*[ncols];
   double **Tlohilo0 = new double*[ncols];
   double **Tlohilo1 = new double*[ncols];
   double **Tlohilo2 = new double*[ncols];
   double **Tlohilo3 = new double*[ncols];
   double **Thilolo0 = new double*[ncols];
   double **Thilolo1 = new double*[ncols];
   double **Thilolo2 = new double*[ncols];
   double **Thilolo3 = new double*[ncols];
   double **Tlololo0 = new double*[ncols];
   double **Tlololo1 = new double*[ncols];
   double **Tlololo2 = new double*[ncols];
   double **Tlololo3 = new double*[ncols];

   for(int i=0; i<ncols; i++)
   {
      Thihihi0[i] = new double[nrc];
      Thihihi1[i] = new double[nrc];
      Thihihi2[i] = new double[nrc];
      Thihihi3[i] = new double[nrc];
      Tlohihi0[i] = new double[nrc];
      Tlohihi1[i] = new double[nrc];
      Tlohihi2[i] = new double[nrc];
      Tlohihi3[i] = new double[nrc];
      Thilohi0[i] = new double[nrc];
      Thilohi1[i] = new double[nrc];
      Thilohi2[i] = new double[nrc];
      Thilohi3[i] = new double[nrc];
      Tlolohi0[i] = new double[nrc];
      Tlolohi1[i] = new double[nrc];
      Tlolohi2[i] = new double[nrc];
      Tlolohi3[i] = new double[nrc];
      Thihilo0[i] = new double[nrc];
      Thihilo1[i] = new double[nrc];
      Thihilo2[i] = new double[nrc];
      Thihilo3[i] = new double[nrc];
      Tlohilo0[i] = new double[nrc];
      Tlohilo1[i] = new double[nrc];
      Tlohilo2[i] = new double[nrc];
      Tlohilo3[i] = new double[nrc];
      Thilolo0[i] = new double[nrc];
      Thilolo1[i] = new double[nrc];
      Thilolo2[i] = new double[nrc];
      Thilolo3[i] = new double[nrc];
      Tlololo0[i] = new double[nrc];
      Tlololo1[i] = new double[nrc];
      Tlololo2[i] = new double[nrc];
      Tlololo3[i] = new double[nrc];
   }
   transpose_od_quarters
      (nrc, ncols,
       Bhihihi0, Bhihihi1, Bhihihi2, Bhihihi3,
       Blohihi0, Blohihi1, Blohihi2, Blohihi3,
       Bhilohi0, Bhilohi1, Bhilohi2, Bhilohi3,
       Blolohi0, Blolohi1, Blolohi2, Blolohi3,
       Bhihilo0, Bhihilo1, Bhihilo2, Bhihilo3,
       Blohilo0, Blohilo1, Blohilo2, Blohilo3,
       Bhilolo0, Bhilolo1, Bhilolo2, Bhilolo3,
       Blololo0, Blololo1, Blololo2, Blololo3,
       Thihihi0, Thihihi1, Thihihi2, Thihihi3,
       Tlohihi0, Tlohihi1, Tlohihi2, Tlohihi3,
       Thilohi0, Thilohi1, Thilohi2, Thilohi3,
       Tlolohi0, Tlolohi1, Tlolohi2, Tlolohi3,
       Thihilo0, Thihilo1, Thihilo2, Thihilo3,
       Tlohilo0, Tlohilo1, Tlohilo2, Tlohilo3,
       Thilolo0, Thilolo1, Thilolo2, Thilolo3,
       Tlololo0, Tlololo1, Tlololo2, Tlololo3);

   double **Chihihi0 = new double*[nrows];
   double **Chihihi1 = new double*[nrows];
   double **Chihihi2 = new double*[nrows];
   double **Chihihi3 = new double*[nrows];
   double **Clohihi0 = new double*[nrows];
   double **Clohihi1 = new double*[nrows];
   double **Clohihi2 = new double*[nrows];
   double **Clohihi3 = new double*[nrows];
   double **Chilohi0 = new double*[nrows];
   double **Chilohi1 = new double*[nrows];
   double **Chilohi2 = new double*[nrows];
   double **Chilohi3 = new double*[nrows];
   double **Clolohi0 = new double*[nrows];
   double **Clolohi1 = new double*[nrows];
   double **Clolohi2 = new double*[nrows];
   double **Clolohi3 = new double*[nrows];
   double **Chihilo0 = new double*[nrows];
   double **Chihilo1 = new double*[nrows];
   double **Chihilo2 = new double*[nrows];
   double **Chihilo3 = new double*[nrows];
   double **Clohilo0 = new double*[nrows];
   double **Clohilo1 = new double*[nrows];
   double **Clohilo2 = new double*[nrows];
   double **Clohilo3 = new double*[nrows];
   double **Chilolo0 = new double*[nrows];
   double **Chilolo1 = new double*[nrows];
   double **Chilolo2 = new double*[nrows];
   double **Chilolo3 = new double*[nrows];
   double **Clololo0 = new double*[nrows];
   double **Clololo1 = new double*[nrows];
   double **Clololo2 = new double*[nrows];
   double **Clololo3 = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Chihihi0[i] = new double[ncols];
      Chihihi1[i] = new double[ncols];
      Chihihi2[i] = new double[ncols];
      Chihihi3[i] = new double[ncols];
      Clohihi0[i] = new double[ncols];
      Clohihi1[i] = new double[ncols];
      Clohihi2[i] = new double[ncols];
      Clohihi3[i] = new double[ncols];
      Chilohi0[i] = new double[ncols];
      Chilohi1[i] = new double[ncols];
      Chilohi2[i] = new double[ncols];
      Chilohi3[i] = new double[ncols];
      Clolohi0[i] = new double[ncols];
      Clolohi1[i] = new double[ncols];
      Clolohi2[i] = new double[ncols];
      Clolohi3[i] = new double[ncols];
      Chihilo0[i] = new double[ncols];
      Chihilo1[i] = new double[ncols];
      Chihilo2[i] = new double[ncols];
      Chihilo3[i] = new double[ncols];
      Clohilo0[i] = new double[ncols];
      Clohilo1[i] = new double[ncols];
      Clohilo2[i] = new double[ncols];
      Clohilo3[i] = new double[ncols];
      Chilolo0[i] = new double[ncols];
      Chilolo1[i] = new double[ncols];
      Chilolo2[i] = new double[ncols];
      Chilolo3[i] = new double[ncols];
      Clololo0[i] = new double[ncols];
      Clololo1[i] = new double[ncols];
      Clololo2[i] = new double[ncols];
      Clololo3[i] = new double[ncols];
   }
   vectored_od_matmatmul
      (nrows, ncols, nrc,
       Ahihihi0, Ahihihi1, Ahihihi2, Ahihihi3,
       Alohihi0, Alohihi1, Alohihi2, Alohihi3,
       Ahilohi0, Ahilohi1, Ahilohi2, Ahilohi3,
       Alolohi0, Alolohi1, Alolohi2, Alolohi3,
       Ahihilo0, Ahihilo1, Ahihilo2, Ahihilo3,
       Alohilo0, Alohilo1, Alohilo2, Alohilo3,
       Ahilolo0, Ahilolo1, Ahilolo2, Ahilolo3,
       Alololo0, Alololo1, Alololo2, Alololo3,
       Thihihi0, Thihihi1, Thihihi2, Thihihi3,
       Tlohihi0, Tlohihi1, Tlohihi2, Tlohihi3,
       Thilohi0, Thilohi1, Thilohi2, Thilohi3,
       Tlolohi0, Tlolohi1, Tlolohi2, Tlolohi3,
       Thihilo0, Thihilo1, Thihilo2, Thihilo3,
       Tlohilo0, Tlohilo1, Tlohilo2, Tlohilo3,
       Thilolo0, Thilolo1, Thilolo2, Thilolo3,
       Tlololo0, Tlololo1, Tlololo2, Tlololo3,
       Chihihi0, Chihihi1, Chihihi2, Chihihi3,
       Clohihi0, Clohihi1, Clohihi2, Clohihi3,
       Chilohi0, Chilohi1, Chilohi2, Chilohi3,
       Clolohi0, Clolohi1, Clolohi2, Clolohi3,
       Chihilo0, Chihilo1, Chihilo2, Chihilo3,
       Clohilo0, Clohilo1, Clohilo2, Clohilo3,
       Chilolo0, Chilolo1, Chilolo2, Chilolo3,
       Clololo0, Clololo1, Clololo2, Clololo3);

   double **Vhihihi = new double*[nrows];
   double **Vlohihi = new double*[nrows];
   double **Vhilohi = new double*[nrows];
   double **Vlolohi = new double*[nrows];
   double **Vhihilo = new double*[nrows];
   double **Vlohilo = new double*[nrows];
   double **Vhilolo = new double*[nrows];
   double **Vlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Vhihihi[i] = new double[ncols];
      Vlohihi[i] = new double[ncols];
      Vhilohi[i] = new double[ncols];
      Vlolohi[i] = new double[ncols];
      Vhihilo[i] = new double[ncols];
      Vlohilo[i] = new double[ncols];
      Vhilolo[i] = new double[ncols];
      Vlololo[i] = new double[ncols];
   }
   to_octo_double_matrix
      (nrows, ncols,
       Chihihi0, Chihihi1, Chihihi2, Chihihi3,
       Clohihi0, Clohihi1, Clohihi2, Clohihi3,
       Chilohi0, Chilohi1, Chilohi2, Chilohi3,
       Clolohi0, Clolohi1, Clolohi2, Clolohi3,
       Chihilo0, Chihilo1, Chihilo2, Chihilo3,
       Clohilo0, Clohilo1, Clohilo2, Clohilo3,
       Chilolo0, Chilolo1, Chilolo2, Chilolo3,
       Clololo0, Clololo1, Clololo2, Clololo3,
       Vhihihi, Vlohihi, Vhilohi, Vlolohi,
       Vhihilo, Vlohilo, Vhilolo, Vlololo);

   cout << "the vectored product A*B :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "V[" << i << "][" << j << "] : "
              << Vhihihi[i][j] << "  " << Vlohihi[i][j] << endl
              << "          "
              << Vhilohi[i][j] << "  " << Vlolohi[i][j] << endl
              << "          "
              << Vhihilo[i][j] << "  " << Vlohilo[i][j] << endl
              << "          "
              << Vhilolo[i][j] << "  " << Vlololo[i][j] << endl;

   double err[8],acc[8];
   err[0] = 0.0; err[1] = 0.0;
   err[2] = 0.0; err[3] = 0.0;
   err[4] = 0.0; err[5] = 0.0;
   err[6] = 0.0; err[7] = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         odf_sub(Chihihi[i][j], Clohihi[i][j], Chilohi[i][j], Clolohi[i][j],
                 Chihilo[i][j], Clohilo[i][j], Chilolo[i][j], Clololo[i][j],
                 Vhihihi[i][j], Vlohihi[i][j], Vhilohi[i][j], Vlolohi[i][j],
                 Vhihilo[i][j], Vlohilo[i][j], Vhilolo[i][j], Vlololo[i][j],
                 &acc[0], &acc[1], &acc[2], &acc[3],
                 &acc[4], &acc[5], &acc[6], &acc[7]);
         odf_inc(&err[0], &err[1], &err[2], &err[3],
                 &err[4], &err[5], &err[6], &err[7],
                 acc[0], acc[1], acc[2], acc[3],
                 acc[4], acc[5], acc[6], acc[7]);
      }

   if(err[0] < 0.0) odf_minus(&err[0], &err[1], &err[2], &err[3],
                              &err[4], &err[5], &err[6], &err[7]);

   cout << "-> error : " << endl; od_write_doubles(err); cout << endl;

   fail = (abs(err[0]) > 1.0E-120);

   if(fail == 1) return fail; // no point to continue

   double **cA = new double*[32*nrows];
   for(int i=0; i<32*nrows; i++) cA[i] = new double[32*nrc];

   od_convolute_quarters
      (nrows, nrc,
       Ahihihi0, Ahihihi1, Ahihihi2, Ahihihi3,
       Alohihi0, Alohihi1, Alohihi2, Alohihi3,
       Ahilohi0, Ahilohi1, Ahilohi2, Ahilohi3,
       Alolohi0, Alolohi1, Alolohi2, Alolohi3,
       Ahihilo0, Ahihilo1, Ahihilo2, Ahihilo3,
       Alohilo0, Alohilo1, Alohilo2, Alohilo3,
       Ahilolo0, Ahilolo1, Ahilolo2, Ahilolo3,
       Alololo0, Alololo1, Alololo2, Alololo3, cA);

   cout << "the convoluted quartered matrix A :" << endl;
   for(int i=0; i<32*nrows; i++)
      for(int j=0; j<32*nrc; j++)
         cout << "cA[" << i << "][" << j << "] : " << cA[i][j] << endl;

   double **sB = new double*[32*nrc];
   for(int i=0; i<32*nrc; i++) sB[i] = new double[ncols];

   od_stack_quarters
      (nrc, ncols,
       Bhihihi0, Bhihihi1, Bhihihi2, Bhihihi3,
       Blohihi0, Blohihi1, Blohihi2, Blohihi3,
       Bhilohi0, Bhilohi1, Bhilohi2, Bhilohi3,
       Blolohi0, Blolohi1, Blolohi2, Blolohi3,
       Bhihilo0, Bhihilo1, Bhihilo2, Bhihilo3,
       Blohilo0, Blohilo1, Blohilo2, Blohilo3,
       Bhilolo0, Bhilolo1, Bhilolo2, Bhilolo3,
       Blololo0, Blololo1, Blololo2, Blololo3, sB);

   cout << "the stacked quartered matrix B :" << endl;
   for(int i=0; i<32*nrc; i++)
      for(int j=0; j<ncols; j++)
         cout << "sB[" << i << "][" << j << "] : " << sB[i][j] << endl;

   return fail;
}
