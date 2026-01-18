/* Tests the collection of functions on vectored hexa double arithmetic.  */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "hexa_double.h"
#include "random16_vectors.h"
#include "hexa_double_functions.h"
#include "vectored_hexa_doubles.h"

int test_quarter_hexa_double ( void );
/*
 * Generates a random positive hexa double, quarters the parts,
 * and then checks if their sum equals the original hexa double.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_hd_product ( int dim );
/*
 * Generates two random vectors of hexa doubles of size dim,
 * and compares their inner product with the vectored inner product.
 * Returns 1 if the test failed, returns 0 otherwise. */

using namespace std;

int main ( void )
{
   int fail;

   srand(time(NULL));

   fail = test_quarter_hexa_double();

   if(fail == 1)
      cout << "\nTest on quarter hexa double failed?!!!\n\n";
   else
      cout << "\nTest on quarter hexa double succeeded.\n\n";

   cout << "Give the dimension : ";
   int dim; cin >> dim;

   fail = test_vectored_hd_product(dim);

   if(fail == 1)
      cout << "\nTest on vectored hexa double product failed?!!!\n\n";
   else
      cout << "\nTest on vectored hexa double product succeeded.\n\n";

   return 0;
}

int test_quarter_hexa_double ( void )
{
   int fail = 0;

   double x[16];
   double xhihihihi0,xhihihihi1,xhihihihi2,xhihihihi3;
   double xlohihihi0,xlohihihi1,xlohihihi2,xlohihihi3;
   double xhilohihi0,xhilohihi1,xhilohihi2,xhilohihi3;
   double xlolohihi0,xlolohihi1,xlolohihi2,xlolohihi3;
   double xhihilohi0,xhihilohi1,xhihilohi2,xhihilohi3;
   double xlohilohi0,xlohilohi1,xlohilohi2,xlohilohi3;
   double xhilolohi0,xhilolohi1,xhilolohi2,xhilolohi3;
   double xlololohi0,xlololohi1,xlololohi2,xlololohi3;
   double xhihihilo0,xhihihilo1,xhihihilo2,xhihihilo3;
   double xlohihilo0,xlohihilo1,xlohihilo2,xlohihilo3;
   double xhilohilo0,xhilohilo1,xhilohilo2,xhilohilo3;
   double xlolohilo0,xlolohilo1,xlolohilo2,xlolohilo3;
   double xhihilolo0,xhihilolo1,xhihilolo2,xhihilolo3;
   double xlohilolo0,xlohilolo1,xlohilolo2,xlohilolo3;
   double xhilololo0,xhilololo1,xhilololo2,xhilololo3;
   double xlolololo0,xlolololo1,xlolololo2,xlolololo3;
   double y[16];
   double e[16];

   random_hexa_double
      (&x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &x[6], &x[7],
       &x[8], &x[9], &x[10], &x[11], &x[12], &x[13], &x[14], &x[15]);

   if(x[0] < 0.0) x[0] = -x[0]; // all parts must be positive
   if(x[1] < 0.0) x[1] = -x[1];
   if(x[2] < 0.0) x[2] = -x[2];
   if(x[3] < 0.0) x[3] = -x[3];
   if(x[4] < 0.0) x[4] = -x[4];
   if(x[5] < 0.0) x[5] = -x[5];
   if(x[6] < 0.0) x[6] = -x[6];
   if(x[7] < 0.0) x[7] = -x[7];
   if(x[8] < 0.0) x[8] = -x[8];
   if(x[9] < 0.0) x[9] = -x[9];
   if(x[10] < 0.0) x[10] = -x[10];
   if(x[11] < 0.0) x[11] = -x[11];
   if(x[12] < 0.0) x[12] = -x[12];
   if(x[13] < 0.0) x[13] = -x[13];
   if(x[14] < 0.0) x[14] = -x[14];
   if(x[15] < 0.0) x[15] = -x[15];

   cout << scientific << setprecision(16);

   cout << "x :"; cout << endl; hd_write_doubles(x); cout << endl;

   quarter_hexa_double
      (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
       x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15],
       &xhihihihi0, &xhihihihi1, &xhihihihi2, &xhihihihi3,
       &xlohihihi0, &xlohihihi1, &xlohihihi2, &xlohihihi3,
       &xhilohihi0, &xhilohihi1, &xhilohihi2, &xhilohihi3,
       &xlolohihi0, &xlolohihi1, &xlolohihi2, &xlolohihi3,
       &xhihilohi0, &xhihilohi1, &xhihilohi2, &xhihilohi3,
       &xlohilohi0, &xlohilohi1, &xlohilohi2, &xlohilohi3,
       &xhilolohi0, &xhilolohi1, &xhilolohi2, &xhilolohi3,
       &xlololohi0, &xlololohi1, &xlololohi2, &xlololohi3,
       &xhihihilo0, &xhihihilo1, &xhihihilo2, &xhihihilo3,
       &xlohihilo0, &xlohihilo1, &xlohihilo2, &xlohihilo3,
       &xhilohilo0, &xhilohilo1, &xhilohilo2, &xhilohilo3,
       &xlolohilo0, &xlolohilo1, &xlolohilo2, &xlolohilo3,
       &xhihilolo0, &xhihilolo1, &xhihilolo2, &xhihilolo3,
       &xlohilolo0, &xlohilolo1, &xlohilolo2, &xlohilolo3,
       &xhilololo0, &xhilololo1, &xhilololo2, &xhilololo3,
       &xlolololo0, &xlolololo1, &xlolololo2, &xlolololo3);

   cout << "xhihihihi0 : " << xhihihihi0 << endl;
   cout << "xhihihihi1 : " << xhihihihi1 << endl;
   cout << "xhihihihi2 : " << xhihihihi2 << endl;
   cout << "xhihihihi3 : " << xhihihihi3 << endl;
   cout << "xlohihihi0 : " << xlohihihi0 << endl;
   cout << "xlohihihi1 : " << xlohihihi1 << endl;
   cout << "xlohihihi2 : " << xlohihihi2 << endl;
   cout << "xlohihihi3 : " << xlohihihi3 << endl;
   cout << "xhilohihi0 : " << xhilohihi0 << endl;
   cout << "xhilohihi1 : " << xhilohihi1 << endl;
   cout << "xhilohihi2 : " << xhilohihi2 << endl;
   cout << "xhilohihi3 : " << xhilohihi3 << endl;
   cout << "xlolohihi0 : " << xlolohihi0 << endl;
   cout << "xlolohihi1 : " << xlolohihi1 << endl;
   cout << "xlolohihi2 : " << xlolohihi2 << endl;
   cout << "xlolohihi3 : " << xlolohihi3 << endl;
   cout << "xhihilohi0 : " << xhihilohi0 << endl;
   cout << "xhihilohi1 : " << xhihilohi1 << endl;
   cout << "xhihilohi2 : " << xhihilohi2 << endl;
   cout << "xhihilohi3 : " << xhihilohi3 << endl;
   cout << "xlohilohi0 : " << xlohilohi0 << endl;
   cout << "xlohilohi1 : " << xlohilohi1 << endl;
   cout << "xlohilohi2 : " << xlohilohi2 << endl;
   cout << "xlohilohi3 : " << xlohilohi3 << endl;
   cout << "xhilolohi0 : " << xhilolohi0 << endl;
   cout << "xhilolohi1 : " << xhilolohi1 << endl;
   cout << "xhilolohi2 : " << xhilolohi2 << endl;
   cout << "xhilolohi3 : " << xhilolohi3 << endl;
   cout << "xlololohi0 : " << xlololohi0 << endl;
   cout << "xlololohi1 : " << xlololohi1 << endl;
   cout << "xlololohi2 : " << xlololohi2 << endl;
   cout << "xlololohi3 : " << xlololohi3 << endl;
   cout << "xhihihilo0 : " << xhihihilo0 << endl;
   cout << "xhihihilo1 : " << xhihihilo1 << endl;
   cout << "xhihihilo2 : " << xhihihilo2 << endl;
   cout << "xhihihilo3 : " << xhihihilo3 << endl;
   cout << "xlohihilo0 : " << xlohihilo0 << endl;
   cout << "xlohihilo1 : " << xlohihilo1 << endl;
   cout << "xlohihilo2 : " << xlohihilo2 << endl;
   cout << "xlohihilo3 : " << xlohihilo3 << endl;
   cout << "xhilohilo0 : " << xhilohilo0 << endl;
   cout << "xhilohilo1 : " << xhilohilo1 << endl;
   cout << "xhilohilo2 : " << xhilohilo2 << endl;
   cout << "xhilohilo3 : " << xhilohilo3 << endl;
   cout << "xlolohilo0 : " << xlolohilo0 << endl;
   cout << "xlolohilo1 : " << xlolohilo1 << endl;
   cout << "xlolohilo2 : " << xlolohilo2 << endl;
   cout << "xlolohilo3 : " << xlolohilo3 << endl;
   cout << "xhihilolo0 : " << xhihilolo0 << endl;
   cout << "xhihilolo1 : " << xhihilolo1 << endl;
   cout << "xhihilolo2 : " << xhihilolo2 << endl;
   cout << "xhihilolo3 : " << xhihilolo3 << endl;
   cout << "xlohilolo0 : " << xlohilolo0 << endl;
   cout << "xlohilolo1 : " << xlohilolo1 << endl;
   cout << "xlohilolo2 : " << xlohilolo2 << endl;
   cout << "xlohilolo3 : " << xlohilolo3 << endl;
   cout << "xhilololo0 : " << xhilololo0 << endl;
   cout << "xhilololo1 : " << xhilololo1 << endl;
   cout << "xhilololo2 : " << xhilololo2 << endl;
   cout << "xhilololo3 : " << xhilololo3 << endl;
   cout << "xlolololo0 : " << xlolololo0 << endl;
   cout << "xlolololo1 : " << xlolololo1 << endl;
   cout << "xlolololo2 : " << xlolololo2 << endl;
   cout << "xlolololo3 : " << xlolololo3 << endl;

   to_hexa_double
      (xhihihihi0, xhihihihi1, xhihihihi2, xhihihihi3,
       xlohihihi0, xlohihihi1, xlohihihi2, xlohihihi3,
       xhilohihi0, xhilohihi1, xhilohihi2, xhilohihi3,
       xlolohihi0, xlolohihi1, xlolohihi2, xlolohihi3,
       xhihilohi0, xhihilohi1, xhihilohi2, xhihilohi3,
       xlohilohi0, xlohilohi1, xlohilohi2, xlohilohi3,
       xhilolohi0, xhilolohi1, xhilolohi2, xhilolohi3,
       xlololohi0, xlololohi1, xlololohi2, xlololohi3,
       xhihihilo0, xhihihilo1, xhihihilo2, xhihihilo3,
       xlohihilo0, xlohihilo1, xlohihilo2, xlohihilo3,
       xhilohilo0, xhilohilo1, xhilohilo2, xhilohilo3,
       xlolohilo0, xlolohilo1, xlolohilo2, xlolohilo3,
       xhihilolo0, xhihilolo1, xhihilolo2, xhihilolo3,
       xlohilolo0, xlohilolo1, xlohilolo2, xlohilolo3,
       xhilololo0, xhilololo1, xhilololo2, xhilololo3,
       xlolololo0, xlolololo1, xlolololo2, xlolololo3,
       &y[0], &y[1], &y[2], &y[3], &y[4], &y[5], &y[6], &y[7],
       &y[8], &y[9], &y[10], &y[11], &y[12], &y[13], &y[14], &y[15]);

   cout << "x :"; cout << endl; hd_write_doubles(x); cout << endl;
   cout << "y :"; cout << endl; hd_write_doubles(y); cout << endl;

   hdf_sub(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
           x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15],
           y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
           y[8], y[9], y[10], y[11], y[12], y[13], y[14], y[15],
           &e[0], &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7],
           &e[8], &e[9], &e[10], &e[11], &e[12], &e[13], &e[14], &e[15]);

   cout << "e :"; cout << endl; hd_write_doubles(e); cout << endl;

   fail = not(e[0] == 0.0) + not(e[1] == 0.0)
        + not(e[2] == 0.0) + not(e[3] == 0.0)
        + not(e[4] == 0.0) + not(e[5] == 0.0)
        + not(e[6] == 0.0) + not(e[7] == 0.0)
        + not(e[8] == 0.0) + not(e[9] == 0.0)
        + not(e[10] == 0.0) + not(e[11] == 0.0)
        + not(e[12] == 0.0) + not(e[13] == 0.0)
        + not(e[14] == 0.0) + not(e[15] == 0.0);

   return fail;
}

int test_vectored_hd_product ( int dim )
{
   int fail = 0;

   double xhihihihi[dim],xlohihihi[dim],xhilohihi[dim],xlolohihi[dim];
   double xhihilohi[dim],xlohilohi[dim],xhilolohi[dim],xlololohi[dim];
   double xhihihilo[dim],xlohihilo[dim],xhilohilo[dim],xlolohilo[dim];
   double xhihilolo[dim],xlohilolo[dim],xhilololo[dim],xlolololo[dim];
   double yhihihihi[dim],ylohihihi[dim],yhilohihi[dim],ylolohihi[dim];
   double yhihilohi[dim],ylohilohi[dim],yhilolohi[dim],ylololohi[dim];
   double yhihihilo[dim],ylohihilo[dim],yhilohilo[dim],ylolohilo[dim];
   double yhihilolo[dim],ylohilolo[dim],yhilololo[dim],ylolololo[dim];

   double x0[dim],x1[dim],x2[dim],x3[dim],x4[dim],x5[dim],x6[dim],x7[dim];
   double x8[dim],x9[dim],x10[dim],x11[dim],x12[dim],x13[dim],x14[dim];
   double x15[dim],x16[dim],x17[dim],x18[dim],x19[dim],x20[dim],x21[dim];
   double x22[dim],x23[dim],x24[dim],x25[dim],x26[dim],x27[dim],x28[dim];
   double x29[dim],x30[dim],x31[dim],x32[dim],x33[dim],x34[dim],x35[dim];
   double x36[dim],x37[dim],x38[dim],x39[dim],x40[dim],x41[dim],x42[dim];
   double x43[dim],x44[dim],x45[dim],x46[dim],x47[dim],x48[dim],x49[dim];
   double x50[dim],x51[dim],x52[dim],x53[dim],x54[dim],x55[dim],x56[dim];
   double x57[dim],x58[dim],x59[dim],x60[dim],x61[dim],x62[dim],x63[dim];

   double y0[dim],y1[dim],y2[dim],y3[dim],y4[dim],y5[dim],y6[dim],y7[dim];
   double y8[dim],y9[dim],y10[dim],y11[dim],y12[dim],y13[dim],y14[dim];
   double y15[dim],y16[dim],y17[dim],y18[dim],y19[dim],y20[dim],y21[dim];
   double y22[dim],y23[dim],y24[dim],y25[dim],y26[dim],y27[dim],y28[dim];
   double y29[dim],y30[dim],y31[dim],y32[dim],y33[dim],y34[dim],y35[dim];
   double y36[dim],y37[dim],y38[dim],y39[dim],y40[dim],y41[dim],y42[dim];
   double y43[dim],y44[dim],y45[dim],y46[dim],y47[dim],y48[dim],y49[dim];
   double y50[dim],y51[dim],y52[dim],y53[dim],y54[dim],y55[dim],y56[dim];
   double y57[dim],y58[dim],y59[dim],y60[dim],y61[dim],y62[dim],y63[dim];

   double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;
   double s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31;
   double s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47;
   double s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63;

   double prd[16],vpd[16],err[16];

   for(int i=0; i<dim; i++)
   {
      random_hexa_double
         (&xhihihihi[i], &xlohihihi[i], &xhilohihi[i], &xlolohihi[i],
          &xhihilohi[i], &xlohilohi[i], &xhilolohi[i], &xlololohi[i],
          &xhihihilo[i], &xlohihilo[i], &xhilohilo[i], &xlolohilo[i],
          &xhihilolo[i], &xlohilolo[i], &xhilololo[i], &xlolololo[i]);

      if(xhihihihi[i] < 0.0) xhihihihi[i] = -xhihihihi[i];
      if(xlohihihi[i] < 0.0) xlohihihi[i] = -xlohihihi[i];
      if(xhilohihi[i] < 0.0) xhilohihi[i] = -xhilohihi[i];
      if(xlolohihi[i] < 0.0) xlolohihi[i] = -xlolohihi[i];
      if(xhihilohi[i] < 0.0) xhihilohi[i] = -xhihilohi[i];
      if(xlohilohi[i] < 0.0) xlohilohi[i] = -xlohilohi[i];
      if(xhilolohi[i] < 0.0) xhilolohi[i] = -xhilolohi[i];
      if(xlololohi[i] < 0.0) xlololohi[i] = -xlololohi[i];
      if(xhihihilo[i] < 0.0) xhihihilo[i] = -xhihihilo[i];
      if(xlohihilo[i] < 0.0) xlohihilo[i] = -xlohihilo[i];
      if(xhilohilo[i] < 0.0) xhilohilo[i] = -xhilohilo[i];
      if(xlolohilo[i] < 0.0) xlolohilo[i] = -xlolohilo[i];
      if(xhihilolo[i] < 0.0) xhihilolo[i] = -xhihilolo[i];
      if(xlohilolo[i] < 0.0) xlohilolo[i] = -xlohilolo[i];
      if(xhilololo[i] < 0.0) xhilololo[i] = -xhilololo[i];
      if(xlolololo[i] < 0.0) xlolololo[i] = -xlolololo[i];

      random_hexa_double
         (&yhihihihi[i], &ylohihihi[i], &yhilohihi[i], &ylolohihi[i],
          &yhihilohi[i], &ylohilohi[i], &yhilolohi[i], &ylololohi[i],
          &yhihihilo[i], &ylohihilo[i], &yhilohilo[i], &ylolohilo[i],
          &yhihilolo[i], &ylohilolo[i], &yhilololo[i], &ylolololo[i]);

      if(yhihihihi[i] < 0.0) yhihihihi[i] = -yhihihihi[i];
      if(ylohihihi[i] < 0.0) ylohihihi[i] = -ylohihihi[i];
      if(yhilohihi[i] < 0.0) yhilohihi[i] = -yhilohihi[i];
      if(ylolohihi[i] < 0.0) ylolohihi[i] = -ylolohihi[i];
      if(yhihilohi[i] < 0.0) yhihilohi[i] = -yhihilohi[i];
      if(ylohilohi[i] < 0.0) ylohilohi[i] = -ylohilohi[i];
      if(yhilolohi[i] < 0.0) yhilolohi[i] = -yhilolohi[i];
      if(ylololohi[i] < 0.0) ylololohi[i] = -ylololohi[i];
      if(yhihihilo[i] < 0.0) yhihihilo[i] = -yhihihilo[i];
      if(ylohihilo[i] < 0.0) ylohihilo[i] = -ylohihilo[i];
      if(yhilohilo[i] < 0.0) yhilohilo[i] = -yhilohilo[i];
      if(ylolohilo[i] < 0.0) ylolohilo[i] = -ylolohilo[i];
      if(yhihilolo[i] < 0.0) yhihilolo[i] = -yhihilolo[i];
      if(ylohilolo[i] < 0.0) ylohilolo[i] = -ylohilolo[i];
      if(yhilololo[i] < 0.0) yhilololo[i] = -yhilololo[i];
      if(ylolololo[i] < 0.0) ylolololo[i] = -ylolololo[i];
   }
   cout << "hexa double vector x :" << endl;
   hd_write_vector
      (dim, xhihihihi, xlohihihi, xhilohihi, xlolohihi,
            xhihilohi, xlohilohi, xhilolohi, xlololohi,
            xhihihilo, xlohihilo, xhilohilo, xlolohilo,
            xhihilolo, xlohilolo, xhilololo, xlolololo);
   cout << "hexa double vector y :" << endl;
   hd_write_vector
      (dim, yhihihihi, ylohihihi, yhilohihi, ylolohihi,
            yhihilohi, ylohilohi, yhilolohi, ylololohi,
            yhihihilo, ylohihilo, yhilohilo, ylolohilo,
            yhihilolo, ylohilolo, yhilololo, ylolololo);

   hexa_double_product
      (dim, xhihihihi, xlohihihi, xhilohihi, xlolohihi,
            xhihilohi, xlohilohi, xhilolohi, xlololohi,
            xhihihilo, xlohihilo, xhilohilo, xlolohilo,
            xhihilolo, xlohilolo, xhilololo, xlolololo,
            yhihihihi, ylohihihi, yhilohihi, ylolohihi,
            yhihilohi, ylohilohi, yhilolohi, ylololohi,
            yhihihilo, ylohihilo, yhilohilo, ylolohilo,
            yhihilolo, ylohilolo, yhilololo, ylolololo,
       &prd[0], &prd[1], &prd[2], &prd[3], &prd[4], &prd[5],
       &prd[6], &prd[7], &prd[8], &prd[9], &prd[10], &prd[11],
       &prd[12], &prd[13], &prd[14], &prd[15]);

   quarter_hd_vector
      (dim, xhihihihi, xlohihihi, xhilohihi, xlolohihi,
            xhihilohi, xlohilohi, xhilolohi, xlololohi,
            xhihihilo, xlohihilo, xhilohilo, xlolohilo,
            xhihilolo, xlohilolo, xhilololo, xlolololo,
       x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15,
       x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29,
       x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43,
       x44, x45, x46, x47, x48, x49, x50, x51, x52, x53, x54, x55, x56, x57,
       x58, x59, x60, x61, x62, x63);
   quarter_hd_vector
      (dim, yhihihihi, ylohihihi, yhilohihi, ylolohihi,
            yhihilohi, ylohilohi, yhilolohi, ylololohi,
            yhihihilo, ylohihilo, yhilohilo, ylolohilo,
            yhihilolo, ylohilolo, yhilololo, ylolololo,
       y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15,
       y16, y17, y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29,
       y30, y31, y32, y33, y34, y35, y36, y37, y38, y39, y40, y41, y42, y43,
       y44, y45, y46, y47, y48, y49, y50, y51, y52, y53, y54, y55, y56, y57,
       y58, y59, y60, y61, y62, y63);

   vectored_hd_product
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,
       x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28,
       x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42,
       x43, x44, x45, x46, x47, x48, x49, x50, x51, x52, x53, x54, x55, x56,
       x57, x58, x59, x60, x61, x62, x63, y0, y1, y2, y3, y4, y5, y6, y7, y8,
       y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22,
       y23, y24, y25, y26, y27, y28, y29, y30, y31, y32, y33, y34, y35, y36,
       y37, y38, y39, y40, y41, y42, y43, y44, y45, y46, y47, y48, y49, y50,
       y51, y52, y53, y54, y55, y56, y57, y58, y59, y60, y61, y62, y63,
       &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &s10, &s11, &s12,
       &s13, &s14, &s15, &s16, &s17, &s18, &s19, &s20, &s21, &s22, &s23, &s24,
       &s25, &s26, &s27, &s28, &s29, &s30, &s31, &s32, &s33, &s34, &s35, &s36,
       &s37, &s38, &s39, &s40, &s41, &s42, &s43, &s44, &s45, &s46, &s47, &s48,
       &s49, &s50, &s51, &s52, &s53, &s54, &s55, &s56, &s57, &s58, &s59, &s60,
       &s61, &s62, &s63);

   to_hexa_double
      (s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15,
       s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29,
       s30, s31, s32, s33, s34, s35, s36, s37, s38, s39, s40, s41, s42, s43,
       s44, s45, s46, s47, s48, s49, s50, s51, s52, s53, s54, s55, s56, s57,
       s58, s59, s60, s61, s62, s63, &vpd[0], &vpd[1], &vpd[2], &vpd[3],
       &vpd[4], &vpd[5], &vpd[6], &vpd[7], &vpd[8], &vpd[9], &vpd[10],
       &vpd[11], &vpd[12], &vpd[13], &vpd[14], &vpd[15]);
 
   cout << "hd x*y :" << endl; hd_write_doubles(prd); cout << endl;
   cout << "vd x*y :" << endl; hd_write_doubles(vpd); cout << endl;

   hdf_sub(prd[0], prd[1], prd[2], prd[3], prd[4], prd[5], prd[6], prd[7],
      prd[8], prd[9], prd[10], prd[11], prd[12], prd[13], prd[14], prd[15],
      vpd[0], vpd[1], vpd[2], vpd[3], vpd[4], vpd[5], vpd[6], vpd[7],
      vpd[8], vpd[9], vpd[10], vpd[11], vpd[12], vpd[13], vpd[14], vpd[15],
      &err[0], &err[1], &err[2], &err[3], &err[4], &err[5], &err[6], &err[7],
      &err[8], &err[9], &err[10], &err[11], &err[12], &err[13], &err[14],
      &err[15]);

   if(err[0] < 0.0)
      hdf_minus(&err[0], &err[1], &err[2], &err[3],
                &err[4], &err[5], &err[6], &err[7],
                &err[8], &err[9], &err[10], &err[11],
                &err[12], &err[13], &err[14], &err[15]);

   cout << " error :" << endl; hd_write_doubles(err); cout << endl;

   fail = (abs(err[0]) > 1.0E-250);

   return fail;
}
