/* Tests the collection of functions on vectored octo double arithmetic.  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "octo_double.h"
#include "random8_vectors.h"
#include "octo_double_functions.h"
#include "vectored_octo_doubles.h"

int test_quarter_octo_double ( void );
/*
 * Generates a random positive octo double, quarters the parts,
 * and then checks if their sum equals the original quad double.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_od_product ( int dim );
/*
 * Generates two random vectors of octo doubles of size dim,
 * and compares their inner producted with the vectored inner product.
 * Returns 1 if the test failed, returns 0 otherwise. */

int main ( void )
{
   int fail,dim;

   fail = test_quarter_octo_double();

   if(fail == 1)
      printf("\nTest on quarter octo double failed?!!!\n\n");
   else
      printf("\nTest on quarter octo double succeeded.\n\n");

   printf("Give the dimension : "); scanf("%d", &dim);

   fail = test_vectored_od_product(dim);

   if(fail == 1)
      printf("\nTest on vectored octo double product failed?!!!\n\n");
   else
      printf("\nTest on vectored octo double product succeeded.\n\n");

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

   srand(time(NULL));

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

   printf("x :\n"); od_write_doubles(x); printf("\n");

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

   printf("xhihihi0 : %.15e\n", xhihihi0);
   printf("xhihihi1 : %.15e\n", xhihihi1);
   printf("xhihihi2 : %.15e\n", xhihihi2);
   printf("xhihihi3 : %.15e\n", xhihihi3);
   printf("xlohihi0 : %.15e\n", xlohihi0);
   printf("xlohihi1 : %.15e\n", xlohihi1);
   printf("xlohihi2 : %.15e\n", xlohihi2);
   printf("xlohihi3 : %.15e\n", xlohihi3);
   printf("xhilohi0 : %.15e\n", xhilohi0);
   printf("xhilohi1 : %.15e\n", xhilohi1);
   printf("xhilohi2 : %.15e\n", xhilohi2);
   printf("xhilohi3 : %.15e\n", xhilohi3);
   printf("xlolohi0 : %.15e\n", xlolohi0);
   printf("xlolohi1 : %.15e\n", xlolohi1);
   printf("xlolohi2 : %.15e\n", xlolohi2);
   printf("xlolohi3 : %.15e\n", xlolohi3);

   printf("xhihilo0 : %.15e\n", xhihilo0);
   printf("xhihilo1 : %.15e\n", xhihilo1);
   printf("xhihilo2 : %.15e\n", xhihilo2);
   printf("xhihilo3 : %.15e\n", xhihilo3);
   printf("xlohilo0 : %.15e\n", xlohilo0);
   printf("xlohilo1 : %.15e\n", xlohilo1);
   printf("xlohilo2 : %.15e\n", xlohilo2);
   printf("xlohilo3 : %.15e\n", xlohilo3);
   printf("xhilolo0 : %.15e\n", xhilolo0);
   printf("xhilolo1 : %.15e\n", xhilolo1);
   printf("xhilolo2 : %.15e\n", xhilolo2);
   printf("xhilolo3 : %.15e\n", xhilolo3);
   printf("xlololo0 : %.15e\n", xlololo0);
   printf("xlololo1 : %.15e\n", xlololo1);
   printf("xlololo2 : %.15e\n", xlololo2);
   printf("xlololo3 : %.15e\n", xlololo3);

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

   printf("x :\n"); od_write_doubles(x); printf("\n");
   printf("y :\n"); od_write_doubles(y); printf("\n");

   odf_sub(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
           y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
           &e[0], &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7]);

   printf("e :\n"); od_write_doubles(e); printf("\n");

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
   printf("octo double vector x :\n");
   od_write_vector
      (dim, xhihihi, xlohihi, xhilohi, xlolohi,
            xhihilo, xlohilo, xhilolo, xlololo);
   printf("octo double vector y :\n");
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

   printf("od x*y :\n"); od_write_doubles(prd); printf("\n");
   printf("vd x*y :\n"); od_write_doubles(vpd); printf("\n");

   odf_sub(prd[0], prd[1], prd[2], prd[3], prd[4], prd[5], prd[6], prd[7],
           vpd[0], vpd[1], vpd[2], vpd[3], vpd[4], vpd[5], vpd[6], vpd[7],
           &err[0], &err[1], &err[2], &err[3],
           &err[4], &err[5], &err[6], &err[7]);

   if(err[0] < 0.0)
      odf_minus(&err[0], &err[1], &err[2], &err[3],
                &err[4], &err[5], &err[6], &err[7]);

   printf(" error :\n"); od_write_doubles(err); printf("\n");

   fail = (abs(err[0]) > 1.0E-120);

   return fail;
}
