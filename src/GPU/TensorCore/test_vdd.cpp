/* Tests the collection of functions to split double doubles. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "double_double.h"
#include "random2_vectors.h"
#include "double_double_functions.h"
#include "vectored_double_doubles.h"

int test_quarter_double_double ( void );
/*
 * Generates a random positive double double, quarters the high and low
 * parts, and then checks if their sum equals the original double double. */

int test_vectored_dd_product ( int dim );
/*
 * Generates two random vectors of double doubles of size dim,
 * and compares their inner producted with the vectored inner product. */

int main ( void )
{
   int fail,dim;

   fail = test_quarter_double_double();

   if(fail == 1)
      printf("\nTest on quarter double double failed?!!!\n\n");
   else
      printf("\nTest on quarter double double succeeded.\n\n");

   printf("Give the dimension : "); scanf("%d", &dim);

   fail = test_vectored_dd_product(dim);

   if(fail == 1)
      printf("\nTest on vectored double double product failed?!!!\n\n");
   else
      printf("\nTest on vectored double double product succeeded.\n\n");

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

   srand(time(NULL));

   random_double_double(&x[0], &x[1]);

   if(x[0] < 0.0) x[0] = -x[0]; // both high and low part
   if(x[1] < 0.0) x[1] = -x[1]; // must be positive

   printf("x : "); dd_write(x, 32); printf("\n");

   quarter_double_double
      (x[0], x[1], &xhi0, &xhi1, &xhi2, &xhi3, &xlo0, &xlo1, &xlo2, &xlo3);

   printf("xhi0 : %.15e\n", xhi0);
   printf("xhi1 : %.15e\n", xhi1);
   printf("xhi2 : %.15e\n", xhi2);
   printf("xhi3 : %.15e\n", xhi3);
   printf("xlo0 : %.15e\n", xlo0);
   printf("xlo1 : %.15e\n", xlo1);
   printf("xlo2 : %.15e\n", xlo2);
   printf("xlo3 : %.15e\n", xlo3);

   to_double_double
      (xhi0, xhi1, xhi2, xhi3, xlo0, xlo1, xlo2, xlo3, &y[0], &y[1]);

   printf("y : "); dd_write(y, 32); printf("\n");

   ddf_sub(x[0], x[1], y[0], y[1], &e[0], &e[1]);

   printf("e : "); dd_write(e, 32); printf("\n");

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
   printf("double double vector x :\n"); dd_write_vector(dim, xhi, xlo);
   printf("double double vector y :\n"); dd_write_vector(dim, yhi, ylo);

   double_double_product(dim, xhi, xlo, yhi, ylo, &prd[0], &prd[1]);

   quarter_dd_vector(dim, xhi, xlo, x0, x1, x2, x3, x4, x5, x6, x7);
   quarter_dd_vector(dim, yhi, ylo, y0, y1, y2, y3, y4, y5, y6, y7);

   vectored_dd_product
      (dim, x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7,
       &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7);

   to_double_double(s0, s1, s2, s3, s4, s5, s6, s7, &vpd[0], &vpd[1]);

   printf("dd x*y : "); dd_write(prd, 32); printf("\n");
   printf("vd x*y : "); dd_write(vpd, 32); printf("\n");

   ddf_sub(prd[0], prd[1], vpd[0], vpd[1], &err[0], &err[1]);

   if(err[0] < 0.0) ddf_minus(&err[0], &err[1]);

   printf(" error : "); dd_write(err, 32); printf("\n");

   fail = (abs(err[0]) > 1.0E-28);

   return fail;
}
