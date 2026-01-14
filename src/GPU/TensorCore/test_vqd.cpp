/* Tests the collection of functions on vectored quad double arithmetic.  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "quad_double.h"
#include "random4_vectors.h"
#include "quad_double_functions.h"
#include "vectored_quad_doubles.h"

int test_quarter_quad_double ( void );
/*
 * Generates a random positive quad double, quarters the parts,
 * and then checks if their sum equals the original quad double. */

int main ( void )
{
   int fail,dim;

   fail = test_quarter_quad_double();

   if(fail == 1)
      printf("\nTest on quarter quad double failed?!!!\n\n");
   else
      printf("\nTest on quarter quad double succeeded.\n\n");

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

   srand(time(NULL));

   random_quad_double(&x[0], &x[1], &x[2], &x[3]);

   if(x[0] < 0.0) x[0] = -x[0]; // all parts must be positive
   if(x[1] < 0.0) x[1] = -x[1];
   if(x[2] < 0.0) x[2] = -x[2];
   if(x[3] < 0.0) x[3] = -x[3];

   printf("x : "); qd_write(x, 64); printf("\n");

   quarter_quad_double
      (x[0], x[1], x[2], x[3],
       &xhihi0, &xhihi1, &xhihi2, &xhihi3,
       &xlohi0, &xlohi1, &xlohi2, &xlohi3,
       &xhilo0, &xhilo1, &xhilo2, &xhilo3,
       &xlolo0, &xlolo1, &xlolo2, &xlolo3);

   printf("xhihi0 : %.15e\n", xhihi0);
   printf("xhihi1 : %.15e\n", xhihi1);
   printf("xhihi2 : %.15e\n", xhihi2);
   printf("xhihi3 : %.15e\n", xhihi3);
   printf("xlohi0 : %.15e\n", xlohi0);
   printf("xlohi1 : %.15e\n", xlohi1);
   printf("xlohi2 : %.15e\n", xlohi2);
   printf("xlohi3 : %.15e\n", xlohi3);
   printf("xhilo0 : %.15e\n", xhilo0);
   printf("xhilo1 : %.15e\n", xhilo1);
   printf("xhilo2 : %.15e\n", xhilo2);
   printf("xhilo3 : %.15e\n", xhilo3);
   printf("xlolo0 : %.15e\n", xlolo0);
   printf("xlolo1 : %.15e\n", xlolo1);
   printf("xlolo2 : %.15e\n", xlolo2);
   printf("xlolo3 : %.15e\n", xlolo3);

   to_quad_double
      (xhihi0, xhihi1, xhihi2, xhihi3, xlohi0, xlohi1, xlohi2, xlohi3,
       xhilo0, xhilo1, xhilo2, xhilo3, xlolo0, xlolo1, xlolo2, xlolo3,
       &y[0], &y[1], &y[2], &y[3]);

   printf("y : "); qd_write(y, 64); printf("\n");

   qdf_sub(x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3],
           &e[0], &e[1], &e[2], &e[3]);

   printf("e : "); qd_write(e, 64); printf("\n");

   fail = not(e[0] == 0.0)
        + not(e[1] == 0.0)
        + not(e[2] == 0.0)
        + not(e[3] == 0.0);

   return fail;
}
