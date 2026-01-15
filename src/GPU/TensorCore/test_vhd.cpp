/* Tests the collection of functions on vectored hexa double arithmetic.  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "hexa_double.h"
#include "random16_vectors.h"
// #include "hexa_double_functions.h"
#include "vectored_hexa_doubles.h"

int test_quarter_hexa_double ( void );
/*
 * Generates a random positive hexa double, quarters the parts,
 * and then checks if their sum equals the original hexa double.
 * Returns 1 if the test failed, returns 0 otherwise. */

int main ( void )
{
   int fail,dim;

   fail = test_quarter_hexa_double();

   if(fail == 1)
      printf("\nTest on quarter hexa double failed?!!!\n\n");
   else
      printf("\nTest on quarter hexa double succeeded.\n\n");

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
   double y[8];
   double e[8];

   srand(time(NULL));

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

   printf("x :\n"); hd_write_doubles(x); printf("\n");

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

   printf("xhihihihi0 : %.15e\n", xhihihihi0);
   printf("xhihihihi1 : %.15e\n", xhihihihi1);
   printf("xhihihihi2 : %.15e\n", xhihihihi2);
   printf("xhihihihi3 : %.15e\n", xhihihihi3);
   printf("xlohihihi0 : %.15e\n", xlohihihi0);
   printf("xlohihihi1 : %.15e\n", xlohihihi1);
   printf("xlohihihi2 : %.15e\n", xlohihihi2);
   printf("xlohihihi3 : %.15e\n", xlohihihi3);
   printf("xhilohihi0 : %.15e\n", xhilohihi0);
   printf("xhilohihi1 : %.15e\n", xhilohihi1);
   printf("xhilohihi2 : %.15e\n", xhilohihi2);
   printf("xhilohihi3 : %.15e\n", xhilohihi3);
   printf("xlolohihi0 : %.15e\n", xlolohihi0);
   printf("xlolohihi1 : %.15e\n", xlolohihi1);
   printf("xlolohihi2 : %.15e\n", xlolohihi2);
   printf("xlolohihi3 : %.15e\n", xlolohihi3);

   printf("xhihilohi0 : %.15e\n", xhihilohi0);
   printf("xhihilohi1 : %.15e\n", xhihilohi1);
   printf("xhihilohi2 : %.15e\n", xhihilohi2);
   printf("xhihilohi3 : %.15e\n", xhihilohi3);
   printf("xlohilohi0 : %.15e\n", xlohilohi0);
   printf("xlohilohi1 : %.15e\n", xlohilohi1);
   printf("xlohilohi2 : %.15e\n", xlohilohi2);
   printf("xlohilohi3 : %.15e\n", xlohilohi3);
   printf("xhilolohi0 : %.15e\n", xhilolohi0);
   printf("xhilolohi1 : %.15e\n", xhilolohi1);
   printf("xhilolohi2 : %.15e\n", xhilolohi2);
   printf("xhilolohi3 : %.15e\n", xhilolohi3);
   printf("xlololohi0 : %.15e\n", xlololohi0);
   printf("xlololohi1 : %.15e\n", xlololohi1);
   printf("xlololohi2 : %.15e\n", xlololohi2);
   printf("xlololohi3 : %.15e\n", xlololohi3);

   printf("xhihihilo0 : %.15e\n", xhihihilo0);
   printf("xhihihilo1 : %.15e\n", xhihihilo1);
   printf("xhihihilo2 : %.15e\n", xhihihilo2);
   printf("xhihihilo3 : %.15e\n", xhihihilo3);
   printf("xlohihilo0 : %.15e\n", xlohihilo0);
   printf("xlohihilo1 : %.15e\n", xlohihilo1);
   printf("xlohihilo2 : %.15e\n", xlohihilo2);
   printf("xlohihilo3 : %.15e\n", xlohihilo3);
   printf("xhilohilo0 : %.15e\n", xhilohilo0);
   printf("xhilohilo1 : %.15e\n", xhilohilo1);
   printf("xhilohilo2 : %.15e\n", xhilohilo2);
   printf("xhilohilo3 : %.15e\n", xhilohilo3);
   printf("xlolohilo0 : %.15e\n", xlolohilo0);
   printf("xlolohilo1 : %.15e\n", xlolohilo1);
   printf("xlolohilo2 : %.15e\n", xlolohilo2);
   printf("xlolohilo3 : %.15e\n", xlolohilo3);

   printf("xhihilolo0 : %.15e\n", xhihilolo0);
   printf("xhihilolo1 : %.15e\n", xhihilolo1);
   printf("xhihilolo2 : %.15e\n", xhihilolo2);
   printf("xhihilolo3 : %.15e\n", xhihilolo3);
   printf("xlohilolo0 : %.15e\n", xlohilolo0);
   printf("xlohilolo1 : %.15e\n", xlohilolo1);
   printf("xlohilolo2 : %.15e\n", xlohilolo2);
   printf("xlohilolo3 : %.15e\n", xlohilolo3);
   printf("xhilololo0 : %.15e\n", xhilololo0);
   printf("xhilololo1 : %.15e\n", xhilololo1);
   printf("xhilololo2 : %.15e\n", xhilololo2);
   printf("xhilololo3 : %.15e\n", xhilololo3);
   printf("xlolololo0 : %.15e\n", xlolololo0);
   printf("xlolololo1 : %.15e\n", xlolololo1);
   printf("xlolololo2 : %.15e\n", xlolololo2);
   printf("xlolololo3 : %.15e\n", xlolololo3);

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

   printf("x :\n"); hd_write_doubles(x); printf("\n");
   printf("y :\n"); hd_write_doubles(y); printf("\n");

   e[0] = x[0] - y[0];
   e[1] = x[1] - y[1];
   e[2] = x[2] - y[2];
   e[3] = x[3] - y[3];
   e[4] = x[4] - y[4];
   e[5] = x[5] - y[5];
   e[6] = x[6] - y[6];
   e[7] = x[7] - y[7];
   e[8] = x[8] - y[8];
   e[9] = x[9] - y[9];
   e[10] = x[10] - y[10];
   e[11] = x[11] - y[11];
   e[12] = x[12] - y[12];
   e[13] = x[13] - y[13];
   e[14] = x[14] - y[14];
   e[15] = x[15] - y[15];

   printf("e :\n"); hd_write_doubles(e); printf("\n");

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
