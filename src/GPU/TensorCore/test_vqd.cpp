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
 * and then checks if their sum equals the original quad double.
 * Returns 1 if the test failed, returns 0 otherwise. */

int test_vectored_qd_product ( int dim );
/*
 * Generates two random vectors of quad doubles of size dim,
 * and compares their inner producted with the vectored inner product.
 * Returns 1 if the test failed, returns 0 otherwise. */

int main ( void )
{
   int fail,dim;

   fail = test_quarter_quad_double();

   if(fail == 1)
      printf("\nTest on quarter quad double failed?!!!\n\n");
   else
      printf("\nTest on quarter quad double succeeded.\n\n");

   printf("Give the dimension : "); scanf("%d", &dim);

   fail = test_vectored_qd_product(dim);

   if(fail == 1)
      printf("\nTest on vectored quad double product failed?!!!\n\n");
   else
      printf("\nTest on vectored quad double product succeeded.\n\n");

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
   printf("quad double vector x :\n");
   qd_write_vector(dim, xhihi, xlohi, xhilo, xlolo);
   printf("quad double vector y :\n");
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

   printf("qd x*y : "); qd_write(prd, 64); printf("\n");
   printf("vd x*y : "); qd_write(vpd, 64); printf("\n");

   qdf_sub(prd[0], prd[1], prd[2], prd[3], vpd[0], vpd[1], vpd[2], vpd[3],
           &err[0], &err[1], &err[2], &err[3]);

   if(err[0] < 0.0) qdf_minus(&err[0], &err[1], &err[2], &err[3]);

   printf(" error : "); qd_write(err, 64); printf("\n");

   fail = (abs(err[0]) > 1.0E-58);

   return fail;
}
