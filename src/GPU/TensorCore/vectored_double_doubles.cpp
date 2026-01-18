/* Collection of functions for vectored double double arithmetic. */

#include <stdio.h>
#include "double_double.h"
#include "double_double_functions.h"
#include "splitting_doubles.h"

void quarter_double_double
 ( double xhi, double xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3 )
{
   quarter_split(xhi, xhi0, xhi1, xhi2, xhi3);
   quarter_split(xlo, xlo0, xlo1, xlo2, xlo3);
}

void quarter_dd_vector
 ( int dim, double *xhi, double *xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3 )
{
   for(int i=0; i<dim; i++)
   {
      quarter_double_double
         (xhi[i], xlo[i],
          &xhi0[i], &xhi1[i], &xhi2[i], &xhi3[i],
          &xlo0[i], &xlo1[i], &xlo2[i], &xlo3[i]);
   }
}

void to_double_double
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3,
   double *xhi, double *xlo )
{
   *xhi = xhi0;
   *xlo = 0.0;

   ddf_inc_d(xhi, xlo, xhi1);
   ddf_inc_d(xhi, xlo, xhi2);
   ddf_inc_d(xhi, xlo, xhi3);
   ddf_inc_d(xhi, xlo, xlo0);
   ddf_inc_d(xhi, xlo, xlo1);
   ddf_inc_d(xhi, xlo, xlo2);
   ddf_inc_d(xhi, xlo, xlo3);
}

void dd_write_vector ( int dim, double *xhi, double *xlo )
{
   double x[2];

   for(int i=0; i<dim; i++)
   {
      x[0] = xhi[i]; x[1] = xlo[i];
      dd_write(x, 32); printf("\n");
   }
}

void double_double_product
 ( int dim, double *xhi, double *xlo, double *yhi, double *ylo,
   double *prdhi, double *prdlo )
{
   double acchi,acclo;

   *prdhi = 0.0;
   *prdlo = 0.0;

   for(int i=0; i<dim; i++)
   {
      ddf_mul(xhi[i], xlo[i], yhi[i], ylo[i], &acchi, &acclo);
      ddf_inc(prdhi, prdlo, acchi, acclo);
   }
}

void double_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo )
{
   double acchi,acclo;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chi[i][j] = 0.0; Clo[i][j] = 0.0;

         for(int k=0; k<dim; k++)
         {
            ddf_mul(Ahi[i][k], Alo[i][k], Bhi[k][j], Blo[k][j],
                    &acchi, &acclo);
            ddf_inc(&Chi[i][j], &Clo[i][j], acchi, acclo);
         }
      }
}

void vectored_dd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7 )
{
   *s0 = 0.0; *s1 = 0.0; *s2 = 0.0; *s3 = 0.0;
   *s4 = 0.0; *s5 = 0.0; *s7 = 0.0; *s7 = 0.0;

   for(int i=0; i<dim; i++)
   {
      *s0 += x0[i]*y0[i];
      *s1 += x0[i]*y1[i] + x1[i]*y0[i];
      *s2 += x0[i]*y2[i] + x1[i]*y1[i] + x2[i]*y0[i];
      *s3 += x0[i]*y3[i] + x1[i]*y2[i] + x2[i]*y1[i] + x3[i]*y0[i];
      *s4 += x0[i]*y4[i] + x1[i]*y3[i] + x2[i]*y2[i] + x3[i]*y1[i]
           + x4[i]*y0[i];
      *s5 += x0[i]*y5[i] + x1[i]*y4[i] + x2[i]*y3[i] + x3[i]*y2[i]
           + x4[i]*y1[i] + x5[i]*y0[i];
      *s6 += x0[i]*y6[i] + x1[i]*y5[i] + x2[i]*y4[i] + x3[i]*y3[i]
           + x4[i]*y2[i] + x5[i]*y1[i] + x6[i]*y0[i];
      *s7 += x0[i]*y7[i] + x1[i]*y6[i] + x2[i]*y5[i] + x3[i]*y4[i]
           + x4[i]*y3[i] + x5[i]*y2[i] + x6[i]*y1[i] + x7[i]*y0[i];
   }
}
