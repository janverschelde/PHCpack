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

void quarter_dd_matrix
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         quarter_double_double
            (Ahi[i][j], Alo[i][j],
             &Ahi0[i][j], &Ahi1[i][j], &Ahi2[i][j], &Ahi3[i][j],
             &Alo0[i][j], &Alo1[i][j], &Alo2[i][j], &Alo3[i][j]);
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

void to_double_double_matrix
 ( int nrows, int ncols,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3,
   double **Ahi, double **Alo )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         to_double_double
            (Ahi0[i][j], Ahi1[i][j], Ahi2[i][j], Ahi3[i][j],
             Alo0[i][j], Alo1[i][j], Alo2[i][j], Alo3[i][j],
             &Ahi[i][j], &Alo[i][j]);
      }
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

void transpose_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         T0[j][i] = A0[i][j]; T1[j][i] = A1[i][j];
         T2[j][i] = A2[i][j]; T3[j][i] = A3[i][j];
         T4[j][i] = A4[i][j]; T5[j][i] = A5[i][j];
         T6[j][i] = A6[i][j]; T7[j][i] = A7[i][j];
      }
}

void vectored_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         vectored_dd_product
            (dim, A0[i], A1[i], A2[i], A3[i], A4[i], A5[i], A6[i], A7[i],
                  B0[j], B1[j], B2[j], B3[j], B4[j], B5[j], B6[j], B7[j],
             &C0[i][j], &C1[i][j], &C2[i][j], &C3[i][j],
             &C4[i][j], &C5[i][j], &C6[i][j], &C7[i][j]);
      }
}
