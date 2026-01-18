/* Collection of functions for vectored quad double arithmetic. */

#include <stdio.h>
#include "quad_double.h"
#include "quad_double_functions.h"
#include "splitting_doubles.h"

void quarter_quad_double
 ( double xhihi, double xlohi, double xhilo, double xlolo,
   double *xhihi0, double *xhihi1, double *xhihi2, double *xhihi3,
   double *xlohi0, double *xlohi1, double *xlohi2, double *xlohi3,
   double *xhilo0, double *xhilo1, double *xhilo2, double *xhilo3,
   double *xlolo0, double *xlolo1, double *xlolo2, double *xlolo3 )
{
   quarter_split(xhihi, xhihi0, xhihi1, xhihi2, xhihi3);
   quarter_split(xlohi, xlohi0, xlohi1, xlohi2, xlohi3);
   quarter_split(xhilo, xhilo0, xhilo1, xhilo2, xhilo3);
   quarter_split(xlolo, xlolo0, xlolo1, xlolo2, xlolo3);
}

void quarter_qd_vector
 ( int dim, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *xhihi0, double *xhihi1, double *xhihi2, double *xhihi3,
   double *xlohi0, double *xlohi1, double *xlohi2, double *xlohi3,
   double *xhilo0, double *xhilo1, double *xhilo2, double *xhilo3,
   double *xlolo0, double *xlolo1, double *xlolo2, double *xlolo3 )
{
   for(int i=0; i<dim; i++)
   {
      quarter_quad_double
         (xhihi[i], xlohi[i], xhilo[i], xlolo[i],
          &xhihi0[i], &xhihi1[i], &xhihi2[i], &xhihi3[i],
          &xlohi0[i], &xlohi1[i], &xlohi2[i], &xlohi3[i],
          &xhilo0[i], &xhilo1[i], &xhilo2[i], &xhilo3[i],
          &xlolo0[i], &xlolo1[i], &xlolo2[i], &xlolo3[i]);
   }
}

void to_quad_double
 ( double xhihi0, double xhihi1, double xhihi2, double xhihi3,
   double xlohi0, double xlohi1, double xlohi2, double xlohi3,
   double xhilo0, double xhilo1, double xhilo2, double xhilo3,
   double xlolo0, double xlolo1, double xlolo2, double xlolo3,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   *xhihi = xhihi0;
   *xlohi = 0.0;
   *xhilo = 0.0;
   *xlolo = 0.0;

   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhihi1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhihi2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhihi3);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi0);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi3);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo0);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo3);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo0);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo3);
}

void qd_write_vector
 ( int dim, double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   double x[4];

   for(int i=0; i<dim; i++)
   {
      x[0] = xhihi[i]; x[1] = xlohi[i];
      x[2] = xhilo[i]; x[3] = xlolo[i];
      qd_write(x, 64); printf("\n");
   }
}

void quad_double_product
 ( int dim, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
            double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *prdhihi, double *prdlohi, double *prdhilo, double *prdlolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   *prdhihi = 0.0; *prdlohi = 0.0;
   *prdhilo = 0.0; *prdlolo = 0.0;

   for(int i=0; i<dim; i++)
   {
      qdf_mul(xhihi[i], xlohi[i], xhilo[i], xlolo[i],
              yhihi[i], ylohi[i], yhilo[i], ylolo[i],
              &acchihi, &acclohi, &acchilo, &acclolo);
      qdf_inc(prdhihi, prdlohi, prdhilo, prdlolo,
              acchihi, acclohi, acchilo, acclolo);
   }
}

void quad_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo,
   double **Chihi, double **Clohi, double **Chilo, double **Clolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chihi[i][j] = 0.0; Clohi[i][j] = 0.0;
         Chilo[i][j] = 0.0; Clolo[i][j] = 0.0;

         for(int k=0; k<dim; k++)
         {
            qdf_mul(Ahihi[i][k], Alohi[i][k], Ahilo[i][k], Alolo[i][k],
                    Bhihi[k][j], Blohi[k][j], Bhilo[k][j], Blolo[k][j],
                    &acchihi, &acclohi, &acchilo, &acclolo);
            qdf_inc(&Chihi[i][j], &Clohi[i][j], &Chilo[i][j], &Clolo[i][j],
                    acchihi, acclohi, acchilo, acclolo);
         }
      }
}

void vectored_qd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *xA, double *xB,
   double *xC, double *xD, double *xE, double *xF,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *yA, double *yB,
   double *yC, double *yD, double *yE, double *yF,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *sA, double *sB,
   double *sC, double *sD, double *sE, double *sF )
{
   *s0 = 0.0; *s1 = 0.0; *s2 = 0.0; *s3 = 0.0;
   *s4 = 0.0; *s5 = 0.0; *s6 = 0.0; *s7 = 0.0;
   *s8 = 0.0; *s9 = 0.0; *sA = 0.0; *sB = 0.0;
   *sC = 0.0; *sD = 0.0; *sE = 0.0; *sF = 0.0;

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
      *s8 += x0[i]*y8[i] + x1[i]*y7[i] + x2[i]*y6[i] + x3[i]*y5[i]
           + x4[i]*y4[i] + x5[i]*y3[i] + x6[i]*y2[i] + x7[i]*y1[i]
           + x8[i]*y0[i];
      *s9 += x0[i]*y9[i] + x1[i]*y8[i] + x2[i]*y7[i] + x3[i]*y6[i]
           + x4[i]*y5[i] + x5[i]*y4[i] + x6[i]*y3[i] + x7[i]*y2[i]
           + x8[i]*y1[i] + x9[i]*y0[i];
      *sA += x0[i]*yA[i] + x1[i]*y9[i] + x2[i]*y8[i] + x3[i]*y7[i]
           + x4[i]*y6[i] + x5[i]*y5[i] + x6[i]*y4[i] + x7[i]*y3[i]
           + x8[i]*y2[i] + x9[i]*y1[i] + xA[i]*y0[i];
      *sB += x0[i]*yB[i] + x1[i]*yA[i] + x2[i]*y9[i] + x3[i]*y8[i]
           + x4[i]*y7[i] + x5[i]*y6[i] + x6[i]*y5[i] + x7[i]*y4[i]
           + x8[i]*y3[i] + x9[i]*y2[i] + xA[i]*y1[i] + xB[i]*y0[i];
      *sC += x0[i]*yC[i] + x1[i]*yB[i] + x2[i]*yA[i] + x3[i]*y9[i]
           + x4[i]*y8[i] + x5[i]*y7[i] + x6[i]*y6[i] + x7[i]*y5[i]
           + x8[i]*y4[i] + x9[i]*y3[i] + xA[i]*y2[i] + xB[i]*y1[i]
           + xC[i]*y0[i];
      *sD += x0[i]*yD[i] + x1[i]*yC[i] + x2[i]*yB[i] + x3[i]*yA[i]
           + x4[i]*y9[i] + x5[i]*y8[i] + x6[i]*y7[i] + x7[i]*y6[i]
           + x8[i]*y5[i] + x9[i]*y4[i] + xA[i]*y3[i] + xB[i]*y2[i]
           + xC[i]*y1[i] + xD[i]*y0[i];
      *sE += x0[i]*yE[i] + x1[i]*yD[i] + x2[i]*yC[i] + x3[i]*yB[i]
           + x4[i]*yA[i] + x5[i]*y9[i] + x6[i]*y8[i] + x7[i]*y7[i]
           + x8[i]*y6[i] + x9[i]*y5[i] + xA[i]*y4[i] + xB[i]*y3[i]
           + xC[i]*y2[i] + xD[i]*y1[i] + xE[i]*y0[i];
      *sF += x0[i]*yF[i] + x1[i]*yE[i] + x2[i]*yD[i] + x3[i]*yC[i]
           + x4[i]*yB[i] + x5[i]*yA[i] + x6[i]*y9[i] + x7[i]*y8[i]
           + x8[i]*y7[i] + x9[i]*y6[i] + xA[i]*y5[i] + xB[i]*y4[i]
           + xC[i]*y3[i] + xD[i]*y2[i] + xE[i]*y1[i] + xF[i]*y0[i];
   }
}
