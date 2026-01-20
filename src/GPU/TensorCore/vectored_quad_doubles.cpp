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

void quarter_qd_matrix
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Ahihi0, double **Ahihi1, double **Ahihi2, double **Ahihi3,
   double **Alohi0, double **Alohi1, double **Alohi2, double **Alohi3,
   double **Ahilo0, double **Ahilo1, double **Ahilo2, double **Ahilo3,
   double **Alolo0, double **Alolo1, double **Alolo2, double **Alolo3 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         quarter_quad_double
            (Ahihi[i][j], Alohi[i][j], Ahilo[i][j], Alolo[i][j],
             &Ahihi0[i][j], &Ahihi1[i][j], &Ahihi2[i][j], &Ahihi3[i][j],
             &Alohi0[i][j], &Alohi1[i][j], &Alohi2[i][j], &Alohi3[i][j],
             &Ahilo0[i][j], &Ahilo1[i][j], &Ahilo2[i][j], &Ahilo3[i][j],
             &Alolo0[i][j], &Alolo1[i][j], &Alolo2[i][j], &Alolo3[i][j]);
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

void to_quad_double_matrix
 ( int nrows, int ncols,
   double **Ahihi0, double **Ahihi1, double **Ahihi2, double **Ahihi3,
   double **Alohi0, double **Alohi1, double **Alohi2, double **Alohi3,
   double **Ahilo0, double **Ahilo1, double **Ahilo2, double **Ahilo3,
   double **Alolo0, double **Alolo1, double **Alolo2, double **Alolo3,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         to_quad_double
            (Ahihi0[i][j], Ahihi1[i][j], Ahihi2[i][j], Ahihi3[i][j],
             Alohi0[i][j], Alohi1[i][j], Alohi2[i][j], Alohi3[i][j],
             Ahilo0[i][j], Ahilo1[i][j], Ahilo2[i][j], Ahilo3[i][j],
             Alolo0[i][j], Alolo1[i][j], Alolo2[i][j], Alolo3[i][j],
             &Ahihi[i][j], &Alohi[i][j], &Ahilo[i][j], &Alolo[i][j]);
      }
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

void transpose_qd_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7,
   double **T8, double **T9, double **T10, double **T11,
   double **T12, double **T13, double **T14, double **T15 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         T0[j][i] = A0[i][j]; T1[j][i] = A1[i][j];
         T2[j][i] = A2[i][j]; T3[j][i] = A3[i][j];
         T4[j][i] = A4[i][j]; T5[j][i] = A5[i][j];
         T6[j][i] = A6[i][j]; T7[j][i] = A7[i][j];
         T8[j][i] = A8[i][j]; T9[j][i] = A9[i][j];
         T10[j][i] = A10[i][j]; T11[j][i] = A11[i][j];
         T12[j][i] = A12[i][j]; T13[j][i] = A13[i][j];
         T14[j][i] = A14[i][j]; T15[j][i] = A15[i][j];
      }
}

void vectored_qd_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **B8, double **B9, double **B10, double **B11,
   double **B12, double **B13, double **B14, double **B15,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7,
   double **C8, double **C9, double **C10, double **C11,
   double **C12, double **C13, double **C14, double **C15 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         vectored_qd_product
            (dim, A0[i], A1[i], A2[i], A3[i], A4[i], A5[i], A6[i], A7[i],
             A8[i], A9[i], A10[i], A11[i], A12[i], A13[i], A14[i], A15[i],
             B0[j], B1[j], B2[j], B3[j], B4[j], B5[j], B6[j], B7[j],
             B8[j], B9[j], B10[j], B11[j], B12[j], B13[j], B14[j], B15[j],
             &C0[i][j], &C1[i][j], &C2[i][j], &C3[i][j],
             &C4[i][j], &C5[i][j], &C6[i][j], &C7[i][j],
             &C8[i][j], &C9[i][j], &C10[i][j], &C11[i][j],
             &C12[i][j], &C13[i][j], &C14[i][j], &C15[i][j]);
      }
}

void qd_convolute_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15, double **cA )
{
   for(int i=0; i<16*nrows; i++)
      for(int j=0; j<16*ncols; j++) cA[i][j] = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         for(int k=0; k<16; k++) cA[16*i+k][16*j+k] = A0[i][j];
         for(int k=0; k<15; k++) cA[16*i+k+1][16*j+k] = A1[i][j];
         for(int k=0; k<14; k++) cA[16*i+k+2][16*j+k] = A2[i][j];
         for(int k=0; k<13; k++) cA[16*i+k+3][16*j+k] = A3[i][j];
         for(int k=0; k<12; k++) cA[16*i+k+4][16*j+k] = A4[i][j];
         for(int k=0; k<11; k++) cA[16*i+k+5][16*j+k] = A5[i][j];
         for(int k=0; k<10; k++) cA[16*i+k+6][16*j+k] = A6[i][j];
         for(int k=0; k<9; k++) cA[16*i+k+7][16*j+k] = A7[i][j];
         for(int k=0; k<8; k++) cA[16*i+k+8][16*j+k] = A8[i][j];
         for(int k=0; k<7; k++) cA[16*i+k+9][16*j+k] = A9[i][j];
         for(int k=0; k<6; k++) cA[16*i+k+10][16*j+k] = A10[i][j];
         for(int k=0; k<5; k++) cA[16*i+k+11][16*j+k] = A11[i][j];
         for(int k=0; k<4; k++) cA[16*i+k+12][16*j+k] = A12[i][j];
         for(int k=0; k<3; k++) cA[16*i+k+13][16*j+k] = A13[i][j];
         for(int k=0; k<2; k++) cA[16*i+k+14][16*j+k] = A14[i][j];
         cA[16*i+15][16*j] = A15[i][j];
      }
}

void qd_stack_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15, double **sA )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         sA[16*i][j] = A0[i][j];
         sA[16*i+1][j] = A1[i][j];
         sA[16*i+2][j] = A2[i][j];
         sA[16*i+3][j] = A3[i][j];
         sA[16*i+4][j] = A4[i][j];
         sA[16*i+5][j] = A5[i][j];
         sA[16*i+6][j] = A6[i][j];
         sA[16*i+7][j] = A7[i][j];
         sA[16*i+8][j] = A8[i][j];
         sA[16*i+9][j] = A9[i][j];
         sA[16*i+10][j] = A10[i][j];
         sA[16*i+11][j] = A11[i][j];
         sA[16*i+12][j] = A12[i][j];
         sA[16*i+13][j] = A13[i][j];
         sA[16*i+14][j] = A14[i][j];
         sA[16*i+15][j] = A15[i][j];
      }
}

void extract_qd_quarters
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7,
   double **D8, double **D9, double **D10, double **D11,
   double **D12, double **D13, double **D14, double **D15 )
{
   for(int i=0; i<16*nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         int k = i % 16;
         int row = i/16;

         if(k == 0) D0[row][j] = qC[i][j];
         if(k == 1) D1[row][j] = qC[i][j];
         if(k == 2) D2[row][j] = qC[i][j];
         if(k == 3) D3[row][j] = qC[i][j];
         if(k == 4) D4[row][j] = qC[i][j];
         if(k == 5) D5[row][j] = qC[i][j];
         if(k == 6) D6[row][j] = qC[i][j];
         if(k == 7) D7[row][j] = qC[i][j];
         if(k == 8) D8[row][j] = qC[i][j];
         if(k == 9) D9[row][j] = qC[i][j];
         if(k == 10) D10[row][j] = qC[i][j];
         if(k == 11) D11[row][j] = qC[i][j];
         if(k == 12) D12[row][j] = qC[i][j];
         if(k == 13) D13[row][j] = qC[i][j];
         if(k == 14) D14[row][j] = qC[i][j];
         if(k == 15) D15[row][j] = qC[i][j];
      }
}
