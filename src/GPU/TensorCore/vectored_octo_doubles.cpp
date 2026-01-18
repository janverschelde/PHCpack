/* Defines the functions with prototypes in vectored_octo_doubles.h. */

#include "octo_double.h"
#include "octo_double_functions.h"
#include "splitting_doubles.h"

void quarter_octo_double
 ( double xhihihi, double xlohihi, double xhilohi, double xlolohi,
   double xhihilo, double xlohilo, double xhilolo, double xlololo,
   double *xhihihi0, double *xhihihi1, double *xhihihi2, double *xhihihi3,
   double *xlohihi0, double *xlohihi1, double *xlohihi2, double *xlohihi3,
   double *xhilohi0, double *xhilohi1, double *xhilohi2, double *xhilohi3,
   double *xlolohi0, double *xlolohi1, double *xlolohi2, double *xlolohi3,
   double *xhihilo0, double *xhihilo1, double *xhihilo2, double *xhihilo3,
   double *xlohilo0, double *xlohilo1, double *xlohilo2, double *xlohilo3,
   double *xhilolo0, double *xhilolo1, double *xhilolo2, double *xhilolo3,
   double *xlololo0, double *xlololo1, double *xlololo2, double *xlololo3 )
{
   quarter_split(xhihihi, xhihihi0, xhihihi1, xhihihi2, xhihihi3);
   quarter_split(xlohihi, xlohihi0, xlohihi1, xlohihi2, xlohihi3);
   quarter_split(xhilohi, xhilohi0, xhilohi1, xhilohi2, xhilohi3);
   quarter_split(xlolohi, xlolohi0, xlolohi1, xlolohi2, xlolohi3);
   quarter_split(xhihilo, xhihilo0, xhihilo1, xhihilo2, xhihilo3);
   quarter_split(xlohilo, xlohilo0, xlohilo1, xlohilo2, xlohilo3);
   quarter_split(xhilolo, xhilolo0, xhilolo1, xhilolo2, xhilolo3);
   quarter_split(xlololo, xlololo0, xlololo1, xlololo2, xlololo3);
}

void quarter_od_vector
 ( int dim,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *xhihihi0, double *xhihihi1, double *xhihihi2, double *xhihihi3,
   double *xlohihi0, double *xlohihi1, double *xlohihi2, double *xlohihi3,
   double *xhilohi0, double *xhilohi1, double *xhilohi2, double *xhilohi3,
   double *xlolohi0, double *xlolohi1, double *xlolohi2, double *xlolohi3,
   double *xhihilo0, double *xhihilo1, double *xhihilo2, double *xhihilo3,
   double *xlohilo0, double *xlohilo1, double *xlohilo2, double *xlohilo3,
   double *xhilolo0, double *xhilolo1, double *xhilolo2, double *xhilolo3,
   double *xlololo0, double *xlololo1, double *xlololo2, double *xlololo3 )
{
   for(int i=0; i<dim; i++)
   {
      quarter_octo_double
         (xhihihi[i], xlohihi[i], xhilohi[i], xlolohi[i],
          xhihilo[i], xlohilo[i], xhilolo[i], xlololo[i],
          &xhihihi0[i], &xhihihi1[i], &xhihihi2[i], &xhihihi3[i],
          &xlohihi0[i], &xlohihi1[i], &xlohihi2[i], &xlohihi3[i],
          &xhilohi0[i], &xhilohi1[i], &xhilohi2[i], &xhilohi3[i],
          &xlolohi0[i], &xlolohi1[i], &xlolohi2[i], &xlolohi3[i],
          &xhihilo0[i], &xhihilo1[i], &xhihilo2[i], &xhihilo3[i],
          &xlohilo0[i], &xlohilo1[i], &xlohilo2[i], &xlohilo3[i],
          &xhilolo0[i], &xhilolo1[i], &xhilolo2[i], &xhilolo3[i],
          &xlololo0[i], &xlololo1[i], &xlololo2[i], &xlololo3[i]);
   }
}

void to_octo_double
 ( double xhihihi0, double xhihihi1, double xhihihi2, double xhihihi3,
   double xlohihi0, double xlohihi1, double xlohihi2, double xlohihi3,
   double xhilohi0, double xhilohi1, double xhilohi2, double xhilohi3,
   double xlolohi0, double xlolohi1, double xlolohi2, double xlolohi3,
   double xhihilo0, double xhihilo1, double xhihilo2, double xhihilo3,
   double xlohilo0, double xlohilo1, double xlohilo2, double xlohilo3,
   double xhilolo0, double xhilolo1, double xhilolo2, double xhilolo3,
   double xlololo0, double xlololo1, double xlololo2, double xlololo3,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo )
{
   *xhihihi = xhihihi0;
   *xlohihi = 0.0;
   *xhilohi = 0.0;
   *xlolohi = 0.0;
   *xhihilo = 0.0;
   *xlohilo = 0.0;
   *xhilolo = 0.0;
   *xlololo = 0.0;

   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihihi1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihihi2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihihi3);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohihi0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohihi1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohihi2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohihi3);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilohi0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilohi1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilohi2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilohi3);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlolohi0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlolohi1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlolohi2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlolohi3);

   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihilo0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihilo1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihilo2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhihilo3);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohilo0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohilo1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohilo2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlohilo3);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilolo0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilolo1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilolo2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xhilolo3);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlololo0);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlololo1);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlololo2);
   odf_inc_d(xhihihi, xlohihi, xhilohi, xlolohi,
             xhihilo, xlohilo, xhilolo, xlololo, xlololo3);
}

void od_write_vector
 ( int dim,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo )
{
   double x[8];

   for(int i=0; i<dim; i++)
   {
      x[0] = xhihihi[i]; x[1] = xlohihi[i];
      x[2] = xhilohi[i]; x[3] = xlolohi[i];
      x[4] = xhihilo[i]; x[5] = xlohilo[i];
      x[6] = xhilolo[i]; x[7] = xlololo[i];
      od_write_doubles(x);
   }
}

void octo_double_product
 ( int dim,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *phihihi, double *plohihi, double *philohi, double *plolohi,
   double *phihilo, double *plohilo, double *philolo, double *plololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   *phihihi = 0.0; *plohihi = 0.0; *philohi = 0.0; *plolohi = 0.0;
   *phihilo = 0.0; *plohilo = 0.0; *philolo = 0.0; *plololo = 0.0;

   for(int i=0; i<dim; i++)
   {
      odf_mul(xhihihi[i], xlohihi[i], xhilohi[i], xlolohi[i],
              xhihilo[i], xlohilo[i], xhilolo[i], xlololo[i],
              yhihihi[i], ylohihi[i], yhilohi[i], ylolohi[i],
              yhihilo[i], ylohilo[i], yhilolo[i], ylololo[i],
              &acchihihi, &acclohihi, &acchilohi, &acclolohi,
              &acchihilo, &acclohilo, &acchilolo, &acclololo);
      odf_inc(phihihi, plohihi, philohi, plolohi,
              phihilo, plohilo, philolo, plololo,
              acchihihi, acclohihi, acchilohi, acclolohi,
              acchihilo, acclohilo, acchilolo, acclololo);
   }
}

void octo_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Bhihihi, double **Blohihi, double **Bhilohi, double **Blolohi,
   double **Bhihilo, double **Blohilo, double **Bhilolo, double **Blololo,
   double **Chihihi, double **Clohihi, double **Chilohi, double **Clolohi,
   double **Chihilo, double **Clohilo, double **Chilolo, double **Clololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chihihi[i][j] = 0.0; Clohihi[i][j] = 0.0;
         Chilohi[i][j] = 0.0; Clolohi[i][j] = 0.0;
         Chihilo[i][j] = 0.0; Clohilo[i][j] = 0.0;
         Chilolo[i][j] = 0.0; Clololo[i][j] = 0.0;

         for(int k=0; k<dim; k++)
         {
            odf_mul(Ahihihi[i][k], Alohihi[i][k], Ahilohi[i][k], Alolohi[i][k],
                    Ahihilo[i][k], Alohilo[i][k], Ahilolo[i][k], Alololo[i][k],
                    Bhihihi[k][j], Blohihi[k][j], Bhilohi[k][j], Blolohi[k][j],
                    Bhihilo[k][j], Blohilo[k][j], Bhilolo[k][j], Blololo[k][j],
                    &acchihihi, &acclohihi, &acchilohi, &acclolohi,
                    &acchihilo, &acclohilo, &acchilolo, &acclololo);
            odf_inc(&Chihihi[i][j], &Clohihi[i][j],
                    &Chilohi[i][j], &Clolohi[i][j],
                    &Chihilo[i][j], &Clohilo[i][j],
                    &Chilolo[i][j], &Clololo[i][j],
                    acchihihi, acclohihi, acchilohi, acclolohi,
                    acchihilo, acclohilo, acchilolo, acclololo);
         }
      }
}

void vectored_od_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *x10, double *x11,
   double *x12, double *x13, double *x14, double *x15,
   double *x16, double *x17, double *x18, double *x19,
   double *x20, double *x21, double *x22, double *x23,
   double *x24, double *x25, double *x26, double *x27,
   double *x28, double *x29, double *x30, double *x31,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *y10, double *y11,
   double *y12, double *y13, double *y14, double *y15,
   double *y16, double *y17, double *y18, double *y19,
   double *y20, double *y21, double *y22, double *y23,
   double *y24, double *y25, double *y26, double *y27,
   double *y28, double *y29, double *y30, double *y31,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *s10, double *s11,
   double *s12, double *s13, double *s14, double *s15,
   double *s16, double *s17, double *s18, double *s19,
   double *s20, double *s21, double *s22, double *s23,
   double *s24, double *s25, double *s26, double *s27,
   double *s28, double *s29, double *s30, double *s31 )
{
   *s0 = 0.0; *s1 = 0.0; *s2 = 0.0; *s3 = 0.0;
   *s4 = 0.0; *s5 = 0.0; *s6 = 0.0; *s7 = 0.0;
   *s8 = 0.0; *s9 = 0.0; *s10 = 0.0; *s11 = 0.0;
   *s12 = 0.0; *s13 = 0.0; *s14 = 0.0; *s15 = 0.0;
   *s16 = 0.0; *s17 = 0.0; *s18 = 0.0; *s19 = 0.0;
   *s20 = 0.0; *s21 = 0.0; *s22 = 0.0; *s23 = 0.0;
   *s24 = 0.0; *s25 = 0.0; *s26 = 0.0; *s27 = 0.0;
   *s28 = 0.0; *s29 = 0.0; *s30 = 0.0; *s31 = 0.0;

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
      *s10 += x0[i]*y10[i] + x1[i]*y9[i] +  x2[i]*y8[i] + x3[i]*y7[i]
            + x4[i]*y6[i]  + x5[i]*y5[i] +  x6[i]*y4[i] + x7[i]*y3[i]
            + x8[i]*y2[i]  + x9[i]*y1[i] + x10[i]*y0[i];
      *s11 += x0[i]*y11[i] + x1[i]*y10[i] + x2[i]*y9[i] + x3[i]*y8[i]
            + x4[i]*y7[i]  + x5[i]*y6[i]  + x6[i]*y5[i] + x7[i]*y4[i]
            + x8[i]*y3[i]  + x9[i]*y2[i] + x10[i]*y1[i] + x11[i]*y0[i];
      *s12 += x0[i]*y12[i] + x1[i]*y11[i] + x2[i]*y10[i] + x3[i]*y9[i]
            + x4[i]*y8[i] + x5[i]*y7[i] + x6[i]*y6[i] + x7[i]*y5[i]
            + x8[i]*y4[i] + x9[i]*y3[i] + x10[i]*y2[i] + x11[i]*y1[i]
            + x12[i]*y0[i];
      *s13 += x0[i]*y13[i] + x1[i]*y12[i] + x2[i]*y11[i] + x3[i]*y10[i]
            + x4[i]*y9[i] + x5[i]*y8[i] + x6[i]*y7[i] + x7[i]*y6[i]
            + x8[i]*y5[i] + x9[i]*y4[i] + x10[i]*y3[i] + x11[i]*y2[i]
            + x12[i]*y1[i] + x13[i]*y0[i];
      *s14 += x0[i]*y14[i] + x1[i]*y13[i] + x2[i]*y12[i] + x3[i]*y11[i]
            + x4[i]*y10[i] + x5[i]*y9[i] + x6[i]*y8[i] + x7[i]*y7[i]
            + x8[i]*y6[i] + x9[i]*y5[i] + x10[i]*y4[i] + x11[i]*y3[i]
            + x12[i]*y2[i] + x13[i]*y1[i] + x14[i]*y0[i];
      *s15 += x0[i]*y15[i] + x1[i]*y14[i] + x2[i]*y13[i] + x3[i]*y12[i]
            + x4[i]*y11[i] + x5[i]*y10[i] + x6[i]*y9[i] + x7[i]*y8[i]
            + x8[i]*y7[i] + x9[i]*y6[i] + x10[i]*y5[i] + x11[i]*y4[i]
            + x12[i]*y3[i] + x13[i]*y2[i] + x14[i]*y1[i] + x15[i]*y0[i];
      *s16 += x0[i]*y16[i] + x1[i]*y15[i] + x2[i]*y14[i] + x3[i]*y13[i]
            + x4[i]*y12[i] + x5[i]*y11[i] + x6[i]*y10[i] + x7[i]*y9[i]
            + x8[i]*y8[i] + x9[i]*y7[i] + x10[i]*y6[i] + x11[i]*y5[i]
            + x12[i]*y4[i] + x13[i]*y3[i] + x14[i]*y2[i] + x15[i]*y1[i]
            + x16[i]*y0[i];
      *s17 += x0[i]*y17[i] + x1[i]*y16[i] + x2[i]*y15[i] + x3[i]*y14[i]
            + x4[i]*y13[i] + x5[i]*y12[i] + x6[i]*y11[i] + x7[i]*y10[i]
            + x8[i]*y9[i] + x9[i]*y8[i] + x10[i]*y7[i] + x11[i]*y6[i]
            + x12[i]*y5[i] + x13[i]*y4[i] + x14[i]*y3[i] + x15[i]*y2[i]
            + x16[i]*y1[i] + x17[i]*y0[i];
      *s18 += x0[i]*y18[i] + x1[i]*y17[i] + x2[i]*y16[i] + x3[i]*y15[i]
            + x4[i]*y14[i] + x5[i]*y13[i] + x6[i]*y12[i] + x7[i]*y11[i]
            + x8[i]*y10[i] + x9[i]*y9[i] + x10[i]*y8[i] + x11[i]*y7[i]
            + x12[i]*y6[i] + x13[i]*y5[i] + x14[i]*y4[i] + x15[i]*y3[i]
            + x16[i]*y2[i] + x17[i]*y1[i] + x18[i]*y0[i];
      *s19 += x0[i]*y19[i] + x1[i]*y18[i] + x2[i]*y17[i] + x3[i]*y16[i]
            + x4[i]*y15[i] + x5[i]*y14[i] + x6[i]*y13[i] + x7[i]*y12[i]
            + x8[i]*y11[i] + x9[i]*y10[i] + x10[i]*y9[i] + x11[i]*y8[i]
            + x12[i]*y7[i] + x13[i]*y6[i] + x14[i]*y5[i] + x15[i]*y4[i]
            + x16[i]*y3[i] + x17[i]*y2[i] + x18[i]*y1[i] + x19[i]*y0[i];
      *s20 += x0[i]*y20[i] + x1[i]*y19[i] + x2[i]*y18[i] + x3[i]*y17[i]
            + x4[i]*y16[i] + x5[i]*y15[i] + x6[i]*y14[i] + x7[i]*y13[i]
            + x8[i]*y12[i] + x9[i]*y11[i] + x10[i]*y10[i] + x11[i]*y9[i]
            + x12[i]*y8[i] + x13[i]*y7[i] + x14[i]*y6[i] + x15[i]*y5[i]
            + x16[i]*y4[i] + x17[i]*y3[i] + x18[i]*y2[i] + x19[i]*y1[i]
            + x20[i]*y0[i];
      *s21 += x0[i]*y21[i] + x1[i]*y20[i] + x2[i]*y19[i] + x3[i]*y18[i]
            + x4[i]*y17[i] + x5[i]*y16[i] + x6[i]*y15[i] + x7[i]*y14[i]
            + x8[i]*y13[i] + x9[i]*y12[i] + x10[i]*y11[i] + x11[i]*y10[i]
            + x12[i]*y9[i] + x13[i]*y8[i] + x14[i]*y7[i] + x15[i]*y6[i]
            + x16[i]*y5[i] + x17[i]*y4[i] + x18[i]*y3[i] + x19[i]*y2[i]
            + x20[i]*y1[i] + x21[i]*y0[i];
      *s22 += x0[i]*y22[i] + x1[i]*y21[i] + x2[i]*y20[i] + x3[i]*y19[i]
            + x4[i]*y18[i] + x5[i]*y17[i] + x6[i]*y16[i] + x7[i]*y15[i]
            + x8[i]*y14[i] + x9[i]*y13[i] + x10[i]*y12[i] + x11[i]*y11[i]
            + x12[i]*y10[i] + x13[i]*y9[i] + x14[i]*y8[i] + x15[i]*y7[i]
            + x16[i]*y6[i] + x17[i]*y5[i] + x18[i]*y4[i] + x19[i]*y3[i]
            + x20[i]*y2[i] + x21[i]*y1[i] + x22[i]*y0[i];
      *s23 += x0[i]*y23[i] + x1[i]*y22[i] + x2[i]*y21[i] + x3[i]*y20[i]
            + x4[i]*y19[i] + x5[i]*y18[i] + x6[i]*y17[i] + x7[i]*y16[i]
            + x8[i]*y15[i] + x9[i]*y14[i] + x10[i]*y13[i] + x11[i]*y12[i]
            + x12[i]*y11[i] + x13[i]*y10[i] + x14[i]*y9[i] + x15[i]*y8[i]
            + x16[i]*y7[i] + x17[i]*y6[i] + x18[i]*y5[i] + x19[i]*y4[i]
            + x20[i]*y3[i] + x21[i]*y2[i] + x22[i]*y1[i] + x23[i]*y0[i];
      *s24 += x0[i]*y24[i] + x1[i]*y23[i] + x2[i]*y22[i] + x3[i]*y21[i]
            + x4[i]*y20[i] + x5[i]*y19[i] + x6[i]*y18[i] + x7[i]*y17[i]
            + x8[i]*y16[i] + x9[i]*y15[i] + x10[i]*y14[i] + x11[i]*y13[i]
            + x12[i]*y12[i] + x13[i]*y11[i] + x14[i]*y10[i] + x15[i]*y9[i]
            + x16[i]*y8[i] + x17[i]*y7[i] + x18[i]*y6[i] + x19[i]*y5[i]
            + x20[i]*y4[i] + x21[i]*y3[i] + x22[i]*y2[i] + x23[i]*y1[i]
            + x24[i]*y0[i];
      *s25 += x0[i]*y25[i] + x1[i]*y24[i] + x2[i]*y23[i] + x3[i]*y22[i]
            + x4[i]*y21[i] + x5[i]*y20[i] + x6[i]*y19[i] + x7[i]*y18[i]
            + x8[i]*y17[i] + x9[i]*y16[i] + x10[i]*y15[i] + x11[i]*y14[i]
            + x12[i]*y13[i] + x13[i]*y12[i] + x14[i]*y11[i] + x15[i]*y10[i]
            + x16[i]*y9[i] + x17[i]*y8[i] + x18[i]*y7[i] + x19[i]*y6[i]
            + x20[i]*y5[i] + x21[i]*y4[i] + x22[i]*y3[i] + x23[i]*y2[i]
            + x24[i]*y1[i] + x25[i]*y0[i];
      *s26 += x0[i]*y26[i] + x1[i]*y25[i] + x2[i]*y24[i] + x3[i]*y23[i]
            + x4[i]*y22[i] + x5[i]*y21[i] + x6[i]*y20[i] + x7[i]*y19[i]
            + x8[i]*y18[i] + x9[i]*y17[i] + x10[i]*y16[i] + x11[i]*y15[i]
            + x12[i]*y14[i] + x13[i]*y13[i] + x14[i]*y12[i] + x15[i]*y11[i]
            + x16[i]*y10[i] + x17[i]*y9[i] + x18[i]*y8[i] + x19[i]*y7[i]
            + x20[i]*y6[i] + x21[i]*y5[i] + x22[i]*y4[i] + x23[i]*y3[i]
            + x24[i]*y2[i] + x25[i]*y1[i] + x26[i]*y0[i];
      *s27 += x0[i]*y27[i] + x1[i]*y26[i] + x2[i]*y25[i] + x3[i]*y24[i]
            + x4[i]*y23[i] + x5[i]*y22[i] + x6[i]*y21[i] + x7[i]*y20[i]
            + x8[i]*y19[i] + x9[i]*y18[i] + x10[i]*y17[i] + x11[i]*y16[i]
            + x12[i]*y15[i] + x13[i]*y14[i] + x14[i]*y13[i] + x15[i]*y12[i]
            + x16[i]*y11[i] + x17[i]*y10[i] + x18[i]*y9[i] + x19[i]*y8[i]
            + x20[i]*y7[i] + x21[i]*y6[i] + x22[i]*y5[i] + x23[i]*y4[i]
            + x24[i]*y3[i] + x25[i]*y2[i] + x26[i]*y1[i] + x27[i]*y0[i];
      *s28 += x0[i]*y28[i] + x1[i]*y27[i] + x2[i]*y26[i] + x3[i]*y25[i]
            + x4[i]*y24[i] + x5[i]*y23[i] + x6[i]*y22[i] + x7[i]*y21[i]
            + x8[i]*y20[i] + x9[i]*y19[i] + x10[i]*y18[i] + x11[i]*y17[i]
            + x12[i]*y16[i] + x13[i]*y15[i] + x14[i]*y14[i] + x15[i]*y13[i]
            + x16[i]*y12[i] + x17[i]*y11[i] + x18[i]*y10[i] + x19[i]*y9[i]
            + x20[i]*y8[i] + x21[i]*y7[i] + x22[i]*y6[i] + x23[i]*y5[i]
            + x24[i]*y4[i] + x25[i]*y3[i] + x26[i]*y2[i] + x27[i]*y1[i]
            + x28[i]*y0[i];
      *s29 += x0[i]*y29[i] + x1[i]*y28[i] + x2[i]*y27[i] + x3[i]*y26[i]
            + x4[i]*y25[i] + x5[i]*y24[i] + x6[i]*y23[i] + x7[i]*y22[i]
            + x8[i]*y21[i] + x9[i]*y20[i] + x10[i]*y19[i] + x11[i]*y18[i]
            + x12[i]*y17[i] + x13[i]*y16[i] + x14[i]*y15[i] + x15[i]*y14[i]
            + x16[i]*y13[i] + x17[i]*y12[i] + x18[i]*y11[i] + x19[i]*y10[i]
            + x20[i]*y9[i] + x21[i]*y8[i] + x22[i]*y7[i] + x23[i]*y6[i]
            + x24[i]*y5[i] + x25[i]*y4[i] + x26[i]*y3[i] + x27[i]*y2[i]
            + x28[i]*y1[i] + x29[i]*y0[i];
      *s30 += x0[i]*y30[i] + x1[i]*y29[i] + x2[i]*y28[i] + x3[i]*y27[i]
            + x4[i]*y26[i] + x5[i]*y25[i] + x6[i]*y24[i] + x7[i]*y23[i]
            + x8[i]*y22[i] + x9[i]*y21[i] + x10[i]*y20[i] + x11[i]*y19[i]
            + x12[i]*y18[i] + x13[i]*y17[i] + x14[i]*y16[i] + x15[i]*y15[i]
            + x16[i]*y14[i] + x17[i]*y13[i] + x18[i]*y12[i] + x19[i]*y11[i]
            + x20[i]*y10[i] + x21[i]*y9[i] + x22[i]*y8[i] + x23[i]*y7[i]
            + x24[i]*y6[i] + x25[i]*y5[i] + x26[i]*y4[i] + x27[i]*y3[i]
            + x28[i]*y2[i] + x29[i]*y1[i] + x30[i]*y0[i];
      *s31 += x0[i]*y31[i] + x1[i]*y30[i] + x2[i]*y29[i] + x3[i]*y28[i]
            + x4[i]*y27[i] + x5[i]*y26[i] + x6[i]*y25[i] + x7[i]*y24[i]
            + x8[i]*y23[i] + x9[i]*y22[i] + x10[i]*y21[i] + x11[i]*y20[i]
            + x12[i]*y19[i] + x13[i]*y18[i] + x14[i]*y17[i] + x15[i]*y16[i]
            + x16[i]*y15[i] + x17[i]*y14[i] + x18[i]*y13[i] + x19[i]*y12[i]
            + x20[i]*y11[i] + x21[i]*y10[i] + x22[i]*y9[i] + x23[i]*y8[i]
            + x24[i]*y7[i] + x25[i]*y6[i] + x26[i]*y5[i] + x27[i]*y4[i]
            + x28[i]*y3[i] + x29[i]*y2[i] + x30[i]*y1[i] + x31[i]*y0[i];
   }
}
