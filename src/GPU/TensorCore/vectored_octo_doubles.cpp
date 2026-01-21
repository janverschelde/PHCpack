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

void quarter_od_matrix
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Ahihihi0, double **Ahihihi1, double **Ahihihi2, double **Ahihihi3,
   double **Alohihi0, double **Alohihi1, double **Alohihi2, double **Alohihi3,
   double **Ahilohi0, double **Ahilohi1, double **Ahilohi2, double **Ahilohi3,
   double **Alolohi0, double **Alolohi1, double **Alolohi2, double **Alolohi3,
   double **Ahihilo0, double **Ahihilo1, double **Ahihilo2, double **Ahihilo3,
   double **Alohilo0, double **Alohilo1, double **Alohilo2, double **Alohilo3,
   double **Ahilolo0, double **Ahilolo1, double **Ahilolo2, double **Ahilolo3,
   double **Alololo0, double **Alololo1, double **Alololo2, double **Alololo3 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         quarter_octo_double
            (Ahihihi[i][j], Alohihi[i][j], Ahilohi[i][j], Alolohi[i][j],
             Ahihilo[i][j], Alohilo[i][j], Ahilolo[i][j], Alololo[i][j],
             &Ahihihi0[i][j], &Ahihihi1[i][j],
             &Ahihihi2[i][j], &Ahihihi3[i][j],
             &Alohihi0[i][j], &Alohihi1[i][j],
             &Alohihi2[i][j], &Alohihi3[i][j],
             &Ahilohi0[i][j], &Ahilohi1[i][j],
             &Ahilohi2[i][j], &Ahilohi3[i][j],
             &Alolohi0[i][j], &Alolohi1[i][j],
             &Alolohi2[i][j], &Alolohi3[i][j],
             &Ahihilo0[i][j], &Ahihilo1[i][j],
             &Ahihilo2[i][j], &Ahihilo3[i][j],
             &Alohilo0[i][j], &Alohilo1[i][j],
             &Alohilo2[i][j], &Alohilo3[i][j],
             &Ahilolo0[i][j], &Ahilolo1[i][j],
             &Ahilolo2[i][j], &Ahilolo3[i][j],
             &Alololo0[i][j], &Alololo1[i][j],
             &Alololo2[i][j], &Alololo3[i][j]);
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

void to_octo_double_matrix
 ( int nrows, int ncols,
   double **Ahihihi0, double **Ahihihi1, double **Ahihihi2, double **Ahihihi3,
   double **Alohihi0, double **Alohihi1, double **Alohihi2, double **Alohihi3,
   double **Ahilohi0, double **Ahilohi1, double **Ahilohi2, double **Ahilohi3,
   double **Alolohi0, double **Alolohi1, double **Alolohi2, double **Alolohi3,
   double **Ahihilo0, double **Ahihilo1, double **Ahihilo2, double **Ahihilo3,
   double **Alohilo0, double **Alohilo1, double **Alohilo2, double **Alohilo3,
   double **Ahilolo0, double **Ahilolo1, double **Ahilolo2, double **Ahilolo3,
   double **Alololo0, double **Alololo1, double **Alololo2, double **Alololo3,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         to_octo_double
            (Ahihihi0[i][j], Ahihihi1[i][j], Ahihihi2[i][j], Ahihihi3[i][j],
             Alohihi0[i][j], Alohihi1[i][j], Alohihi2[i][j], Alohihi3[i][j],
             Ahilohi0[i][j], Ahilohi1[i][j], Ahilohi2[i][j], Ahilohi3[i][j],
             Alolohi0[i][j], Alolohi1[i][j], Alolohi2[i][j], Alolohi3[i][j],
             Ahihilo0[i][j], Ahihilo1[i][j], Ahihilo2[i][j], Ahihilo3[i][j],
             Alohilo0[i][j], Alohilo1[i][j], Alohilo2[i][j], Alohilo3[i][j],
             Ahilolo0[i][j], Ahilolo1[i][j], Ahilolo2[i][j], Ahilolo3[i][j],
             Alololo0[i][j], Alololo1[i][j], Alololo2[i][j], Alololo3[i][j],
             &Ahihihi[i][j], &Alohihi[i][j], &Ahilohi[i][j], &Alolohi[i][j],
             &Ahihilo[i][j], &Alohilo[i][j], &Ahilolo[i][j], &Alololo[i][j]);
      }
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

void transpose_od_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **A16, double **A17, double **A18, double **A19,
   double **A20, double **A21, double **A22, double **A23,
   double **A24, double **A25, double **A26, double **A27,
   double **A28, double **A29, double **A30, double **A31,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7,
   double **T8, double **T9, double **T10, double **T11,
   double **T12, double **T13, double **T14, double **T15,
   double **T16, double **T17, double **T18, double **T19,
   double **T20, double **T21, double **T22, double **T23,
   double **T24, double **T25, double **T26, double **T27,
   double **T28, double **T29, double **T30, double **T31 )
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
         T16[j][i] = A16[i][j]; T17[j][i] = A17[i][j];
         T18[j][i] = A18[i][j]; T19[j][i] = A19[i][j];
         T20[j][i] = A20[i][j]; T21[j][i] = A21[i][j];
         T22[j][i] = A22[i][j]; T23[j][i] = A23[i][j];
         T24[j][i] = A24[i][j]; T25[j][i] = A25[i][j];
         T26[j][i] = A26[i][j]; T27[j][i] = A27[i][j];
         T28[j][i] = A28[i][j]; T29[j][i] = A29[i][j];
         T30[j][i] = A30[i][j]; T31[j][i] = A31[i][j];
      }
}

void vectored_od_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **A16, double **A17, double **A18, double **A19,
   double **A20, double **A21, double **A22, double **A23,
   double **A24, double **A25, double **A26, double **A27,
   double **A28, double **A29, double **A30, double **A31,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **B8, double **B9, double **B10, double **B11,
   double **B12, double **B13, double **B14, double **B15,
   double **B16, double **B17, double **B18, double **B19,
   double **B20, double **B21, double **B22, double **B23,
   double **B24, double **B25, double **B26, double **B27,
   double **B28, double **B29, double **B30, double **B31,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7,
   double **C8, double **C9, double **C10, double **C11,
   double **C12, double **C13, double **C14, double **C15,
   double **C16, double **C17, double **C18, double **C19,
   double **C20, double **C21, double **C22, double **C23,
   double **C24, double **C25, double **C26, double **C27,
   double **C28, double **C29, double **C30, double **C31 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         vectored_od_product
            (dim, A0[i], A1[i], A2[i], A3[i], A4[i], A5[i], A6[i], A7[i],
             A8[i], A9[i], A10[i], A11[i], A12[i], A13[i], A14[i], A15[i],
             A16[i], A17[i], A18[i], A19[i], A20[i], A21[i], A22[i], A23[i],
             A24[i], A25[i], A26[i], A27[i], A28[i], A29[i], A30[i], A31[i],
             B0[j], B1[j], B2[j], B3[j], B4[j], B5[j], B6[j], B7[j],
             B8[j], B9[j], B10[j], B11[j], B12[j], B13[j], B14[j], B15[j],
             B16[j], B17[j], B18[j], B19[j], B20[j], B21[j], B22[j], B23[j],
             B24[j], B25[j], B26[j], B27[j], B28[j], B29[j], B30[j], B31[j],
             &C0[i][j], &C1[i][j], &C2[i][j], &C3[i][j],
             &C4[i][j], &C5[i][j], &C6[i][j], &C7[i][j],
             &C8[i][j], &C9[i][j], &C10[i][j], &C11[i][j],
             &C12[i][j], &C13[i][j], &C14[i][j], &C15[i][j],
             &C16[i][j], &C17[i][j], &C18[i][j], &C19[i][j],
             &C20[i][j], &C21[i][j], &C22[i][j], &C23[i][j],
             &C24[i][j], &C25[i][j], &C26[i][j], &C27[i][j],
             &C28[i][j], &C29[i][j], &C30[i][j], &C31[i][j]);
      }
}

void od_convolute_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **A16, double **A17, double **A18, double **A19,
   double **A20, double **A21, double **A22, double **A23,
   double **A24, double **A25, double **A26, double **A27,
   double **A28, double **A29, double **A30, double **A31, double **cA )
{
   for(int i=0; i<32*nrows; i++)
      for(int j=0; j<32*ncols; j++) cA[i][j] = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         for(int k=0; k<32; k++) cA[32*i+k][32*j+k] = A0[i][j];
         for(int k=0; k<31; k++) cA[32*i+k+1][32*j+k] = A1[i][j];
         for(int k=0; k<30; k++) cA[32*i+k+2][32*j+k] = A2[i][j];
         for(int k=0; k<29; k++) cA[32*i+k+3][32*j+k] = A3[i][j];
         for(int k=0; k<28; k++) cA[32*i+k+4][32*j+k] = A4[i][j];
         for(int k=0; k<27; k++) cA[32*i+k+5][32*j+k] = A5[i][j];
         for(int k=0; k<26; k++) cA[32*i+k+6][32*j+k] = A6[i][j];
         for(int k=0; k<25; k++) cA[32*i+k+7][32*j+k] = A7[i][j];
         for(int k=0; k<24; k++) cA[32*i+k+8][32*j+k] = A8[i][j];
         for(int k=0; k<23; k++) cA[32*i+k+9][32*j+k] = A9[i][j];
         for(int k=0; k<22; k++) cA[32*i+k+10][32*j+k] = A10[i][j];
         for(int k=0; k<21; k++) cA[32*i+k+11][32*j+k] = A11[i][j];
         for(int k=0; k<20; k++) cA[32*i+k+12][32*j+k] = A12[i][j];
         for(int k=0; k<19; k++) cA[32*i+k+13][32*j+k] = A13[i][j];
         for(int k=0; k<18; k++) cA[32*i+k+14][32*j+k] = A14[i][j];
         for(int k=0; k<17; k++) cA[32*i+k+15][32*j+k] = A15[i][j];
         for(int k=0; k<16; k++) cA[32*i+k+16][32*j+k] = A16[i][j];
         for(int k=0; k<15; k++) cA[32*i+k+17][32*j+k] = A17[i][j];
         for(int k=0; k<14; k++) cA[32*i+k+18][32*j+k] = A18[i][j];
         for(int k=0; k<13; k++) cA[32*i+k+19][32*j+k] = A19[i][j];
         for(int k=0; k<12; k++) cA[32*i+k+20][32*j+k] = A20[i][j];
         for(int k=0; k<11; k++) cA[32*i+k+21][32*j+k] = A21[i][j];
         for(int k=0; k<10; k++) cA[32*i+k+22][32*j+k] = A22[i][j];
         for(int k=0; k<9; k++) cA[32*i+k+23][32*j+k] = A23[i][j];
         for(int k=0; k<8; k++) cA[32*i+k+24][32*j+k] = A24[i][j];
         for(int k=0; k<7; k++) cA[32*i+k+25][32*j+k] = A25[i][j];
         for(int k=0; k<6; k++) cA[32*i+k+26][32*j+k] = A26[i][j];
         for(int k=0; k<5; k++) cA[32*i+k+27][32*j+k] = A27[i][j];
         for(int k=0; k<4; k++) cA[32*i+k+28][32*j+k] = A28[i][j];
         for(int k=0; k<3; k++) cA[32*i+k+29][32*j+k] = A29[i][j];
         for(int k=0; k<2; k++) cA[32*i+k+30][32*j+k] = A30[i][j];
         cA[32*i+31][32*j] = A31[i][j];
      }
}

void od_stack_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **A16, double **A17, double **A18, double **A19,
   double **A20, double **A21, double **A22, double **A23,
   double **A24, double **A25, double **A26, double **A27,
   double **A28, double **A29, double **A30, double **A31, double **sA )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         sA[32*i][j] = A0[i][j];
         sA[32*i+1][j] = A1[i][j];
         sA[32*i+2][j] = A2[i][j];
         sA[32*i+3][j] = A3[i][j];
         sA[32*i+4][j] = A4[i][j];
         sA[32*i+5][j] = A5[i][j];
         sA[32*i+6][j] = A6[i][j];
         sA[32*i+7][j] = A7[i][j];
         sA[32*i+8][j] = A8[i][j];
         sA[32*i+9][j] = A9[i][j];
         sA[32*i+10][j] = A10[i][j];
         sA[32*i+11][j] = A11[i][j];
         sA[32*i+12][j] = A12[i][j];
         sA[32*i+13][j] = A13[i][j];
         sA[32*i+14][j] = A14[i][j];
         sA[32*i+15][j] = A15[i][j];
         sA[32*i+16][j] = A16[i][j];
         sA[32*i+17][j] = A17[i][j];
         sA[32*i+18][j] = A18[i][j];
         sA[32*i+19][j] = A19[i][j];
         sA[32*i+20][j] = A20[i][j];
         sA[32*i+21][j] = A21[i][j];
         sA[32*i+22][j] = A22[i][j];
         sA[32*i+23][j] = A23[i][j];
         sA[32*i+24][j] = A24[i][j];
         sA[32*i+25][j] = A25[i][j];
         sA[32*i+26][j] = A26[i][j];
         sA[32*i+27][j] = A27[i][j];
         sA[32*i+28][j] = A28[i][j];
         sA[32*i+29][j] = A29[i][j];
         sA[32*i+30][j] = A30[i][j];
         sA[32*i+31][j] = A31[i][j];
      }
}

void extract_od_quarters
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7,
   double **D8, double **D9, double **D10, double **D11,
   double **D12, double **D13, double **D14, double **D15,
   double **D16, double **D17, double **D18, double **D19,
   double **D20, double **D21, double **D22, double **D23,
   double **D24, double **D25, double **D26, double **D27,
   double **D28, double **D29, double **D30, double **D31 )
{
   for(int i=0; i<32*nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         int k = i % 32;
         int row = i/32;

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
         if(k == 16) D16[row][j] = qC[i][j];
         if(k == 17) D17[row][j] = qC[i][j];
         if(k == 18) D18[row][j] = qC[i][j];
         if(k == 19) D19[row][j] = qC[i][j];
         if(k == 20) D20[row][j] = qC[i][j];
         if(k == 21) D21[row][j] = qC[i][j];
         if(k == 22) D22[row][j] = qC[i][j];
         if(k == 23) D23[row][j] = qC[i][j];
         if(k == 24) D24[row][j] = qC[i][j];
         if(k == 25) D25[row][j] = qC[i][j];
         if(k == 26) D26[row][j] = qC[i][j];
         if(k == 27) D27[row][j] = qC[i][j];
         if(k == 28) D28[row][j] = qC[i][j];
         if(k == 29) D29[row][j] = qC[i][j];
         if(k == 30) D30[row][j] = qC[i][j];
         if(k == 31) D31[row][j] = qC[i][j];
      }
}
