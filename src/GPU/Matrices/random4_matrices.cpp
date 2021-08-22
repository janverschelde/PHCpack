// The file random4_matrices.cpp defines the functions specified in
// the file random4_matrices.h.

#include "random4_vectors.h"
#include "random4_matrices.h"
#include "quad_double_functions.h"

void random_dbl4_upper_matrix
 ( int rows, int cols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Ahihi[i][j] = 0.0; Alohi[i][j] = 0.0;
         Ahilo[i][j] = 0.0; Alolo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
         random_quad_double
            (&Ahihi[i][j],&Alohi[i][j],&Ahilo[i][j],&Alolo[i][j]);
   }
}

void random_cmplx4_upper_matrix
 ( int rows, int cols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo )
{
   double rnd_hihi,rnd_lohi,rnd_hilo,rnd_lolo;
   double cosrnd_hihi,cosrnd_lohi,cosrnd_hilo,cosrnd_lolo;
   double sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Arehihi[i][j] = 0.0; Arelohi[i][j] = 0.0;
         Arehilo[i][j] = 0.0; Arelolo[i][j] = 0.0;
         Aimhihi[i][j] = 0.0; Aimlohi[i][j] = 0.0;
         Aimhilo[i][j] = 0.0; Aimlolo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
      {
         random_quad_double(&rnd_hihi,&rnd_lohi,&rnd_hilo,&rnd_lolo);
         sinrnd_hihi = rnd_hihi; sinrnd_lohi = rnd_lohi;
         sinrnd_hilo = rnd_hilo; sinrnd_lolo = rnd_lolo;

         double y_hihi,y_lohi,y_hilo,y_lolo;  // work around to compute cos

         qdf_sqr(sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo,
                     &y_hihi,    &y_lohi,    &y_hilo,    &y_lolo);
         qdf_minus(&y_hihi,&y_lohi,&y_hilo,&y_lolo);       // y = -sin^2
         qdf_inc_d(&y_hihi,&y_lohi,&y_hilo,&y_lolo,1.0);   // y = 1 - sin^2
                                                      // cos is sqrt(1-sin^2)
         qdf_sqrt(y_hihi,y_lohi,y_hilo,y_lolo,
                  &cosrnd_hihi,&cosrnd_lohi,&cosrnd_hilo,&cosrnd_lolo);

         Arehihi[i][j] = cosrnd_hihi; Arelohi[i][j] = cosrnd_lohi;
         Arehilo[i][j] = cosrnd_hilo; Arelolo[i][j] = cosrnd_lolo;
         Aimhihi[i][j] = sinrnd_hihi; Aimlohi[i][j] = sinrnd_lohi;
         Aimhilo[i][j] = sinrnd_hilo; Aimlolo[i][j] = sinrnd_lolo;
      }
   }
}

void random_dbl4_matrix
 ( int rows, int cols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
         random_quad_double
            (&Ahihi[i][j],&Alohi[i][j],&Ahilo[i][j],&Alolo[i][j]);
   }
}

void random_cmplx4_matrix
 ( int rows, int cols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo )
{
   double rnd_hihi,rnd_lohi,rnd_hilo,rnd_lolo;
   double cosrnd_hihi,cosrnd_lohi,cosrnd_hilo,cosrnd_lolo;
   double sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
      {
         random_quad_double(&rnd_hihi,&rnd_lohi,&rnd_hilo,&rnd_lolo);
         sinrnd_hihi = rnd_hihi; sinrnd_lohi = rnd_lohi;
         sinrnd_hilo = rnd_hilo; sinrnd_lolo = rnd_lolo;

         double y_hihi,y_lohi,y_hilo,y_lolo;  // work around to compute cos

         qdf_sqr(sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo,
                     &y_hihi,    &y_lohi,    &y_hilo,    &y_lolo);
         qdf_minus(&y_hihi,&y_lohi,&y_hilo,&y_lolo);       // y = -sin^2
         qdf_inc_d(&y_hihi,&y_lohi,&y_hilo,&y_lolo,1.0);   // y = 1 - sin^2
                                                      // cos is sqrt(1-sin^2)
         qdf_sqrt(y_hihi,y_lohi,y_hilo,y_lolo,
                  &cosrnd_hihi,&cosrnd_lohi,&cosrnd_hilo,&cosrnd_lolo);

         Arehihi[i][j] = cosrnd_hihi; Arelohi[i][j] = cosrnd_lohi;
         Arehilo[i][j] = cosrnd_hilo; Arelolo[i][j] = cosrnd_lolo;
         Aimhihi[i][j] = sinrnd_hihi; Aimlohi[i][j] = sinrnd_lohi;
         Aimhilo[i][j] = sinrnd_hilo; Aimlolo[i][j] = sinrnd_lolo;
      }
   }
}
