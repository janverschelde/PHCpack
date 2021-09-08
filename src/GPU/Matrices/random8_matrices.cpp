// The file random8_matrices.cpp defines the functions specified in
// the file random8_matrices.h.

#include "random8_vectors.h"
#include "random8_matrices.h"
#include "octo_double_functions.h"

void random_dbl8_upper_matrix
 ( int rows, int cols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Ahihihi[i][j] = 0.0; Alohihi[i][j] = 0.0;
         Ahilohi[i][j] = 0.0; Alolohi[i][j] = 0.0;
         Ahihilo[i][j] = 0.0; Alohilo[i][j] = 0.0;
         Ahilolo[i][j] = 0.0; Alololo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
         random_octo_double
            (&Ahihihi[i][j],&Alohihi[i][j],&Ahilohi[i][j],&Alolohi[i][j],
             &Ahihilo[i][j],&Alohilo[i][j],&Ahilolo[i][j],&Alololo[i][j]);
   }
}

void random_cmplx8_upper_matrix
 ( int rows, int cols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo )
{
   double rnd_hihihi,rnd_lohihi,rnd_hilohi,rnd_lolohi;
   double rnd_hihilo,rnd_lohilo,rnd_hilolo,rnd_lololo;
   double cosrnd_hihihi,cosrnd_lohihi,cosrnd_hilohi,cosrnd_lolohi;
   double cosrnd_hihilo,cosrnd_lohilo,cosrnd_hilolo,cosrnd_lololo;
   double sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi;
   double sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Arehihihi[i][j] = 0.0; Arelohihi[i][j] = 0.0;
         Arehilohi[i][j] = 0.0; Arelolohi[i][j] = 0.0;
         Arehihilo[i][j] = 0.0; Arelohilo[i][j] = 0.0;
         Arehilolo[i][j] = 0.0; Arelololo[i][j] = 0.0;
         Aimhihihi[i][j] = 0.0; Aimlohihi[i][j] = 0.0;
         Aimhilohi[i][j] = 0.0; Aimlolohi[i][j] = 0.0;
         Aimhihilo[i][j] = 0.0; Aimlohilo[i][j] = 0.0;
         Aimhilolo[i][j] = 0.0; Aimlololo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
      {
         random_octo_double
            (&rnd_hihihi,&rnd_lohihi,&rnd_hilohi,&rnd_lolohi,
             &rnd_hihilo,&rnd_lohilo,&rnd_hilolo,&rnd_lololo);
         sinrnd_hihihi = rnd_hihihi; sinrnd_lohihi = rnd_lohihi;
         sinrnd_hilohi = rnd_hilohi; sinrnd_lolohi = rnd_lolohi;
         sinrnd_hihilo = rnd_hihilo; sinrnd_lohilo = rnd_lohilo;
         sinrnd_hilolo = rnd_hilolo; sinrnd_lololo = rnd_lololo;

                                               // work around to compute cos
         double y_hihihi,y_lohihi,y_hilohi,y_lolohi; 
         double y_hihilo,y_lohilo,y_hilolo,y_lololo;

         odf_sqr(sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi,
                 sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo,
                     &y_hihihi,    &y_lohihi,    &y_hilohi,    &y_lolohi,
                     &y_hihilo,    &y_lohilo,    &y_hilolo,    &y_lololo);
         odf_minus(&y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
                   &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo);  // y = -sin^2
                                                           // y = 1 - sin^2
         odf_inc_d(&y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
                   &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo,1.0);
                                                      // cos is sqrt(1-sin^2)
         odf_sqrt(y_hihihi,y_lohihi,y_hilohi,y_lolohi,
                  y_hihilo,y_lohilo,y_hilolo,y_lololo,
                  &cosrnd_hihihi,&cosrnd_lohihi,&cosrnd_hilohi,&cosrnd_lolohi,
                  &cosrnd_hihilo,&cosrnd_lohilo,&cosrnd_hilolo,&cosrnd_lololo);

         Arehihihi[i][j] = cosrnd_hihihi; Arelohihi[i][j] = cosrnd_lohihi;
         Arehilohi[i][j] = cosrnd_hilohi; Arelolohi[i][j] = cosrnd_lolohi;
         Arehihilo[i][j] = cosrnd_hihilo; Arelohilo[i][j] = cosrnd_lohilo;
         Arehilolo[i][j] = cosrnd_hilolo; Arelololo[i][j] = cosrnd_lololo;
         Aimhihihi[i][j] = sinrnd_hihihi; Aimlohihi[i][j] = sinrnd_lohihi;
         Aimhilohi[i][j] = sinrnd_hilohi; Aimlolohi[i][j] = sinrnd_lolohi;
         Aimhihilo[i][j] = sinrnd_hihilo; Aimlohilo[i][j] = sinrnd_lohilo;
         Aimhilolo[i][j] = sinrnd_hilolo; Aimlololo[i][j] = sinrnd_lololo;
      }
   }
}

void random_dbl8_matrix
 ( int rows, int cols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
         random_octo_double
            (&Ahihihi[i][j],&Alohihi[i][j],&Ahilohi[i][j],&Alolohi[i][j],
             &Ahihilo[i][j],&Alohilo[i][j],&Ahilolo[i][j],&Alololo[i][j]);
   }
}

void random_cmplx8_matrix
 ( int rows, int cols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo )
{
   double rnd_hihihi,rnd_lohihi,rnd_hilohi,rnd_lolohi;
   double rnd_hihilo,rnd_lohilo,rnd_hilolo,rnd_lololo;
   double cosrnd_hihihi,cosrnd_lohihi,cosrnd_hilohi,cosrnd_lolohi;
   double cosrnd_hihilo,cosrnd_lohilo,cosrnd_hilolo,cosrnd_lololo;
   double sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi;
   double sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
      {
         random_octo_double
            (&rnd_hihihi,&rnd_lohihi,&rnd_hilohi,&rnd_lolohi,
             &rnd_hihilo,&rnd_lohilo,&rnd_hilolo,&rnd_lololo);
         sinrnd_hihihi = rnd_hihihi; sinrnd_lohihi = rnd_lohihi;
         sinrnd_hilohi = rnd_hilohi; sinrnd_lolohi = rnd_lolohi;
         sinrnd_hihilo = rnd_hihilo; sinrnd_lohilo = rnd_lohilo;
         sinrnd_hilolo = rnd_hilolo; sinrnd_lololo = rnd_lololo;

         // work around to compute cos
         double y_hihihi,y_lohihi,y_hilohi,y_lolohi; 
         double y_hihilo,y_lohilo,y_hilolo,y_lololo;

         odf_sqr(sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi,
                 sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo,
                     &y_hihihi,    &y_lohihi,    &y_hilohi,    &y_lolohi,
                     &y_hihilo,    &y_lohilo,    &y_hilolo,    &y_lololo);
         odf_minus(&y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
                   &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo);  // y = -sin^2
                                                           // y = 1 - sin^2
         odf_inc_d(&y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
                   &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo,1.0); 
                                                      // cos is sqrt(1-sin^2)
         odf_sqrt(y_hihihi,y_lohihi,y_hilohi,y_lolohi,
                  y_hihilo,y_lohilo,y_hilolo,y_lololo,
                  &cosrnd_hihihi,&cosrnd_lohihi,&cosrnd_hilohi,&cosrnd_lolohi,
                  &cosrnd_hihilo,&cosrnd_lohilo,&cosrnd_hilolo,&cosrnd_lololo);

         Arehihihi[i][j] = cosrnd_hihihi; Arelohihi[i][j] = cosrnd_lohihi;
         Arehilohi[i][j] = cosrnd_hilohi; Arelolohi[i][j] = cosrnd_lolohi;
         Arehihilo[i][j] = cosrnd_hihilo; Arelohilo[i][j] = cosrnd_lohilo;
         Arehilolo[i][j] = cosrnd_hilolo; Arelololo[i][j] = cosrnd_lololo;
         Aimhihihi[i][j] = sinrnd_hihihi; Aimlohihi[i][j] = sinrnd_lohihi;
         Aimhilohi[i][j] = sinrnd_hilohi; Aimlolohi[i][j] = sinrnd_lolohi;
         Aimhihilo[i][j] = sinrnd_hihilo; Aimlohilo[i][j] = sinrnd_lohilo;
         Aimhilolo[i][j] = sinrnd_hilolo; Aimlololo[i][j] = sinrnd_lololo;
      }
   }
}
