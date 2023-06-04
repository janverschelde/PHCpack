// The file dbl8_indexed_coefficients.cpp defines the functions specified in
// the file dbl8_indexed_coefficients.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random8_monomials.h"
#include "random8_vectors.h"
#include "random8_series.h"
#include "random8_polynomials.h"

using namespace std;

int dbl8_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihihi, double **cstlohihi,
   double **csthilohi, double **cstlolohi,
   double **csthihilo, double **cstlohilo,
   double **csthilolo, double **cstlololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo, int vrblvl )
{
   const int degp1 = deg+1;
   double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
   double rndhihilo,rndlohilo,rndhilolo,rndlololo;

   for(int i=0; i<dim; i++)
   {
      csthihihi[i] = new double[degp1];
      cstlohihi[i] = new double[degp1];
      csthilohi[i] = new double[degp1];
      cstlolohi[i] = new double[degp1];
      csthihilo[i] = new double[degp1];
      cstlohilo[i] = new double[degp1];
      csthilolo[i] = new double[degp1];
      cstlololo[i] = new double[degp1];

      random_dbl8_exponential
         (deg,&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
              &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo,
          csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
          csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;

         for(int j=0; j<=deg; j++)
            cout << csthihihi[i][j] << "  " << cstlohihi[i][j] << endl
                 << csthilohi[i][j] << "  " << cstlolohi[i][j] << endl
                 << csthihilo[i][j] << "  " << cstlohilo[i][j] << endl
                 << csthilolo[i][j] << "  " << cstlololo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double*[nbr[i]];
      cfflohihi[i] = new double*[nbr[i]];
      cffhilohi[i] = new double*[nbr[i]];
      cfflolohi[i] = new double*[nbr[i]];
      cffhihilo[i] = new double*[nbr[i]];
      cfflohilo[i] = new double*[nbr[i]];
      cffhilolo[i] = new double*[nbr[i]];
      cfflololo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffhihihi[i][j] = new double[degp1];
         cfflohihi[i][j] = new double[degp1];
         cffhilohi[i][j] = new double[degp1];
         cfflolohi[i][j] = new double[degp1];
         cffhihilo[i][j] = new double[degp1];
         cfflohilo[i][j] = new double[degp1];
         cffhilolo[i][j] = new double[degp1];
         cfflololo[i][j] = new double[degp1];

         random_dbl8_exponential
            (deg,&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                 &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo,
             cffhihihi[i][j],cfflohihi[i][j],cffhilohi[i][j],cfflolohi[i][j],
             cffhihilo[i][j],cfflohilo[i][j],cffhilolo[i][j],cfflololo[i][j]);
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<nbr[i]; j++)
         {
            cout << "-> coefficients of monomial " << j
                 << " of polynomial " << i << " :" << endl;

            for(int k=0; k<=deg; k++)
               cout << cffhihihi[i][j][k] << "  "
                    << cfflohihi[i][j][k] << endl
                    << cffhilohi[i][j][k] << "  "
                    << cfflolohi[i][j][k] << endl
                    << cffhihilo[i][j][k] << "  "
                    << cfflohilo[i][j][k] << endl
                    << cffhilolo[i][j][k] << "  "
                    << cfflololo[i][j][k] << endl;
         }
      }
   }
   return 0;
}

int cmplx8_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehihihi, double **cstrelohihi,
   double **cstrehilohi, double **cstrelolohi,
   double **cstrehihilo, double **cstrelohilo,
   double **cstrehilolo, double **cstrelololo,
   double **cstimhihihi, double **cstimlohihi,
   double **cstimhilohi, double **cstimlolohi,
   double **cstimhihilo, double **cstimlohilo,
   double **cstimhilolo, double **cstimlololo,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo, int vrblvl )
{
   const int degp1 = deg+1;
   double rndrehihihi,rndrelohihi,rndrehilohi,rndrelolohi;
   double rndrehihilo,rndrelohilo,rndrehilolo,rndrelololo;
   double rndimhihihi,rndimlohihi,rndimhilohi,rndimlolohi;
   double rndimhihilo,rndimlohilo,rndimhilolo,rndimlololo;

   for(int i=0; i<dim; i++)
   {
      cstrehihihi[i] = new double[degp1];
      cstrelohihi[i] = new double[degp1];
      cstrehilohi[i] = new double[degp1];
      cstrelolohi[i] = new double[degp1];
      cstrehihilo[i] = new double[degp1];
      cstrelohilo[i] = new double[degp1];
      cstrehilolo[i] = new double[degp1];
      cstrelololo[i] = new double[degp1];
      cstimhihihi[i] = new double[degp1];
      cstimlohihi[i] = new double[degp1];
      cstimhilohi[i] = new double[degp1];
      cstimlolohi[i] = new double[degp1];
      cstimhihilo[i] = new double[degp1];
      cstimlohilo[i] = new double[degp1];
      cstimhilolo[i] = new double[degp1];
      cstimlololo[i] = new double[degp1];
/*
      random_cmplx8_exponential
         (deg,&rndrehihihi,&rndrelohihi,&rndrehilohi,&rndrelolohi,
              &rndrehihilo,&rndrelohilo,&rndrehilolo,&rndrelololo,
              &rndimhihihi,&rndimlohihi,&rndimhilohi,&rndimlolohi,
              &rndimhihilo,&rndimlohilo,&rndimhilolo,&rndimlololo,
          cstrehihihi[i],cstrelohihi[i],cstrehilohi[i],cstrelolohi[i],
          cstrehihilo[i],cstrelohilo[i],cstrehilolo[i],cstrelololo[i],
          cstimhihihi[i],cstimlohihi[i],cstimhilohi[i],cstimlolohi[i],
          cstimhihilo[i],cstimlohilo[i],cstimhilolo[i],cstimlololo[i]);
 */
      for(int k=0; k<=deg; k++)
         random_octo_complex
            (&cstrehihihi[i][k],&cstrelohihi[i][k],
             &cstrehilohi[i][k],&cstrelolohi[i][k],
             &cstrehihilo[i][k],&cstrelohilo[i][k],
             &cstrehilolo[i][k],&cstrelololo[i][k],
             &cstimhihihi[i][k],&cstimlohihi[i][k],
             &cstimhilohi[i][k],&cstimlolohi[i][k],
             &cstimhihilo[i][k],&cstimlohilo[i][k],
             &cstimhilolo[i][k],&cstimlololo[i][k]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cstrehihihi[i][j] << "  " << cstrelohihi[i][j] << endl
                 << cstrehilohi[i][j] << "  " << cstrelolohi[i][j] << endl
                 << cstrehihilo[i][j] << "  " << cstrelohilo[i][j] << endl
                 << cstrehilolo[i][j] << "  " << cstrelololo[i][j] << endl
                 << cstimhihihi[i][j] << "  " << cstimlohihi[i][j] << endl
                 << cstimhilohi[i][j] << "  " << cstimlolohi[i][j] << endl
                 << cstimhihilo[i][j] << "  " << cstimlohilo[i][j] << endl
                 << cstimhilolo[i][j] << "  " << cstimlololo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffrehihihi[i] = new double*[nbr[i]];
      cffrelohihi[i] = new double*[nbr[i]];
      cffrehilohi[i] = new double*[nbr[i]];
      cffrelolohi[i] = new double*[nbr[i]];
      cffrehihilo[i] = new double*[nbr[i]];
      cffrelohilo[i] = new double*[nbr[i]];
      cffrehilolo[i] = new double*[nbr[i]];
      cffrelololo[i] = new double*[nbr[i]];
      cffimhihihi[i] = new double*[nbr[i]];
      cffimlohihi[i] = new double*[nbr[i]];
      cffimhilohi[i] = new double*[nbr[i]];
      cffimlolohi[i] = new double*[nbr[i]];
      cffimhihilo[i] = new double*[nbr[i]];
      cffimlohilo[i] = new double*[nbr[i]];
      cffimhilolo[i] = new double*[nbr[i]];
      cffimlololo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffrehihihi[i][j] = new double[degp1];
         cffrelohihi[i][j] = new double[degp1];
         cffrehilohi[i][j] = new double[degp1];
         cffrelolohi[i][j] = new double[degp1];
         cffrehihilo[i][j] = new double[degp1];
         cffrelohilo[i][j] = new double[degp1];
         cffrehilolo[i][j] = new double[degp1];
         cffrelololo[i][j] = new double[degp1];
         cffimhihihi[i][j] = new double[degp1];
         cffimlohihi[i][j] = new double[degp1];
         cffimhilohi[i][j] = new double[degp1];
         cffimlolohi[i][j] = new double[degp1];
         cffimhihilo[i][j] = new double[degp1];
         cffimlohilo[i][j] = new double[degp1];
         cffimhilolo[i][j] = new double[degp1];
         cffimlololo[i][j] = new double[degp1];
/*
         random_cmplx8_exponential
            (deg,&rndrehihihi,&rndrelohihi,&rndrehilohi,&rndrelolohi,
                 &rndrehihilo,&rndrelohilo,&rndrehilolo,&rndrelololo,
                 &rndimhihihi,&rndimlohihi,&rndimhilohi,&rndimlolohi,
                 &rndimhihilo,&rndimlohilo,&rndimhilolo,&rndimlololo,
             cffrehihihi[i][j],cffrelohihi[i][j],
             cffrehilohi[i][j],cffrelolohi[i][j],
             cffrehihilo[i][j],cffrelohilo[i][j],
             cffrehilolo[i][j],cffrelololo[i][j],
             cffimhihihi[i][j],cffimlohihi[i][j],
             cffimhilohi[i][j],cffimlolohi[i][j],
             cffimhihilo[i][j],cffimlohilo[i][j],
             cffimhilolo[i][j],cffimlololo[i][j]);
 */
         for(int k=0; k<=deg; k++)
            random_octo_complex
               (&cffrehihihi[i][j][k],&cffrelohihi[i][j][k],
                &cffrehilohi[i][j][k],&cffrelolohi[i][j][k],
                &cffrehihilo[i][j][k],&cffrelohilo[i][j][k],
                &cffrehilolo[i][j][k],&cffrelololo[i][j][k],
                &cffimhihihi[i][j][k],&cffimlohihi[i][j][k],
                &cffimhilohi[i][j][k],&cffimlolohi[i][j][k],
                &cffimhihilo[i][j][k],&cffimlohilo[i][j][k],
                &cffimhilolo[i][j][k],&cffimlololo[i][j][k]);
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<nbr[i]; j++)
         {
            cout << "-> coefficients of monomial " << j
                 << " of polynomial " << i << " :" << endl;

            for(int k=0; k<=deg; k++)
               cout << cffrehihihi[i][j][k] << "  "
                    << cffrelohihi[i][j][k] << endl
                    << cffrehilohi[i][j][k] << "  "
                    << cffrelolohi[i][j][k] << endl
                    << cffrehihilo[i][j][k] << "  "
                    << cffrelohilo[i][j][k] << endl
                    << cffrehilolo[i][j][k] << "  "
                    << cffrelololo[i][j][k] << endl
                    << cffimhihihi[i][j][k] << "  "
                    << cffimlohihi[i][j][k] << endl
                    << cffimhilohi[i][j][k] << "  "
                    << cffimlolohi[i][j][k] << endl
                    << cffimhihilo[i][j][k] << "  "
                    << cffimlohilo[i][j][k] << endl
                    << cffimhilolo[i][j][k] << "  "
                    << cffimlololo[i][j][k] << endl;
         }
      }
   }
   return 0;
}
