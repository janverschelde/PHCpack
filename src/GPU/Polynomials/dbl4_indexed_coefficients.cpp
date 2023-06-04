// The file dbl4_indexed_coefficients.cpp defines the functions specified in
// the file dbl4_indexed_coefficients.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random4_monomials.h"
#include "random4_series.h"
#include "random4_vectors.h"
#include "random4_polynomials.h"

using namespace std;

int dbl4_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihi, double **cstlohi, double **csthilo, double **cstlolo,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   int vrblvl )
{
   const int degp1 = deg+1;
   double rndhihi,rndlohi,rndhilo,rndlolo;

   for(int i=0; i<dim; i++)
   {
      csthihi[i] = new double[deg+1];
      cstlohi[i] = new double[deg+1];
      csthilo[i] = new double[deg+1];
      cstlolo[i] = new double[deg+1];

      random_dbl4_exponential
         (deg,&rndhihi,&rndlohi,&rndhilo,&rndlolo,
          csthihi[i],cstlohi[i],csthilo[i],cstlolo[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;

         for(int j=0; j<=deg; j++)
            cout << csthihi[i][j] << "  " << cstlohi[i][j] << endl
                 << csthilo[i][j] << "  " << cstlolo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffhihi[i] = new double*[nbr[i]];
      cfflohi[i] = new double*[nbr[i]];
      cffhilo[i] = new double*[nbr[i]];
      cfflolo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffhihi[i][j] = new double[deg+1];
         cfflohi[i][j] = new double[deg+1];
         cffhilo[i][j] = new double[deg+1];
         cfflolo[i][j] = new double[deg+1];

         random_dbl4_exponential
            (deg,&rndhihi,&rndlohi,&rndhilo,&rndlolo,
             cffhihi[i][j],cfflohi[i][j],cffhilo[i][j],cfflolo[i][j]);
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
               cout << cffhihi[i][j][k] << "  " << cfflohi[i][j][k] << endl
                    << cffhilo[i][j][k] << "  " << cfflolo[i][j][k] << endl;
         }
      }
   }
   return 0;
}

int cmplx4_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehihi, double **cstrelohi,
   double **cstrehilo, double **cstrelolo,
   double **cstimhihi, double **cstimlohi,
   double **cstimhilo, double **cstimlolo,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo, int vrblvl )
{
   const int degp1 = deg+1;
   double rndrehihi,rndrelohi,rndrehilo,rndrelolo;
   double rndimhihi,rndimlohi,rndimhilo,rndimlolo;

   for(int i=0; i<dim; i++)
   {
      cstrehihi[i] = new double[degp1];
      cstrelohi[i] = new double[degp1];
      cstrehilo[i] = new double[degp1];
      cstrelolo[i] = new double[degp1];
      cstimhihi[i] = new double[degp1];
      cstimlohi[i] = new double[degp1];
      cstimhilo[i] = new double[degp1];
      cstimlolo[i] = new double[degp1];
/*
      random_cmplx4_exponential
         (deg,&rndrehihi,&rndrelohi,&rndrehilo,&rndrelolo,
              &rndimhihi,&rndimlohi,&rndimhilo,&rndimlolo,
          cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
          cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i]);
 */
      for(int k=0; k<=deg; k++)
         random_quad_double_complex
            (&cstrehihi[i][k],&cstrelohi[i][k],
             &cstrehilo[i][k],&cstrelolo[i][k],
             &cstimhihi[i][k],&cstimlohi[i][k],
             &cstimhilo[i][k],&cstimlolo[i][k]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;

         for(int j=0; j<=deg; j++)
            cout << cstrehihi[i][j] << "  " << cstrelohi[i][j] << endl
                 << cstrehilo[i][j] << "  " << cstrelolo[i][j] << endl
                 << cstimhihi[i][j] << "  " << cstimlohi[i][j] << endl
                 << cstimhilo[i][j] << "  " << cstimlolo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffrehihi[i] = new double*[nbr[i]];
      cffrelohi[i] = new double*[nbr[i]];
      cffrehilo[i] = new double*[nbr[i]];
      cffrelolo[i] = new double*[nbr[i]];
      cffimhihi[i] = new double*[nbr[i]];
      cffimlohi[i] = new double*[nbr[i]];
      cffimhilo[i] = new double*[nbr[i]];
      cffimlolo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffrehihi[i][j] = new double[degp1];
         cffrelohi[i][j] = new double[degp1];
         cffrehilo[i][j] = new double[degp1];
         cffrelolo[i][j] = new double[degp1];
         cffimhihi[i][j] = new double[degp1];
         cffimlohi[i][j] = new double[degp1];
         cffimhilo[i][j] = new double[degp1];
         cffimlolo[i][j] = new double[degp1];
/*
         random_cmplx4_exponential
            (deg,&rndrehihi,&rndrelohi,&rndrehilo,&rndrelolo,
                 &rndimhihi,&rndimlohi,&rndimhilo,&rndimlolo,
             cffrehihi[i][j],cffrelohi[i][j],cffrehilo[i][j],cffrelolo[i][j],
             cffimhihi[i][j],cffimlohi[i][j],cffimhilo[i][j],cffimlolo[i][j]);
 */
         for(int k=0; k<=deg; k++)
            random_quad_double_complex
               (&cffrehihi[i][j][k],&cffrelohi[i][j][k],
                &cffrehilo[i][j][k],&cffrelolo[i][j][k],
                &cffimhihi[i][j][k],&cffimlohi[i][j][k],
                &cffimhilo[i][j][k],&cffimlolo[i][j][k]);
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
               cout << cffrehihi[i][j][k] << "  "
                    << cffrelohi[i][j][k] << endl
                    << cffrehilo[i][j][k] << "  "
                    << cffrelolo[i][j][k] << endl
                    << cffimhihi[i][j][k] << "  "
                    << cffimlohi[i][j][k] << endl
                    << cffimhilo[i][j][k] << "  "
                    << cffimlolo[i][j][k] << endl;
         }
      }
   }
   return 0;
}
