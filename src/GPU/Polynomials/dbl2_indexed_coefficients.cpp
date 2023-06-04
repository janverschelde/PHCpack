// The file dbl2_indexed_coefficients.cpp defines the functions specified in
// the file dbl2_indexed_coefficients.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random2_vectors.h"
#include "random2_monomials.h"
#include "random2_polynomials.h"
#include "random2_series.h"

using namespace std;

int dbl2_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthi, double **cstlo, double ***cffhi, double ***cfflo,
   int vrblvl )
{
   const int degp1 = deg+1;
   double rndhi,rndlo;

   for(int i=0; i<dim; i++)
   {
      csthi[i] = new double[deg+1];
      cstlo[i] = new double[deg+1];

      random_dbl2_exponential(deg,&rndhi,&rndlo,csthi[i],cstlo[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << csthi[i][j] << "  " << cstlo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffhi[i] = new double*[nbr[i]];
      cfflo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffhi[i][j] = new double[deg+1];
         cfflo[i][j] = new double[deg+1];

         random_dbl2_exponential(deg,&rndhi,&rndlo,cffhi[i][j],cfflo[i][j]);
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
               cout << cffhi[i][j][k] << "  " << cfflo[i][j][k] << endl;
         }
      }
   }
   return 0;
}

int cmplx2_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehi, double **cstrelo, double **cstimhi, double **cstimlo,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
   int vrblvl )
{
   const int degp1 = deg+1;
   double rndrehi,rndrelo,rndimhi,rndimlo;

   for(int i=0; i<dim; i++)
   {
      cstrehi[i] = new double[degp1];
      cstrelo[i] = new double[degp1];
      cstimhi[i] = new double[degp1];
      cstimlo[i] = new double[degp1];
/*
      random_cmplx2_exponential
         (deg,&rndrehi,&rndrelo,&rndimhi,&rndimlo,
          cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i]);
 */
      for(int k=0; k<=deg; k++)
         random_double_double_complex
            (&cstrehi[i][k],&cstrelo[i][k],&cstimhi[i][k],&cstimlo[i][k]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Constant coefficient series :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;

         for(int j=0; j<=deg; j++)
            cout << cstrehi[i][j] << "  " << cstrelo[i][j] << endl
                 << cstimhi[i][j] << "  " << cstimlo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffrehi[i] = new double*[nbr[i]];
      cffrelo[i] = new double*[nbr[i]];
      cffimhi[i] = new double*[nbr[i]];
      cffimlo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffrehi[i][j] = new double[degp1];
         cffrelo[i][j] = new double[degp1];
         cffimhi[i][j] = new double[degp1];
         cffimlo[i][j] = new double[degp1];
 /*
         random_cmplx2_exponential
            (deg,&rndrehi,&rndrelo,&rndimhi,&rndimlo,
             cffrehi[i][j],cffrelo[i][j],cffimhi[i][j],cffimlo[i][j]);
  */
         for(int k=0; k<=deg; k++)
            random_double_double_complex
               (&cffrehi[i][j][k],&cffrelo[i][j][k],
                &cffimhi[i][j][k],&cffimlo[i][j][k]);
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
               cout << cffrehi[i][j][k] << "  " << cffrelo[i][j][k] << endl
                    << cffimhi[i][j][k] << "  " << cffimlo[i][j][k] << endl;
         }
      }
   }
   return 0;
}
