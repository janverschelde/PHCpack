// The file dbl_indexed_coefficients.cpp defines the functions specified in
// the file dbl_indexed_coefficients.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random_monomials.h"
#include "random_series.h"

using namespace std;

int dbl_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cst, double ***cff, int vrblvl )
{
   const int degp1 = deg+1;
   double rnd;

   for(int i=0; i<dim; i++)
   {
      cst[i] = new double[degp1];
      random_dbl_exponential(deg,&rnd,cst[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << cst[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cff[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cff[i][j] = new double[degp1];
         random_dbl_exponential(deg,&rnd,cff[i][j]);
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
            for(int k=0; k<=deg; k++) cout << cff[i][j][k] << endl;
         }
      }
   }
   return 0;
}

int cmplx_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstre, double **cstim, double ***cffre, double ***cffim,
   int vrblvl )
{
   const int degp1 = deg+1;
   double rndre,rndim;

   for(int i=0; i<dim; i++)
   {
      cstre[i] = new double[degp1];
      cstim[i] = new double[degp1];

      // random_cmplx_exponential(deg,&rndre,&rndim,cstre[i],cstim[i]);

      for(int k=0; k<=deg; k++)
      {
         rndre = random_angle();        
         cstre[i][k] = cos(rndre);
         cstim[i][k] = sin(rndre);
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "constant coefficient series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cstre[i][j] << "  " << cstim[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
   {
      cffre[i] = new double*[nbr[i]];
      cffim[i] = new double*[nbr[i]];
      for(int j=0; j<nbr[i]; j++)
      {
         cffre[i][j] = new double[degp1];
         cffim[i][j] = new double[degp1];

         // random_cmplx_exponential(deg,&rndre,&rndim,cffre[i][j],cffim[i][j]);

         for(int k=0; k<=deg; k++)
         {
            rndre = random_angle();
            cffre[i][j][k] = cos(rndre);
            cffim[i][j][k] = sin(rndre);
         }
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
               cout << cffre[i][j][k] << "  " << cffim[i][j][k] << endl;
         }
      }
   }
   return 0;
}
