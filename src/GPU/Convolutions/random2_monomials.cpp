// The file random2_monomials.cpp defines functions specified 
// in random2_monomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_monomials.h"
#include "double_double_functions.h"
#include "random2_vectors.h"
#include "random2_series.h"
#include "random2_monomials.h"

bool make_real2_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhi, double *cfflo )
{
   bool fail;

   if(nvr > dim)
   {
      std::cout << "ERROR: nvr = " << nvr << " > " << dim << " dim"
                << std::endl;

      return true;
   }
   else
   {
      for(int i=0; i<=deg; i++) random_double_double(&cffhi[i],&cfflo[i]);

      for(int i=0; i<nvr; i++)
      {
         exp[i] = 1 + (rand() % pwr);
         do
         {
            idx[i] = rand() % dim;
            fail = sorted_insert(i,idx);
         }
         while(fail);
      }
      return false;
   }
}

bool make_complex2_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo )
{
   bool fail;

   if(nvr > dim)
   {
      std::cout << "ERROR: nvr = " << nvr << " > " << dim << " dim"
                << std::endl;

      return true;
   }
   else
   {
      double rndhi,rndlo,sinhi,sinlo;

      for(int i=0; i<=deg; i++)
      {
         random_double_double(&rndhi,&rndlo);           // random cos

         cffrehi[i] = rndhi; cffrelo[i] = rndlo;        // cos(angle)
         ddf_sqr(rndhi,rndlo,&sinhi,&sinlo);            // cos^2(angle)
         ddf_minus(&sinhi,&sinlo);                      // -cos^2(angle)
         ddf_inc_d(&sinhi,&sinlo,1.0);                  // 1-cos^2(angle)
         ddf_sqrt(sinhi,sinlo,&cffimhi[i],&cffimlo[i]); // sin is sqrt
      }
      for(int i=0; i<nvr; i++)
      {
         exp[i] = 1 + (rand() % pwr);
         do
         {
            idx[i] = rand() % dim;
            fail = sorted_insert(i,idx);
         }
         while(fail);
      }
      return false;
   }
}

void make_real2_input ( int dim, int deg, double **datahi, double **datalo )
{
   double rndhi,rndlo;
   double* pluxhi = new double[deg+1];
   double* pluxlo = new double[deg+1];
   double* minxhi = new double[deg+1];
   double* minxlo = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_dbl2_exponentials(deg,&rndhi,&rndlo,pluxhi,pluxlo,minxhi,minxlo);
      for(int j=0; j<=deg; j++)
      {
         datahi[i][j]   = pluxhi[j]; datalo[i][j]   = pluxlo[j];
         datahi[i+1][j] = minxhi[j]; datalo[i+1][j] = minxlo[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datahi[dim-1][0] = 1.0; datalo[dim-1][0] = 0.0;
      for(int j=1; j<=deg; j++)
      {
         datahi[dim-1][j] = 0.0; datalo[dim-1][j] = 0.0;
      }
   }
   free(pluxhi); free(pluxlo); free(minxhi); free(minxlo); 
}

void make_complex2_input
 ( int dim, int deg, double **datarehi, double **datarelo,
   double **dataimhi, double **dataimlo )
{
   double rndrehi,rndrelo,rndimhi,rndimlo;
   double* pluxrehi = new double[deg+1];
   double* pluxrelo = new double[deg+1];
   double* pluximhi = new double[deg+1];
   double* pluximlo = new double[deg+1];
   double* minxrehi = new double[deg+1];
   double* minxrelo = new double[deg+1];
   double* minximhi = new double[deg+1];
   double* minximlo = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_cmplx2_exponentials(deg,
         &rndrehi,&rndrelo,&rndimhi,&rndimlo,
         pluxrehi,pluxrelo,pluximhi,pluximlo,
         minxrehi,minxrelo,minximhi,minximlo);

      for(int j=0; j<=deg; j++)
      {
         datarehi[i][j] = pluxrehi[j]; datarelo[i][j] = pluxrelo[j];
         dataimhi[i][j] = pluximhi[j]; dataimlo[i][j] = pluximlo[j];
         datarehi[i+1][j] = minxrehi[j]; datarelo[i+1][j] = minxrelo[j];
         dataimhi[i+1][j] = minximhi[j]; dataimlo[i+1][j] = minximlo[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datarehi[dim-1][0] = 1.0; datarelo[dim-1][0] = 0.0;
      dataimhi[dim-1][0] = 0.0; dataimlo[dim-1][0] = 0.0;
      for(int j=1; j<=deg; j++)
      {
         datarehi[dim-1][j] = 0.0; datarelo[dim-1][j] = 0.0;
         dataimhi[dim-1][j] = 0.0; dataimlo[dim-1][j] = 0.0;
      }
   }
   free(pluxrehi); free(pluxrelo); free(pluximhi); free(pluximlo);
   free(minxrehi); free(minxrelo); free(minximhi); free(minximlo); 
}
