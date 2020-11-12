// The file random_monomials.cpp defines functions specified 
// in random_monomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random_series.h"
#include "random_monomials.h"

bool sorted_insert ( int n, int *data )
{
   if(n == 0)
      return false;
   else
   {
      int nbr = data[n];
      int idx = n;

      for(int i=0; i<n; i++)
         if(data[i] >= nbr)
         {
            idx = i; break;
         }

      if(idx == n)
         return false;                  // sequence is already sorted
      else
      {
         if(data[idx] == nbr)           // found duplicate number
            return true;
         else
         {
            for(int i=n; i>idx; i--)
               data[i] = data[i-1];     // shift the numbers

            data[idx] = nbr;            // insert number
            return false;
         }
      }
   }
}

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp, double *cff )
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
      for(int i=0; i<=deg; i++) cff[i] = random_double();

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

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffre, double *cffim )
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
      double rnd;

      for(int i=0; i<=deg; i++)
      {
         rnd = random_angle();        
         cffre[i] = cos(rnd);
         cffim[i] = sin(rnd);
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

void common_factors ( int nvr, int *exp, int *nbrfac, int *expfac )
{
   *nbrfac = 0;

   for(int i=0; i<nvr; i++)
   {
      if(exp[i] <= 1)
         expfac[i] = 0;
      else
      {
         expfac[i] = exp[i] - 1;
         *nbrfac = *nbrfac + 1;
      }
   }
}

void make_real_input ( int dim, int deg, double **data )
{
   double rnd;
   double* plux = new double[deg+1];
   double* minx = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_dbl_exponentials(deg,&rnd,plux,minx);

      for(int j=0; j<=deg; j++)
      {
         data[i][j] = plux[j]; data[i+1][j] = minx[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      data[dim-1][0] = 1.0;
      for(int j=1; j<=deg; j++) data[dim-1][j] = 0.0;
   }
   free(plux); free(minx);
}

void make_complex_input
 ( int dim, int deg, double **datare, double **dataim )
{
   double rndre,rndim;
   double* pluxre = new double[deg+1];
   double* pluxim = new double[deg+1];
   double* minxre = new double[deg+1];
   double* minxim = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_cmplx_exponentials
         (deg,&rndre,&rndim,pluxre,pluxim,minxre,minxim);

      for(int j=0; j<=deg; j++)
      {
         datare[i][j]   = pluxre[j]; dataim[i][j]   = pluxim[j];
         datare[i+1][j] = minxre[j]; dataim[i+1][j] = minxim[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datare[dim-1][0] = 1.0; dataim[dim-1][0] = 0.0;
      for(int j=1; j<=deg; j++)
      {
         datare[dim-1][j] = 0.0;
         dataim[dim-1][j] = 0.0;
      }
   }
   free(pluxre); free(pluxim); free(minxre); free(minxim);
}
