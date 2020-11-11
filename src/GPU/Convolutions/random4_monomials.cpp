// The file random4_monomials.cpp defines functions specified 
// in random4_monomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_monomials.h"
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "random4_series.h"
#include "random4_monomials.h"

bool make_real4_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo )
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
      for(int i=0; i<=deg; i++)
         random_quad_double(&cffhihi[i],&cfflohi[i],&cffhilo[i],&cfflolo[i]);

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

bool make_complex4_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo )
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
      double rndhihi,rndlohi,rndhilo,rndlolo;
      double sinhihi,sinlohi,sinhilo,sinlolo;

      for(int i=0; i<=deg; i++)
      {
         random_quad_double
            (&rndhihi,&rndlohi,&rndhilo,&rndlolo);           // random cos

         cffrehihi[i] = rndhihi; cffrelohi[i] = rndlohi;     // cos(angle)
         cffrehilo[i] = rndhilo; cffrelolo[i] = rndlolo;
         qdf_sqrt(rndhihi,rndlohi,rndhilo,rndlolo,
                  &sinhihi,&sinlohi,&sinhilo,&sinlolo);      // cos^(angle)
         qdf_minus(&sinhihi,&sinlohi,&sinhilo,&sinlolo);     // -cos^(angle)
         qdf_inc_d(&sinhihi,&sinlohi,&sinhilo,&sinlolo,1.0); // 1-cos^2(angle)
         qdf_sqrt(sinhihi,sinlohi,sinhilo,sinlolo,
                  &cffimhihi[i],&cffimlohi[i],
                  &cffimhilo[i],&cffimlolo[i]);              // sin is sqrt
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

void make_real4_input
 ( int dim, int deg, double **datahihi, double **datalohi,
                     double **datahilo, double **datalolo )
{
   double rndhihi,rndlohi,rndhilo,rndlolo;
   double* pluxhihi = new double[deg+1];
   double* pluxlohi = new double[deg+1];
   double* pluxhilo = new double[deg+1];
   double* pluxlolo = new double[deg+1];
   double* minxhihi = new double[deg+1];
   double* minxlohi = new double[deg+1];
   double* minxhilo = new double[deg+1];
   double* minxlolo = new double[deg+1];

   for(int i=0; i<dim; i++)
   {
      random_dbl4_exponentials
         (deg,&rndhihi,&rndlohi,&rndhilo,&rndlolo,
              pluxhihi,pluxlohi,pluxhilo,pluxlolo,
              minxhihi,minxlohi,minxhilo,minxlolo);

      for(int j=0; j<=deg; j++)
      {
         datahihi[i][j] = pluxhihi[j];
         datalohi[i][j] = pluxlohi[j];
         datalolo[i][j] = pluxhilo[j];
         datalolo[i][j] = pluxlolo[j];
      }
   }
   free(pluxhihi); free(pluxlohi); free(pluxhilo); free(pluxlolo);
   free(minxhihi); free(minxlohi); free(minxhilo); free(minxlolo); 
}

void make_complex4_input
 ( int dim, int deg,
   double **datarehihi, double **datarelohi,
   double **datarehilo, double **datarelolo,
   double **dataimhihi, double **dataimlohi,
   double **dataimhilo, double **dataimlolo )
{
   double rndrehihi,rndrelohi,rndrehilo,rndrelolo;
   double rndimhihi,rndimlohi,rndimhilo,rndimlolo;
   double* pluxrehihi = new double[deg+1];
   double* pluxrelohi = new double[deg+1];
   double* pluxrehilo = new double[deg+1];
   double* pluxrelolo = new double[deg+1];
   double* pluximhihi = new double[deg+1];
   double* pluximlohi = new double[deg+1];
   double* pluximhilo = new double[deg+1];
   double* pluximlolo = new double[deg+1];
   double* minxrehihi = new double[deg+1];
   double* minxrelohi = new double[deg+1];
   double* minxrehilo = new double[deg+1];
   double* minxrelolo = new double[deg+1];
   double* minximhihi = new double[deg+1];
   double* minximlohi = new double[deg+1];
   double* minximhilo = new double[deg+1];
   double* minximlolo = new double[deg+1];

   for(int i=0; i<dim; i++)
   {
      random_cmplx4_exponentials(deg,
         &rndrehihi,&rndrelohi,&rndrehilo,&rndrelolo,
         &rndimhihi,&rndimlohi,&rndimhilo,&rndimlolo,
         pluxrehihi,pluxrelohi,pluxrehilo,pluxrelolo,
         pluximhihi,pluximlohi,pluximhilo,pluximlolo,
         minxrehihi,minxrelohi,minxrehilo,minxrelolo,
         minximhihi,minximlohi,minximhilo,minximlolo);

      for(int j=0; j<=deg; j++)
      {
         datarehihi[i][j] = pluxrehihi[j]; datarelohi[i][j] = pluxrelohi[j];
         datarehilo[i][j] = pluxrehilo[j]; datarelolo[i][j] = pluxrelolo[j];
         dataimhihi[i][j] = pluximhihi[j]; dataimlohi[i][j] = pluximlohi[j];
         dataimhilo[i][j] = pluximhilo[j]; dataimlolo[i][j] = pluximlolo[j];
      }
   }
   free(pluxrehihi); free(pluxrelohi); free(pluxrehilo); free(pluxrelolo);
   free(pluximhihi); free(pluximlohi); free(pluximhilo); free(pluximlolo);
   free(minxrehihi); free(minxrelohi); free(minxrehilo); free(minxrelolo);
   free(minximhihi); free(minximlohi); free(minximhilo); free(minximlolo); 
}
