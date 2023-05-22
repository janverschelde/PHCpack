// The file random8_monomials.cpp defines functions specified 
// in random8_monomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_monomials.h"
#include "octo_double_functions.h"
#include "random8_vectors.h"
#include "random8_series.h"
#include "random8_monomials.h"

bool make_real8_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo )
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
         random_octo_double
            (&cffhihihi[i],&cfflohihi[i],&cffhilohi[i],&cfflolohi[i],
             &cffhihilo[i],&cfflohilo[i],&cffhilolo[i],&cfflololo[i]);

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

void random_octo_complex
 ( double *rehihihi, double *relohihi, double *rehilohi, double *relolohi,
   double *rehihilo, double *relohilo, double *rehilolo, double *relololo,
   double *imhihihi, double *imlohihi, double *imhilohi, double *imlolohi,
   double *imhihilo, double *imlohilo, double *imhilolo, double *imlololo )
{
   double sinhihihi,sinlohihi,sinhilohi,sinlolohi;
   double sinhihilo,sinlohilo,sinhilolo,sinlololo;

   random_octo_double                                           // random cos
      (rehihihi,relohihi,rehilohi,relolohi,
       rehihilo,relohilo,rehilolo,relololo);
                                                                // cos(angle)
   odf_sqrt(*rehihihi,*relohihi,*rehilohi,*relolohi,
            *rehihilo,*relohilo,*rehilolo,*relololo,
            &sinhihihi,&sinlohihi,&sinhilohi,&sinlolohi,
            &sinhihilo,&sinlohilo,&sinhilolo,&sinlololo);      // cos^(angle)
   odf_minus(&sinhihihi,&sinlohihi,&sinhilohi,&sinlolohi,
             &sinhihilo,&sinlohilo,&sinhilolo,&sinlololo);
                                                              // -cos^(angle)
   odf_inc_d(&sinhihihi,&sinlohihi,&sinhilohi,&sinlolohi,
             &sinhihilo,&sinlohilo,&sinhilolo,&sinlololo,1.0);
                                                            // 1-cos^2(angle)
   odf_sqrt(sinhihihi,sinlohihi,sinhilohi,sinlolohi,
            sinhihilo,sinlohilo,sinhilolo,sinlololo,
            imhihihi,imlohihi,imhilohi,imlolohi,
            imhihilo,imlohilo,imhilolo,imlololo);              // sin is sqrt
}

bool make_complex8_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo )
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
      double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
      double rndhihilo,rndlohilo,rndhilolo,rndlololo;
      double sinhihihi,sinlohihi,sinhilohi,sinlolohi;
      double sinhihilo,sinlohilo,sinhilolo,sinlololo;

      for(int i=0; i<=deg; i++)
      {
         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);      // random cos

         cffrehihihi[i] = rndhihihi; cffrelohihi[i] = rndlohihi;
         cffrehilohi[i] = rndhilohi; cffrelolohi[i] = rndlolohi;
         cffrehihilo[i] = rndhihilo; cffrelohilo[i] = rndlohilo;
         cffrehilolo[i] = rndhilolo; cffrelololo[i] = rndlololo;
                                                                 // cos(angle)
         odf_sqrt(rndhihihi,rndlohihi,rndhilohi,rndlolohi,
                  rndhihilo,rndlohilo,rndhilolo,rndlololo,
                  &sinhihihi,&sinlohihi,&sinhilohi,&sinlolohi,
                  &sinhihilo,&sinlohilo,&sinhilolo,&sinlololo); // cos^(angle)
         odf_minus(&sinhihihi,&sinlohihi,&sinhilohi,&sinlolohi,
                   &sinhihilo,&sinlohilo,&sinhilolo,&sinlololo);
                                                               // -cos^(angle)
         odf_inc_d(&sinhihihi,&sinlohihi,&sinhilohi,&sinlolohi,
                   &sinhihilo,&sinlohilo,&sinhilolo,&sinlololo,1.0);
                                                             // 1-cos^2(angle)
         odf_sqrt(sinhihihi,sinlohihi,sinhilohi,sinlolohi,
                  sinhihilo,sinlohilo,sinhilolo,sinlololo,
                  &cffimhihihi[i],&cffimlohihi[i],
                  &cffimhilohi[i],&cffimlolohi[i],
                  &cffimhihilo[i],&cffimlohilo[i],
                  &cffimhilolo[i],&cffimlololo[i]);            // sin is sqrt
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

void make_real8_input
 ( int dim, int deg,
   double **datahihihi, double **datalohihi,
   double **datahilohi, double **datalolohi,
   double **datahihilo, double **datalohilo,
   double **datahilolo, double **datalololo )
{
   double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
   double rndhihilo,rndlohilo,rndhilolo,rndlololo;
   double* pluxhihihi = new double[deg+1];
   double* pluxlohihi = new double[deg+1];
   double* pluxhilohi = new double[deg+1];
   double* pluxlolohi = new double[deg+1];
   double* pluxhihilo = new double[deg+1];
   double* pluxlohilo = new double[deg+1];
   double* pluxhilolo = new double[deg+1];
   double* pluxlololo = new double[deg+1];
   double* minxhihihi = new double[deg+1];
   double* minxlohihi = new double[deg+1];
   double* minxhilohi = new double[deg+1];
   double* minxlolohi = new double[deg+1];
   double* minxhihilo = new double[deg+1];
   double* minxlohilo = new double[deg+1];
   double* minxhilolo = new double[deg+1];
   double* minxlololo = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_dbl8_exponentials
         (deg,&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
              &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo,
              pluxhihihi,pluxlohihi,pluxhilohi,pluxlolohi,
              pluxhihilo,pluxlohilo,pluxhilolo,pluxlololo,
              minxhihihi,minxlohihi,minxhilohi,minxlolohi,
              minxhihilo,minxlohilo,minxhilolo,minxlololo);

      for(int j=0; j<=deg; j++)
      {
         datahihihi[i][j] = pluxhihihi[j]; datalohihi[i][j] = pluxlohihi[j];
         datahilohi[i][j] = pluxhilohi[j]; datalolohi[i][j] = pluxlolohi[j];
         datahihilo[i][j] = pluxhihilo[j]; datalohilo[i][j] = pluxlohilo[j];
         datahilolo[i][j] = pluxhilolo[j]; datalololo[i][j] = pluxlololo[j];
         datahihihi[i+1][j] = minxhihihi[j];
         datalohihi[i+1][j] = minxlohihi[j];
         datahilohi[i+1][j] = minxhilohi[j];
         datalolohi[i+1][j] = minxlolohi[j];
         datahihilo[i+1][j] = minxhihilo[j];
         datalohilo[i+1][j] = minxlohilo[j];
         datahilolo[i+1][j] = minxhilolo[j];
         datalololo[i+1][j] = minxlololo[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datahihihi[dim-1][0] = 1.0; datalohihi[dim-1][0] = 0.0;
      datahilohi[dim-1][0] = 0.0; datalolohi[dim-1][0] = 0.0;
      datahihilo[dim-1][0] = 0.0; datalohilo[dim-1][0] = 0.0;
      datahilolo[dim-1][0] = 0.0; datalololo[dim-1][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         datahihihi[dim-1][j] = 0.0; datalohihi[dim-1][j] = 0.0;
         datahilohi[dim-1][j] = 0.0; datalolohi[dim-1][j] = 0.0;
         datahihilo[dim-1][j] = 0.0; datalohilo[dim-1][j] = 0.0;
         datahilolo[dim-1][j] = 0.0; datalololo[dim-1][j] = 0.0;
      }
   }
   free(pluxhihihi); free(pluxlohihi); free(pluxhilohi); free(pluxlolohi);
   free(pluxhihilo); free(pluxlohilo); free(pluxhilolo); free(pluxlololo);
   free(minxhihihi); free(minxlohihi); free(minxhilohi); free(minxlolohi); 
   free(minxhihilo); free(minxlohilo); free(minxhilolo); free(minxlololo); 
}

void make_complex8_input
 ( int dim, int deg,
   double **datarehihihi, double **datarelohihi,
   double **datarehilohi, double **datarelolohi,
   double **datarehihilo, double **datarelohilo,
   double **datarehilolo, double **datarelololo,
   double **dataimhihihi, double **dataimlohihi,
   double **dataimhilohi, double **dataimlolohi,
   double **dataimhihilo, double **dataimlohilo,
   double **dataimhilolo, double **dataimlololo )
{
   double rndrehihihi,rndrelohihi,rndrehilohi,rndrelolohi;
   double rndrehihilo,rndrelohilo,rndrehilolo,rndrelololo;
   double rndimhihihi,rndimlohihi,rndimhilohi,rndimlolohi;
   double rndimhihilo,rndimlohilo,rndimhilolo,rndimlololo;
   double* pluxrehihihi = new double[deg+1];
   double* pluxrelohihi = new double[deg+1];
   double* pluxrehilohi = new double[deg+1];
   double* pluxrelolohi = new double[deg+1];
   double* pluxrehihilo = new double[deg+1];
   double* pluxrelohilo = new double[deg+1];
   double* pluxrehilolo = new double[deg+1];
   double* pluxrelololo = new double[deg+1];
   double* pluximhihihi = new double[deg+1];
   double* pluximlohihi = new double[deg+1];
   double* pluximhilohi = new double[deg+1];
   double* pluximlolohi = new double[deg+1];
   double* pluximhihilo = new double[deg+1];
   double* pluximlohilo = new double[deg+1];
   double* pluximhilolo = new double[deg+1];
   double* pluximlololo = new double[deg+1];
   double* minxrehihihi = new double[deg+1];
   double* minxrelohihi = new double[deg+1];
   double* minxrehilohi = new double[deg+1];
   double* minxrelolohi = new double[deg+1];
   double* minxrehihilo = new double[deg+1];
   double* minxrelohilo = new double[deg+1];
   double* minxrehilolo = new double[deg+1];
   double* minxrelololo = new double[deg+1];
   double* minximhihihi = new double[deg+1];
   double* minximlohihi = new double[deg+1];
   double* minximhilohi = new double[deg+1];
   double* minximlolohi = new double[deg+1];
   double* minximhihilo = new double[deg+1];
   double* minximlohilo = new double[deg+1];
   double* minximhilolo = new double[deg+1];
   double* minximlololo = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_cmplx8_exponentials(deg,
         &rndrehihihi,&rndrelohihi,&rndrehilohi,&rndrelolohi,
         &rndrehihilo,&rndrelohilo,&rndrehilolo,&rndrelololo,
         &rndimhihihi,&rndimlohihi,&rndimhilohi,&rndimlolohi,
         &rndimhihilo,&rndimlohilo,&rndimhilolo,&rndimlololo,
         pluxrehihihi,pluxrelohihi,pluxrehilohi,pluxrelolohi,
         pluxrehihilo,pluxrelohilo,pluxrehilolo,pluxrelololo,
         pluximhihihi,pluximlohihi,pluximhilohi,pluximlolohi,
         pluximhihilo,pluximlohilo,pluximhilolo,pluximlololo,
         minxrehihihi,minxrelohihi,minxrehilohi,minxrelolohi,
         minxrehihilo,minxrelohilo,minxrehilolo,minxrelololo,
         minximhihihi,minximlohihi,minximhilohi,minximlolohi,
         minximhihilo,minximlohilo,minximhilolo,minximlololo);

      for(int j=0; j<=deg; j++)
      {
         datarehihihi[i][j] = pluxrehihihi[j];
         datarelohihi[i][j] = pluxrelohihi[j];
         datarehilohi[i][j] = pluxrehilohi[j];
         datarelolohi[i][j] = pluxrelolohi[j];
         datarehihilo[i][j] = pluxrehihilo[j];
         datarelohilo[i][j] = pluxrelohilo[j];
         datarehilolo[i][j] = pluxrehilolo[j];
         datarelololo[i][j] = pluxrelololo[j];
         dataimhihihi[i][j] = pluximhihihi[j];
         dataimlohihi[i][j] = pluximlohihi[j];
         dataimhilohi[i][j] = pluximhilohi[j];
         dataimlolohi[i][j] = pluximlolohi[j];
         dataimhihilo[i][j] = pluximhihilo[j];
         dataimlohilo[i][j] = pluximlohilo[j];
         dataimhilolo[i][j] = pluximhilolo[j];
         dataimlololo[i][j] = pluximlololo[j];
         datarehihihi[i+1][j] = minxrehihihi[j];
         datarelohihi[i+1][j] = minxrelohihi[j];
         datarehilohi[i+1][j] = minxrehilohi[j];
         datarelolohi[i+1][j] = minxrelolohi[j];
         datarehihilo[i+1][j] = minxrehihilo[j];
         datarelohilo[i+1][j] = minxrelohilo[j];
         datarehilolo[i+1][j] = minxrehilolo[j];
         datarelololo[i+1][j] = minxrelololo[j];
         dataimhihihi[i+1][j] = minximhihihi[j];
         dataimlohihi[i+1][j] = minximlohihi[j];
         dataimhilohi[i+1][j] = minximhilohi[j];
         dataimlolohi[i+1][j] = minximlolohi[j];
         dataimhihilo[i+1][j] = minximhihilo[j];
         dataimlohilo[i+1][j] = minximlohilo[j];
         dataimhilolo[i+1][j] = minximhilolo[j];
         dataimlololo[i+1][j] = minximlololo[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datarehihihi[dim-1][0] = 1.0; dataimhihihi[dim-1][0] = 0.0;
      datarelohihi[dim-1][0] = 0.0; dataimlohihi[dim-1][0] = 0.0;
      datarehilohi[dim-1][0] = 0.0; dataimhilohi[dim-1][0] = 0.0;
      datarelolohi[dim-1][0] = 0.0; dataimlolohi[dim-1][0] = 0.0;
      datarehihilo[dim-1][0] = 0.0; dataimhihilo[dim-1][0] = 0.0;
      datarelohilo[dim-1][0] = 0.0; dataimlohilo[dim-1][0] = 0.0;
      datarehilolo[dim-1][0] = 0.0; dataimhilolo[dim-1][0] = 0.0;
      datarelololo[dim-1][0] = 0.0; dataimlololo[dim-1][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         datarehihihi[dim-1][j] = 0.0; dataimhihihi[dim-1][j] = 0.0;
         datarelohihi[dim-1][j] = 0.0; dataimlohihi[dim-1][j] = 0.0;
         datarehilohi[dim-1][j] = 0.0; dataimhilohi[dim-1][j] = 0.0;
         datarelolohi[dim-1][j] = 0.0; dataimlolohi[dim-1][j] = 0.0;
         datarehihilo[dim-1][j] = 0.0; dataimhihilo[dim-1][j] = 0.0;
         datarelohilo[dim-1][j] = 0.0; dataimlohilo[dim-1][j] = 0.0;
         datarehilolo[dim-1][j] = 0.0; dataimhilolo[dim-1][j] = 0.0;
         datarelololo[dim-1][j] = 0.0; dataimlololo[dim-1][j] = 0.0;
      }
   }
   free(pluxrehihihi); free(pluxrelohihi);
   free(pluxrehilohi); free(pluxrelolohi);
   free(pluxrehihilo); free(pluxrelohilo);
   free(pluxrehilolo); free(pluxrelololo);
   free(pluximhihihi); free(pluximlohihi);
   free(pluximhilohi); free(pluximlolohi);
   free(pluximhihilo); free(pluximlohilo);
   free(pluximhilolo); free(pluximlololo);
   free(minxrehihihi); free(minxrelohihi);
   free(minxrehilohi); free(minxrelolohi);
   free(minxrehihilo); free(minxrelohilo);
   free(minxrehilolo); free(minxrelololo);
   free(minximhihihi); free(minximlohihi);
   free(minximhilohi); free(minximlolohi); 
   free(minximhihilo); free(minximlohilo);
   free(minximhilolo); free(minximlololo); 
}
