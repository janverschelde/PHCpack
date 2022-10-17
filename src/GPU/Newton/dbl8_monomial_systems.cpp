// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"
#include "random8_vectors.h"
#include "random8_series.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomial_systems.h"

using namespace std;

void make_complex8_exponentials
 ( int dim, int deg,
   double **srehihihi, double **srelohihi,
   double **srehilohi, double **srelolohi,
   double **srehihilo, double **srelohilo,
   double **srehilolo, double **srelololo,
   double **simhihihi, double **simlohihi,
   double **simhilohi, double **simlolohi,
   double **simhihilo, double **simlohilo,
   double **simhilolo, double **simlololo )
{
   double xrehihihi,xrelohihi,xrehilohi,xrelolohi;
   double xrehihilo,xrelohilo,xrehilolo,xrelololo;
   double ximhihihi,ximlohihi,ximhilohi,ximlolohi;
   double ximhihilo,ximlohilo,ximhilolo,ximlololo;
   double yhihihi,ylohihi,yhilohi,ylolohi;
   double yhihilo,ylohilo,yhilolo,ylololo;

   for(int i=0; i<dim; i++)
   {                                            // cosine of some angle
      random_octo_double
         (&xrehihihi,&xrelohihi,&xrehilohi,&xrelolohi,
          &xrehihilo,&xrelohilo,&xrehilolo,&xrelololo);
 
      odf_sqr(xrehihihi,xrelohihi,xrehilohi,xrelolohi,
              xrehihilo,xrelohilo,xrehilolo,xrelololo,
              &yhihihi,&ylohihi,&yhilohi,&ylolohi,
              &yhihilo,&ylohilo,&yhilolo,&ylololo);        // y = cos^2
      odf_minus(&yhihihi,&ylohihi,&yhilohi,&ylolohi,
                &yhihilo,&ylohilo,&yhilolo,&ylololo);      // y = -cos^2
      odf_inc_d(&yhihihi,&ylohihi,&yhilohi,&ylolohi,
                &yhihilo,&ylohilo,&yhilolo,&ylololo,1.0);  // y = 1 - cos^2
      odf_sqrt(yhihihi,ylohihi,yhilohi,ylolohi,
               yhihilo,ylohilo,yhilolo,ylololo,
               &ximhihihi,&ximlohihi,&ximhilohi,&ximlolohi,
               &ximhihilo,&ximlohilo,&ximhilolo,&ximlololo);
     // sin is sqrt(1-cos^2)

      cmplx8_exponential
         (deg,xrehihihi,xrelohihi,xrehilohi,xrelolohi,
              xrehihilo,xrelohilo,xrehilolo,xrelololo,
              ximhihihi,ximlohihi,ximhilohi,ximlolohi,
              ximhihilo,ximlohilo,ximhilolo,ximlololo,
          srehihihi[i],srelohihi[i],srehilohi[i],srelolohi[i],
          srehihilo[i],srelohilo[i],srehilolo[i],srelololo[i],
          simhihihi[i],simlohihi[i],simhilohi[i],simlolohi[i],
          simhihilo[i],simlohilo[i],simhilolo[i],simlololo[i]);
   }
}

void evaluate_complex8_monomials
 ( int dim, int deg, int **rowsA,
   double **srehihihi, double **srelohihi,
   double **srehilohi, double **srelolohi,
   double **srehihilo, double **srelohilo,
   double **srehilolo, double **srelololo,
   double **simhihihi, double **simlohihi,
   double **simhilohi, double **simlolohi,
   double **simhihilo, double **simlohilo,
   double **simhilolo, double **simlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo )
{
   const int degp1 = deg+1;

   double *accrehihihi = new double[degp1]; // accumulates product
   double *accrelohihi = new double[degp1];
   double *accrehilohi = new double[degp1];
   double *accrelolohi = new double[degp1];
   double *accrehihilo = new double[degp1];
   double *accrelohilo = new double[degp1];
   double *accrehilolo = new double[degp1];
   double *accrelololo = new double[degp1];
   double *accimhihihi = new double[degp1];
   double *accimlohihi = new double[degp1];
   double *accimhilohi = new double[degp1];
   double *accimlolohi = new double[degp1];
   double *accimhihilo = new double[degp1];
   double *accimlohilo = new double[degp1];
   double *accimhilolo = new double[degp1];
   double *accimlololo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx8_product
                  (deg,srehihihi[j],srelohihi[j],srehilohi[j],srelolohi[j],
                       srehihilo[j],srelohilo[j],srehilolo[j],srelololo[j],
                       simhihihi[j],simlohihi[j],simhilohi[j],simlolohi[j],
                       simhihilo[j],simlohilo[j],simhilolo[j],simlololo[j],
                   rhsrehihihi[i],rhsrelohihi[i],rhsrehilohi[i],rhsrelolohi[i],
                   rhsrehihilo[i],rhsrelohilo[i],rhsrehilolo[i],rhsrelololo[i],
                   rhsimhihihi[i],rhsimlohihi[i],rhsimhilohi[i],rhsimlolohi[i],
                   rhsimhihilo[i],rhsimlohilo[i],rhsimhilolo[i],rhsimlololo[i],
                   accrehihihi,accrelohihi,accrehilohi,accrelolohi,
                   accrehihilo,accrelohilo,accrehilolo,accrelololo,
                   accimhihihi,accimlohihi,accimhilohi,accimlolohi,
                   accimhihilo,accimlohilo,accimhilolo,accimlololo);

               for(int L=0; L<degp1; L++)
               {
                  rhsrehihihi[i][L] = accrehihihi[L];
                  rhsrelohihi[i][L] = accrelohihi[L];
                  rhsrehilohi[i][L] = accrehilohi[L];
                  rhsrelolohi[i][L] = accrelolohi[L];
                  rhsrehihilo[i][L] = accrehihilo[L];
                  rhsrelohilo[i][L] = accrelohilo[L];
                  rhsrehilolo[i][L] = accrehilolo[L];
                  rhsrelololo[i][L] = accrelololo[L];
                  rhsimhihihi[i][L] = accimhihihi[L];
                  rhsimlohihi[i][L] = accimlohihi[L];
                  rhsimhilohi[i][L] = accimhilohi[L];
                  rhsimlolohi[i][L] = accimlolohi[L];
                  rhsimhihilo[i][L] = accimhihilo[L];
                  rhsimlohilo[i][L] = accimlohilo[L];
                  rhsimhilolo[i][L] = accimhilolo[L];
                  rhsimlololo[i][L] = accimlololo[L];
               }
            }
         }
      }
   }
   free(accrehihihi); free(accrelohihi); free(accrehilohi); free(accrelolohi);
   free(accrehihilo); free(accrelohilo); free(accrehilolo); free(accrelololo);
   free(accimhihihi); free(accimlohihi); free(accimhilohi); free(accimlolohi);
   free(accimhihilo); free(accimlohilo); free(accimhilolo); free(accimlololo);
}
