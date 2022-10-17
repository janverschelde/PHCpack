// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "random4_series.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomial_systems.h"

using namespace std;

void make_complex4_exponentials
 ( int dim, int deg,
   double **srehihi, double **srelohi, double **srehilo, double **srelolo,
   double **simhihi, double **simlohi, double **simhilo, double **simlolo )
{
   double xrehihi,xrelohi,xrehilo,xrelolo;
   double ximhihi,ximlohi,ximhilo,ximlolo;
   double yhihi,ylohi,yhilo,ylolo;

   for(int i=0; i<dim; i++)
   {                                            // cosine of some angle
      random_quad_double(&xrehihi,&xrelohi,&xrehilo,&xrelolo);
 
      qdf_sqr(xrehihi,xrelohi,xrehilo,xrelolo,
              &yhihi,&ylohi,&yhilo,&ylolo);        // y = cos^2
      qdf_minus(&yhihi,&ylohi,&yhilo,&ylolo);      // y = -cos^2
      qdf_inc_d(&yhihi,&ylohi,&yhilo,&ylolo,1.0);  // y = 1 - cos^2
      qdf_sqrt(yhihi,ylohi,yhilo,ylolo,
               &ximhihi,&ximlohi,&ximhilo,&ximlolo); // sin is sqrt(1-cos^2)

      cmplx4_exponential
         (deg,xrehihi,xrelohi,xrehilo,xrelolo,
              ximhihi,ximlohi,ximhilo,ximlolo,
          srehihi[i],srelohi[i],srehilo[i],srelolo[i],
          simhihi[i],simlohi[i],simhilo[i],simlolo[i]);
   }
}

void evaluate_complex4_monomials
 ( int dim, int deg, int **rowsA,
   double **srehihi, double **srelohi, double **srehilo, double **srelolo,
   double **simhihi, double **simlohi, double **simhilo, double **simlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo )
{
   const int degp1 = deg+1;

   double *accrehihi = new double[degp1]; // accumulates product
   double *accrelohi = new double[degp1];
   double *accrehilo = new double[degp1];
   double *accrelolo = new double[degp1];
   double *accimhihi = new double[degp1];
   double *accimlohi = new double[degp1];
   double *accimhilo = new double[degp1];
   double *accimlolo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx4_product
                  (deg,srehihi[j],srelohi[j],srehilo[j],srelolo[j],
                       simhihi[j],simlohi[j],simhilo[j],simlolo[j],
                   rhsrehihi[i],rhsrelohi[i],rhsrehilo[i],rhsrelolo[i],
                   rhsimhihi[i],rhsimlohi[i],rhsimhilo[i],rhsimlolo[i],
                   accrehihi,accrelohi,accrehilo,accrelolo,
                   accimhihi,accimlohi,accimhilo,accimlolo);

               for(int L=0; L<degp1; L++)
               {
                  rhsrehihi[i][L] = accrehihi[L];
                  rhsrelohi[i][L] = accrelohi[L];
                  rhsrehilo[i][L] = accrehilo[L];
                  rhsrelolo[i][L] = accrelolo[L];
                  rhsimhihi[i][L] = accimhihi[L];
                  rhsimlohi[i][L] = accimlohi[L];
                  rhsimhilo[i][L] = accimhilo[L];
                  rhsimlolo[i][L] = accimlolo[L];
               }
            }
         }
      }
   }
   free(accrehihi); free(accrelohi); free(accrehilo); free(accrelolo);
   free(accimhihi); free(accimlohi); free(accimhilo); free(accimlolo);
}
