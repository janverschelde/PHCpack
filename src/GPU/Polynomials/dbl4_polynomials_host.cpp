/* The file dbl4_polynomials_host.cpp defines functions specified
 * in dbl4_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "quad_double_functions.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"
#include "dbl4_polynomials_host.h"

void CPU_dbl4_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl4_product(deg,
            inputhihi[ix1],inputlohi[ix1],inputhilo[ix1],inputlolo[ix1],
              cffhihi[i],    cfflohi[i],    cffhilo[i],    cfflolo[i],
            forwardhihi[0],forwardlohi[0],forwardhilo[0],forwardlolo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            qdf_inc(&outputhihi[dim][j],&outputlohi[dim][j],
                    &outputhilo[dim][j],&outputlolo[dim][j],
                    forwardhihi[0][j],  forwardlohi[0][j],
                    forwardhilo[0][j],  forwardlolo[0][j]);
            // output[ix1][j] += cff[i][j];
            qdf_inc(&outputhihi[ix1][j],&outputlohi[ix1][j],
                    &outputhilo[ix1][j],&outputlolo[ix1][j],
                        cffhihi[i][j],      cfflohi[i][j],
                        cffhilo[i][j],      cfflolo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl4_product(deg,
                cffhihi[i],    cfflohi[i],    cffhilo[i],    cfflolo[i],
              inputhihi[ix1],inputlohi[ix1],inputhilo[ix1],inputlolo[ix1],
            forwardhihi[0],forwardlohi[0],forwardhilo[0],forwardlolo[0]);
         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
            qdf_inc(&outputhihi[ix2][j],&outputlohi[ix2][j],
                    &outputhilo[ix2][j],&outputlolo[ix2][j],
                    forwardhihi[0][j],   forwardlohi[0][j],
                    forwardhilo[0][j],   forwardlolo[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl4_product(deg,
                 cffhihi[i],     cfflohi[i],     cffhilo[i],     cfflolo[i],
               inputhihi[ix2], inputlohi[ix2], inputhilo[ix2], inputlolo[ix2],
            backwardhihi[0],backwardlohi[0],backwardhilo[0],backwardlolo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
            qdf_inc( &outputhihi[ix1][j],&outputlohi[ix1][j],
                     &outputhilo[ix1][j],&outputlolo[ix1][j],
                    backwardhihi[0][j], backwardlohi[0][j],
                    backwardhilo[0][j], backwardlolo[0][j]);

         CPU_dbl4_product(deg,
            forwardhihi[0],forwardlohi[0],forwardhilo[0],forwardlolo[0],
              inputhihi[ix2],inputlohi[ix2],inputhilo[ix2],inputlolo[ix2],
            forwardhihi[1],forwardlohi[1],forwardhilo[1],forwardlolo[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
            qdf_inc(&outputhihi[dim][j],&outputlohi[dim][j],
                    &outputhilo[dim][j],&outputlolo[dim][j],
                    forwardhihi[1][j],  forwardlohi[1][j],
                    forwardhilo[1][j],  forwardlolo[1][j]);
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl4_speel(nvr[i],deg,idx[i],
                 cffhihi[i],  cfflohi[i],  cffhilo[i],  cfflolo[i],
               inputhihi,   inputlohi,   inputhilo,   inputlolo,
             forwardhihi, forwardlohi, forwardhilo, forwardlolo,
            backwardhihi,backwardlohi,backwardhilo,backwardlolo,
               crosshihi,   crosslohi,   crosshilo,   crosslolo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
            qdf_inc(&outputhihi[dim][j],&outputlohi[dim][j],
                    &outputhilo[dim][j],&outputlolo[dim][j],
                    forwardhihi[ix1][j],forwardlohi[ix1][j],
                    forwardhilo[ix1][j],forwardlolo[ix1][j]);

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
            qdf_inc(&outputhihi[ix2][j],&outputlohi[ix2][j],
                    &outputhilo[ix2][j],&outputlolo[ix2][j],
                    forwardhihi[ix1][j],forwardlohi[ix1][j],
                    forwardhilo[ix1][j],forwardlolo[ix1][j]);

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
            qdf_inc( &outputhihi[ix2][j], &outputlohi[ix2][j],
                     &outputhilo[ix2][j], &outputlolo[ix2][j],
                    backwardhihi[ix1][j],backwardlohi[ix1][j],
                    backwardhilo[ix1][j],backwardlolo[ix1][j]);

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
               qdf_inc(&outputhihi[ix2][j],&outputlohi[ix2][j],
                       &outputhilo[ix2][j],&outputlolo[ix2][j],
                         crosshihi[k-1][j],  crosslohi[k-1][j],
                         crosshilo[k-1][j],  crosslolo[k-1][j]);
         }
      }
   }
}

void CPU_cmplx4_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffrehihi, double **cffrelohi, 
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **backwardrehihi, double **backwardrelohi,
   double **backwardrehilo, double **backwardrelolo,
   double **backwardimhihi, double **backwardimlohi,
   double **backwardimhilo, double **backwardimlolo,
   double **crossrehihi, double **crossrelohi,
   double **crossrehilo, double **crossrelolo,
   double **crossimhihi, double **crossimlohi,
   double **crossimhilo, double **crossimlolo, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_cmplx4_product(deg,
            inputrehihi[ix1],inputrelohi[ix1],
            inputrehilo[ix1],inputrelolo[ix1],
            inputimhihi[ix1],inputimlohi[ix1],
            inputimhilo[ix1],inputimlolo[ix1],
            cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
            cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
            forwardrehihi[0],forwardrelohi[0],
            forwardrehilo[0],forwardrelolo[0],
            forwardimhihi[0],forwardimlohi[0],
            forwardimhilo[0],forwardimlolo[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            qdf_inc
               (&outputrehihi[dim][j],&outputrelohi[dim][j],
                &outputrehilo[dim][j],&outputrelolo[dim][j],
                forwardrehihi[0][j],forwardrelohi[0][j],
                forwardrehilo[0][j],forwardrelolo[0][j]);
            qdf_inc
               (&outputimhihi[dim][j],&outputimlohi[dim][j],
                &outputimhilo[dim][j],&outputimlolo[dim][j],
                forwardimhihi[0][j],forwardimlohi[0][j],
                forwardimhilo[0][j],forwardimlolo[0][j]);
            // output[ix1][j] += cff[i][j];
            qdf_inc
               (&outputrehihi[ix1][j],&outputrelohi[ix1][j],
                &outputrehilo[ix1][j],&outputrelolo[ix1][j],
                cffrehihi[i][j],cffrelohi[i][j],
                cffrehilo[i][j],cffrelolo[i][j]);
            qdf_inc
               (&outputimhihi[ix1][j],&outputimlohi[ix1][j],
                &outputimhilo[ix1][j],&outputimlolo[ix1][j],
                cffimhihi[i][j],cffimlohi[i][j],
                cffimhilo[i][j],cffimlolo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_cmplx4_product(deg,
            cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
            cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
            inputrehihi[ix1],inputrelohi[ix1],
            inputrehilo[ix1],inputrelolo[ix1],
            inputimhihi[ix1],inputimlohi[ix1],
            inputimhilo[ix1],inputimlolo[ix1],
            forwardrehihi[0],forwardrelohi[0],
            forwardrehilo[0],forwardrelolo[0],
            forwardimhihi[0],forwardimlohi[0],
            forwardimhilo[0],forwardimlolo[0]);

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
         {
            qdf_inc
               (&outputrehihi[ix2][j],&outputrelohi[ix2][j],
                &outputrehilo[ix2][j],&outputrelolo[ix2][j],
                forwardrehihi[0][j],forwardrelohi[0][j],
                forwardrehilo[0][j],forwardrelolo[0][j]);
            qdf_inc
               (&outputimhihi[ix2][j],&outputimlohi[ix2][j],
                &outputimhilo[ix2][j],&outputimlolo[ix2][j],
                forwardimhihi[0][j],forwardimlohi[0][j],
                forwardimhilo[0][j],forwardimlolo[0][j]);
         }
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_cmplx4_product(deg,
            cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
            cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
            inputrehihi[ix2],inputrelohi[ix2],
            inputrehilo[ix2],inputrelolo[ix2],
            inputimhihi[ix2],inputimlohi[ix2],
            inputimhilo[ix2],inputimlolo[ix2],
            backwardrehihi[0],backwardrelohi[0],
            backwardrehilo[0],backwardrelolo[0],
            backwardimhihi[0],backwardimlohi[0],
            backwardimhilo[0],backwardimlolo[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
         {
            qdf_inc
               (&outputrehihi[ix1][j],&outputrelohi[ix1][j],
                &outputrehilo[ix1][j],&outputrelolo[ix1][j],
                backwardrehihi[0][j],backwardrelohi[0][j],
                backwardrehilo[0][j],backwardrelolo[0][j]);
            qdf_inc
               (&outputimhihi[ix1][j],&outputimlohi[ix1][j],
                &outputimhilo[ix1][j],&outputimlolo[ix1][j],
                backwardimhihi[0][j],backwardimlohi[0][j],
                backwardimhilo[0][j],backwardimlolo[0][j]);
         }
         CPU_cmplx4_product(deg,
            forwardrehihi[0],forwardrelohi[0],
            forwardrehilo[0],forwardrelolo[0],
            forwardimhihi[0],forwardimlohi[0],
            forwardimhilo[0],forwardimlolo[0],
            inputrehihi[ix2],inputrelohi[ix2],
            inputrehilo[ix2],inputrelolo[ix2],
            inputimhihi[ix2],inputimlohi[ix2],
            inputimhilo[ix2],inputimlolo[ix2],
            forwardrehihi[1],forwardrelohi[1],
            forwardrehilo[1],forwardrelolo[1],
            forwardimhihi[1],forwardimlohi[1],
            forwardimhilo[1],forwardimlolo[1]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
         {
            qdf_inc
               (&outputrehihi[dim][j],&outputrelohi[dim][j],
                &outputrehilo[dim][j],&outputrelolo[dim][j],
                forwardrehihi[1][j],forwardrelohi[1][j],
                forwardrehilo[1][j],forwardrelolo[1][j]);
            qdf_inc
               (&outputimhihi[dim][j],&outputimlohi[dim][j],
                &outputimhilo[dim][j],&outputimlolo[dim][j],
                forwardimhihi[1][j],forwardimlohi[1][j],
                forwardimhilo[1][j],forwardimlolo[1][j]);
         }
      }
      else if(nvr[i] > 2)
      {
         CPU_cmplx4_speel
            (nvr[i],deg,idx[i],
             cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
             cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
             inputrehihi,inputrelohi,inputrehilo,inputrelolo,
             inputimhihi,inputimlohi,inputimhilo,inputimlolo,
             forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
             forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
             backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
             backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,
             crossrehihi,crossrelohi,crossrehilo,crossrelolo,
             crossimhihi,crossimlohi,crossimhilo,crossimlolo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
         {
            qdf_inc
               (&outputrehihi[dim][j],&outputrelohi[dim][j],
                &outputrehilo[dim][j],&outputrelolo[dim][j],
                forwardrehihi[ix1][j],forwardrelohi[ix1][j],
                forwardrehilo[ix1][j],forwardrelolo[ix1][j]);
            qdf_inc
               (&outputimhihi[dim][j],&outputimlohi[dim][j],
                &outputimhilo[dim][j],&outputimlolo[dim][j],
                forwardimhihi[ix1][j],forwardimlohi[ix1][j],
                forwardimhilo[ix1][j],forwardimlolo[ix1][j]);
         }
         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
         {
            qdf_inc
               (&outputrehihi[ix2][j],&outputrelohi[ix2][j],
                &outputrehilo[ix2][j],&outputrelolo[ix2][j],
                forwardrehihi[ix1][j],forwardrelohi[ix1][j],
                forwardrehilo[ix1][j],forwardrelolo[ix1][j]);
            qdf_inc
               (&outputimhihi[ix2][j],&outputimlohi[ix2][j],
                &outputimhilo[ix2][j],&outputimlolo[ix2][j],
                forwardimhihi[ix1][j],forwardimlohi[ix1][j],
                forwardimhilo[ix1][j],forwardimlolo[ix1][j]);
         }
         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
         {
            qdf_inc
               (&outputrehihi[ix2][j],&outputrelohi[ix2][j],
                &outputrehilo[ix2][j],&outputrelolo[ix2][j],
                backwardrehihi[ix1][j],backwardrelohi[ix1][j],
                backwardrehilo[ix1][j],backwardrelolo[ix1][j]);
            qdf_inc
               (&outputimhihi[ix2][j],&outputimlohi[ix2][j],
                &outputimhilo[ix2][j],&outputimlolo[ix2][j],
                backwardimhihi[ix1][j],backwardimlohi[ix1][j],
                backwardimhilo[ix1][j],backwardimlolo[ix1][j]);
         }
         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
            {
               qdf_inc
                  (&outputrehihi[ix2][j],&outputrelohi[ix2][j],
                   &outputrehilo[ix2][j],&outputrelolo[ix2][j],
                   crossrehihi[k-1][j],crossrelohi[k-1][j],
                   crossrehilo[k-1][j],crossrelolo[k-1][j]);
               qdf_inc
                  (&outputimhihi[ix2][j],&outputimlohi[ix2][j],
                   &outputimhilo[ix2][j],&outputimlolo[ix2][j],
                   crossimhihi[k-1][j],crossimlohi[k-1][j],
                   crossimhilo[k-1][j],crossimlolo[k-1][j]);
            }
         }
      }
   }
}

void CPU_dbl4_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double **forwardhihi = new double*[dim];
   double **forwardlohi = new double*[dim];
   double **forwardhilo = new double*[dim];
   double **forwardlolo = new double*[dim];
   double **backwardhihi = new double*[dim-1]; // in case dim = 2
   double **backwardlohi = new double*[dim-1];
   double **backwardhilo = new double*[dim-1];
   double **backwardlolo = new double*[dim-1];
   double **crosshihi = new double*[dim-1];    // in case dim = 2
   double **crosslohi = new double*[dim-1];
   double **crosshilo = new double*[dim-1];
   double **crosslolo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardhihi[i] = new double[deg+1];
      forwardlohi[i] = new double[deg+1];
      forwardhilo[i] = new double[deg+1];
      forwardlolo[i] = new double[deg+1];
      backwardhihi[i] = new double[deg+1];
      backwardlohi[i] = new double[deg+1];
      backwardhilo[i] = new double[deg+1];
      backwardlolo[i] = new double[deg+1];
      crosshihi[i] = new double[deg+1];
      crosslohi[i] = new double[deg+1];
      crosshilo[i] = new double[deg+1];
      crosslolo[i] = new double[deg+1];
   }
   forwardhihi[dim-1] = new double[deg+1];
   forwardlohi[dim-1] = new double[deg+1];
   forwardhilo[dim-1] = new double[deg+1];
   forwardlolo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputhihi[dim][i] = csthihi[i];
      outputlohi[dim][i] = cstlohi[i];
      outputhilo[dim][i] = csthilo[i];
      outputlolo[dim][i] = cstlolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhihi[i][j] = 0.0;
         outputlohi[i][j] = 0.0;
         outputhilo[i][j] = 0.0;
         outputlolo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl4_poly_speel
      (dim,nbr,deg,nvr,idx,
            cffhihi,     cfflohi,     cffhilo,     cfflolo,
          inputhihi,   inputlohi,   inputhilo,   inputlolo,
         outputhihi,  outputlohi,  outputhilo,  outputlolo,
        forwardhihi, forwardlohi, forwardhilo, forwardlolo,
       backwardhihi,backwardlohi,backwardhilo,backwardlolo,
          crosshihi,   crosslohi,   crosshilo,   crosslolo,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
   }
   for(int i=0; i<dim-1; i++)
   {
      free(forwardhihi[i]); free(backwardhihi[i]); free(crosshihi[i]);
      free(forwardlohi[i]); free(backwardlohi[i]); free(crosslohi[i]);
      free(forwardhilo[i]); free(backwardhilo[i]); free(crosshilo[i]);
      free(forwardlolo[i]); free(backwardlolo[i]); free(crosslolo[i]);
   }
   free(forwardhihi[dim-1]); free(forwardlohi[dim-1]);
   free(forwardhilo[dim-1]); free(forwardlolo[dim-1]);
   free(forwardhihi); free(backwardhihi); free(crosshihi);
   free(forwardlohi); free(backwardlohi); free(crosslohi);
   free(forwardhilo); free(backwardhilo); free(crosshilo);
   free(forwardlolo); free(backwardlolo); free(crosslolo);
}

void CPU_cmplx4_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo, 
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo, 
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double **forwardrehihi = new double*[dim];
   double **forwardrelohi = new double*[dim];
   double **forwardrehilo = new double*[dim];
   double **forwardrelolo = new double*[dim];
   double **forwardimhihi = new double*[dim];
   double **forwardimlohi = new double*[dim];
   double **forwardimhilo = new double*[dim];
   double **forwardimlolo = new double*[dim];
   double **backwardrehihi = new double*[dim-1]; // in case dim = 2
   double **backwardrelohi = new double*[dim-1];
   double **backwardrehilo = new double*[dim-1];
   double **backwardrelolo = new double*[dim-1];
   double **backwardimhihi = new double*[dim-1];
   double **backwardimlohi = new double*[dim-1];
   double **backwardimhilo = new double*[dim-1];
   double **backwardimlolo = new double*[dim-1];
   double **crossrehihi = new double*[dim-1];    // in case dim = 2
   double **crossrelohi = new double*[dim-1];
   double **crossrehilo = new double*[dim-1];
   double **crossrelolo = new double*[dim-1];
   double **crossimhihi = new double*[dim-1];
   double **crossimlohi = new double*[dim-1];
   double **crossimhilo = new double*[dim-1];
   double **crossimlolo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardrehihi[i] = new double[deg+1];
      forwardrelohi[i] = new double[deg+1];
      forwardrehilo[i] = new double[deg+1];
      forwardrelolo[i] = new double[deg+1];
      forwardimhihi[i] = new double[deg+1];
      forwardimlohi[i] = new double[deg+1];
      forwardimhilo[i] = new double[deg+1];
      forwardimlolo[i] = new double[deg+1];
      backwardrehihi[i] = new double[deg+1];
      backwardrelohi[i] = new double[deg+1];
      backwardrehilo[i] = new double[deg+1];
      backwardrelolo[i] = new double[deg+1];
      backwardimhihi[i] = new double[deg+1];
      backwardimlohi[i] = new double[deg+1];
      backwardimhilo[i] = new double[deg+1];
      backwardimlolo[i] = new double[deg+1];
      crossrehihi[i] = new double[deg+1];
      crossrelohi[i] = new double[deg+1];
      crossrehilo[i] = new double[deg+1];
      crossrelolo[i] = new double[deg+1];
      crossimhihi[i] = new double[deg+1];
      crossimlohi[i] = new double[deg+1];
      crossimhilo[i] = new double[deg+1];
      crossimlolo[i] = new double[deg+1];
   }
   forwardrehihi[dim-1] = new double[deg+1];
   forwardrelohi[dim-1] = new double[deg+1];
   forwardrehilo[dim-1] = new double[deg+1];
   forwardrelolo[dim-1] = new double[deg+1];
   forwardimhihi[dim-1] = new double[deg+1];
   forwardimlohi[dim-1] = new double[deg+1];
   forwardimhilo[dim-1] = new double[deg+1];
   forwardimlolo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputrehihi[dim][i] = cstrehihi[i];
      outputrelohi[dim][i] = cstrelohi[i];
      outputrehilo[dim][i] = cstrehilo[i];
      outputrelolo[dim][i] = cstrelolo[i];
      outputimhihi[dim][i] = cstimhihi[i];
      outputimlohi[dim][i] = cstimlohi[i];
      outputimhilo[dim][i] = cstimhilo[i];
      outputimlolo[dim][i] = cstimlolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrehihi[i][j] = 0.0;
         outputrelohi[i][j] = 0.0;
         outputrehilo[i][j] = 0.0;
         outputrelolo[i][j] = 0.0;
         outputimhihi[i][j] = 0.0;
         outputimlohi[i][j] = 0.0;
         outputimhilo[i][j] = 0.0;
         outputimlolo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_cmplx4_poly_speel
      (dim,nbr,deg,nvr,idx,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       inputrehihi,inputrelohi,inputrehilo,inputrelolo,
       inputimhihi,inputimlohi,inputimhilo,inputimlolo,
       outputrehihi,outputrelohi,outputrehilo,outputrelolo,
       outputimhihi,outputimlohi,outputimhilo,outputimlolo,
       forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
       forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
       backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
       backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,
       crossrehihi,crossrelohi,crossrehilo,crossrelolo,
       crossimhihi,crossimlohi,crossimhilo,crossimlolo,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
   }
   for(int i=0; i<dim-1; i++)
   {
      free(forwardrehihi[i]); free(backwardrehihi[i]); free(crossrehihi[i]);
      free(forwardrelohi[i]); free(backwardrelohi[i]); free(crossrelohi[i]);
      free(forwardrehilo[i]); free(backwardrehilo[i]); free(crossrehilo[i]);
      free(forwardrelolo[i]); free(backwardrelolo[i]); free(crossrelolo[i]);
      free(forwardimhihi[i]); free(backwardimhihi[i]); free(crossimhihi[i]);
      free(forwardimlohi[i]); free(backwardimlohi[i]); free(crossimlohi[i]);
      free(forwardimhilo[i]); free(backwardimhilo[i]); free(crossimhilo[i]);
      free(forwardimlolo[i]); free(backwardimlolo[i]); free(crossimlolo[i]);
   }
   free(forwardrehihi[dim-1]); free(forwardrelohi[dim-1]);
   free(forwardrehilo[dim-1]); free(forwardrelolo[dim-1]);
   free(forwardimhihi[dim-1]); free(forwardimlohi[dim-1]);
   free(forwardimhilo[dim-1]); free(forwardimlolo[dim-1]);
   free(forwardrehihi); free(backwardrehihi); free(crossrehihi);
   free(forwardrelohi); free(backwardrelohi); free(crossrelohi);
   free(forwardrehilo); free(backwardrehilo); free(crossrehilo);
   free(forwardrelolo); free(backwardrelolo); free(crossrelolo);
   free(forwardimhihi); free(backwardimhihi); free(crossimhihi);
   free(forwardimlohi); free(backwardimlohi); free(crossimlohi);
   free(forwardimhilo); free(backwardimhilo); free(crossimhilo);
   free(forwardimlolo); free(backwardimlolo); free(crossimlolo);
}

void CPU_dbl4_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo,
   ConvolutionJob job, bool verbose )
{
   const int inp1tp = job.get_first_type();
   const int inp1ix = job.get_first_input();
   const int inp2tp = job.get_second_type();
   const int inp2ix = job.get_second_input();
   const int outptp = job.get_output_type();
   const int outidx = job.get_output_index();

   if(outptp == 1) // forward product either initializes or accumulates
   {
      if(verbose) cout << "-> computing f[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl4_product(deg,
                cffhihi,cfflohi,cffhilo,cfflolo,
              inputhihi[inp2ix],  inputlohi[inp2ix],
              inputhilo[inp2ix],  inputlolo[inp2ix],
            forwardhihi[outidx],forwardlohi[outidx],
            forwardhilo[outidx],forwardlolo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl4_product(deg,
               inputhihi[inp1ix],inputlohi[inp1ix],
               inputhilo[inp1ix],inputlolo[inp1ix],
               cffhihi,cfflohi,cffhilo,cfflolo,
               forwardhihi[outidx],forwardlohi[outidx],
               forwardhilo[outidx],forwardlolo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl4_product(deg,
                 inputhihi[inp1ix],  inputlohi[inp1ix],
                 inputhilo[inp1ix],  inputlolo[inp1ix],
               forwardhihi[inp2ix],forwardlohi[inp2ix],
               forwardhilo[inp2ix],forwardlolo[inp2ix],
               forwardhihi[outidx],forwardlohi[outidx],
               forwardhilo[outidx],forwardlolo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl4_product(deg,
              crosshihi[inp1ix],  crosslohi[inp1ix],
              crosshilo[inp1ix],  crosslolo[inp1ix],
              inputhihi[inp2ix],  inputlohi[inp2ix],
              inputhilo[inp2ix],  inputlolo[inp2ix],
            forwardhihi[outidx],forwardlohi[outidx],
            forwardhilo[outidx],forwardlolo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl4_product(deg,
               inputhihi[inp1ix],inputlohi[inp1ix],
               inputhilo[inp1ix],inputlolo[inp1ix],
               cffhihi,cfflohi,cffhilo,cfflolo,
               forwardhihi[outidx],forwardlohi[outidx],
               forwardhilo[outidx],forwardlolo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl4_product(deg,
               forwardhihi[inp1ix],forwardlohi[inp1ix],
               forwardhilo[inp1ix],forwardlolo[inp1ix],
                 inputhihi[inp2ix],  inputlohi[inp2ix],
                 inputhilo[inp2ix],  inputlolo[inp2ix],
               forwardhihi[outidx],forwardlohi[outidx],
               forwardhilo[outidx],forwardlolo[outidx]);
         }
      }
   }
   else if(outptp == 2) // backward product either initializes or accumulates
   {
      if(verbose) cout << "-> computing b[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
            CPU_dbl4_product(deg,
                    cffhihi,cfflohi,cffhilo,cfflolo,
                  inputhihi[inp2ix],   inputlohi[inp2ix],
                  inputhilo[inp2ix],   inputlolo[inp2ix],
               backwardhihi[outidx],backwardlohi[outidx],
               backwardhilo[outidx],backwardlolo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl4_product(deg,
               cffhihi,cfflohi,cffhilo,cfflolo,
               backwardhihi[inp2ix],backwardlohi[inp2ix],
               backwardhilo[inp2ix],backwardlolo[inp2ix],
               backwardhihi[outidx],backwardlohi[outidx],
               backwardhilo[outidx],backwardlolo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl4_product(deg,
                  inputhihi[inp1ix],   inputlohi[inp1ix],
                  inputhilo[inp1ix],   inputlolo[inp1ix],
                  inputhihi[inp2ix],   inputlohi[inp2ix],
                  inputhilo[inp2ix],   inputlolo[inp2ix],
               backwardhihi[outidx],backwardlohi[outidx],
               backwardhilo[outidx],backwardlolo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl4_product(deg,
                  inputhihi[inp1ix],   inputlohi[inp1ix],
                  inputhilo[inp1ix],   inputlolo[inp1ix],
               backwardhihi[inp2ix],backwardlohi[inp2ix],
               backwardhilo[inp2ix],backwardlolo[inp2ix],
               backwardhihi[outidx],backwardlohi[outidx],
               backwardhilo[outidx],backwardlolo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl4_product(deg,
               backwardhihi[inp1ix],backwardlohi[inp1ix],
               backwardhilo[inp1ix],backwardlolo[inp1ix],
               cffhihi,cfflohi,cffhilo,cfflolo,
               backwardhihi[outidx],backwardlohi[outidx],
               backwardhilo[outidx],backwardlolo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl4_product(deg,
               backwardhihi[inp1ix],backwardlohi[inp1ix],
               backwardhilo[inp1ix],backwardlolo[inp1ix],
                  inputhihi[inp2ix],   inputlohi[inp2ix],
                  inputhilo[inp2ix],   inputlolo[inp2ix],
               backwardhihi[outidx],backwardlohi[outidx],
               backwardhilo[outidx],backwardlolo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl4_product(deg,
            cffhihi,cfflohi,cffhilo,cfflolo,
            inputhihi[inp2ix],inputlohi[inp2ix],
            inputhilo[inp2ix],inputlolo[inp2ix],
            crosshihi[outidx],crosslohi[outidx],
            crosshilo[outidx],crosslolo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl4_product(deg,
              inputhihi[inp1ix],  inputlohi[inp1ix],
              inputhilo[inp1ix],  inputlolo[inp1ix],
            forwardhihi[inp2ix],forwardlohi[inp2ix],
            forwardhilo[inp2ix],forwardlolo[inp2ix],
              crosshihi[outidx],  crosslohi[outidx],
              crosshilo[outidx],  crosslolo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl4_product(deg,
              forwardhihi[inp1ix],forwardlohi[inp1ix],
              forwardhilo[inp1ix],forwardlolo[inp1ix],
                inputhihi[inp2ix],  inputlohi[inp2ix],
                inputhilo[inp2ix],  inputlolo[inp2ix],
                crosshihi[outidx],  crosslohi[outidx],
                crosshilo[outidx],  crosslolo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl4_product(deg,
               forwardhihi[inp1ix], forwardlohi[inp1ix],
               forwardhilo[inp1ix], forwardlolo[inp1ix],
              backwardhihi[inp2ix],backwardlohi[inp2ix],
              backwardhilo[inp2ix],backwardlolo[inp2ix],
                 crosshihi[outidx],   crosslohi[outidx],
                 crosshilo[outidx],   crosslolo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl4_product(deg,
            backwardhihi[inp1ix],backwardlohi[inp1ix],
            backwardhilo[inp1ix],backwardlolo[inp1ix],
             forwardhihi[inp2ix], forwardlohi[inp2ix],
             forwardhilo[inp2ix], forwardlolo[inp2ix],
               crosshihi[outidx],   crosslohi[outidx],
               crosshilo[outidx],   crosslolo[outidx]);
      }
   }
}

void CPU_cmplx4_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrehihi, double *cffrelohi,
   double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi,
   double *cffimhilo, double *cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **backwardrehihi, double **backwardrelohi,
   double **backwardrehilo, double **backwardrelolo,
   double **backwardimhihi, double **backwardimlohi,
   double **backwardimhilo, double **backwardimlolo,
   double **crossrehihi, double **crossrelohi,
   double **crossrehilo, double **crossrelolo,
   double **crossimhihi, double **crossimlohi,
   double **crossimhilo, double **crossimlolo,
   ConvolutionJob job, bool verbose )
{
   const int inp1tp = job.get_first_type();
   const int inp1ix = job.get_first_input();
   const int inp2tp = job.get_second_type();
   const int inp2ix = job.get_second_input();
   const int outptp = job.get_output_type();
   const int outidx = job.get_output_index();

   if(outptp == 1) // forward product either initializes or accumulates
   {
      if(verbose) cout << "-> computing f[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_cmplx4_product(deg,
            cffrehihi,cffrelohi,cffrehilo,cffrelolo,
            cffimhihi,cffimlohi,cffimhilo,cffimlolo,
            inputrehihi[inp2ix],inputrelohi[inp2ix],
            inputrehilo[inp2ix],inputrelolo[inp2ix],
            inputimhihi[inp2ix],inputimlohi[inp2ix],
            inputimhilo[inp2ix],inputimlolo[inp2ix],
            forwardrehihi[outidx],forwardrelohi[outidx],
            forwardrehilo[outidx],forwardrelolo[outidx],
            forwardimhihi[outidx],forwardimlohi[outidx],
            forwardimhilo[outidx],forwardimlolo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx4_product(deg,
               inputrehihi[inp1ix],inputrelohi[inp1ix],
               inputrehilo[inp1ix],inputrelolo[inp1ix],
               inputimhihi[inp1ix],inputimlohi[inp1ix],
               inputimhilo[inp1ix],inputimlolo[inp1ix],
               cffrehihi,cffrelohi,cffrehilo,cffrelolo,
               cffimhihi,cffimlohi,cffimhilo,cffimlolo,
               forwardrehihi[outidx],forwardrelohi[outidx],
               forwardrehilo[outidx],forwardrelolo[outidx],
               forwardimhihi[outidx],forwardimlohi[outidx],
               forwardimhilo[outidx],forwardimlolo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_cmplx4_product(deg,
               inputrehihi[inp1ix],inputrelohi[inp1ix],
               inputrehilo[inp1ix],inputrelolo[inp1ix],
               inputimhihi[inp1ix],inputimlohi[inp1ix],
               inputimhilo[inp1ix],inputimlolo[inp1ix],
               forwardrehihi[inp2ix],forwardrelohi[inp2ix],
               forwardrehilo[inp2ix],forwardrelolo[inp2ix],
               forwardimhihi[inp2ix],forwardimlohi[inp2ix],
               forwardimhilo[inp2ix],forwardimlolo[inp2ix],
               forwardrehihi[outidx],forwardrelohi[outidx],
               forwardrehilo[outidx],forwardrelolo[outidx],
               forwardimhihi[outidx],forwardimlohi[outidx],
               forwardimhilo[outidx],forwardimlolo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_cmplx4_product(deg,
            crossrehihi[inp1ix],crossrelohi[inp1ix],
            crossrehilo[inp1ix],crossrelolo[inp1ix],
            crossimhihi[inp1ix],crossimlohi[inp1ix],
            crossimhilo[inp1ix],crossimlolo[inp1ix],
            inputrehihi[inp2ix],inputrelohi[inp2ix],
            inputrehilo[inp2ix],inputrelolo[inp2ix],
            inputimhihi[inp2ix],inputimlohi[inp2ix],
            inputimhilo[inp2ix],inputimlolo[inp2ix],
            forwardrehihi[outidx],forwardrelohi[outidx],
            forwardrehilo[outidx],forwardrelolo[outidx],
            forwardimhihi[outidx],forwardimlohi[outidx],
            forwardimhilo[outidx],forwardimlolo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx4_product(deg,
               inputrehihi[inp1ix],inputrelohi[inp1ix],
               inputrehilo[inp1ix],inputrelolo[inp1ix],
               inputimhihi[inp1ix],inputimlohi[inp1ix],
               inputimhilo[inp1ix],inputimlolo[inp1ix],
               cffrehihi,cffrelohi,cffrehilo,cffrelolo,
               cffimhihi,cffimlohi,cffimhilo,cffimlolo,
               forwardrehihi[outidx],forwardrelohi[outidx],
               forwardrehilo[outidx],forwardrelolo[outidx],
               forwardimhihi[outidx],forwardimlohi[outidx],
               forwardimhilo[outidx],forwardimlolo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx4_product(deg,
               forwardrehihi[inp1ix],forwardrelohi[inp1ix],
               forwardrehilo[inp1ix],forwardrelolo[inp1ix],
               forwardimhihi[inp1ix],forwardimlohi[inp1ix],
               forwardimhilo[inp1ix],forwardimlolo[inp1ix],
               inputrehihi[inp2ix],inputrelohi[inp2ix],
               inputrehilo[inp2ix],inputrelolo[inp2ix],
               inputimhihi[inp2ix],inputimlohi[inp2ix],
               inputimhilo[inp2ix],inputimlolo[inp2ix],
               forwardrehihi[outidx],forwardrelohi[outidx],
               forwardrehilo[outidx],forwardrelolo[outidx],
               forwardimhihi[outidx],forwardimlohi[outidx],
               forwardimhilo[outidx],forwardimlolo[outidx]);
         }
      }
   }
   else if(outptp == 2) // backward product either initializes or accumulates
   {
      if(verbose) cout << "-> computing b[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
            CPU_cmplx4_product(deg,
               cffrehihi,cffrelohi,cffrehilo,cffrelolo,
               cffimhihi,cffimlohi,cffimhilo,cffimlolo,
               inputrehihi[inp2ix],inputrelohi[inp2ix],
               inputrehilo[inp2ix],inputrelolo[inp2ix],
               inputimhihi[inp2ix],inputimlohi[inp2ix],
               inputimhilo[inp2ix],inputimlolo[inp2ix],
               backwardrehihi[outidx],backwardrelohi[outidx],
               backwardrehilo[outidx],backwardrelolo[outidx],
               backwardimhihi[outidx],backwardimlohi[outidx],
               backwardimhilo[outidx],backwardimlolo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_cmplx4_product(deg,
               cffrehihi,cffrelohi,cffrehilo,cffrelolo,
               cffimhihi,cffimlohi,cffimhilo,cffimlolo,
               backwardrehihi[inp2ix],backwardrelohi[inp2ix],
               backwardrehilo[inp2ix],backwardrelolo[inp2ix],
               backwardimhihi[inp2ix],backwardimlohi[inp2ix],
               backwardimhilo[inp2ix],backwardimlolo[inp2ix],
               backwardrehihi[outidx],backwardrelohi[outidx],
               backwardrehilo[outidx],backwardrelolo[outidx],
               backwardimhihi[outidx],backwardimlohi[outidx],
               backwardimhilo[outidx],backwardimlolo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_cmplx4_product(deg,
               inputrehihi[inp1ix],inputrelohi[inp1ix],
               inputrehilo[inp1ix],inputrelolo[inp1ix],
               inputimhihi[inp1ix],inputimlohi[inp1ix],
               inputimhilo[inp1ix],inputimlolo[inp1ix],
               inputrehihi[inp2ix],inputrelohi[inp2ix],
               inputrehilo[inp2ix],inputrelolo[inp2ix],
               inputimhihi[inp2ix],inputimlohi[inp2ix],
               inputimhilo[inp2ix],inputimlolo[inp2ix],
               backwardrehihi[outidx],backwardrelohi[outidx],
               backwardrehilo[outidx],backwardrelolo[outidx],
               backwardimhihi[outidx],backwardimlohi[outidx],
               backwardimhilo[outidx],backwardimlolo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_cmplx4_product(deg,
               inputrehihi[inp1ix],inputrelohi[inp1ix],
               inputrehilo[inp1ix],inputrelolo[inp1ix],
               inputimhihi[inp1ix],inputimlohi[inp1ix],
               inputimhilo[inp1ix],inputimlolo[inp1ix],
               backwardrehihi[inp2ix],backwardrelohi[inp2ix],
               backwardrehilo[inp2ix],backwardrelolo[inp2ix],
               backwardimhihi[inp2ix],backwardimlohi[inp2ix],
               backwardimhilo[inp2ix],backwardimlolo[inp2ix],
               backwardrehihi[outidx],backwardrelohi[outidx],
               backwardrehilo[outidx],backwardrelolo[outidx],
               backwardimhihi[outidx],backwardimlohi[outidx],
               backwardimhilo[outidx],backwardimlolo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_cmplx4_product(deg,
               backwardrehihi[inp1ix],backwardrelohi[inp1ix],
               backwardrehilo[inp1ix],backwardrelolo[inp1ix],
               backwardimhihi[inp1ix],backwardimlohi[inp1ix],
               backwardimhilo[inp1ix],backwardimlolo[inp1ix],
               cffrehihi,cffrelohi,cffrehilo,cffrelolo,
               cffimhihi,cffimlohi,cffimhilo,cffimlolo,
               backwardrehihi[outidx],backwardrelohi[outidx],
               backwardrehilo[outidx],backwardrelolo[outidx],
               backwardimhihi[outidx],backwardimlohi[outidx],
               backwardimhilo[outidx],backwardimlolo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx4_product(deg,
               backwardrehihi[inp1ix],backwardrelohi[inp1ix],
               backwardrehilo[inp1ix],backwardrelolo[inp1ix],
               backwardimhihi[inp1ix],backwardimlohi[inp1ix],
               backwardimhilo[inp1ix],backwardimlolo[inp1ix],
               inputrehihi[inp2ix],inputrelohi[inp2ix],
               inputrehilo[inp2ix],inputrelolo[inp2ix],
               inputimhihi[inp2ix],inputimlohi[inp2ix],
               inputimhilo[inp2ix],inputimlolo[inp2ix],
               backwardrehihi[outidx],backwardrelohi[outidx],
               backwardrehilo[outidx],backwardrelolo[outidx],
               backwardimhihi[outidx],backwardimlohi[outidx],
               backwardimhilo[outidx],backwardimlolo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_cmplx4_product(deg,
            cffrehihi,cffrelohi,cffrehilo,cffrelolo,
            cffimhihi,cffimlohi,cffimhilo,cffimlolo,
            inputrehihi[inp2ix],inputrelohi[inp2ix],
            inputrehilo[inp2ix],inputrelolo[inp2ix],
            inputimhihi[inp2ix],inputimlohi[inp2ix],
            inputimhilo[inp2ix],inputimlolo[inp2ix],
            crossrehihi[outidx],crossrelohi[outidx],
            crossrehilo[outidx],crossrelolo[outidx],
            crossimhihi[outidx],crossimlohi[outidx],
            crossimhilo[outidx],crossimlolo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx4_product(deg,
            inputrehihi[inp1ix],inputrelohi[inp1ix],
            inputrehilo[inp1ix],inputrelolo[inp1ix],
            inputimhihi[inp1ix],inputimlohi[inp1ix],
            inputimhilo[inp1ix],inputimlolo[inp1ix],
            forwardrehihi[inp2ix],forwardrelohi[inp2ix],
            forwardrehilo[inp2ix],forwardrelolo[inp2ix],
            forwardimhihi[inp2ix],forwardimlohi[inp2ix],
            forwardimhilo[inp2ix],forwardimlolo[inp2ix],
            crossrehihi[outidx],crossrelohi[outidx],
            crossrehilo[outidx],crossrelolo[outidx],
            crossimhihi[outidx],crossimlohi[outidx],
            crossimhilo[outidx],crossimlolo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_cmplx4_product(deg,
              forwardrehihi[inp1ix],forwardrelohi[inp1ix],
              forwardrehilo[inp1ix],forwardrelolo[inp1ix],
              forwardimhihi[inp1ix],forwardimlohi[inp1ix],
              forwardimhilo[inp1ix],forwardimlolo[inp1ix],
              inputrehihi[inp2ix],inputrelohi[inp2ix],
              inputrehilo[inp2ix],inputrelolo[inp2ix],
              inputimhihi[inp2ix],inputimlohi[inp2ix],
              inputimhilo[inp2ix],inputimlolo[inp2ix],
              crossrehihi[outidx],crossrelohi[outidx],
              crossrehilo[outidx],crossrelolo[outidx],
              crossimhihi[outidx],crossimlohi[outidx],
              crossimhilo[outidx],crossimlolo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_cmplx4_product(deg,
              forwardrehihi[inp1ix],forwardrelohi[inp1ix],
              forwardrehilo[inp1ix],forwardrelolo[inp1ix],
              forwardimhihi[inp1ix],forwardimlohi[inp1ix],
              forwardimhilo[inp1ix],forwardimlolo[inp1ix],
              backwardrehihi[inp2ix],backwardrelohi[inp2ix],
              backwardrehilo[inp2ix],backwardrelolo[inp2ix],
              backwardimhihi[inp2ix],backwardimlohi[inp2ix],
              backwardimhilo[inp2ix],backwardimlolo[inp2ix],
              crossrehihi[outidx],crossrelohi[outidx],
              crossrehilo[outidx],crossrelolo[outidx],
              crossimhihi[outidx],crossimlohi[outidx],
              crossimhilo[outidx],crossimlolo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx4_product(deg,
            backwardrehihi[inp1ix],backwardrelohi[inp1ix],
            backwardrehilo[inp1ix],backwardrelolo[inp1ix],
            backwardimhihi[inp1ix],backwardimlohi[inp1ix],
            backwardimhilo[inp1ix],backwardimlolo[inp1ix],
            forwardrehihi[inp2ix],forwardrelohi[inp2ix],
            forwardrehilo[inp2ix],forwardrelolo[inp2ix],
            forwardimhihi[inp2ix],forwardimlohi[inp2ix],
            forwardimhilo[inp2ix],forwardimlolo[inp2ix],
            crossrehihi[outidx],crossrelohi[outidx],
            crossrehilo[outidx],crossrelolo[outidx],
            crossimhihi[outidx],crossimlohi[outidx],
            crossimhilo[outidx],crossimlolo[outidx]);
      }
   }
}

void CPU_dbl4_add_job
 ( int deg,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo,
   AdditionJob job, bool verbose )
{
   const int adtype = job.get_addition_type();
   const int intype = job.get_increment_type();
   const int updmon = job.get_update_monomial();
   const int updidx = job.get_update_index();
   const int incmon = job.get_increment_monomial();
   const int incidx = job.get_increment_index();

   if(adtype == 1)
   {
      if(incmon < 0)
      {
         if(incidx < 0)
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cst[i];
               qdf_inc(&forwardhihi[updmon][updidx][i],
                       &forwardlohi[updmon][updidx][i],
                       &forwardhilo[updmon][updidx][i],
                       &forwardlolo[updmon][updidx][i],
                       csthihi[i],cstlohi[i],csthilo[i],cstlolo[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               qdf_inc(&forwardhihi[updmon][updidx][i],
                       &forwardlohi[updmon][updidx][i],
                       &forwardhilo[updmon][updidx][i],
                       &forwardlolo[updmon][updidx][i],
                       cffhihi[incidx][i],cfflohi[incidx][i],
                       cffhilo[incidx][i],cfflolo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            qdf_inc(&forwardhihi[updmon][updidx][i],
                    &forwardlohi[updmon][updidx][i],
                    &forwardhilo[updmon][updidx][i],
                    &forwardlolo[updmon][updidx][i],
                    forwardhihi[incmon][incidx][i],
                    forwardlohi[incmon][incidx][i],
                    forwardhilo[incmon][incidx][i],
                    forwardlolo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            qdf_inc(&forwardhihi[updmon][updidx][i],
                    &forwardlohi[updmon][updidx][i],
                    &forwardhilo[updmon][updidx][i],
                    &forwardlolo[updmon][updidx][i],
                    backwardhihi[incmon][incidx][i],
                    backwardlohi[incmon][incidx][i],
                    backwardhilo[incmon][incidx][i],
                    backwardlolo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            qdf_inc(&forwardhihi[updmon][updidx][i],
                    &forwardlohi[updmon][updidx][i],
                    &forwardhilo[updmon][updidx][i],
                    &forwardlolo[updmon][updidx][i],
                    crosshihi[incmon][incidx][i],
                    crosslohi[incmon][incidx][i],
                    crosshilo[incmon][incidx][i],
                    crosslolo[incmon][incidx][i]);
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
            qdf_inc(&backwardhihi[updmon][updidx][i],
                    &backwardlohi[updmon][updidx][i],
                    &backwardhilo[updmon][updidx][i],
                    &backwardlolo[updmon][updidx][i],
                    cffhihi[incidx][i],cfflohi[incidx][i],
                    cffhilo[incidx][i],cfflolo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            qdf_inc(&backwardhihi[updmon][updidx][i],
                    &backwardlohi[updmon][updidx][i],
                    &backwardhilo[updmon][updidx][i],
                    &backwardlolo[updmon][updidx][i],
                    forwardhihi[incmon][incidx][i],
                    forwardlohi[incmon][incidx][i],
                    forwardhilo[incmon][incidx][i],
                    forwardlolo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            qdf_inc(&backwardhihi[updmon][updidx][i],
                    &backwardlohi[updmon][updidx][i],
                    &backwardhilo[updmon][updidx][i],
                    &backwardlolo[updmon][updidx][i],
                    backwardhihi[incmon][incidx][i],
                    backwardlohi[incmon][incidx][i],
                    backwardhilo[incmon][incidx][i],
                    backwardlolo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            qdf_inc(&backwardhihi[updmon][updidx][i],
                    &backwardlohi[updmon][updidx][i],
                    &backwardhilo[updmon][updidx][i],
                    &backwardlolo[updmon][updidx][i],
                    crosshihi[incmon][incidx][i],
                    crosslohi[incmon][incidx][i],
                    crosshilo[incmon][incidx][i],
                    crosslolo[incmon][incidx][i]);
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
            qdf_inc(&crosshihi[updmon][updidx][i],
                    &crosslohi[updmon][updidx][i],
                    &crosshilo[updmon][updidx][i],
                    &crosslolo[updmon][updidx][i],
                    cffhihi[incidx][i],cfflohi[incidx][i],
                    cffhilo[incidx][i],cfflolo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            qdf_inc(&crosshihi[updmon][updidx][i],
                    &crosslohi[updmon][updidx][i],
                    &crosshilo[updmon][updidx][i],
                    &crosslolo[updmon][updidx][i],
                    forwardhihi[incmon][incidx][i],
                    forwardlohi[incmon][incidx][i],
                    forwardhilo[incmon][incidx][i],
                    forwardlolo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            qdf_inc(&crosshihi[updmon][updidx][i],
                    &crosslohi[updmon][updidx][i],
                    &crosshilo[updmon][updidx][i],
                    &crosslolo[updmon][updidx][i],
                    backwardhihi[incmon][incidx][i],
                    backwardlohi[incmon][incidx][i],
                    backwardhilo[incmon][incidx][i],
                    backwardlolo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            qdf_inc(&crosshihi[updmon][updidx][i],
                    &crosslohi[updmon][updidx][i],
                    &crosshilo[updmon][updidx][i],
                    &crosslolo[updmon][updidx][i],
                    crosshihi[incmon][incidx][i],
                    crosslohi[incmon][incidx][i],
                    crosshilo[incmon][incidx][i],
                    crosslolo[incmon][incidx][i]);
      }
   }
}

void CPU_cmplx4_add_job
 ( int deg, double *cstrehihi, double *cstrelohi,
            double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double ***forwardrehihi, double ***forwardrelohi,
   double ***forwardrehilo, double ***forwardrelolo,
   double ***forwardimhihi, double ***forwardimlohi,
   double ***forwardimhilo, double ***forwardimlolo,
   double ***backwardrehihi, double ***backwardrelohi,
   double ***backwardrehilo, double ***backwardrelolo, 
   double ***backwardimhihi, double ***backwardimlohi,
   double ***backwardimhilo, double ***backwardimlolo, 
   double ***crossrehihi, double ***crossrelohi,
   double ***crossrehilo, double ***crossrelolo,
   double ***crossimhihi, double ***crossimlohi,
   double ***crossimhilo, double ***crossimlolo,
   AdditionJob job, bool verbose )
{
   const int adtype = job.get_addition_type();
   const int intype = job.get_increment_type();
   const int updmon = job.get_update_monomial();
   const int updidx = job.get_update_index();
   const int incmon = job.get_increment_monomial();
   const int incidx = job.get_increment_index();

   if(adtype == 1)
   {
      if(incmon < 0)
      {
         if(incidx < 0)
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cst[i];
            {
               qdf_inc(&forwardrehihi[updmon][updidx][i],
                       &forwardrelohi[updmon][updidx][i],
                       &forwardrehilo[updmon][updidx][i],
                       &forwardrelolo[updmon][updidx][i],
                       cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i]);
               qdf_inc(&forwardimhihi[updmon][updidx][i],
                       &forwardimlohi[updmon][updidx][i],
                       &forwardimhilo[updmon][updidx][i],
                       &forwardimlolo[updmon][updidx][i],
                       cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i]);
            }
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
            {
               qdf_inc(&forwardrehihi[updmon][updidx][i],
                       &forwardrelohi[updmon][updidx][i],
                       &forwardrehilo[updmon][updidx][i],
                       &forwardrelolo[updmon][updidx][i],
                       cffrehihi[incidx][i],cffrelohi[incidx][i],
                       cffrehilo[incidx][i],cffrelolo[incidx][i]);
               qdf_inc(&forwardimhihi[updmon][updidx][i],
                       &forwardimlohi[updmon][updidx][i],
                       &forwardimhilo[updmon][updidx][i],
                       &forwardimlolo[updmon][updidx][i],
                       cffimhihi[incidx][i],cffimlohi[incidx][i],
                       cffimhilo[incidx][i],cffimlolo[incidx][i]);
            }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            qdf_inc(&forwardrehihi[updmon][updidx][i],
                    &forwardrelohi[updmon][updidx][i],
                    &forwardrehilo[updmon][updidx][i],
                    &forwardrelolo[updmon][updidx][i],
                    forwardrehihi[incmon][incidx][i],
                    forwardrelohi[incmon][incidx][i],
                    forwardrehilo[incmon][incidx][i],
                    forwardrelolo[incmon][incidx][i]);
            qdf_inc(&forwardimhihi[updmon][updidx][i],
                    &forwardimlohi[updmon][updidx][i],
                    &forwardimhilo[updmon][updidx][i],
                    &forwardimlolo[updmon][updidx][i],
                    forwardimhihi[incmon][incidx][i],
                    forwardimlohi[incmon][incidx][i],
                    forwardimhilo[incmon][incidx][i],
                    forwardimlolo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            qdf_inc(&forwardrehihi[updmon][updidx][i],
                    &forwardrelohi[updmon][updidx][i],
                    &forwardrehilo[updmon][updidx][i],
                    &forwardrelolo[updmon][updidx][i],
                    backwardrehihi[incmon][incidx][i],
                    backwardrelohi[incmon][incidx][i],
                    backwardrehilo[incmon][incidx][i],
                    backwardrelolo[incmon][incidx][i]);
            qdf_inc(&forwardimhihi[updmon][updidx][i],
                    &forwardimlohi[updmon][updidx][i],
                    &forwardimhilo[updmon][updidx][i],
                    &forwardimlolo[updmon][updidx][i],
                    backwardimhihi[incmon][incidx][i],
                    backwardimlohi[incmon][incidx][i],
                    backwardimhilo[incmon][incidx][i],
                    backwardimlolo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            qdf_inc(&forwardrehihi[updmon][updidx][i],
                    &forwardrelohi[updmon][updidx][i],
                    &forwardrehilo[updmon][updidx][i],
                    &forwardrelolo[updmon][updidx][i],
                    crossrehihi[incmon][incidx][i],
                    crossrelohi[incmon][incidx][i],
                    crossrehilo[incmon][incidx][i],
                    crossrelolo[incmon][incidx][i]);
            qdf_inc(&forwardimhihi[updmon][updidx][i],
                    &forwardimlohi[updmon][updidx][i],
                    &forwardimhilo[updmon][updidx][i],
                    &forwardimlolo[updmon][updidx][i],
                    crossimhihi[incmon][incidx][i],
                    crossimlohi[incmon][incidx][i],
                    crossimhilo[incmon][incidx][i],
                    crossimlolo[incmon][incidx][i]);
         }
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
         {
            qdf_inc(&backwardrehihi[updmon][updidx][i],
                    &backwardrelohi[updmon][updidx][i],
                    &backwardrehilo[updmon][updidx][i],
                    &backwardrelolo[updmon][updidx][i],
                    cffrehihi[incidx][i],cffrelohi[incidx][i],
                    cffrehilo[incidx][i],cffrelolo[incidx][i]);
            qdf_inc(&backwardimhihi[updmon][updidx][i],
                    &backwardimlohi[updmon][updidx][i],
                    &backwardimhilo[updmon][updidx][i],
                    &backwardimlolo[updmon][updidx][i],
                    cffimhihi[incidx][i],cffimlohi[incidx][i],
                    cffimhilo[incidx][i],cffimlolo[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            qdf_inc(&backwardrehihi[updmon][updidx][i],
                    &backwardrelohi[updmon][updidx][i],
                    &backwardrehilo[updmon][updidx][i],
                    &backwardrelolo[updmon][updidx][i],
                    forwardrehihi[incmon][incidx][i],
                    forwardrelohi[incmon][incidx][i],
                    forwardrehilo[incmon][incidx][i],
                    forwardrelolo[incmon][incidx][i]);
            qdf_inc(&backwardimhihi[updmon][updidx][i],
                    &backwardimlohi[updmon][updidx][i],
                    &backwardimhilo[updmon][updidx][i],
                    &backwardimlolo[updmon][updidx][i],
                    forwardimhihi[incmon][incidx][i],
                    forwardimlohi[incmon][incidx][i],
                    forwardimhilo[incmon][incidx][i],
                    forwardimlolo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            qdf_inc(&backwardrehihi[updmon][updidx][i],
                    &backwardrelohi[updmon][updidx][i],
                    &backwardrehilo[updmon][updidx][i],
                    &backwardrelolo[updmon][updidx][i],
                    backwardrehihi[incmon][incidx][i],
                    backwardrelohi[incmon][incidx][i],
                    backwardrehilo[incmon][incidx][i],
                    backwardrelolo[incmon][incidx][i]);
            qdf_inc(&backwardimhihi[updmon][updidx][i],
                    &backwardimlohi[updmon][updidx][i],
                    &backwardimhilo[updmon][updidx][i],
                    &backwardimlolo[updmon][updidx][i],
                    backwardimhihi[incmon][incidx][i],
                    backwardimlohi[incmon][incidx][i],
                    backwardimhilo[incmon][incidx][i],
                    backwardimlolo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            qdf_inc(&backwardrehihi[updmon][updidx][i],
                    &backwardrelohi[updmon][updidx][i],
                    &backwardrehilo[updmon][updidx][i],
                    &backwardrelolo[updmon][updidx][i],
                    crossrehihi[incmon][incidx][i],
                    crossrelohi[incmon][incidx][i],
                    crossrehilo[incmon][incidx][i],
                    crossrelolo[incmon][incidx][i]);
            qdf_inc(&backwardimhihi[updmon][updidx][i],
                    &backwardimlohi[updmon][updidx][i],
                    &backwardimhilo[updmon][updidx][i],
                    &backwardimlolo[updmon][updidx][i],
                    crossimhihi[incmon][incidx][i],
                    crossimlohi[incmon][incidx][i],
                    crossimhilo[incmon][incidx][i],
                    crossimlolo[incmon][incidx][i]);
         }
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
         {
            qdf_inc(&crossrehihi[updmon][updidx][i],
                    &crossrelohi[updmon][updidx][i],
                    &crossrehilo[updmon][updidx][i],
                    &crossrelolo[updmon][updidx][i],
                    cffrehihi[incidx][i],cffrelohi[incidx][i],
                    cffrehilo[incidx][i],cffrelolo[incidx][i]);
            qdf_inc(&crossimhihi[updmon][updidx][i],
                    &crossimlohi[updmon][updidx][i],
                    &crossimhilo[updmon][updidx][i],
                    &crossimlolo[updmon][updidx][i],
                    cffimhihi[incidx][i],cffimlohi[incidx][i],
                    cffimhilo[incidx][i],cffimlolo[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            qdf_inc(&crossrehihi[updmon][updidx][i],
                    &crossrelohi[updmon][updidx][i],
                    &crossrehilo[updmon][updidx][i],
                    &crossrelolo[updmon][updidx][i],
                    forwardrehihi[incmon][incidx][i],
                    forwardrelohi[incmon][incidx][i],
                    forwardrehilo[incmon][incidx][i],
                    forwardrelolo[incmon][incidx][i]);
            qdf_inc(&crossimhihi[updmon][updidx][i],
                    &crossimlohi[updmon][updidx][i],
                    &crossimhilo[updmon][updidx][i],
                    &crossimlolo[updmon][updidx][i],
                    forwardimhihi[incmon][incidx][i],
                    forwardimlohi[incmon][incidx][i],
                    forwardimhilo[incmon][incidx][i],
                    forwardimlolo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            qdf_inc(&crossrehihi[updmon][updidx][i],
                    &crossrelohi[updmon][updidx][i],
                    &crossrehilo[updmon][updidx][i],
                    &crossrelolo[updmon][updidx][i],
                    backwardrehihi[incmon][incidx][i],
                    backwardrelohi[incmon][incidx][i],
                    backwardrehilo[incmon][incidx][i],
                    backwardrelolo[incmon][incidx][i]);
            qdf_inc(&crossimhihi[updmon][updidx][i],
                    &crossimlohi[updmon][updidx][i],
                    &crossimhilo[updmon][updidx][i],
                    &crossimlolo[updmon][updidx][i],
                    backwardimhihi[incmon][incidx][i],
                    backwardimlohi[incmon][incidx][i],
                    backwardimhilo[incmon][incidx][i],
                    backwardimlolo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            qdf_inc(&crossrehihi[updmon][updidx][i],
                    &crossrelohi[updmon][updidx][i],
                    &crossrehilo[updmon][updidx][i],
                    &crossrelolo[updmon][updidx][i],
                    crossrehihi[incmon][incidx][i],
                    crossrelohi[incmon][incidx][i],
                    crossrehilo[incmon][incidx][i],
                    crossrelolo[incmon][incidx][i]);
            qdf_inc(&crossimhihi[updmon][updidx][i],
                    &crossimlohi[updmon][updidx][i],
                    &crossimhilo[updmon][updidx][i],
                    &crossimlolo[updmon][updidx][i],
                    crossimhihi[incmon][incidx][i],
                    crossimlohi[incmon][incidx][i],
                    crossimhilo[incmon][incidx][i],
                    crossimlolo[incmon][incidx][i]);
         }
      }
   }
}

void CPU_dbl4_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo )
{
   for(int i=0; i<=deg; i++)
   {
      outputhihi[dim][i] = csthihi[i];
      outputlohi[dim][i] = cstlohi[i];
      outputhilo[dim][i] = csthilo[i];
      outputlolo[dim][i] = cstlolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhihi[i][j] = 0.0;
         outputlohi[i][j] = 0.0;
         outputhilo[i][j] = 0.0;
         outputlolo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      int ix0 = idx[k][0];   // first variable in monomial k
      int ix1 = nvr[k]-1;    // last forward has the value
      int ix2 = nvr[k]-2;    // next to last forward has last derivative
                             // last backward has the first derivative
      int ixn = idx[k][ix1]; // index of the last variable in monomial k

      for(int i=0; i<=deg; i++) // value is last forward location
         // output[dim][i] = output[dim][i] + forward[k][ix1][i];
         qdf_inc(&outputhihi[dim][i],   &outputlohi[dim][i],
                 &outputhilo[dim][i],   &outputlolo[dim][i],
                 forwardhihi[k][ix1][i],forwardlohi[k][ix1][i],
                 forwardhilo[k][ix1][i],forwardlolo[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                        cffhihi[k][i],      cfflohi[k][i],
                        cffhilo[k][i],      cfflolo[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            qdf_inc(&outputhihi[ixn][i],   &outputlohi[ixn][i],
                    &outputhilo[ixn][i],   &outputlolo[ixn][i],
                    forwardhihi[k][ix2][i],forwardlohi[k][ix2][i],
                    forwardhilo[k][ix2][i],forwardlolo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            qdf_inc( &outputhihi[ix0][i],    &outputlohi[ix0][i],
                     &outputhilo[ix0][i],    &outputlolo[ix0][i],
                    backwardhihi[k][ix2][i],backwardlohi[k][ix2][i],
                    backwardhilo[k][ix2][i],backwardlolo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  qdf_inc(&outputhihi[ix0][i], &outputlohi[ix0][i],
                          &outputhilo[ix0][i], &outputlolo[ix0][i],
                            crosshihi[k][j-1][i],crosslohi[k][j-1][i],
                            crosshilo[k][j-1][i],crosslolo[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_cmplx4_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   double ***forwardrehihi, double ***forwardrelohi,
   double ***forwardrehilo, double ***forwardrelolo,
   double ***forwardimhihi, double ***forwardimlohi,
   double ***forwardimhilo, double ***forwardimlolo,
   double ***backwardrehihi, double ***backwardrelohi,
   double ***backwardrehilo, double ***backwardrelolo,
   double ***backwardimhihi, double ***backwardimlohi,
   double ***backwardimhilo, double ***backwardimlolo,
   double ***crossrehihi, double ***crossrelohi,
   double ***crossrehilo, double ***crossrelolo,
   double ***crossimhihi, double ***crossimlohi,
   double ***crossimhilo, double ***crossimlolo )
{
   for(int i=0; i<=deg; i++)
   {
      outputrehihi[dim][i] = cstrehihi[i];
      outputrelohi[dim][i] = cstrelohi[i];
      outputrehilo[dim][i] = cstrehilo[i];
      outputrelolo[dim][i] = cstrelolo[i];
      outputimhihi[dim][i] = cstimhihi[i];
      outputimlohi[dim][i] = cstimlohi[i];
      outputimhilo[dim][i] = cstimhilo[i];
      outputimlolo[dim][i] = cstimlolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrehihi[i][j] = 0.0; outputrelohi[i][j] = 0.0;
         outputrehilo[i][j] = 0.0; outputrelolo[i][j] = 0.0;
         outputimhihi[i][j] = 0.0; outputimlohi[i][j] = 0.0;
         outputimhilo[i][j] = 0.0; outputimlolo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      int ix0 = idx[k][0];   // first variable in monomial k
      int ix1 = nvr[k]-1;    // last forward has the value
      int ix2 = nvr[k]-2;    // next to last forward has last derivative
                             // last backward has the first derivative
      int ixn = idx[k][ix1]; // index of the last variable in monomial k

      for(int i=0; i<=deg; i++) // value is last forward location
         // output[dim][i] = output[dim][i] + forward[k][ix1][i];
      {
         qdf_inc(&outputrehihi[dim][i],&outputrelohi[dim][i],
                 &outputrehilo[dim][i],&outputrelolo[dim][i],
                 forwardrehihi[k][ix1][i],forwardrelohi[k][ix1][i],
                 forwardrehilo[k][ix1][i],forwardrelolo[k][ix1][i]);
         qdf_inc(&outputimhihi[dim][i],&outputimlohi[dim][i],
                 &outputimhilo[dim][i],&outputimlolo[dim][i],
                 forwardimhihi[k][ix1][i],forwardimlohi[k][ix1][i],
                 forwardimhilo[k][ix1][i],forwardimlolo[k][ix1][i]);
      }
      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
         {
            qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                    &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                    cffrehihi[k][i],cffrelohi[k][i],
                    cffrehilo[k][i],cffrelolo[k][i]);
            qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                    &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                    cffimhihi[k][i],cffimlohi[k][i],
                    cffimhilo[k][i],cffimlolo[k][i]);
         }
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            qdf_inc(&outputrehihi[ixn][i],&outputrelohi[ixn][i],
                    &outputrehilo[ixn][i],&outputrelolo[ixn][i],
                    forwardrehihi[k][ix2][i],forwardrelohi[k][ix2][i],
                    forwardrehilo[k][ix2][i],forwardrelolo[k][ix2][i]);
            qdf_inc(&outputimhihi[ixn][i],&outputimlohi[ixn][i],
                    &outputimhilo[ixn][i],&outputimlolo[ixn][i],
                    forwardimhihi[k][ix2][i],forwardimlohi[k][ix2][i],
                    forwardimhilo[k][ix2][i],forwardimlolo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                    &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                    backwardrehihi[k][ix2][i],backwardrelohi[k][ix2][i],
                    backwardrehilo[k][ix2][i],backwardrelolo[k][ix2][i]);
            qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                    &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                    backwardimhihi[k][ix2][i],backwardimlohi[k][ix2][i],
                    backwardimhilo[k][ix2][i],backwardimlolo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
               {
                  qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                          &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                          crossrehihi[k][j-1][i],crossrelohi[k][j-1][i],
                          crossrehilo[k][j-1][i],crossrelolo[k][j-1][i]);
                  qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                          &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                          crossimhihi[k][j-1][i],crossimlohi[k][j-1][i],
                          crossimhilo[k][j-1][i],crossimlolo[k][j-1][i]);
               }
            }
         }
      }
   }
}

void CPU_dbl4_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo,
   AdditionJobs jobs, bool verbose )
{
   for(int k=0; k<jobs.get_depth(); k++)
   {
      if(verbose) cout << "executing addition jobs at layer "
                       << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
      {
         AdditionJob job = jobs.get_job(k,i);
         if(verbose) cout << "job " << i << " : " << job << endl;

         CPU_dbl4_add_job(deg,
                 csthihi,     cstlohi,     csthilo,     cstlolo,
                 cffhihi,     cfflohi,     cffhilo,     cfflolo,
             forwardhihi, forwardlohi, forwardhilo, forwardlolo,
            backwardhihi,backwardlohi,backwardhilo,backwardlolo,
               crosshihi,   crosslohi,   crosshilo,   crosslolo,
            job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputhihi[dim][i] = forwardhihi[lastmon][lastidx][i];
      outputlohi[dim][i] = forwardlohi[lastmon][lastidx][i];
      outputhilo[dim][i] = forwardhilo[lastmon][lastidx][i];
      outputlolo[dim][i] = forwardlolo[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(verbose)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++)
         {
            outputhihi[0][i] = 0.0;
            outputlohi[0][i] = 0.0;
            outputhilo[0][i] = 0.0;
            outputlolo[0][i] = 0.0;
         }
      }
      else
      {
         if(verbose)
            cout << "updating derivative 0 with coefficient "
                 << difidx << endl;

         for(int i=0; i<=deg; i++)
         {
            outputhihi[0][i] = cffhihi[difidx][i];
            outputlohi[0][i] = cfflohi[difidx][i];
            outputhilo[0][i] = cffhilo[difidx][i];
            outputlolo[0][i] = cfflolo[difidx][i];
         }
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0] - 2;
      
      if(verbose)
         cout << "Updating derivative 0, ix0 = " << ix0
              << ", ix2 = " << ix2
              << " : b[" << ix0 << "," << ix2 << "]" << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = backward[ix0][ix2][i];
      {
         outputhihi[0][i] = backwardhihi[ix0][ix2][i];
         outputlohi[0][i] = backwardlohi[ix0][ix2][i];
         outputhilo[0][i] = backwardhilo[ix0][ix2][i];
         outputlolo[0][i] = backwardlolo[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(verbose)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) 
            {
               outputhihi[k][i] = 0.0;
               outputlohi[k][i] = 0.0;
               outputhilo[k][i] = 0.0;
               outputlolo[k][i] = 0.0;
            }
         }
         else
         {
            if(verbose)
               cout << "updating derivative " << k 
                    << " with coefficient " << difidx << endl;

            for(int i=0; i<=deg; i++)
            {
               outputhihi[k][i] = cffhihi[difidx][i];
               outputlohi[k][i] = cfflohi[difidx][i];
               outputhilo[k][i] = cffhilo[difidx][i];
               outputlolo[k][i] = cfflolo[difidx][i];
            }
         }
      }
      else
      {
         int ix0 = jobs.get_differential_index(k,cnt);

         if(idx[ix0][0] == k) // k is first variable of monomial
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : b[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = backward[ix0][ix2][i];
            {
               outputhihi[k][i] = backwardhihi[ix0][ix2][i];
               outputlohi[k][i] = backwardlohi[ix0][ix2][i];
               outputhilo[k][i] = backwardhilo[ix0][ix2][i];
               outputlolo[k][i] = backwardlolo[ix0][ix2][i];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : f[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = forward[ix0][ix2][i];
            {
               outputhihi[k][i] = forwardhihi[ix0][ix2][i];
               outputlohi[k][i] = forwardlohi[ix0][ix2][i];
               outputhilo[k][i] = forwardhilo[ix0][ix2][i];
               outputlolo[k][i] = forwardlolo[ix0][ix2][i];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : c[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = cross[ix0][ix2][i];
            {
               outputhihi[k][i] = crosshihi[ix0][ix2][i];
               outputlohi[k][i] = crosslohi[ix0][ix2][i];
               outputhilo[k][i] = crosshilo[ix0][ix2][i];
               outputlolo[k][i] = crosslolo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_cmplx4_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   double ***forwardrehihi, double ***forwardrelohi,
   double ***forwardrehilo, double ***forwardrelolo,
   double ***forwardimhihi, double ***forwardimlohi,
   double ***forwardimhilo, double ***forwardimlolo,
   double ***backwardrehihi, double ***backwardrelohi,
   double ***backwardrehilo, double ***backwardrelolo,
   double ***backwardimhihi, double ***backwardimlohi,
   double ***backwardimhilo, double ***backwardimlolo,
   double ***crossrehihi, double ***crossrelohi,
   double ***crossrehilo, double ***crossrelolo,
   double ***crossimhihi, double ***crossimlohi,
   double ***crossimhilo, double ***crossimlolo,
   AdditionJobs jobs, bool verbose )
{
   for(int k=0; k<jobs.get_depth(); k++)
   {
      if(verbose) cout << "executing addition jobs at layer "
                       << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
      {
         AdditionJob job = jobs.get_job(k,i);
         if(verbose) cout << "job " << i << " : " << job << endl;

         CPU_cmplx4_add_job(deg,
            cstrehihi,cstrelohi,cstrehilo,cstrelolo,
            cstimhihi,cstimlohi,cstimhilo,cstimlolo,
            cffrehihi,cffrelohi,cffrehilo,cffrelolo,
            cffimhihi,cffimlohi,cffimhilo,cffimlolo,
            forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
            forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
            backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
            backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,
            crossrehihi,crossrelohi,crossrehilo,crossrelolo,
            crossimhihi,crossimlohi,crossimhilo,crossimlolo,job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputrehihi[dim][i] = forwardrehihi[lastmon][lastidx][i];
      outputrelohi[dim][i] = forwardrelohi[lastmon][lastidx][i];
      outputrehilo[dim][i] = forwardrehilo[lastmon][lastidx][i];
      outputrelolo[dim][i] = forwardrelolo[lastmon][lastidx][i];
      outputimhihi[dim][i] = forwardimhihi[lastmon][lastidx][i];
      outputimlohi[dim][i] = forwardimlohi[lastmon][lastidx][i];
      outputimhilo[dim][i] = forwardimhilo[lastmon][lastidx][i];
      outputimlolo[dim][i] = forwardimlolo[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputrehihi[0][i] = 0.0; outputrelohi[0][i] = 0.0;
         outputrehilo[0][i] = 0.0; outputrelolo[0][i] = 0.0;
         outputimhihi[0][i] = 0.0; outputimlohi[0][i] = 0.0; 
         outputimhilo[0][i] = 0.0; outputimlolo[0][i] = 0.0;
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0] - 2;
      
      if(verbose)
         cout << "Updating derivative 0, ix0 = " << ix0
              << ", ix2 = " << ix2
              << " : b[" << ix0 << "," << ix2 << "]" << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = backward[ix0][ix2][i];
      {
         outputrehihi[0][i] = backwardrehihi[ix0][ix2][i];
         outputrelohi[0][i] = backwardrelohi[ix0][ix2][i];
         outputrehilo[0][i] = backwardrehilo[ix0][ix2][i];
         outputrelolo[0][i] = backwardrelolo[ix0][ix2][i];
         outputimhihi[0][i] = backwardimhihi[ix0][ix2][i];
         outputimlohi[0][i] = backwardimlohi[ix0][ix2][i];
         outputimhilo[0][i] = backwardimhilo[ix0][ix2][i];
         outputimlolo[0][i] = backwardimlolo[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputrehihi[k][i] = 0.0; outputrelohi[k][i] = 0.0;
            outputrehilo[k][i] = 0.0; outputrelolo[k][i] = 0.0;
            outputimhihi[k][i] = 0.0; outputimlohi[k][i] = 0.0; 
            outputimhilo[k][i] = 0.0; outputimlolo[k][i] = 0.0;
         }
      }
      else
      {
         int ix0 = jobs.get_differential_index(k,cnt);

         if(idx[ix0][0] == k) // k is first variable of monomial
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : b[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = backward[ix0][ix2][i];
            {
               outputrehihi[k][i] = backwardrehihi[ix0][ix2][i];
               outputrelohi[k][i] = backwardrelohi[ix0][ix2][i];
               outputrehilo[k][i] = backwardrehilo[ix0][ix2][i];
               outputrelolo[k][i] = backwardrelolo[ix0][ix2][i];
               outputimhihi[k][i] = backwardimhihi[ix0][ix2][i];
               outputimlohi[k][i] = backwardimlohi[ix0][ix2][i];
               outputimhilo[k][i] = backwardimhilo[ix0][ix2][i];
               outputimlolo[k][i] = backwardimlolo[ix0][ix2][i];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : f[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = forward[ix0][ix2][i];
            {
               outputrehihi[k][i] = forwardrehihi[ix0][ix2][i];
               outputrelohi[k][i] = forwardrelohi[ix0][ix2][i];
               outputrehilo[k][i] = forwardrehilo[ix0][ix2][i];
               outputrelolo[k][i] = forwardrelolo[ix0][ix2][i];
               outputimhihi[k][i] = forwardimhihi[ix0][ix2][i];
               outputimlohi[k][i] = forwardimlohi[ix0][ix2][i];
               outputimhilo[k][i] = forwardimhilo[ix0][ix2][i];
               outputimlolo[k][i] = forwardimlolo[ix0][ix2][i];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : c[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = cross[ix0][ix2][i];
            {
               outputrehihi[k][i] = crossrehihi[ix0][ix2][i];
               outputrelohi[k][i] = crossrelohi[ix0][ix2][i];
               outputrehilo[k][i] = crossrehilo[ix0][ix2][i];
               outputrelolo[k][i] = crossrelolo[ix0][ix2][i];
               outputimhihi[k][i] = crossimhihi[ix0][ix2][i];
               outputimlohi[k][i] = crossimlohi[ix0][ix2][i];
               outputimhilo[k][i] = crossimhilo[ix0][ix2][i];
               outputimlolo[k][i] = crossimlolo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl4_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double ***forwardhihi = new double**[nbr];
   double ***forwardlohi = new double**[nbr];
   double ***forwardhilo = new double**[nbr];
   double ***forwardlolo = new double**[nbr];
   double ***backwardhihi = new double**[nbr];
   double ***backwardlohi = new double**[nbr];
   double ***backwardhilo = new double**[nbr];
   double ***backwardlolo = new double**[nbr];
   double ***crosshihi = new double**[nbr];
   double ***crosslohi = new double**[nbr];
   double ***crosshilo = new double**[nbr];
   double ***crosslolo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardhihi[k] = new double*[nvrk];
      forwardlohi[k] = new double*[nvrk];
      forwardhilo[k] = new double*[nvrk];
      forwardlolo[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardhihi[k][i] = new double[deg+1];
         forwardlohi[k][i] = new double[deg+1];
         forwardhilo[k][i] = new double[deg+1];
         forwardlolo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardhihi[k] = new double*[nvrk-1];
         backwardlohi[k] = new double*[nvrk-1];
         backwardhilo[k] = new double*[nvrk-1];
         backwardlolo[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardhihi[k][i] = new double[deg+1];
            backwardlohi[k][i] = new double[deg+1];
            backwardhilo[k][i] = new double[deg+1];
            backwardlolo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crosshihi[k] = new double*[nvrk-2];
         crosslohi[k] = new double*[nvrk-2];
         crosshilo[k] = new double*[nvrk-2];
         crosslolo[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crosshihi[k][i] = new double[deg+1];
            crosslohi[k][i] = new double[deg+1];
            crosshilo[k][i] = new double[deg+1];
            crosslolo[k][i] = new double[deg+1];
         }
      }
   }
   clock_t start = clock();
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(vrb) cout << "executing convolution jobs at layer "
                   << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = cnvjobs.get_job(k,i);
         if(vrb) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_dbl4_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffhihi[monidx],cfflohi[monidx],cffhilo[monidx],cfflolo[monidx],
             inputhihi,inputlohi,inputhilo,inputlolo,
              forwardhihi[monidx], forwardlohi[monidx],
              forwardhilo[monidx], forwardlolo[monidx],
             backwardhihi[monidx],backwardlohi[monidx],
             backwardhilo[monidx],backwardlolo[monidx],
                crosshihi[monidx],   crosslohi[monidx],
                crosshilo[monidx],   crosslolo[monidx],job,vrb);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl4_poly_addjobs
      (dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
            cffhihi,     cfflohi,     cffhilo,     cfflolo,
          inputhihi,   inputlohi,   inputhilo,   inputlolo,
         outputhihi,  outputlohi,  outputhilo,  outputlolo,
        forwardhihi, forwardlohi, forwardhilo, forwardlolo,
       backwardhihi,backwardlohi,backwardhilo,backwardlolo,
          crosshihi,   crosslohi,   crosshilo,   crosslolo,addjobs,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
   }
   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k];

      for(int i=0; i<nvrk; i++)
      {
         free(forwardhihi[k][i]);
         free(forwardlohi[k][i]);
         free(forwardhilo[k][i]);
         free(forwardlolo[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardhihi[k][i]);
                      free(backwardlohi[k][i]);
                      free(backwardhilo[k][i]);
                      free(backwardlolo[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crosshihi[k][i]);
                      free(crosslohi[k][i]);
                      free(crosshilo[k][i]);
                      free(crosslolo[k][i]);
                   }
   }
   free(forwardhihi); free(backwardhihi); free(crosshihi);
   free(forwardlohi); free(backwardlohi); free(crosslohi);
   free(forwardhilo); free(backwardhilo); free(crosshilo);
   free(forwardlolo); free(backwardlolo); free(crosslolo);
}

void CPU_cmplx4_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double ***forwardrehihi = new double**[nbr];
   double ***forwardrelohi = new double**[nbr];
   double ***forwardrehilo = new double**[nbr];
   double ***forwardrelolo = new double**[nbr];
   double ***forwardimhihi = new double**[nbr];
   double ***forwardimlohi = new double**[nbr];
   double ***forwardimhilo = new double**[nbr];
   double ***forwardimlolo = new double**[nbr];
   double ***backwardrehihi = new double**[nbr];
   double ***backwardrelohi = new double**[nbr];
   double ***backwardrehilo = new double**[nbr];
   double ***backwardrelolo = new double**[nbr];
   double ***backwardimhihi = new double**[nbr];
   double ***backwardimlohi = new double**[nbr];
   double ***backwardimhilo = new double**[nbr];
   double ***backwardimlolo = new double**[nbr];
   double ***crossrehihi = new double**[nbr];
   double ***crossrelohi = new double**[nbr];
   double ***crossrehilo = new double**[nbr];
   double ***crossrelolo = new double**[nbr];
   double ***crossimhihi = new double**[nbr];
   double ***crossimlohi = new double**[nbr];
   double ***crossimhilo = new double**[nbr];
   double ***crossimlolo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardrehihi[k] = new double*[nvrk];
      forwardrelohi[k] = new double*[nvrk];
      forwardrehilo[k] = new double*[nvrk];
      forwardrelolo[k] = new double*[nvrk];
      forwardimhihi[k] = new double*[nvrk];
      forwardimlohi[k] = new double*[nvrk];
      forwardimhilo[k] = new double*[nvrk];
      forwardimlolo[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardrehihi[k][i] = new double[deg+1];
         forwardrelohi[k][i] = new double[deg+1];
         forwardrehilo[k][i] = new double[deg+1];
         forwardrelolo[k][i] = new double[deg+1];
         forwardimhihi[k][i] = new double[deg+1];
         forwardimlohi[k][i] = new double[deg+1];
         forwardimhilo[k][i] = new double[deg+1];
         forwardimlolo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardrehihi[k] = new double*[nvrk-1];
         backwardrelohi[k] = new double*[nvrk-1];
         backwardrehilo[k] = new double*[nvrk-1];
         backwardrelolo[k] = new double*[nvrk-1];
         backwardimhihi[k] = new double*[nvrk-1];
         backwardimlohi[k] = new double*[nvrk-1];
         backwardimhilo[k] = new double*[nvrk-1];
         backwardimlolo[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardrehihi[k][i] = new double[deg+1];
            backwardrelohi[k][i] = new double[deg+1];
            backwardrehilo[k][i] = new double[deg+1];
            backwardrelolo[k][i] = new double[deg+1];
            backwardimhihi[k][i] = new double[deg+1];
            backwardimlohi[k][i] = new double[deg+1];
            backwardimhilo[k][i] = new double[deg+1];
            backwardimlolo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crossrehihi[k] = new double*[nvrk-2];
         crossrelohi[k] = new double*[nvrk-2];
         crossrehilo[k] = new double*[nvrk-2];
         crossrelolo[k] = new double*[nvrk-2];
         crossimhihi[k] = new double*[nvrk-2];
         crossimlohi[k] = new double*[nvrk-2];
         crossimhilo[k] = new double*[nvrk-2];
         crossimlolo[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crossrehihi[k][i] = new double[deg+1];
            crossrelohi[k][i] = new double[deg+1];
            crossrehilo[k][i] = new double[deg+1];
            crossrelolo[k][i] = new double[deg+1];
            crossimhihi[k][i] = new double[deg+1];
            crossimlohi[k][i] = new double[deg+1];
            crossimhilo[k][i] = new double[deg+1];
            crossimlolo[k][i] = new double[deg+1];
         }
      }
   }
   clock_t start = clock();
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(vrb) cout << "executing convolution jobs at layer "
                   << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = cnvjobs.get_job(k,i);
         if(vrb) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_cmplx4_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffrehihi[monidx],cffrelohi[monidx],
             cffrehilo[monidx],cffrelolo[monidx],
             cffimhihi[monidx],cffimlohi[monidx],
             cffimhilo[monidx],cffimlolo[monidx],
             inputrehihi,inputrelohi,inputrehilo,inputrelolo,
             inputimhihi,inputimlohi,inputimhilo,inputimlolo,
             forwardrehihi[monidx],forwardrelohi[monidx],
             forwardrehilo[monidx],forwardrelolo[monidx],
             forwardimhihi[monidx],forwardimlohi[monidx],
             forwardimhilo[monidx],forwardimlolo[monidx],
             backwardrehihi[monidx],backwardrelohi[monidx],
             backwardrehilo[monidx],backwardrelolo[monidx],
             backwardimhihi[monidx],backwardimlohi[monidx],
             backwardimhilo[monidx],backwardimlolo[monidx],
             crossrehihi[monidx],crossrelohi[monidx],
             crossrehilo[monidx],crossrelolo[monidx],
             crossimhihi[monidx],crossimlohi[monidx],
             crossimhilo[monidx],crossimlolo[monidx],job,vrb);
      }
   }
 /*
   CPU_cmplx4_poly_updates
      (dim,nbr,deg,nvr,idx,
       cstrehihi,cstrelohi,cstrehilo,cstrelolo,
       cstimhihi,cstimlohi,cstimhilo,cstimlolo,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       inputrehihi,inputrelohi,inputrehilo,inputrelolo,
       inputimhihi,inputimlohi,inputimhilo,inputimlolo,
       outputrehihi,outputrelohi,outputrehilo,outputrelolo,
       outputimhihi,outputimlohi,outputimhilo,outputimlolo,
       forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
       forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
       backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
       backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,
       crossrehihi,crossrelohi,crossrehilo,crossrelolo,
       crossimhihi,crossimlohi,crossimhilo,crossimlolo);
  */
   CPU_cmplx4_poly_addjobs
      (dim,nbr,deg,nvr,idx,
       cstrehihi,cstrelohi,cstrehilo,cstrelolo,
       cstimhihi,cstimlohi,cstimhilo,cstimlolo,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       inputrehihi,inputrelohi,inputrehilo,inputrelolo,
       inputimhihi,inputimlohi,inputimhilo,inputimlolo,
       outputrehihi,outputrelohi,outputrehilo,outputrelolo,
       outputimhihi,outputimlohi,outputimhilo,outputimlolo,
       forwardrehihi,forwardrelohi,forwardrehilo,forwardrelolo,
       forwardimhihi,forwardimlohi,forwardimhilo,forwardimlolo,
       backwardrehihi,backwardrelohi,backwardrehilo,backwardrelolo,
       backwardimhihi,backwardimlohi,backwardimhilo,backwardimlolo,
       crossrehihi,crossrelohi,crossrehilo,crossrelolo,
       crossimhihi,crossimlohi,crossimhilo,crossimlolo,addjobs,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
      cout << scientific << setprecision(16);
   }
   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k];

      for(int i=0; i<nvrk; i++)
      {
         free(forwardrehihi[k][i]); free(forwardrelohi[k][i]);
         free(forwardrehilo[k][i]); free(forwardrelolo[k][i]);
         free(forwardimhihi[k][i]); free(forwardimlohi[k][i]);
         free(forwardimhilo[k][i]); free(forwardimlolo[k][i]);
      }
      if(nvrk > 1)
         for(int i=0; i<nvrk-1; i++)
         {
            free(backwardrehihi[k][i]); free(backwardrelohi[k][i]);
            free(backwardrehilo[k][i]); free(backwardrelolo[k][i]);
            free(backwardimhihi[k][i]); free(backwardimlohi[k][i]);
            free(backwardimhilo[k][i]); free(backwardimlolo[k][i]);
         }
      if(nvrk > 2)
         for(int i=0; i<nvrk-2; i++)
         {
            free(crossrehihi[k][i]); free(crossrelohi[k][i]);
            free(crossrehilo[k][i]); free(crossrelolo[k][i]);
            free(crossimhihi[k][i]); free(crossimlohi[k][i]);
            free(crossimhilo[k][i]); free(crossimlolo[k][i]);
         }
   }
   free(forwardrehihi); free(backwardrehihi); free(crossrehihi);
   free(forwardrelohi); free(backwardrelohi); free(crossrelohi);
   free(forwardrehilo); free(backwardrehilo); free(crossrehilo);
   free(forwardrelolo); free(backwardrelolo); free(crossrelolo);
   free(forwardimhihi); free(backwardimhihi); free(crossimhihi);
   free(forwardimlohi); free(backwardimlohi); free(crossimlohi);
   free(forwardimhilo); free(backwardimhilo); free(crossimhilo);
   free(forwardimlolo); free(backwardimlolo); free(crossimlolo);
}
