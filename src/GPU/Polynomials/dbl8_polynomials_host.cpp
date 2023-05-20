/* The file dbl8_polynomials_host.cpp defines functions specified
 * in dbl8_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "octo_double_functions.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"
#include "dbl8_polynomials_host.h"

void CPU_dbl8_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
   double **forwardhihihi, double **forwardlohihi,
   double **forwardhilohi, double **forwardlolohi,
   double **forwardhihilo, double **forwardlohilo,
   double **forwardhilolo, double **forwardlololo,
   double **backwardhihihi, double **backwardlohihi,
   double **backwardhilohi, double **backwardlolohi,
   double **backwardhihilo, double **backwardlohilo,
   double **backwardhilolo, double **backwardlololo,
   double **crosshihihi, double **crosslohihi,
   double **crosshilohi, double **crosslolohi,
   double **crosshihilo, double **crosslohilo,
   double **crosshilolo, double **crosslololo, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl8_product(deg,inputhihihi[ix1],inputlohihi[ix1],
                              inputhilohi[ix1],inputlolohi[ix1],
                              inputhihilo[ix1],inputlohilo[ix1],
                              inputhilolo[ix1],inputlololo[ix1],
                                cffhihihi[i],    cfflohihi[i],
                                cffhilohi[i],    cfflolohi[i],
                                cffhihilo[i],    cfflohilo[i],
                                cffhilolo[i],    cfflololo[i],
                              forwardhihihi[0],forwardlohihi[0],
                              forwardhilohi[0],forwardlolohi[0],
                              forwardhihilo[0],forwardlohilo[0],
                              forwardhilolo[0],forwardlololo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            odf_inc(&outputhihihi[dim][j],&outputlohihi[dim][j],
                    &outputhilohi[dim][j],&outputlolohi[dim][j],
                    &outputhihilo[dim][j],&outputlohilo[dim][j],
                    &outputhilolo[dim][j],&outputlololo[dim][j],
                    forwardhihihi[0][j],  forwardlohihi[0][j],
                    forwardhilohi[0][j],  forwardlolohi[0][j],
                    forwardhihilo[0][j],  forwardlohilo[0][j],
                    forwardhilolo[0][j],  forwardlololo[0][j]);
            // output[ix1][j] += cff[i][j];
            odf_inc(&outputhihihi[ix1][j],&outputlohihi[ix1][j],
                    &outputhilohi[ix1][j],&outputlolohi[ix1][j],
                    &outputhihilo[ix1][j],&outputlohilo[ix1][j],
                    &outputhilolo[ix1][j],&outputlololo[ix1][j],
                        cffhihihi[i][j],      cfflohihi[i][j],
                        cffhilohi[i][j],      cfflolohi[i][j],
                        cffhihilo[i][j],      cfflohilo[i][j],
                        cffhilolo[i][j],      cfflololo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl8_product(deg, cffhihihi[i],    cfflohihi[i],
                               cffhilohi[i],    cfflolohi[i],
                               cffhihilo[i],    cfflohilo[i],
                               cffhilolo[i],    cfflololo[i],
                             inputhihihi[ix1],inputlohihi[ix1],
                             inputhilohi[ix1],inputlolohi[ix1],
                             inputhihilo[ix1],inputlohilo[ix1],
                             inputhilolo[ix1],inputlololo[ix1],
                           forwardhihihi[0],forwardlohihi[0],
                           forwardhilohi[0],forwardlolohi[0],
                           forwardhihilo[0],forwardlohilo[0],
                           forwardhilolo[0],forwardlololo[0]);
         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
            odf_inc(&outputhihihi[ix2][j],&outputlohihi[ix2][j],
                    &outputhilohi[ix2][j],&outputlolohi[ix2][j],
                    &outputhihilo[ix2][j],&outputlohilo[ix2][j],
                    &outputhilolo[ix2][j],&outputlololo[ix2][j],
                    forwardhihihi[0][j],   forwardlohihi[0][j],
                    forwardhilolo[0][j],   forwardlolohi[0][j],
                    forwardhihilo[0][j],   forwardlohilo[0][j],
                    forwardhilolo[0][j],   forwardlololo[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl8_product(deg, cffhihihi[i],     cfflohihi[i],
                               cfflohihi[i],     cfflolohi[i],
                               cffhihilo[i],     cfflohilo[i],
                               cffhilolo[i],     cfflololo[i],
                             inputhihihi[ix2], inputlohihi[ix2],
                             inputhilohi[ix2], inputlolohi[ix2],
                             inputhihilo[ix2], inputlohilo[ix2],
                             inputhilolo[ix2], inputlololo[ix2],
                          backwardhihihi[0],backwardlohihi[0],
                          backwardhilohi[0],backwardlolohi[0],
                          backwardhihilo[0],backwardlohilo[0],
                          backwardhilolo[0],backwardlololo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
            odf_inc( &outputhihihi[ix1][j],&outputlohihi[ix1][j],
                     &outputhilohi[ix1][j],&outputlolohi[ix1][j],
                     &outputhihilo[ix1][j],&outputlohilo[ix1][j],
                     &outputhilolo[ix1][j],&outputlololo[ix1][j],
                    backwardhihihi[0][j], backwardlohihi[0][j],
                    backwardhilohi[0][j], backwardlolohi[0][j],
                    backwardhihilo[0][j], backwardlohilo[0][j],
                    backwardhilolo[0][j], backwardlololo[0][j]);

         CPU_dbl8_product(deg,forwardhihihi[0],forwardlohihi[0],
                              forwardhilohi[0],forwardlolohi[0],
                              forwardhihilo[0],forwardlohilo[0],
                              forwardhilolo[0],forwardlololo[0],
                                inputhihihi[ix2],inputlohihi[ix2],
                                inputhilohi[ix2],inputlolohi[ix2],
                                inputhihilo[ix2],inputlohilo[ix2],
                                inputhilolo[ix2],inputlololo[ix2],
                              forwardhihihi[1],forwardlohihi[1],
                              forwardhilohi[1],forwardlolohi[1],
                              forwardhihilo[1],forwardlohilo[1],
                              forwardhilolo[1],forwardlololo[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
            odf_inc(&outputhihihi[dim][j],&outputlohihi[dim][j],
                    &outputhilohi[dim][j],&outputlolohi[dim][j],
                    &outputhihilo[dim][j],&outputlohilo[dim][j],
                    &outputhilolo[dim][j],&outputlololo[dim][j],
                    forwardhihihi[1][j],  forwardlohihi[1][j],
                    forwardhilohi[1][j],  forwardlolohi[1][j],
                    forwardhihilo[1][j],  forwardlohilo[1][j],
                    forwardhilolo[1][j],  forwardlololo[1][j]);
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl8_speel(nvr[i],deg,idx[i],
                 cffhihihi[i],  cfflohihi[i],  cffhilohi[i],  cfflolohi[i],
                 cffhihilo[i],  cfflohilo[i],  cffhilolo[i],  cfflololo[i],
               inputhihihi,   inputlohihi,   inputhilohi,   inputlolohi,
               inputhihilo,   inputlohilo,   inputhilolo,   inputlololo,
             forwardhihihi, forwardlohihi, forwardhilohi, forwardlolohi,
             forwardhihilo, forwardlohilo, forwardhilolo, forwardlololo,
            backwardhihihi,backwardlohihi,backwardhilohi,backwardlolohi,
            backwardhihilo,backwardlohilo,backwardhilolo,backwardlololo,
               crosshihihi,   crosslohihi,   crosshilohi,   crosslolohi,
               crosshihilo,   crosslohilo,   crosshilolo,   crosslololo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
            odf_inc(&outputhihihi[dim][j],&outputlohihi[dim][j],
                    &outputhilohi[dim][j],&outputlolohi[dim][j],
                    &outputhihilo[dim][j],&outputlohilo[dim][j],
                    &outputhilolo[dim][j],&outputlololo[dim][j],
                    forwardhihihi[ix1][j],forwardlohihi[ix1][j],
                    forwardhilohi[ix1][j],forwardlolohi[ix1][j],
                    forwardhihilo[ix1][j],forwardlohilo[ix1][j],
                    forwardhilolo[ix1][j],forwardlololo[ix1][j]);

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
            odf_inc(&outputhihihi[ix2][j],&outputlohihi[ix2][j],
                    &outputhilohi[ix2][j],&outputlolohi[ix2][j],
                    &outputhihilo[ix2][j],&outputlohilo[ix2][j],
                    &outputhilolo[ix2][j],&outputlololo[ix2][j],
                    forwardhihihi[ix1][j],forwardlohihi[ix1][j],
                    forwardhilohi[ix1][j],forwardlolohi[ix1][j],
                    forwardhihilo[ix1][j],forwardlohilo[ix1][j],
                    forwardhilolo[ix1][j],forwardlololo[ix1][j]);

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
            odf_inc( &outputhihihi[ix2][j], &outputlohihi[ix2][j],
                     &outputhilohi[ix2][j], &outputlolohi[ix2][j],
                     &outputhihilo[ix2][j], &outputlohilo[ix2][j],
                     &outputhilolo[ix2][j], &outputlololo[ix2][j],
                    backwardhihihi[ix1][j],backwardlohihi[ix1][j],
                    backwardhilohi[ix1][j],backwardlolohi[ix1][j],
                    backwardhihilo[ix1][j],backwardlohilo[ix1][j],
                    backwardhilolo[ix1][j],backwardlololo[ix1][j]);

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
               odf_inc(&outputhihihi[ix2][j],&outputlohihi[ix2][j],
                       &outputhilohi[ix2][j],&outputlolohi[ix2][j],
                       &outputhihilo[ix2][j],&outputlohilo[ix2][j],
                       &outputhilolo[ix2][j],&outputlololo[ix2][j],
                         crosshihihi[k-1][j],  crosslohihi[k-1][j],
                         crosshilohi[k-1][j],  crosslolohi[k-1][j],
                         crosshihilo[k-1][j],  crosslohilo[k-1][j],
                         crosshilolo[k-1][j],  crosslololo[k-1][j]);
         }
      }
   }
}

void CPU_cmplx8_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo,
   double **forwardrehihihi, double **forwardrelohihi,
   double **forwardrehilohi, double **forwardrelolohi,
   double **forwardrehihilo, double **forwardrelohilo,
   double **forwardrehilolo, double **forwardrelololo,
   double **forwardimhihihi, double **forwardimlohihi,
   double **forwardimhilohi, double **forwardimlolohi,
   double **forwardimhihilo, double **forwardimlohilo,
   double **forwardimhilolo, double **forwardimlololo,
   double **backwardrehihihi, double **backwardrelohihi,
   double **backwardrehilohi, double **backwardrelolohi,
   double **backwardrehihilo, double **backwardrelohilo,
   double **backwardrehilolo, double **backwardrelololo,
   double **backwardimhihihi, double **backwardimlohihi,
   double **backwardimhilohi, double **backwardimlolohi,
   double **backwardimhihilo, double **backwardimlohilo,
   double **backwardimhilolo, double **backwardimlololo,
   double **crossrehihihi, double **crossrelohihi,
   double **crossrehilohi, double **crossrelolohi,
   double **crossrehihilo, double **crossrelohilo,
   double **crossrehilolo, double **crossrelololo,
   double **crossimhihihi, double **crossimlohihi,
   double **crossimhilohi, double **crossimlolohi,
   double **crossimhihilo, double **crossimlohilo,
   double **crossimhilolo, double **crossimlololo,
   bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_cmplx8_product(deg,
            inputrehihihi[ix1],inputrelohihi[ix1],
            inputrehilohi[ix1],inputrelolohi[ix1],
            inputrehihilo[ix1],inputrelohilo[ix1],
            inputrehilolo[ix1],inputrelololo[ix1],
            inputimhihihi[ix1],inputimlohihi[ix1],
            inputimhilohi[ix1],inputimlolohi[ix1],
            inputimhihilo[ix1],inputimlohilo[ix1],
            inputimhilolo[ix1],inputimlololo[ix1],
            cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
            cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
            cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
            cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
            forwardrehihihi[0],forwardrelohihi[0],
            forwardrehilohi[0],forwardrelolohi[0],
            forwardrehihilo[0],forwardrelohilo[0],
            forwardrehilolo[0],forwardrelololo[0],
            forwardimhihihi[0],forwardimlohihi[0],
            forwardimhilohi[0],forwardimlolohi[0],
            forwardimhihilo[0],forwardimlohilo[0],
            forwardimhilolo[0],forwardimlololo[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            odf_inc
               (&outputrehihihi[dim][j],&outputrelohihi[dim][j],
                &outputrehilohi[dim][j],&outputrelolohi[dim][j],
                &outputrehihilo[dim][j],&outputrelohilo[dim][j],
                &outputrehilolo[dim][j],&outputrelololo[dim][j],
                forwardrehihihi[0][j],forwardrelohihi[0][j],
                forwardrehilohi[0][j],forwardrelolohi[0][j],
                forwardrehihilo[0][j],forwardrelohilo[0][j],
                forwardrehilolo[0][j],forwardrelololo[0][j]);
            odf_inc
               (&outputimhihihi[dim][j],&outputimlohihi[dim][j],
                &outputimhilohi[dim][j],&outputimlolohi[dim][j],
                &outputimhihilo[dim][j],&outputimlohilo[dim][j],
                &outputimhilolo[dim][j],&outputimlololo[dim][j],
                forwardimhihihi[0][j],forwardimlohihi[0][j],
                forwardimhilohi[0][j],forwardimlolohi[0][j],
                forwardimhihilo[0][j],forwardimlohilo[0][j],
                forwardimhilolo[0][j],forwardimlololo[0][j]);
            // output[ix1][j] += cff[i][j];
            odf_inc
               (&outputrehihihi[ix1][j],&outputrelohihi[ix1][j],
                &outputrehilohi[ix1][j],&outputrelolohi[ix1][j],
                &outputrehihilo[ix1][j],&outputrelohilo[ix1][j],
                &outputrehilolo[ix1][j],&outputrelololo[ix1][j],
                cffrehihihi[i][j],cffrelohihi[i][j],
                cffrehilohi[i][j],cffrelolohi[i][j],
                cffrehihilo[i][j],cffrelohilo[i][j],
                cffrehilolo[i][j],cffrelololo[i][j]);
            odf_inc
               (&outputimhihihi[ix1][j],&outputimlohihi[ix1][j],
                &outputimhilohi[ix1][j],&outputimlolohi[ix1][j],
                &outputimhihilo[ix1][j],&outputimlohilo[ix1][j],
                &outputimhilolo[ix1][j],&outputimlololo[ix1][j],
                cffimhihihi[i][j],cffimlohihi[i][j],
                cffimhilohi[i][j],cffimlolohi[i][j],
                cffimhihilo[i][j],cffimlohilo[i][j],
                cffimhilolo[i][j],cffimlololo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_cmplx8_product(deg,
            cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
            cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
            cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
            cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
            inputrehihihi[ix1],inputrelohihi[ix1],
            inputrehilohi[ix1],inputrelolohi[ix1],
            inputrehihilo[ix1],inputrelohilo[ix1],
            inputrehilolo[ix1],inputrelololo[ix1],
            inputimhihihi[ix1],inputimlohihi[ix1],
            inputimhilohi[ix1],inputimlolohi[ix1],
            inputimhihilo[ix1],inputimlohilo[ix1],
            inputimhilolo[ix1],inputimlololo[ix1],
            forwardrehihihi[0],forwardrelohihi[0],
            forwardrehilohi[0],forwardrelolohi[0],
            forwardrehihilo[0],forwardrelohilo[0],
            forwardrehilolo[0],forwardrelololo[0],
            forwardimhihihi[0],forwardimlohihi[0],
            forwardimhilohi[0],forwardimlolohi[0],
            forwardimhihilo[0],forwardimlohilo[0],
            forwardimhilolo[0],forwardimlololo[0]);

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
         {
            odf_inc
               (&outputrehihihi[ix2][j],&outputrelohihi[ix2][j],
                &outputrehilohi[ix2][j],&outputrelolohi[ix2][j],
                &outputrehihilo[ix2][j],&outputrelohilo[ix2][j],
                &outputrehilolo[ix2][j],&outputrelololo[ix2][j],
                forwardrehihihi[0][j],forwardrelohihi[0][j],
                forwardrehilohi[0][j],forwardrelolohi[0][j],
                forwardrehihilo[0][j],forwardrelohilo[0][j],
                forwardrehilolo[0][j],forwardrelololo[0][j]);
            odf_inc
               (&outputimhihihi[ix2][j],&outputimlohihi[ix2][j],
                &outputimhilohi[ix2][j],&outputimlolohi[ix2][j],
                &outputimhihilo[ix2][j],&outputimlohilo[ix2][j],
                &outputimhilolo[ix2][j],&outputimlololo[ix2][j],
                forwardimhihihi[0][j],forwardimlohihi[0][j],
                forwardimhilohi[0][j],forwardimlolohi[0][j],
                forwardimhihilo[0][j],forwardimlohilo[0][j],
                forwardimhilolo[0][j],forwardimlololo[0][j]);
         }
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_cmplx8_product(deg,
            cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
            cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
            cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
            cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
            inputrehihihi[ix2],inputrelohihi[ix2],
            inputrehilohi[ix2],inputrelolohi[ix2],
            inputrehihilo[ix2],inputrelohilo[ix2],
            inputrehilolo[ix2],inputrelololo[ix2],
            inputimhihihi[ix2],inputimlohihi[ix2],
            inputimhilohi[ix2],inputimlolohi[ix2],
            inputimhihilo[ix2],inputimlohilo[ix2],
            inputimhilolo[ix2],inputimlololo[ix2],
            backwardrehihihi[0],backwardrelohihi[0],
            backwardrehilohi[0],backwardrelolohi[0],
            backwardrehihilo[0],backwardrelohilo[0],
            backwardrehilolo[0],backwardrelololo[0],
            backwardimhihihi[0],backwardimlohihi[0],
            backwardimhilohi[0],backwardimlolohi[0],
            backwardimhihilo[0],backwardimlohilo[0],
            backwardimhilolo[0],backwardimlololo[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
         {
            odf_inc
               (&outputrehihihi[ix1][j],&outputrelohihi[ix1][j],
                &outputrehilohi[ix1][j],&outputrelolohi[ix1][j],
                &outputrehihilo[ix1][j],&outputrelohilo[ix1][j],
                &outputrehilolo[ix1][j],&outputrelololo[ix1][j],
                backwardrehihihi[0][j],backwardrelohihi[0][j],
                backwardrehilohi[0][j],backwardrelolohi[0][j],
                backwardrehihilo[0][j],backwardrelohilo[0][j],
                backwardrehilolo[0][j],backwardrelololo[0][j]);
            odf_inc
               (&outputimhihihi[ix1][j],&outputimlohihi[ix1][j],
                &outputimhilohi[ix1][j],&outputimlolohi[ix1][j],
                &outputimhihilo[ix1][j],&outputimlohilo[ix1][j],
                &outputimhilolo[ix1][j],&outputimlololo[ix1][j],
                backwardimhihihi[0][j],backwardimlohihi[0][j],
                backwardimhilohi[0][j],backwardimlolohi[0][j],
                backwardimhihilo[0][j],backwardimlohilo[0][j],
                backwardimhilolo[0][j],backwardimlololo[0][j]);
         }
         CPU_cmplx8_product(deg,
            forwardrehihihi[0],forwardrelohihi[0],
            forwardrehilohi[0],forwardrelolohi[0],
            forwardrehihilo[0],forwardrelohilo[0],
            forwardrehilolo[0],forwardrelololo[0],
            forwardimhihihi[0],forwardimlohihi[0],
            forwardimhilohi[0],forwardimlolohi[0],
            forwardimhihilo[0],forwardimlohilo[0],
            forwardimhilolo[0],forwardimlololo[0],
            inputrehihihi[ix2],inputrelohihi[ix2],
            inputrehilohi[ix2],inputrelolohi[ix2],
            inputrehihilo[ix2],inputrelohilo[ix2],
            inputrehilolo[ix2],inputrelololo[ix2],
            inputimhihihi[ix2],inputimlohihi[ix2],
            inputimhilohi[ix2],inputimlolohi[ix2],
            inputimhihilo[ix2],inputimlohilo[ix2],
            inputimhilolo[ix2],inputimlololo[ix2],
            forwardrehihihi[1],forwardrelohihi[1],
            forwardrehilohi[1],forwardrelolohi[1],
            forwardrehihilo[1],forwardrelohilo[1],
            forwardrehilolo[1],forwardrelololo[1],
            forwardimhihihi[1],forwardimlohihi[1],
            forwardimhilohi[1],forwardimlolohi[1],
            forwardimhihilo[1],forwardimlohilo[1],
            forwardimhilolo[1],forwardimlololo[1]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
         {
            odf_inc
               (&outputrehihihi[dim][j],&outputrelohihi[dim][j],
                &outputrehilohi[dim][j],&outputrelolohi[dim][j],
                &outputrehihilo[dim][j],&outputrelohilo[dim][j],
                &outputrehilolo[dim][j],&outputrelololo[dim][j],
                forwardrehihihi[1][j],forwardrelohihi[1][j],
                forwardrehilohi[1][j],forwardrelolohi[1][j],
                forwardrehihilo[1][j],forwardrelohilo[1][j],
                forwardrehilolo[1][j],forwardrelololo[1][j]);
            odf_inc
               (&outputimhihihi[dim][j],&outputimlohihi[dim][j],
                &outputimhilohi[dim][j],&outputimlolohi[dim][j],
                &outputimhihilo[dim][j],&outputimlohilo[dim][j],
                &outputimhilolo[dim][j],&outputimlololo[dim][j],
                forwardimhihihi[1][j],forwardimlohihi[1][j],
                forwardimhilohi[1][j],forwardimlolohi[1][j],
                forwardimhihilo[1][j],forwardimlohilo[1][j],
                forwardimhilolo[1][j],forwardimlololo[1][j]);
         }
      }
      else if(nvr[i] > 2)
      {
         CPU_cmplx8_speel
            (nvr[i],deg,idx[i],
             cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
             cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
             cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
             cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
             inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
             inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
             inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
             inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
             forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
             forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
             forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
             forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
             backwardrehihihi,backwardrelohihi,
             backwardrehilohi,backwardrelolohi,
             backwardrehihilo,backwardrelohilo,
             backwardrehilolo,backwardrelololo,
             backwardimhihihi,backwardimlohihi,
             backwardimhilohi,backwardimlolohi,
             backwardimhihilo,backwardimlohilo,
             backwardimhilolo,backwardimlololo,
             crossrehihihi,crossrelohihi,crossrehilohi,crossrelolohi,
             crossrehihilo,crossrelohilo,crossrehilolo,crossrelololo,
             crossimhihihi,crossimlohihi,crossimhilohi,crossimlolohi,
             crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
         {
            odf_inc
               (&outputrehihihi[dim][j],&outputrelohihi[dim][j],
                &outputrehilohi[dim][j],&outputrelolohi[dim][j],
                &outputrehihilo[dim][j],&outputrelohilo[dim][j],
                &outputrehilolo[dim][j],&outputrelololo[dim][j],
                forwardrehihihi[ix1][j],forwardrelohihi[ix1][j],
                forwardrehilohi[ix1][j],forwardrelolohi[ix1][j],
                forwardrehihilo[ix1][j],forwardrelohilo[ix1][j],
                forwardrehilolo[ix1][j],forwardrelololo[ix1][j]);
            odf_inc
               (&outputimhihihi[dim][j],&outputimlohihi[dim][j],
                &outputimhilohi[dim][j],&outputimlolohi[dim][j],
                &outputimhihilo[dim][j],&outputimlohilo[dim][j],
                &outputimhilolo[dim][j],&outputimlololo[dim][j],
                forwardimhihihi[ix1][j],forwardimlohihi[ix1][j],
                forwardimhilohi[ix1][j],forwardimlolohi[ix1][j],
                forwardimhihilo[ix1][j],forwardimlohilo[ix1][j],
                forwardimhilolo[ix1][j],forwardimlololo[ix1][j]);
         }
         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
         {
            odf_inc
               (&outputrehihihi[ix2][j],&outputrelohihi[ix2][j],
                &outputrehilohi[ix2][j],&outputrelolohi[ix2][j],
                &outputrehihilo[ix2][j],&outputrelohilo[ix2][j],
                &outputrehilolo[ix2][j],&outputrelololo[ix2][j],
                forwardrehihihi[ix1][j],forwardrelohihi[ix1][j],
                forwardrehilohi[ix1][j],forwardrelolohi[ix1][j],
                forwardrehihilo[ix1][j],forwardrelohilo[ix1][j],
                forwardrehilolo[ix1][j],forwardrelololo[ix1][j]);
            odf_inc
               (&outputimhihihi[ix2][j],&outputimlohihi[ix2][j],
                &outputimhilohi[ix2][j],&outputimlolohi[ix2][j],
                &outputimhihilo[ix2][j],&outputimlohilo[ix2][j],
                &outputimhilolo[ix2][j],&outputimlololo[ix2][j],
                forwardimhihihi[ix1][j],forwardimlohihi[ix1][j],
                forwardimhilohi[ix1][j],forwardimlolohi[ix1][j],
                forwardimhihilo[ix1][j],forwardimlohilo[ix1][j],
                forwardimhilolo[ix1][j],forwardimlololo[ix1][j]);
         }
         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
         {
            odf_inc
               (&outputrehihihi[ix2][j],&outputrelohihi[ix2][j],
                &outputrehilohi[ix2][j],&outputrelolohi[ix2][j],
                &outputrehihilo[ix2][j],&outputrelohilo[ix2][j],
                &outputrehilolo[ix2][j],&outputrelololo[ix2][j],
                backwardrehihihi[ix1][j],backwardrelohihi[ix1][j],
                backwardrehilohi[ix1][j],backwardrelolohi[ix1][j],
                backwardrehihilo[ix1][j],backwardrelohilo[ix1][j],
                backwardrehilolo[ix1][j],backwardrelololo[ix1][j]);
            odf_inc
               (&outputimhihihi[ix2][j],&outputimlohihi[ix2][j],
                &outputimhilohi[ix2][j],&outputimlolohi[ix2][j],
                &outputimhihilo[ix2][j],&outputimlohilo[ix2][j],
                &outputimhilolo[ix2][j],&outputimlololo[ix2][j],
                backwardimhihihi[ix1][j],backwardimlohihi[ix1][j],
                backwardimhilohi[ix1][j],backwardimlolohi[ix1][j],
                backwardimhihilo[ix1][j],backwardimlohilo[ix1][j],
                backwardimhilolo[ix1][j],backwardimlololo[ix1][j]);
         }
         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
            {
               odf_inc
                  (&outputrehihihi[ix2][j],&outputrelohihi[ix2][j],
                   &outputrehilohi[ix2][j],&outputrelolohi[ix2][j],
                   &outputrehihilo[ix2][j],&outputrelohilo[ix2][j],
                   &outputrehilolo[ix2][j],&outputrelololo[ix2][j],
                   crossrehihihi[k-1][j],crossrelohihi[k-1][j],
                   crossrehilohi[k-1][j],crossrelolohi[k-1][j],
                   crossrehihilo[k-1][j],crossrelohilo[k-1][j],
                   crossrehilolo[k-1][j],crossrelololo[k-1][j]);
               odf_inc
                  (&outputimhihihi[ix2][j],&outputimlohihi[ix2][j],
                   &outputimhilohi[ix2][j],&outputimlolohi[ix2][j],
                   &outputimhihilo[ix2][j],&outputimlohilo[ix2][j],
                   &outputimhilolo[ix2][j],&outputimlololo[ix2][j],
                   crossimhihihi[k-1][j],crossimlohihi[k-1][j],
                   crossimhilohi[k-1][j],crossimlolohi[k-1][j],
                   crossimhihilo[k-1][j],crossimlohilo[k-1][j],
                   crossimhilolo[k-1][j],crossimlololo[k-1][j]);
            }
         }
      }
   }
}

void CPU_dbl8_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo, 
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo, 
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   double *elapsedsec, bool verbose )
{
   double **forwardhihihi = new double*[dim];
   double **forwardhilohi = new double*[dim];
   double **forwardhihilo = new double*[dim];
   double **forwardhilolo = new double*[dim];
   double **forwardlohihi = new double*[dim];
   double **forwardlolohi = new double*[dim];
   double **forwardlohilo = new double*[dim];
   double **forwardlololo = new double*[dim];
   double **backwardhihihi = new double*[dim-1]; // in case dim = 2
   double **backwardhilohi = new double*[dim-1];
   double **backwardhihilo = new double*[dim-1];
   double **backwardhilolo = new double*[dim-1];
   double **backwardlohihi = new double*[dim-1]; 
   double **backwardlolohi = new double*[dim-1];
   double **backwardlohilo = new double*[dim-1];
   double **backwardlololo = new double*[dim-1];
   double **crosshihihi = new double*[dim-1];    // in case dim = 2
   double **crosshilohi = new double*[dim-1];
   double **crosshihilo = new double*[dim-1];
   double **crosshilolo = new double*[dim-1];
   double **crosslohihi = new double*[dim-1];
   double **crosslolohi = new double*[dim-1];
   double **crosslohilo = new double*[dim-1];
   double **crosslololo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardhihihi[i] = new double[deg+1];
      forwardhilohi[i] = new double[deg+1];
      forwardhihilo[i] = new double[deg+1];
      forwardhilolo[i] = new double[deg+1];
      forwardlohihi[i] = new double[deg+1];
      forwardlolohi[i] = new double[deg+1];
      forwardlohilo[i] = new double[deg+1];
      forwardlololo[i] = new double[deg+1];
      backwardhihihi[i] = new double[deg+1];
      backwardhilohi[i] = new double[deg+1];
      backwardhihilo[i] = new double[deg+1];
      backwardhilolo[i] = new double[deg+1];
      backwardlohihi[i] = new double[deg+1];
      backwardlolohi[i] = new double[deg+1];
      backwardlohilo[i] = new double[deg+1];
      backwardlololo[i] = new double[deg+1];
      crosshihihi[i] = new double[deg+1];
      crosshilohi[i] = new double[deg+1];
      crosshihilo[i] = new double[deg+1];
      crosshilolo[i] = new double[deg+1];
      crosslohihi[i] = new double[deg+1];
      crosslolohi[i] = new double[deg+1];
      crosslohilo[i] = new double[deg+1];
      crosslololo[i] = new double[deg+1];
   }
   forwardhihihi[dim-1] = new double[deg+1];
   forwardhilohi[dim-1] = new double[deg+1];
   forwardhihilo[dim-1] = new double[deg+1];
   forwardhilolo[dim-1] = new double[deg+1];
   forwardlohihi[dim-1] = new double[deg+1];
   forwardlolohi[dim-1] = new double[deg+1];
   forwardlohilo[dim-1] = new double[deg+1];
   forwardlololo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputhihihi[dim][i] = csthihihi[i];
      outputhilohi[dim][i] = csthilohi[i];
      outputhihilo[dim][i] = csthihilo[i];
      outputhilolo[dim][i] = csthilolo[i];
      outputlohihi[dim][i] = cstlohihi[i];
      outputlolohi[dim][i] = cstlolohi[i];
      outputlohilo[dim][i] = cstlohilo[i];
      outputlololo[dim][i] = cstlololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhihihi[i][j] = 0.0;
         outputhilohi[i][j] = 0.0;
         outputhihilo[i][j] = 0.0;
         outputhilolo[i][j] = 0.0;
         outputlohihi[i][j] = 0.0;
         outputlolohi[i][j] = 0.0;
         outputlohilo[i][j] = 0.0;
         outputlololo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl8_poly_speel
      (dim,nbr,deg,nvr,idx,
            cffhihihi,     cffhilohi,     cffhihilo,     cffhilolo,
            cfflohihi,     cfflolohi,     cfflohilo,     cfflololo,
          inputhihihi,   inputhilohi,   inputhihilo,   inputhilolo,
          inputlohihi,   inputlolohi,   inputlohilo,   inputlololo,
         outputhihihi,  outputhilohi,  outputhihilo,  outputhilolo,
         outputlohihi,  outputlolohi,  outputlohilo,  outputlololo,
        forwardhihihi, forwardhilohi, forwardhihilo, forwardhilolo,
        forwardlohihi, forwardlolohi, forwardlohilo, forwardlololo,
       backwardhihihi,backwardhilohi,backwardhihilo,backwardhilolo,
       backwardlohihi,backwardlolohi,backwardlohilo,backwardlololo,
          crosshihihi,   crosshilohi,   crosshihilo,   crosshilolo,
          crosslohihi,   crosslolohi,   crosslohilo,   crosslololo,verbose);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(verbose)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
   }
   for(int i=0; i<dim-1; i++)
   {
      free(forwardhihihi[i]); free(backwardhihihi[i]); free(crosshihihi[i]);
      free(forwardhilohi[i]); free(backwardhilohi[i]); free(crosshilohi[i]);
      free(forwardhihilo[i]); free(backwardhihilo[i]); free(crosshihilo[i]);
      free(forwardhilolo[i]); free(backwardhilolo[i]); free(crosshilolo[i]);
      free(forwardlohihi[i]); free(backwardlohihi[i]); free(crosslohihi[i]);
      free(forwardlolohi[i]); free(backwardlolohi[i]); free(crosslolohi[i]);
      free(forwardlohilo[i]); free(backwardlohilo[i]); free(crosslohilo[i]);
      free(forwardlololo[i]); free(backwardlololo[i]); free(crosslololo[i]);
   }
   free(forwardhihihi[dim-1]); free(forwardhilohi[dim-1]);
   free(forwardhihilo[dim-1]); free(forwardhilolo[dim-1]);
   free(forwardlohihi[dim-1]); free(forwardlolohi[dim-1]);
   free(forwardlohilo[dim-1]); free(forwardlololo[dim-1]);
   free(forwardhihihi); free(backwardhihihi); free(crosshihihi);
   free(forwardhilohi); free(backwardhilohi); free(crosshilohi);
   free(forwardhihilo); free(backwardhihilo); free(crosshihilo);
   free(forwardhilolo); free(backwardhilolo); free(crosshilolo);
   free(forwardlohihi); free(backwardlohihi); free(crosslohihi);
   free(forwardlolohi); free(backwardlolohi); free(crosslolohi);
   free(forwardlohilo); free(backwardlohilo); free(crosslohilo);
   free(forwardlololo); free(backwardlololo); free(crosslololo);
}

void CPU_dbl8_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhihihi, double *cffhilohi,
   double *cffhihilo, double *cffhilolo,
   double *cfflohihi, double *cfflolohi,
   double *cfflohilo, double *cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo,
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo,
   double **forwardhihihi, double **forwardhilohi,
   double **forwardhihilo, double **forwardhilolo,
   double **forwardlohihi, double **forwardlolohi,
   double **forwardlohilo, double **forwardlololo,
   double **backwardhihihi, double **backwardhilohi,
   double **backwardhihilo, double **backwardhilolo,
   double **backwardlohihi, double **backwardlolohi,
   double **backwardlohilo, double **backwardlololo,
   double **crosshihihi, double **crosshilohi,
   double **crosshihilo, double **crosshilolo,
   double **crosslohihi, double **crosslolohi,
   double **crosslohilo, double **crosslololo,
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
         CPU_dbl8_product(deg,
                cffhihihi,cffhilohi,cffhihilo,cffhilolo,
                cfflohihi,cfflolohi,cfflohilo,cfflololo,
              inputhihihi[inp2ix],  inputhilohi[inp2ix],
              inputhihilo[inp2ix],  inputhilolo[inp2ix],
              inputlohihi[inp2ix],  inputlolohi[inp2ix],
              inputlohilo[inp2ix],  inputlololo[inp2ix],
            forwardhihihi[outidx],forwardhilohi[outidx],
            forwardhihilo[outidx],forwardhilolo[outidx],
            forwardlohihi[outidx],forwardlolohi[outidx],
            forwardlohilo[outidx],forwardlololo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl8_product(deg,
               inputhihihi[inp1ix],inputhilohi[inp1ix],
               inputhihilo[inp1ix],inputhilolo[inp1ix],
               inputlohihi[inp1ix],inputlolohi[inp1ix],
               inputlohilo[inp1ix],inputlololo[inp1ix],
               cffhihihi,cffhilohi,cffhihilo,cffhilolo,
               cfflohihi,cfflolohi,cfflohilo,cfflololo,
               forwardhihihi[outidx],forwardhilohi[outidx],
               forwardhihilo[outidx],forwardhilolo[outidx],
               forwardlohihi[outidx],forwardlolohi[outidx],
               forwardlohilo[outidx],forwardlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
                 inputhihihi[inp1ix],  inputhilohi[inp1ix],
                 inputhihilo[inp1ix],  inputhilolo[inp1ix],
                 inputlohihi[inp1ix],  inputlolohi[inp1ix],
                 inputlohilo[inp1ix],  inputlololo[inp1ix],
               forwardhihihi[inp2ix],forwardhilohi[inp2ix],
               forwardhihilo[inp2ix],forwardhilolo[inp2ix],
               forwardlohihi[inp2ix],forwardlolohi[inp2ix],
               forwardlohilo[inp2ix],forwardlololo[inp2ix],
               forwardhihihi[outidx],forwardhilohi[outidx],
               forwardhihilo[outidx],forwardhilolo[outidx],
               forwardlohihi[outidx],forwardlolohi[outidx],
               forwardlohilo[outidx],forwardlololo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
              crosshihihi[inp1ix],  crosshilohi[inp1ix],
              crosshihilo[inp1ix],  crosshilolo[inp1ix],
              crosslohihi[inp1ix],  crosslolohi[inp1ix],
              crosslohilo[inp1ix],  crosslololo[inp1ix],
              inputhihihi[inp2ix],  inputhilohi[inp2ix],
              inputhihilo[inp2ix],  inputhilolo[inp2ix],
              inputlohihi[inp2ix],  inputlolohi[inp2ix],
              inputlohilo[inp2ix],  inputlololo[inp2ix],
            forwardhihihi[outidx],forwardhilohi[outidx],
            forwardhihilo[outidx],forwardhilolo[outidx],
            forwardlohihi[outidx],forwardlolohi[outidx],
            forwardlohilo[outidx],forwardlololo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl8_product(deg,
               inputhihihi[inp1ix],inputhilohi[inp1ix],
               inputhihilo[inp1ix],inputhilolo[inp1ix],
               inputlohihi[inp1ix],inputlolohi[inp1ix],
               inputlohilo[inp1ix],inputlololo[inp1ix],
               cffhihihi,cffhilohi,cffhihilo,cffhilolo,
               cfflohihi,cfflolohi,cfflohilo,cfflololo,
               forwardhihihi[outidx],forwardhilohi[outidx],
               forwardhihilo[outidx],forwardhilolo[outidx],
               forwardlohihi[outidx],forwardlolohi[outidx],
               forwardlohilo[outidx],forwardlololo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
               forwardhihihi[inp1ix],forwardhilohi[inp1ix],
               forwardhihilo[inp1ix],forwardhilolo[inp1ix],
               forwardlohihi[inp1ix],forwardlolohi[inp1ix],
               forwardlohilo[inp1ix],forwardlololo[inp1ix],
                 inputhihihi[inp2ix],  inputhilohi[inp2ix],
                 inputhihilo[inp2ix],  inputhilolo[inp2ix],
                 inputlohihi[inp2ix],  inputlolohi[inp2ix],
                 inputlohilo[inp2ix],  inputlololo[inp2ix],
               forwardhihihi[outidx],forwardhilohi[outidx],
               forwardhihilo[outidx],forwardhilolo[outidx],
               forwardlohihi[outidx],forwardlolohi[outidx],
               forwardlohilo[outidx],forwardlololo[outidx]);
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
            CPU_dbl8_product(deg,
                    cffhihihi,cffhilohi,cffhihilo,cffhilolo,
                    cfflohihi,cfflolohi,cfflohilo,cfflololo,
                  inputhihihi[inp2ix],   inputhilohi[inp2ix],
                  inputhihilo[inp2ix],   inputhilolo[inp2ix],
                  inputlohihi[inp2ix],   inputlolohi[inp2ix],
                  inputlohilo[inp2ix],   inputlololo[inp2ix],
               backwardhihihi[outidx],backwardhilohi[outidx],
               backwardhihilo[outidx],backwardhilolo[outidx],
               backwardlohihi[outidx],backwardlolohi[outidx],
               backwardlohilo[outidx],backwardlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
               cffhihihi,cffhilohi,cffhihilo,cffhilolo,
               cfflohihi,cfflolohi,cfflohilo,cfflololo,
               backwardhihihi[inp2ix],backwardhilohi[inp2ix],
               backwardhihilo[inp2ix],backwardhilolo[inp2ix],
               backwardlohihi[inp2ix],backwardlolohi[inp2ix],
               backwardlohilo[inp2ix],backwardlololo[inp2ix],
               backwardhihihi[outidx],backwardhilohi[outidx],
               backwardhihilo[outidx],backwardhilolo[outidx],
               backwardlohihi[outidx],backwardlolohi[outidx],
               backwardlohilo[outidx],backwardlololo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl8_product(deg,
                  inputhihihi[inp1ix],   inputhilohi[inp1ix],
                  inputhihilo[inp1ix],   inputhilolo[inp1ix],
                  inputlohihi[inp1ix],   inputlolohi[inp1ix],
                  inputlohilo[inp1ix],   inputlololo[inp1ix],
                  inputhihihi[inp2ix],   inputhilohi[inp2ix],
                  inputhihilo[inp2ix],   inputhilolo[inp2ix],
                  inputlohihi[inp2ix],   inputlolohi[inp2ix],
                  inputlohilo[inp2ix],   inputlololo[inp2ix],
               backwardhihihi[outidx],backwardhilohi[outidx],
               backwardhihilo[outidx],backwardhilolo[outidx],
               backwardlohihi[outidx],backwardlolohi[outidx],
               backwardlohilo[outidx],backwardlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
                  inputhihihi[inp1ix],   inputhilohi[inp1ix],
                  inputhihilo[inp1ix],   inputhilolo[inp1ix],
                  inputlohihi[inp1ix],   inputlolohi[inp1ix],
                  inputlohilo[inp1ix],   inputlololo[inp1ix],
               backwardhihihi[inp2ix],backwardhilohi[inp2ix],
               backwardhihilo[inp2ix],backwardhilolo[inp2ix],
               backwardlohihi[inp2ix],backwardlolohi[inp2ix],
               backwardlohilo[inp2ix],backwardlololo[inp2ix],
               backwardhihihi[outidx],backwardhilohi[outidx],
               backwardhihilo[outidx],backwardhilolo[outidx],
               backwardlohihi[outidx],backwardlolohi[outidx],
               backwardlohilo[outidx],backwardlololo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl8_product(deg,
               backwardhihihi[inp1ix],backwardhilohi[inp1ix],
               backwardhihilo[inp1ix],backwardhilolo[inp1ix],
               backwardlohihi[inp1ix],backwardlolohi[inp1ix],
               backwardlohilo[inp1ix],backwardlololo[inp1ix],
               cffhihihi,cffhilohi,cffhihilo,cffhilolo,
               cfflohihi,cfflolohi,cfflohilo,cfflololo,
               backwardhihihi[outidx],backwardhilohi[outidx],
               backwardhihilo[outidx],backwardhilolo[outidx],
               backwardlohihi[outidx],backwardlolohi[outidx],
               backwardlohilo[outidx],backwardlololo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
               backwardhihihi[inp1ix],backwardhilohi[inp1ix],
               backwardhihilo[inp1ix],backwardhilolo[inp1ix],
               backwardlohihi[inp1ix],backwardlolohi[inp1ix],
               backwardlohilo[inp1ix],backwardlololo[inp1ix],
                  inputhihihi[inp2ix],   inputhilohi[inp2ix],
                  inputhihilo[inp2ix],   inputhilolo[inp2ix],
                  inputlohihi[inp2ix],   inputlolohi[inp2ix],
                  inputlohilo[inp2ix],   inputlololo[inp2ix],
               backwardhihihi[outidx],backwardhilohi[outidx],
               backwardhihilo[outidx],backwardhilolo[outidx],
               backwardlohihi[outidx],backwardlolohi[outidx],
               backwardlohilo[outidx],backwardlololo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
            cffhihihi,cffhilohi,cffhihilo,cffhilolo,
            cfflohihi,cfflolohi,cfflohilo,cfflololo,
            inputhihihi[inp2ix],inputhilohi[inp2ix],
            inputhihilo[inp2ix],inputhilolo[inp2ix],
            inputlohihi[inp2ix],inputlolohi[inp2ix],
            inputlohilo[inp2ix],inputlololo[inp2ix],
            crosshihihi[outidx],crosshilohi[outidx],
            crosshihilo[outidx],crosshilolo[outidx],
            crosslohihi[outidx],crosslolohi[outidx],
            crosslohilo[outidx],crosslololo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
              inputhihihi[inp1ix],  inputhilohi[inp1ix],
              inputhihilo[inp1ix],  inputhilolo[inp1ix],
              inputlohihi[inp1ix],  inputlolohi[inp1ix],
              inputlohilo[inp1ix],  inputlololo[inp1ix],
            forwardhihihi[inp2ix],forwardhilohi[inp2ix],
            forwardhihilo[inp2ix],forwardhilolo[inp2ix],
            forwardlohihi[inp2ix],forwardlolohi[inp2ix],
            forwardlohilo[inp2ix],forwardlololo[inp2ix],
              crosshihihi[outidx],  crosshilohi[outidx],
              crosshihilo[outidx],  crosshilolo[outidx],
              crosslohihi[outidx],  crosslolohi[outidx],
              crosslohilo[outidx],  crosslololo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl8_product(deg,
              forwardhihihi[inp1ix],forwardhilohi[inp1ix],
              forwardhihilo[inp1ix],forwardhilolo[inp1ix],
              forwardlohihi[inp1ix],forwardlolohi[inp1ix],
              forwardlohilo[inp1ix],forwardlololo[inp1ix],
                inputhihihi[inp2ix],  inputhilohi[inp2ix],
                inputhihilo[inp2ix],  inputhilolo[inp2ix],
                inputlohihi[inp2ix],  inputlolohi[inp2ix],
                inputlohilo[inp2ix],  inputlololo[inp2ix],
                crosshihihi[outidx],  crosshilohi[outidx],
                crosshihilo[outidx],  crosshilolo[outidx],
                crosslohihi[outidx],  crosslolohi[outidx],
                crosslohilo[outidx],  crosslololo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl8_product(deg,
               forwardhihihi[inp1ix], forwardhilohi[inp1ix],
               forwardhihilo[inp1ix], forwardhilolo[inp1ix],
               forwardlohihi[inp1ix], forwardlolohi[inp1ix],
               forwardlohilo[inp1ix], forwardlololo[inp1ix],
              backwardhihihi[inp2ix],backwardhilohi[inp2ix],
              backwardhihilo[inp2ix],backwardhilolo[inp2ix],
              backwardlohihi[inp2ix],backwardlolohi[inp2ix],
              backwardlohilo[inp2ix],backwardlololo[inp2ix],
                 crosshihihi[outidx],   crosshilohi[outidx],
                 crosshihilo[outidx],   crosshilolo[outidx],
                 crosslohihi[outidx],   crosslolohi[outidx],
                 crosslohilo[outidx],   crosslololo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
            backwardhihihi[inp1ix],backwardhilohi[inp1ix],
            backwardhihilo[inp1ix],backwardhilolo[inp1ix],
            backwardlohihi[inp1ix],backwardlolohi[inp1ix],
            backwardlohilo[inp1ix],backwardlololo[inp1ix],
             forwardhihihi[inp2ix], forwardhilohi[inp2ix],
             forwardhihilo[inp2ix], forwardhilolo[inp2ix],
             forwardlohihi[inp2ix], forwardlolohi[inp2ix],
             forwardlohilo[inp2ix], forwardlololo[inp2ix],
               crosshihihi[outidx],   crosshilohi[outidx],
               crosshihilo[outidx],   crosshilolo[outidx],
               crosslohihi[outidx],   crosslolohi[outidx],
               crosslohilo[outidx],   crosslololo[outidx]);
      }
   }
}

void CPU_dbl8_add_job
 ( int deg,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double ***forwardhihihi, double ***forwardhilohi,
   double ***forwardhihilo, double ***forwardhilolo,
   double ***forwardlohihi, double ***forwardlolohi,
   double ***forwardlohilo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardhilohi,
   double ***backwardhihilo, double ***backwardhilolo, 
   double ***backwardlohihi, double ***backwardlolohi,
   double ***backwardlohilo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosshilohi,
   double ***crosshihilo, double ***crosshilolo,
   double ***crosslohihi, double ***crosslolohi,
   double ***crosslohilo, double ***crosslololo,
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
               odf_inc(&forwardhihihi[updmon][updidx][i],
                       &forwardhilohi[updmon][updidx][i],
                       &forwardhihilo[updmon][updidx][i],
                       &forwardhilolo[updmon][updidx][i],
                       &forwardlohihi[updmon][updidx][i],
                       &forwardlolohi[updmon][updidx][i],
                       &forwardlohilo[updmon][updidx][i],
                       &forwardlololo[updmon][updidx][i],
                       csthihihi[i],csthilohi[i],csthihilo[i],csthilolo[i],
                       cstlohihi[i],cstlolohi[i],cstlohilo[i],cstlololo[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               odf_inc(&forwardhihihi[updmon][updidx][i],
                       &forwardhilohi[updmon][updidx][i],
                       &forwardhihilo[updmon][updidx][i],
                       &forwardhilolo[updmon][updidx][i],
                       &forwardlohihi[updmon][updidx][i],
                       &forwardlolohi[updmon][updidx][i],
                       &forwardlohilo[updmon][updidx][i],
                       &forwardlololo[updmon][updidx][i],
                       cffhihihi[incidx][i],cffhilohi[incidx][i],
                       cffhihilo[incidx][i],cffhilolo[incidx][i],
                       cfflohihi[incidx][i],cfflolohi[incidx][i],
                       cfflohilo[incidx][i],cfflololo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            odf_inc(&forwardhihihi[updmon][updidx][i],
                    &forwardhilohi[updmon][updidx][i],
                    &forwardhihilo[updmon][updidx][i],
                    &forwardhilolo[updmon][updidx][i],
                    &forwardlohihi[updmon][updidx][i],
                    &forwardlolohi[updmon][updidx][i],
                    &forwardlohilo[updmon][updidx][i],
                    &forwardlololo[updmon][updidx][i],
                    forwardhihihi[incmon][incidx][i],
                    forwardhilohi[incmon][incidx][i],
                    forwardhihilo[incmon][incidx][i],
                    forwardhilolo[incmon][incidx][i],
                    forwardlohihi[incmon][incidx][i],
                    forwardlolohi[incmon][incidx][i],
                    forwardlohilo[incmon][incidx][i],
                    forwardlololo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            odf_inc(&forwardhihihi[updmon][updidx][i],
                    &forwardhilohi[updmon][updidx][i],
                    &forwardhihilo[updmon][updidx][i],
                    &forwardhilolo[updmon][updidx][i],
                    &forwardlohihi[updmon][updidx][i],
                    &forwardlolohi[updmon][updidx][i],
                    &forwardlohilo[updmon][updidx][i],
                    &forwardlololo[updmon][updidx][i],
                    backwardhihihi[incmon][incidx][i],
                    backwardhilohi[incmon][incidx][i],
                    backwardhihilo[incmon][incidx][i],
                    backwardhilolo[incmon][incidx][i],
                    backwardlohihi[incmon][incidx][i],
                    backwardlolohi[incmon][incidx][i],
                    backwardlohilo[incmon][incidx][i],
                    backwardlololo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            odf_inc(&forwardhihihi[updmon][updidx][i],
                    &forwardhilohi[updmon][updidx][i],
                    &forwardhihilo[updmon][updidx][i],
                    &forwardhilolo[updmon][updidx][i],
                    &forwardlohihi[updmon][updidx][i],
                    &forwardlolohi[updmon][updidx][i],
                    &forwardlohilo[updmon][updidx][i],
                    &forwardlololo[updmon][updidx][i],
                    crosshihihi[incmon][incidx][i],
                    crosshilohi[incmon][incidx][i],
                    crosshihilo[incmon][incidx][i],
                    crosshilolo[incmon][incidx][i],
                    crosslohihi[incmon][incidx][i],
                    crosslolohi[incmon][incidx][i],
                    crosslohilo[incmon][incidx][i],
                    crosslololo[incmon][incidx][i]);
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    cffhihihi[incidx][i],cffhilohi[incidx][i],
                    cffhihilo[incidx][i],cffhilolo[incidx][i],
                    cfflohihi[incidx][i],cfflolohi[incidx][i],
                    cfflohilo[incidx][i],cfflololo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    forwardhihihi[incmon][incidx][i],
                    forwardhilohi[incmon][incidx][i],
                    forwardhihilo[incmon][incidx][i],
                    forwardhilolo[incmon][incidx][i],
                    forwardlohihi[incmon][incidx][i],
                    forwardlolohi[incmon][incidx][i],
                    forwardlohilo[incmon][incidx][i],
                    forwardlololo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    backwardhihihi[incmon][incidx][i],
                    backwardhilohi[incmon][incidx][i],
                    backwardhihilo[incmon][incidx][i],
                    backwardhilolo[incmon][incidx][i],
                    backwardlohihi[incmon][incidx][i],
                    backwardlolohi[incmon][incidx][i],
                    backwardlohilo[incmon][incidx][i],
                    backwardlololo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    crosshihihi[incmon][incidx][i],
                    crosshilohi[incmon][incidx][i],
                    crosshihilo[incmon][incidx][i],
                    crosshilolo[incmon][incidx][i],
                    crosslohihi[incmon][incidx][i],
                    crosslolohi[incmon][incidx][i],
                    crosslohilo[incmon][incidx][i],
                    crosslololo[incmon][incidx][i]);
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    cffhihihi[incidx][i],cffhilohi[incidx][i],
                    cffhihilo[incidx][i],cffhilolo[incidx][i],
                    cfflohihi[incidx][i],cfflolohi[incidx][i],
                    cfflohilo[incidx][i],cfflololo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    forwardhihihi[incmon][incidx][i],
                    forwardhilohi[incmon][incidx][i],
                    forwardhihilo[incmon][incidx][i],
                    forwardhilolo[incmon][incidx][i],
                    forwardlohihi[incmon][incidx][i],
                    forwardlolohi[incmon][incidx][i],
                    forwardlohilo[incmon][incidx][i],
                    forwardlololo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    backwardhihihi[incmon][incidx][i],
                    backwardhilohi[incmon][incidx][i],
                    backwardhihilo[incmon][incidx][i],
                    backwardhilolo[incmon][incidx][i],
                    backwardlohihi[incmon][incidx][i],
                    backwardlolohi[incmon][incidx][i],
                    backwardlohilo[incmon][incidx][i],
                    backwardlololo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    crosshihihi[incmon][incidx][i],
                    crosshilohi[incmon][incidx][i],
                    crosshihilo[incmon][incidx][i],
                    crosshilolo[incmon][incidx][i],
                    crosslohihi[incmon][incidx][i],
                    crosslolohi[incmon][incidx][i],
                    crosslohilo[incmon][incidx][i],
                    crosslololo[incmon][incidx][i]);
      }
   }
}

void CPU_dbl8_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo, 
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo, 
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   double ***forwardhihihi, double ***forwardhilohi,
   double ***forwardhihilo, double ***forwardhilolo,
   double ***forwardlohihi, double ***forwardlolohi,
   double ***forwardlohilo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardhilohi,
   double ***backwardhihilo, double ***backwardhilolo, 
   double ***backwardlohihi, double ***backwardlolohi,
   double ***backwardlohilo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosshilohi,
   double ***crosshihilo, double ***crosshilolo,
   double ***crosslohihi, double ***crosslolohi,
   double ***crosslohilo, double ***crosslololo )
{
   for(int i=0; i<=deg; i++)
   {
      outputhihihi[dim][i] = csthihihi[i];
      outputhilohi[dim][i] = csthilohi[i];
      outputhihilo[dim][i] = csthihilo[i];
      outputhilolo[dim][i] = csthilolo[i];
      outputlohihi[dim][i] = cstlohihi[i];
      outputlolohi[dim][i] = cstlolohi[i];
      outputlohilo[dim][i] = cstlohilo[i];
      outputlololo[dim][i] = cstlololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhihihi[i][j] = 0.0;
         outputhilohi[i][j] = 0.0;
         outputhihilo[i][j] = 0.0;
         outputhilolo[i][j] = 0.0;
         outputlohihi[i][j] = 0.0;
         outputlolohi[i][j] = 0.0;
         outputlohilo[i][j] = 0.0;
         outputlololo[i][j] = 0.0;
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
         odf_inc(&outputhihihi[dim][i],   &outputhilohi[dim][i],
                 &outputhihilo[dim][i],   &outputhilolo[dim][i],
                 &outputlohihi[dim][i],   &outputlolohi[dim][i],
                 &outputlohilo[dim][i],   &outputlololo[dim][i],
                 forwardhihihi[k][ix1][i],forwardhilohi[k][ix1][i],
                 forwardhihilo[k][ix1][i],forwardhilolo[k][ix1][i],
                 forwardlohihi[k][ix1][i],forwardlolohi[k][ix1][i],
                 forwardlohilo[k][ix1][i],forwardlololo[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            odf_inc(&outputhihihi[ix0][i],&outputhilohi[ix0][i],
                    &outputhihilo[ix0][i],&outputhilolo[ix0][i],
                    &outputlohihi[ix0][i],&outputlolohi[ix0][i],
                    &outputlohilo[ix0][i],&outputlololo[ix0][i],
                        cffhihihi[k][i],      cffhilohi[k][i],
                        cffhihilo[k][i],      cffhilolo[k][i],
                        cfflohihi[k][i],      cfflolohi[k][i],
                        cfflohilo[k][i],      cfflololo[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            odf_inc(&outputhihihi[ixn][i],   &outputhilohi[ixn][i],
                    &outputhihilo[ixn][i],   &outputhilolo[ixn][i],
                    &outputlohihi[ixn][i],   &outputlolohi[ixn][i],
                    &outputlohilo[ixn][i],   &outputlololo[ixn][i],
                    forwardhihihi[k][ix2][i],forwardhilohi[k][ix2][i],
                    forwardhihilo[k][ix2][i],forwardhilolo[k][ix2][i],
                    forwardlohihi[k][ix2][i],forwardlolohi[k][ix2][i],
                    forwardlohilo[k][ix2][i],forwardlololo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            odf_inc( &outputhihihi[ix0][i],    &outputhilohi[ix0][i],
                     &outputhihilo[ix0][i],    &outputhilolo[ix0][i],
                     &outputlohihi[ix0][i],    &outputlolohi[ix0][i],
                     &outputlohilo[ix0][i],    &outputlololo[ix0][i],
                    backwardhihihi[k][ix2][i],backwardhilohi[k][ix2][i],
                    backwardhihilo[k][ix2][i],backwardhilolo[k][ix2][i],
                    backwardlohihi[k][ix2][i],backwardlolohi[k][ix2][i],
                    backwardlohilo[k][ix2][i],backwardlololo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  odf_inc(&outputhihihi[ix0][i], &outputhilohi[ix0][i],
                          &outputhihilo[ix0][i], &outputhilolo[ix0][i],
                          &outputlohihi[ix0][i], &outputlolohi[ix0][i],
                          &outputlohilo[ix0][i], &outputlololo[ix0][i],
                            crosshihihi[k][j-1][i],crosshilohi[k][j-1][i],
                            crosshihilo[k][j-1][i],crosshilolo[k][j-1][i],
                            crosslohihi[k][j-1][i],crosslolohi[k][j-1][i],
                            crosslohilo[k][j-1][i],crosslololo[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_dbl8_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo, 
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo, 
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   double ***forwardhihihi, double ***forwardhilohi,
   double ***forwardhihilo, double ***forwardhilolo,
   double ***forwardlohihi, double ***forwardlolohi,
   double ***forwardlohilo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardhilohi,
   double ***backwardhihilo, double ***backwardhilolo, 
   double ***backwardlohihi, double ***backwardlolohi,
   double ***backwardlohilo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosshilohi,
   double ***crosshihilo, double ***crosshilolo,
   double ***crosslohihi, double ***crosslolohi,
   double ***crosslohilo, double ***crosslololo,
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

         CPU_dbl8_add_job(deg,
                 csthihihi,     csthilohi,     csthihilo,     csthilolo,
                 cstlohihi,     cstlolohi,     cstlohilo,     cstlololo,
                 cffhihihi,     cffhilohi,     cffhihilo,     cffhilolo,
                 cfflohihi,     cfflolohi,     cfflohilo,     cfflololo,
             forwardhihihi, forwardhilohi, forwardhihilo, forwardhilolo,
             forwardlohihi, forwardlolohi, forwardlohilo, forwardlololo,
            backwardhihihi,backwardhilohi,backwardhihilo,backwardhilolo,
            backwardlohihi,backwardlolohi,backwardlohilo,backwardlololo,
               crosshihihi,   crosshilohi,   crosshihilo,   crosshilolo,
               crosslohihi,   crosslolohi,   crosslohilo,   crosslololo,
            job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputhihihi[dim][i] = forwardhihihi[lastmon][lastidx][i];
      outputhilohi[dim][i] = forwardhilohi[lastmon][lastidx][i];
      outputhihilo[dim][i] = forwardhihilo[lastmon][lastidx][i];
      outputhilolo[dim][i] = forwardhilolo[lastmon][lastidx][i];
      outputlohihi[dim][i] = forwardlohihi[lastmon][lastidx][i];
      outputlolohi[dim][i] = forwardlolohi[lastmon][lastidx][i];
      outputlohilo[dim][i] = forwardlohilo[lastmon][lastidx][i];
      outputlololo[dim][i] = forwardlololo[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputhihihi[0][i] = 0.0;
         outputhilohi[0][i] = 0.0;
         outputhihilo[0][i] = 0.0;
         outputhilolo[0][i] = 0.0;
         outputlohihi[0][i] = 0.0;
         outputlolohi[0][i] = 0.0;
         outputlohilo[0][i] = 0.0;
         outputlololo[0][i] = 0.0;
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
         outputhihihi[0][i] = backwardhihihi[ix0][ix2][i];
         outputhilohi[0][i] = backwardhilohi[ix0][ix2][i];
         outputhihilo[0][i] = backwardhihilo[ix0][ix2][i];
         outputhilolo[0][i] = backwardhilolo[ix0][ix2][i];
         outputlohihi[0][i] = backwardlohihi[ix0][ix2][i];
         outputlolohi[0][i] = backwardlolohi[ix0][ix2][i];
         outputlohilo[0][i] = backwardlohilo[ix0][ix2][i];
         outputlololo[0][i] = backwardlololo[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputhihihi[k][i] = 0.0;
            outputhilohi[k][i] = 0.0;
            outputhihilo[k][i] = 0.0;
            outputhilolo[k][i] = 0.0;
            outputlohihi[k][i] = 0.0;
            outputlolohi[k][i] = 0.0;
            outputlohilo[k][i] = 0.0;
            outputlololo[k][i] = 0.0;
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
               outputhihihi[k][i] = backwardhihihi[ix0][ix2][i];
               outputhilohi[k][i] = backwardhilohi[ix0][ix2][i];
               outputhihilo[k][i] = backwardhihilo[ix0][ix2][i];
               outputhilolo[k][i] = backwardhilolo[ix0][ix2][i];
               outputlohihi[k][i] = backwardlohihi[ix0][ix2][i];
               outputlolohi[k][i] = backwardlolohi[ix0][ix2][i];
               outputlohilo[k][i] = backwardlohilo[ix0][ix2][i];
               outputlololo[k][i] = backwardlololo[ix0][ix2][i];
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
               outputhihihi[k][i] = forwardhihihi[ix0][ix2][i];
               outputhilohi[k][i] = forwardhilohi[ix0][ix2][i];
               outputhihilo[k][i] = forwardhihilo[ix0][ix2][i];
               outputhilolo[k][i] = forwardhilolo[ix0][ix2][i];
               outputlohihi[k][i] = forwardlohihi[ix0][ix2][i];
               outputlolohi[k][i] = forwardlolohi[ix0][ix2][i];
               outputlohilo[k][i] = forwardlohilo[ix0][ix2][i];
               outputlololo[k][i] = forwardlololo[ix0][ix2][i];
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
               outputhihihi[k][i] = crosshihihi[ix0][ix2][i];
               outputhilohi[k][i] = crosshilohi[ix0][ix2][i];
               outputhihilo[k][i] = crosshihilo[ix0][ix2][i];
               outputhilolo[k][i] = crosshilolo[ix0][ix2][i];
               outputlohihi[k][i] = crosslohihi[ix0][ix2][i];
               outputlolohi[k][i] = crosslolohi[ix0][ix2][i];
               outputlohilo[k][i] = crosslohilo[ix0][ix2][i];
               outputlololo[k][i] = crosslololo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl8_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo, 
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo, 
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose )
{
   double ***forwardhihihi = new double**[nbr];
   double ***forwardhilohi = new double**[nbr];
   double ***forwardhihilo = new double**[nbr];
   double ***forwardhilolo = new double**[nbr];
   double ***forwardlohihi = new double**[nbr];
   double ***forwardlolohi = new double**[nbr];
   double ***forwardlohilo = new double**[nbr];
   double ***forwardlololo = new double**[nbr];
   double ***backwardhihihi = new double**[nbr];
   double ***backwardhilohi = new double**[nbr];
   double ***backwardhihilo = new double**[nbr];
   double ***backwardhilolo = new double**[nbr];
   double ***backwardlohihi = new double**[nbr];
   double ***backwardlolohi = new double**[nbr];
   double ***backwardlohilo = new double**[nbr];
   double ***backwardlololo = new double**[nbr];
   double ***crosshihihi = new double**[nbr];
   double ***crosshilohi = new double**[nbr];
   double ***crosshihilo = new double**[nbr];
   double ***crosshilolo = new double**[nbr];
   double ***crosslohihi = new double**[nbr];
   double ***crosslolohi = new double**[nbr];
   double ***crosslohilo = new double**[nbr];
   double ***crosslololo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardhihihi[k] = new double*[nvrk];
      forwardhilohi[k] = new double*[nvrk];
      forwardhihilo[k] = new double*[nvrk];
      forwardhilolo[k] = new double*[nvrk];
      forwardlohihi[k] = new double*[nvrk];
      forwardlolohi[k] = new double*[nvrk];
      forwardlohilo[k] = new double*[nvrk];
      forwardlololo[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardhihihi[k][i] = new double[deg+1];
         forwardhilohi[k][i] = new double[deg+1];
         forwardhihilo[k][i] = new double[deg+1];
         forwardhilolo[k][i] = new double[deg+1];
         forwardlohihi[k][i] = new double[deg+1];
         forwardlolohi[k][i] = new double[deg+1];
         forwardlohilo[k][i] = new double[deg+1];
         forwardlololo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardhihihi[k] = new double*[nvrk-1];
         backwardhilohi[k] = new double*[nvrk-1];
         backwardhihilo[k] = new double*[nvrk-1];
         backwardhilolo[k] = new double*[nvrk-1];
         backwardlohihi[k] = new double*[nvrk-1];
         backwardlolohi[k] = new double*[nvrk-1];
         backwardlohilo[k] = new double*[nvrk-1];
         backwardlololo[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardhihihi[k][i] = new double[deg+1];
            backwardhilohi[k][i] = new double[deg+1];
            backwardhihilo[k][i] = new double[deg+1];
            backwardhilolo[k][i] = new double[deg+1];
            backwardlohihi[k][i] = new double[deg+1];
            backwardlolohi[k][i] = new double[deg+1];
            backwardlohilo[k][i] = new double[deg+1];
            backwardlololo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crosshihihi[k] = new double*[nvrk-2];
         crosshilohi[k] = new double*[nvrk-2];
         crosshihilo[k] = new double*[nvrk-2];
         crosshilolo[k] = new double*[nvrk-2];
         crosslohihi[k] = new double*[nvrk-2];
         crosslolohi[k] = new double*[nvrk-2];
         crosslohilo[k] = new double*[nvrk-2];
         crosslololo[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crosshihihi[k][i] = new double[deg+1];
            crosshilohi[k][i] = new double[deg+1];
            crosshihilo[k][i] = new double[deg+1];
            crosshilolo[k][i] = new double[deg+1];
            crosslohihi[k][i] = new double[deg+1];
            crosslolohi[k][i] = new double[deg+1];
            crosslohilo[k][i] = new double[deg+1];
            crosslololo[k][i] = new double[deg+1];
         }
      }
   }
   clock_t start = clock();
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(verbose) cout << "executing convolution jobs at layer "
                       << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = cnvjobs.get_job(k,i);
         if(verbose) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_dbl8_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffhihihi[monidx],cffhilohi[monidx],
             cffhihilo[monidx],cffhilolo[monidx],
             cfflohihi[monidx],cfflolohi[monidx],
             cfflohilo[monidx],cfflololo[monidx],
             inputhihihi,inputhilohi,inputhihilo,inputhilolo,
             inputlohihi,inputlolohi,inputlohilo,inputlololo,
              forwardhihihi[monidx], forwardhilohi[monidx],
              forwardhihilo[monidx], forwardhilolo[monidx],
              forwardlohihi[monidx], forwardlolohi[monidx],
              forwardlohilo[monidx], forwardlololo[monidx],
             backwardhihihi[monidx],backwardhilohi[monidx],
             backwardhihilo[monidx],backwardhilolo[monidx],
             backwardlohihi[monidx],backwardlolohi[monidx],
             backwardlohilo[monidx],backwardlololo[monidx],
                crosshihihi[monidx],   crosshilohi[monidx],
                crosshihilo[monidx],   crosshilolo[monidx],
                crosslohihi[monidx],   crosslolohi[monidx],
                crosslohilo[monidx],   crosslololo[monidx],job,verbose);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl8_poly_addjobs
      (dim,nbr,deg,nvr,idx,
            csthihihi,     csthilohi,     csthihilo,     csthilolo,
            cstlohihi,     cstlolohi,     cstlohilo,     cstlololo,
            cffhihihi,     cffhilohi,     cffhihilo,     cffhilolo,
            cfflohihi,     cfflolohi,     cfflohilo,     cfflololo,
          inputhihihi,   inputhilohi,   inputhihilo,   inputhilolo,
          inputlohihi,   inputlolohi,   inputlohilo,   inputlololo,
         outputhihihi,  outputhilohi,  outputhihilo,  outputhilolo,
         outputlohihi,  outputlolohi,  outputlohilo,  outputlololo,
        forwardhihihi, forwardhilohi, forwardhihilo, forwardhilolo,
        forwardlohihi, forwardlolohi, forwardlohilo, forwardlololo,
       backwardhihihi,backwardhilohi,backwardhihilo,backwardhilolo,
       backwardlohihi,backwardlolohi,backwardlohilo,backwardlololo,
          crosshihihi,   crosshilohi,   crosshihilo,   crosshilolo,
          crosslohihi,   crosslolohi,   crosslohilo,   crosslololo,
       addjobs,verbose);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(verbose)
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
         free(forwardhihihi[k][i]);
         free(forwardhilohi[k][i]);
         free(forwardhihilo[k][i]);
         free(forwardhilolo[k][i]);
         free(forwardlohihi[k][i]);
         free(forwardlolohi[k][i]);
         free(forwardlohilo[k][i]);
         free(forwardlololo[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardhihihi[k][i]);
                      free(backwardhilohi[k][i]);
                      free(backwardhihilo[k][i]);
                      free(backwardhilolo[k][i]);
                      free(backwardlohihi[k][i]);
                      free(backwardlolohi[k][i]);
                      free(backwardlohilo[k][i]);
                      free(backwardlololo[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crosshihihi[k][i]);
                      free(crosshilohi[k][i]);
                      free(crosshihilo[k][i]);
                      free(crosshilolo[k][i]);
                      free(crosslohihi[k][i]);
                      free(crosslolohi[k][i]);
                      free(crosslohilo[k][i]);
                      free(crosslololo[k][i]);
                   }
   }
   free(forwardhihihi); free(backwardhihihi); free(crosshihihi);
   free(forwardhilohi); free(backwardhilohi); free(crosshilohi);
   free(forwardhihilo); free(backwardhihilo); free(crosshihilo);
   free(forwardhilolo); free(backwardhilolo); free(crosshilolo);
   free(forwardlohihi); free(backwardlohihi); free(crosslohihi);
   free(forwardlolohi); free(backwardlolohi); free(crosslolohi);
   free(forwardlohilo); free(backwardlohilo); free(crosslohilo);
   free(forwardlololo); free(backwardlololo); free(crosslololo);
}
