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
                    forwardhilohi[0][j],   forwardlolohi[0][j],
                    forwardhihilo[0][j],   forwardlohilo[0][j],
                    forwardhilolo[0][j],   forwardlololo[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl8_product(deg, cffhihihi[i],     cfflohihi[i],
                               cffhilohi[i],     cfflolohi[i],
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
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
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
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double **forwardhihihi = new double*[dim];
   double **forwardlohihi = new double*[dim];
   double **forwardhilohi = new double*[dim];
   double **forwardlolohi = new double*[dim];
   double **forwardhihilo = new double*[dim];
   double **forwardlohilo = new double*[dim];
   double **forwardhilolo = new double*[dim];
   double **forwardlololo = new double*[dim];
   double **backwardhihihi = new double*[dim-1]; // in case dim = 2
   double **backwardlohihi = new double*[dim-1];
   double **backwardhilohi = new double*[dim-1];
   double **backwardlolohi = new double*[dim-1];
   double **backwardhihilo = new double*[dim-1]; 
   double **backwardlohilo = new double*[dim-1];
   double **backwardhilolo = new double*[dim-1];
   double **backwardlololo = new double*[dim-1];
   double **crosshihihi = new double*[dim-1];    // in case dim = 2
   double **crosslohihi = new double*[dim-1];
   double **crosshilohi = new double*[dim-1];
   double **crosslolohi = new double*[dim-1];
   double **crosshihilo = new double*[dim-1];
   double **crosslohilo = new double*[dim-1];
   double **crosshilolo = new double*[dim-1];
   double **crosslololo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardhihihi[i] = new double[deg+1];
      forwardlohihi[i] = new double[deg+1];
      forwardhilohi[i] = new double[deg+1];
      forwardlolohi[i] = new double[deg+1];
      forwardhihilo[i] = new double[deg+1];
      forwardlohilo[i] = new double[deg+1];
      forwardhilolo[i] = new double[deg+1];
      forwardlololo[i] = new double[deg+1];
      backwardhihihi[i] = new double[deg+1];
      backwardlohihi[i] = new double[deg+1];
      backwardhilohi[i] = new double[deg+1];
      backwardlolohi[i] = new double[deg+1];
      backwardhihilo[i] = new double[deg+1];
      backwardlohilo[i] = new double[deg+1];
      backwardhilolo[i] = new double[deg+1];
      backwardlololo[i] = new double[deg+1];
      crosshihihi[i] = new double[deg+1];
      crosslohihi[i] = new double[deg+1];
      crosshilohi[i] = new double[deg+1];
      crosslolohi[i] = new double[deg+1];
      crosshihilo[i] = new double[deg+1];
      crosslohilo[i] = new double[deg+1];
      crosshilolo[i] = new double[deg+1];
      crosslololo[i] = new double[deg+1];
   }
   forwardhihihi[dim-1] = new double[deg+1];
   forwardlohihi[dim-1] = new double[deg+1];
   forwardhilohi[dim-1] = new double[deg+1];
   forwardlolohi[dim-1] = new double[deg+1];
   forwardhihilo[dim-1] = new double[deg+1];
   forwardlohilo[dim-1] = new double[deg+1];
   forwardhilolo[dim-1] = new double[deg+1];
   forwardlololo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputhihihi[dim][i] = csthihihi[i];
      outputlohihi[dim][i] = cstlohihi[i];
      outputhilohi[dim][i] = csthilohi[i];
      outputlolohi[dim][i] = cstlolohi[i];
      outputhihilo[dim][i] = csthihilo[i];
      outputlohilo[dim][i] = cstlohilo[i];
      outputhilolo[dim][i] = csthilolo[i];
      outputlololo[dim][i] = cstlololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhihihi[i][j] = 0.0;
         outputlohihi[i][j] = 0.0;
         outputhilohi[i][j] = 0.0;
         outputlolohi[i][j] = 0.0;
         outputhihilo[i][j] = 0.0;
         outputlohilo[i][j] = 0.0;
         outputhilolo[i][j] = 0.0;
         outputlololo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl8_poly_speel
      (dim,nbr,deg,nvr,idx,
            cffhihihi,     cfflohihi,     cffhilohi,     cfflolohi,
            cffhihilo,     cfflohilo,     cffhilolo,     cfflololo,
          inputhihihi,   inputlohihi,   inputhilohi,   inputlolohi,
          inputhihilo,   inputlohilo,   inputhilolo,   inputlololo,
         outputhihihi,  outputlohihi,  outputhilohi,  outputlolohi,
         outputhihilo,  outputlohilo,  outputhilolo,  outputlololo,
        forwardhihihi, forwardlohihi, forwardhilohi, forwardlolohi,
        forwardhihilo, forwardlohilo, forwardhilolo, forwardlololo,
       backwardhihihi,backwardlohihi,backwardhilohi,backwardlolohi,
       backwardhihilo,backwardlohilo,backwardhilolo,backwardlololo,
          crosshihihi,   crosslohihi,   crosshilohi,   crosslolohi,
          crosshihilo,   crosslohilo,   crosshilolo,   crosslololo,vrb);
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
      free(forwardhihihi[i]); free(backwardhihihi[i]); free(crosshihihi[i]);
      free(forwardlohihi[i]); free(backwardlohihi[i]); free(crosslohihi[i]);
      free(forwardhilohi[i]); free(backwardhilohi[i]); free(crosshilohi[i]);
      free(forwardlolohi[i]); free(backwardlolohi[i]); free(crosslolohi[i]);
      free(forwardhihilo[i]); free(backwardhihilo[i]); free(crosshihilo[i]);
      free(forwardlohilo[i]); free(backwardlohilo[i]); free(crosslohilo[i]);
      free(forwardhilolo[i]); free(backwardhilolo[i]); free(crosshilolo[i]);
      free(forwardlololo[i]); free(backwardlololo[i]); free(crosslololo[i]);
   }
   free(forwardhihihi[dim-1]); free(forwardlohihi[dim-1]);
   free(forwardhilohi[dim-1]); free(forwardlolohi[dim-1]);
   free(forwardhihilo[dim-1]); free(forwardlohilo[dim-1]);
   free(forwardhilolo[dim-1]); free(forwardlololo[dim-1]);
   free(forwardhihihi); free(backwardhihihi); free(crosshihihi);
   free(forwardlohihi); free(backwardlohihi); free(crosslohihi);
   free(forwardhilohi); free(backwardhilohi); free(crosshilohi);
   free(forwardlolohi); free(backwardlolohi); free(crosslolohi);
   free(forwardhihilo); free(backwardhihilo); free(crosshihilo);
   free(forwardlohilo); free(backwardlohilo); free(crosslohilo);
   free(forwardhilolo); free(backwardhilolo); free(crosshilolo);
   free(forwardlololo); free(backwardlololo); free(crosslololo);
}

void CPU_cmplx8_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
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
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double **forwardrehihihi = new double*[dim];
   double **forwardrelohihi = new double*[dim];
   double **forwardrehilohi = new double*[dim];
   double **forwardrelolohi = new double*[dim];
   double **forwardrehihilo = new double*[dim];
   double **forwardrelohilo = new double*[dim];
   double **forwardrehilolo = new double*[dim];
   double **forwardrelololo = new double*[dim];
   double **forwardimhihihi = new double*[dim];
   double **forwardimlohihi = new double*[dim];
   double **forwardimhilohi = new double*[dim];
   double **forwardimlolohi = new double*[dim];
   double **forwardimhihilo = new double*[dim];
   double **forwardimlohilo = new double*[dim];
   double **forwardimhilolo = new double*[dim];
   double **forwardimlololo = new double*[dim];
   double **backwardrehihihi = new double*[dim-1]; // in case dim = 2
   double **backwardrelohihi = new double*[dim-1];
   double **backwardrehilohi = new double*[dim-1];
   double **backwardrelolohi = new double*[dim-1];
   double **backwardrehihilo = new double*[dim-1];
   double **backwardrelohilo = new double*[dim-1];
   double **backwardrehilolo = new double*[dim-1];
   double **backwardrelololo = new double*[dim-1];
   double **backwardimhihihi = new double*[dim-1];
   double **backwardimlohihi = new double*[dim-1];
   double **backwardimhilohi = new double*[dim-1];
   double **backwardimlolohi = new double*[dim-1];
   double **backwardimhihilo = new double*[dim-1];
   double **backwardimlohilo = new double*[dim-1];
   double **backwardimhilolo = new double*[dim-1];
   double **backwardimlololo = new double*[dim-1];
   double **crossrehihihi = new double*[dim-1];    // in case dim = 2
   double **crossrelohihi = new double*[dim-1];
   double **crossrehilohi = new double*[dim-1];
   double **crossrelolohi = new double*[dim-1];
   double **crossrehihilo = new double*[dim-1];
   double **crossrelohilo = new double*[dim-1];
   double **crossrehilolo = new double*[dim-1];
   double **crossrelololo = new double*[dim-1];
   double **crossimhihihi = new double*[dim-1];
   double **crossimlohihi = new double*[dim-1];
   double **crossimhilohi = new double*[dim-1];
   double **crossimlolohi = new double*[dim-1];
   double **crossimhihilo = new double*[dim-1];
   double **crossimlohilo = new double*[dim-1];
   double **crossimhilolo = new double*[dim-1];
   double **crossimlololo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardrehihihi[i] = new double[deg+1];
      forwardrelohihi[i] = new double[deg+1];
      forwardrehilohi[i] = new double[deg+1];
      forwardrelolohi[i] = new double[deg+1];
      forwardrehihilo[i] = new double[deg+1];
      forwardrelohilo[i] = new double[deg+1];
      forwardrehilolo[i] = new double[deg+1];
      forwardrelololo[i] = new double[deg+1];
      forwardimhihihi[i] = new double[deg+1];
      forwardimlohihi[i] = new double[deg+1];
      forwardimhilohi[i] = new double[deg+1];
      forwardimlolohi[i] = new double[deg+1];
      forwardimhihilo[i] = new double[deg+1];
      forwardimlohilo[i] = new double[deg+1];
      forwardimhilolo[i] = new double[deg+1];
      forwardimlololo[i] = new double[deg+1];
      backwardrehihihi[i] = new double[deg+1];
      backwardrelohihi[i] = new double[deg+1];
      backwardrehilohi[i] = new double[deg+1];
      backwardrelolohi[i] = new double[deg+1];
      backwardrehihilo[i] = new double[deg+1];
      backwardrelohilo[i] = new double[deg+1];
      backwardrehilolo[i] = new double[deg+1];
      backwardrelololo[i] = new double[deg+1];
      backwardimhihihi[i] = new double[deg+1];
      backwardimlohihi[i] = new double[deg+1];
      backwardimhilohi[i] = new double[deg+1];
      backwardimlolohi[i] = new double[deg+1];
      backwardimhihilo[i] = new double[deg+1];
      backwardimlohilo[i] = new double[deg+1];
      backwardimhilolo[i] = new double[deg+1];
      backwardimlololo[i] = new double[deg+1];
      crossrehihihi[i] = new double[deg+1];
      crossrelohihi[i] = new double[deg+1];
      crossrehilohi[i] = new double[deg+1];
      crossrelolohi[i] = new double[deg+1];
      crossrehihilo[i] = new double[deg+1];
      crossrelohilo[i] = new double[deg+1];
      crossrehilolo[i] = new double[deg+1];
      crossrelololo[i] = new double[deg+1];
      crossimhihihi[i] = new double[deg+1];
      crossimlohihi[i] = new double[deg+1];
      crossimhilohi[i] = new double[deg+1];
      crossimlolohi[i] = new double[deg+1];
      crossimhihilo[i] = new double[deg+1];
      crossimlohilo[i] = new double[deg+1];
      crossimhilolo[i] = new double[deg+1];
      crossimlololo[i] = new double[deg+1];
   }
   forwardrehihihi[dim-1] = new double[deg+1];
   forwardrelohihi[dim-1] = new double[deg+1];
   forwardrehilohi[dim-1] = new double[deg+1];
   forwardrelolohi[dim-1] = new double[deg+1];
   forwardrehihilo[dim-1] = new double[deg+1];
   forwardrelohilo[dim-1] = new double[deg+1];
   forwardrehilolo[dim-1] = new double[deg+1];
   forwardrelololo[dim-1] = new double[deg+1];
   forwardimhihihi[dim-1] = new double[deg+1];
   forwardimlohihi[dim-1] = new double[deg+1];
   forwardimhilohi[dim-1] = new double[deg+1];
   forwardimlolohi[dim-1] = new double[deg+1];
   forwardimhihilo[dim-1] = new double[deg+1];
   forwardimlohilo[dim-1] = new double[deg+1];
   forwardimhilolo[dim-1] = new double[deg+1];
   forwardimlololo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputrehihihi[dim][i] = cstrehihihi[i];
      outputrelohihi[dim][i] = cstrelohihi[i];
      outputrehilohi[dim][i] = cstrehilohi[i];
      outputrelolohi[dim][i] = cstrelolohi[i];
      outputrehihilo[dim][i] = cstrehihilo[i];
      outputrelohilo[dim][i] = cstrelohilo[i];
      outputrehilolo[dim][i] = cstrehilolo[i];
      outputrelololo[dim][i] = cstrelololo[i];
      outputimhihihi[dim][i] = cstimhihihi[i];
      outputimlohihi[dim][i] = cstimlohihi[i];
      outputimhilohi[dim][i] = cstimhilohi[i];
      outputimlolohi[dim][i] = cstimlolohi[i];
      outputimhihilo[dim][i] = cstimhihilo[i];
      outputimlohilo[dim][i] = cstimlohilo[i];
      outputimhilolo[dim][i] = cstimhilolo[i];
      outputimlololo[dim][i] = cstimlololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrehihihi[i][j] = 0.0;
         outputrelohihi[i][j] = 0.0;
         outputrehilohi[i][j] = 0.0;
         outputrelolohi[i][j] = 0.0;
         outputrehihilo[i][j] = 0.0;
         outputrelohilo[i][j] = 0.0;
         outputrehilolo[i][j] = 0.0;
         outputrelololo[i][j] = 0.0;
         outputimhihihi[i][j] = 0.0;
         outputimlohihi[i][j] = 0.0;
         outputimhilohi[i][j] = 0.0;
         outputimlolohi[i][j] = 0.0;
         outputimhihilo[i][j] = 0.0;
         outputimlohilo[i][j] = 0.0;
         outputimhilolo[i][j] = 0.0;
         outputimlololo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_cmplx8_poly_speel
      (dim,nbr,deg,nvr,idx,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
       inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
       inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
       inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
       outputrehihihi,outputrelohihi,outputrehilohi,outputrelolohi,
       outputrehihilo,outputrelohilo,outputrehilolo,outputrelololo,
       outputimhihihi,outputimlohihi,outputimhilohi,outputimlolohi,
       outputimhihilo,outputimlohilo,outputimhilolo,outputimlololo,
       forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
       forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
       forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
       forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
       backwardrehihihi,backwardrelohihi,backwardrehilohi,backwardrelolohi,
       backwardrehihilo,backwardrelohilo,backwardrehilolo,backwardrelololo,
       backwardimhihihi,backwardimlohihi,backwardimhilohi,backwardimlolohi,
       backwardimhihilo,backwardimlohilo,backwardimhilolo,backwardimlololo,
       crossrehihihi,crossrelohihi,crossrehilohi,crossrelolohi,
       crossrehihilo,crossrelohilo,crossrehilolo,crossrelololo,
       crossimhihihi,crossimlohihi,crossimhilohi,crossimlolohi,
       crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo,vrb);
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
      free(forwardrehihihi[i]); free(forwardrelohihi[i]);
      free(forwardrehilohi[i]); free(forwardrelolohi[i]);
      free(forwardrehihilo[i]); free(forwardrelohilo[i]);
      free(forwardrehilolo[i]); free(forwardrelololo[i]);
      free(forwardimhihihi[i]); free(forwardimlohihi[i]);
      free(forwardimhilohi[i]); free(forwardimlolohi[i]);
      free(forwardimhihilo[i]); free(forwardimlohilo[i]);
      free(forwardimhilolo[i]); free(forwardimlololo[i]);
      free(backwardrehihihi[i]); free(backwardrelohihi[i]);
      free(backwardrehilohi[i]); free(backwardrelolohi[i]);
      free(backwardrehihilo[i]); free(backwardrelohilo[i]);
      free(backwardrehilolo[i]); free(backwardrelololo[i]);
      free(backwardimhihihi[i]); free(backwardimlohihi[i]);
      free(backwardimhilohi[i]); free(backwardimlolohi[i]);
      free(backwardimhihilo[i]); free(backwardimlohilo[i]);
      free(backwardimhilolo[i]); free(backwardimlololo[i]);
      free(crossrehihihi[i]); free(crossrelohihi[i]);
      free(crossrehilohi[i]); free(crossrelolohi[i]);
      free(crossrehihilo[i]); free(crossrelohilo[i]);
      free(crossrehilolo[i]); free(crossrelololo[i]);
      free(crossimhihihi[i]); free(crossimlohihi[i]);
      free(crossimhilohi[i]); free(crossimlolohi[i]);
      free(crossimhihilo[i]); free(crossimlohilo[i]);
      free(crossimhilolo[i]); free(crossimlololo[i]);
   }
   free(forwardrehihihi[dim-1]); free(forwardrelohihi[dim-1]);
   free(forwardrehilohi[dim-1]); free(forwardrelolohi[dim-1]);
   free(forwardrehihilo[dim-1]); free(forwardrelohilo[dim-1]);
   free(forwardrehilolo[dim-1]); free(forwardrelololo[dim-1]);
   free(forwardimhihihi[dim-1]); free(forwardimlohihi[dim-1]);
   free(forwardimhilohi[dim-1]); free(forwardimlolohi[dim-1]);
   free(forwardimhihilo[dim-1]); free(forwardimlohilo[dim-1]);
   free(forwardimhilolo[dim-1]); free(forwardimlololo[dim-1]);
   free(forwardrehihihi); free(forwardrelohihi);
   free(forwardrehilohi); free(forwardrelolohi);
   free(forwardrehihilo); free(forwardrelohilo);
   free(forwardrehilolo); free(forwardrelololo);
   free(forwardimhihihi); free(forwardimlohihi);
   free(forwardimhilohi); free(forwardimlolohi);
   free(forwardimhihilo); free(forwardimlohilo);
   free(forwardimhilolo); free(forwardimlololo);
   free(backwardrehihihi); free(backwardrelohihi);
   free(backwardrehilohi); free(backwardrelolohi);
   free(backwardrehihilo); free(backwardrelohilo);
   free(backwardrehilolo); free(backwardrelololo);
   free(backwardimhihihi); free(backwardimlohihi);
   free(backwardimhilohi); free(backwardimlolohi);
   free(backwardimhihilo); free(backwardimlohilo);
   free(backwardimhilolo); free(backwardimlololo);
   free(crossrehihihi); free(crossrelohihi);
   free(crossrehilohi); free(crossrelolohi);
   free(crossrehihilo); free(crossrelohilo);
   free(crossrehilolo); free(crossrelololo);
   free(crossimhihihi); free(crossimlohihi);
   free(crossimhilohi); free(crossimlolohi);
   free(crossimhihilo); free(crossimlohilo);
   free(crossimhilolo); free(crossimlololo);
}

void CPU_dbl8_conv_job
 ( int deg, int nvr, int *idx,

   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
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
   double **crosshilolo, double **crosslololo,
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
                cffhihihi,cfflohihi,cffhilohi,cfflolohi,
                cffhihilo,cfflohilo,cffhilolo,cfflololo,
              inputhihihi[inp2ix],  inputlohihi[inp2ix],
              inputhilohi[inp2ix],  inputlolohi[inp2ix],
              inputhihilo[inp2ix],  inputlohilo[inp2ix],
              inputhilolo[inp2ix],  inputlololo[inp2ix],
            forwardhihihi[outidx],forwardlohihi[outidx],
            forwardhilohi[outidx],forwardlolohi[outidx],
            forwardhihilo[outidx],forwardlohilo[outidx],
            forwardhilolo[outidx],forwardlololo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl8_product(deg,
               inputhihihi[inp1ix],inputlohihi[inp1ix],
               inputhilohi[inp1ix],inputlolohi[inp1ix],
               inputhihilo[inp1ix],inputlohilo[inp1ix],
               inputhilolo[inp1ix],inputlololo[inp1ix],
               cffhihihi,cfflohihi,cffhilohi,cfflolohi,
               cffhihilo,cfflohilo,cffhilolo,cfflololo,
               forwardhihihi[outidx],forwardlohihi[outidx],
               forwardhilohi[outidx],forwardlolohi[outidx],
               forwardhihilo[outidx],forwardlohilo[outidx],
               forwardhilolo[outidx],forwardlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
                 inputhihihi[inp1ix],  inputlohihi[inp1ix],
                 inputhilohi[inp1ix],  inputlolohi[inp1ix],
                 inputhihilo[inp1ix],  inputlohilo[inp1ix],
                 inputhilolo[inp1ix],  inputlololo[inp1ix],
               forwardhihihi[inp2ix],forwardlohihi[inp2ix],
               forwardhilohi[inp2ix],forwardlolohi[inp2ix],
               forwardhihilo[inp2ix],forwardlohilo[inp2ix],
               forwardhilolo[inp2ix],forwardlololo[inp2ix],
               forwardhihihi[outidx],forwardlohihi[outidx],
               forwardhilohi[outidx],forwardlolohi[outidx],
               forwardhihilo[outidx],forwardlohilo[outidx],
               forwardhilolo[outidx],forwardlololo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
              crosshihihi[inp1ix],  crosslohihi[inp1ix],
              crosshilohi[inp1ix],  crosslolohi[inp1ix],
              crosshihilo[inp1ix],  crosslohilo[inp1ix],
              crosshilolo[inp1ix],  crosslololo[inp1ix],
              inputhihihi[inp2ix],  inputlohihi[inp2ix],
              inputhilohi[inp2ix],  inputlolohi[inp2ix],
              inputhihilo[inp2ix],  inputlohilo[inp2ix],
              inputhilolo[inp2ix],  inputlololo[inp2ix],
            forwardhihihi[outidx],forwardlohihi[outidx],
            forwardhilohi[outidx],forwardlolohi[outidx],
            forwardhihilo[outidx],forwardlohilo[outidx],
            forwardhilolo[outidx],forwardlololo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl8_product(deg,
               inputhihihi[inp1ix],inputlohihi[inp1ix],
               inputhilohi[inp1ix],inputlolohi[inp1ix],
               inputhihilo[inp1ix],inputlohilo[inp1ix],
               inputhilolo[inp1ix],inputlololo[inp1ix],
               cffhihihi,cfflohihi,cffhilohi,cfflolohi,
               cffhihilo,cfflohilo,cffhilolo,cfflololo,
               forwardhihihi[outidx],forwardlohihi[outidx],
               forwardhilohi[outidx],forwardlolohi[outidx],
               forwardhihilo[outidx],forwardlohilo[outidx],
               forwardhilolo[outidx],forwardlololo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
               forwardhihihi[inp1ix],forwardlohihi[inp1ix],
               forwardhilohi[inp1ix],forwardlolohi[inp1ix],
               forwardhihilo[inp1ix],forwardlohilo[inp1ix],
               forwardhilolo[inp1ix],forwardlololo[inp1ix],
                 inputhihihi[inp2ix],  inputlohihi[inp2ix],
                 inputhilohi[inp2ix],  inputlolohi[inp2ix],
                 inputhihilo[inp2ix],  inputlohilo[inp2ix],
                 inputhilolo[inp2ix],  inputlololo[inp2ix],
               forwardhihihi[outidx],forwardlohihi[outidx],
               forwardhilohi[outidx],forwardlolohi[outidx],
               forwardhihilo[outidx],forwardlohilo[outidx],
               forwardhilolo[outidx],forwardlololo[outidx]);
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
                    cffhihihi,cfflohihi,cffhilohi,cfflolohi,
                    cffhihilo,cfflohilo,cffhilolo,cfflololo,
                  inputhihihi[inp2ix],   inputlohihi[inp2ix],
                  inputhilohi[inp2ix],   inputlolohi[inp2ix],
                  inputhihilo[inp2ix],   inputlohilo[inp2ix],
                  inputhilolo[inp2ix],   inputlololo[inp2ix],
               backwardhihihi[outidx],backwardlohihi[outidx],
               backwardhilohi[outidx],backwardlolohi[outidx],
               backwardhihilo[outidx],backwardlohilo[outidx],
               backwardhilolo[outidx],backwardlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
               cffhihihi,cfflohihi,cffhilohi,cfflolohi,
               cffhihilo,cfflohilo,cffhilolo,cfflololo,
               backwardhihihi[inp2ix],backwardlohihi[inp2ix],
               backwardhilohi[inp2ix],backwardlolohi[inp2ix],
               backwardhihilo[inp2ix],backwardlohilo[inp2ix],
               backwardhilolo[inp2ix],backwardlololo[inp2ix],
               backwardhihihi[outidx],backwardlohihi[outidx],
               backwardhilohi[outidx],backwardlolohi[outidx],
               backwardhihilo[outidx],backwardlohilo[outidx],
               backwardhilolo[outidx],backwardlololo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl8_product(deg,
                  inputhihihi[inp1ix],   inputlohihi[inp1ix],
                  inputhilohi[inp1ix],   inputlolohi[inp1ix],
                  inputhihilo[inp1ix],   inputlohilo[inp1ix],
                  inputhilolo[inp1ix],   inputlololo[inp1ix],
                  inputhihihi[inp2ix],   inputlohihi[inp2ix],
                  inputhilohi[inp2ix],   inputlolohi[inp2ix],
                  inputhihilo[inp2ix],   inputlohilo[inp2ix],
                  inputhilolo[inp2ix],   inputlololo[inp2ix],
               backwardhihihi[outidx],backwardlohihi[outidx],
               backwardhilohi[outidx],backwardlolohi[outidx],
               backwardhihilo[outidx],backwardlohilo[outidx],
               backwardhilolo[outidx],backwardlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
                  inputhihihi[inp1ix],   inputlohihi[inp1ix],
                  inputhilohi[inp1ix],   inputlolohi[inp1ix],
                  inputhihilo[inp1ix],   inputlohilo[inp1ix],
                  inputhilolo[inp1ix],   inputlololo[inp1ix],
               backwardhihihi[inp2ix],backwardlohihi[inp2ix],
               backwardhilohi[inp2ix],backwardlolohi[inp2ix],
               backwardhihilo[inp2ix],backwardlohilo[inp2ix],
               backwardhilolo[inp2ix],backwardlololo[inp2ix],
               backwardhihihi[outidx],backwardlohihi[outidx],
               backwardhilohi[outidx],backwardlolohi[outidx],
               backwardhihilo[outidx],backwardlohilo[outidx],
               backwardhilolo[outidx],backwardlololo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl8_product(deg,
               backwardhihihi[inp1ix],backwardlohihi[inp1ix],
               backwardhilohi[inp1ix],backwardlolohi[inp1ix],
               backwardhihilo[inp1ix],backwardlohilo[inp1ix],
               backwardhilolo[inp1ix],backwardlololo[inp1ix],
               cffhihihi,cfflohihi,cffhilohi,cfflolohi,
               cffhihilo,cfflohilo,cffhilolo,cfflololo,
               backwardhihihi[outidx],backwardlohihi[outidx],
               backwardhilohi[outidx],backwardlolohi[outidx],
               backwardhihilo[outidx],backwardlohilo[outidx],
               backwardhilolo[outidx],backwardlololo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl8_product(deg,
               backwardhihihi[inp1ix],backwardlohihi[inp1ix],
               backwardhilohi[inp1ix],backwardlolohi[inp1ix],
               backwardhihilo[inp1ix],backwardlohilo[inp1ix],
               backwardhilolo[inp1ix],backwardlololo[inp1ix],
                  inputhihihi[inp2ix],   inputlohihi[inp2ix],
                  inputhilohi[inp2ix],   inputlolohi[inp2ix],
                  inputhihilo[inp2ix],   inputlohilo[inp2ix],
                  inputhilolo[inp2ix],   inputlololo[inp2ix],
               backwardhihihi[outidx],backwardlohihi[outidx],
               backwardhilohi[outidx],backwardlolohi[outidx],
               backwardhihilo[outidx],backwardlohilo[outidx],
               backwardhilolo[outidx],backwardlololo[outidx]);
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
            cffhihihi,cfflohihi,cffhilohi,cfflolohi,
            cffhihilo,cfflohilo,cffhilolo,cfflololo,
            inputhihihi[inp2ix],inputlohihi[inp2ix],
            inputhilohi[inp2ix],inputlolohi[inp2ix],
            inputhihilo[inp2ix],inputlohilo[inp2ix],
            inputhilolo[inp2ix],inputlololo[inp2ix],
            crosshihihi[outidx],crosslohihi[outidx],
            crosshilohi[outidx],crosslolohi[outidx],
            crosshihilo[outidx],crosslohilo[outidx],
            crosshilolo[outidx],crosslololo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
              inputhihihi[inp1ix],  inputlohihi[inp1ix],
              inputhilohi[inp1ix],  inputlolohi[inp1ix],
              inputhihilo[inp1ix],  inputlohilo[inp1ix],
              inputhilolo[inp1ix],  inputlololo[inp1ix],
            forwardhihihi[inp2ix],forwardlohihi[inp2ix],
            forwardhilohi[inp2ix],forwardlolohi[inp2ix],
            forwardhihilo[inp2ix],forwardlohilo[inp2ix],
            forwardhilolo[inp2ix],forwardlololo[inp2ix],
              crosshihihi[outidx],  crosslohihi[outidx],
              crosshilohi[outidx],  crosslolohi[outidx],
              crosshihilo[outidx],  crosslohilo[outidx],
              crosshilolo[outidx],  crosslololo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl8_product(deg,
              forwardhihihi[inp1ix],forwardlohihi[inp1ix],
              forwardhilohi[inp1ix],forwardlolohi[inp1ix],
              forwardhihilo[inp1ix],forwardlohilo[inp1ix],
              forwardhilolo[inp1ix],forwardlololo[inp1ix],
                inputhihihi[inp2ix],  inputlohihi[inp2ix],
                inputhilohi[inp2ix],  inputlolohi[inp2ix],
                inputhihilo[inp2ix],  inputlohilo[inp2ix],
                inputhilolo[inp2ix],  inputlololo[inp2ix],
                crosshihihi[outidx],  crosslohihi[outidx],
                crosshilohi[outidx],  crosslolohi[outidx],
                crosshihilo[outidx],  crosslohilo[outidx],
                crosshilolo[outidx],  crosslololo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl8_product(deg,
               forwardhihihi[inp1ix], forwardlohihi[inp1ix],
               forwardhilohi[inp1ix], forwardlolohi[inp1ix],
               forwardhihilo[inp1ix], forwardlohilo[inp1ix],
               forwardhilolo[inp1ix], forwardlololo[inp1ix],
              backwardhihihi[inp2ix],backwardlohihi[inp2ix],
              backwardhilohi[inp2ix],backwardlolohi[inp2ix],
              backwardhihilo[inp2ix],backwardlohilo[inp2ix],
              backwardhilolo[inp2ix],backwardlololo[inp2ix],
                 crosshihihi[outidx],   crosslohihi[outidx],
                 crosshilohi[outidx],   crosslolohi[outidx],
                 crosshihilo[outidx],   crosslohilo[outidx],
                 crosshilolo[outidx],   crosslololo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl8_product(deg,
            backwardhihihi[inp1ix],backwardlohihi[inp1ix],
            backwardhilohi[inp1ix],backwardlolohi[inp1ix],
            backwardhihilo[inp1ix],backwardlohilo[inp1ix],
            backwardhilolo[inp1ix],backwardlololo[inp1ix],
             forwardhihihi[inp2ix], forwardlohihi[inp2ix],
             forwardhilohi[inp2ix], forwardlolohi[inp2ix],
             forwardhihilo[inp2ix], forwardlohilo[inp2ix],
             forwardhilolo[inp2ix], forwardlololo[inp2ix],
               crosshihihi[outidx],   crosslohihi[outidx],
               crosshilohi[outidx],   crosslolohi[outidx],
               crosshihilo[outidx],   crosslohilo[outidx],
               crosshilolo[outidx],   crosslololo[outidx]);
      }
   }
}

void CPU_cmplx8_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
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
         CPU_cmplx8_product(deg,
            cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
            cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
            cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
            cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
            inputrehihihi[inp2ix],inputrelohihi[inp2ix],
            inputrehilohi[inp2ix],inputrelolohi[inp2ix],
            inputrehihilo[inp2ix],inputrelohilo[inp2ix],
            inputrehilolo[inp2ix],inputrelololo[inp2ix],
            inputimhihihi[inp2ix],inputimlohihi[inp2ix],
            inputimhilohi[inp2ix],inputimlolohi[inp2ix],
            inputimhihilo[inp2ix],inputimlohilo[inp2ix],
            inputimhilolo[inp2ix],inputimlololo[inp2ix],
            forwardrehihihi[outidx],forwardrelohihi[outidx],
            forwardrehilohi[outidx],forwardrelolohi[outidx],
            forwardrehihilo[outidx],forwardrelohilo[outidx],
            forwardrehilolo[outidx],forwardrelololo[outidx],
            forwardimhihihi[outidx],forwardimlohihi[outidx],
            forwardimhilohi[outidx],forwardimlolohi[outidx],
            forwardimhihilo[outidx],forwardimlohilo[outidx],
            forwardimhilolo[outidx],forwardimlololo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx8_product(deg,
               inputrehihihi[inp1ix],inputrelohihi[inp1ix],
               inputrehilohi[inp1ix],inputrelolohi[inp1ix],
               inputrehihilo[inp1ix],inputrelohilo[inp1ix],
               inputrehilolo[inp1ix],inputrelololo[inp1ix],
               inputimhihihi[inp1ix],inputimlohihi[inp1ix],
               inputimhilohi[inp1ix],inputimlolohi[inp1ix],
               inputimhihilo[inp1ix],inputimlohilo[inp1ix],
               inputimhilolo[inp1ix],inputimlololo[inp1ix],
               cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
               forwardrehihihi[outidx],forwardrelohihi[outidx],
               forwardrehilohi[outidx],forwardrelolohi[outidx],
               forwardrehihilo[outidx],forwardrelohilo[outidx],
               forwardrehilolo[outidx],forwardrelololo[outidx],
               forwardimhihihi[outidx],forwardimlohihi[outidx],
               forwardimhilohi[outidx],forwardimlolohi[outidx],
               forwardimhihilo[outidx],forwardimlohilo[outidx],
               forwardimhilolo[outidx],forwardimlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_cmplx8_product(deg,
               inputrehihihi[inp1ix],inputrelohihi[inp1ix],
               inputrehilohi[inp1ix],inputrelolohi[inp1ix],
               inputrehihilo[inp1ix],inputrelohilo[inp1ix],
               inputrehilolo[inp1ix],inputrelololo[inp1ix],
               inputimhihihi[inp1ix],inputimlohihi[inp1ix],
               inputimhilohi[inp1ix],inputimlolohi[inp1ix],
               inputimhihilo[inp1ix],inputimlohilo[inp1ix],
               inputimhilolo[inp1ix],inputimlololo[inp1ix],
               forwardrehihihi[inp2ix],forwardrelohihi[inp2ix],
               forwardrehilohi[inp2ix],forwardrelolohi[inp2ix],
               forwardrehihilo[inp2ix],forwardrelohilo[inp2ix],
               forwardrehilolo[inp2ix],forwardrelololo[inp2ix],
               forwardimhihihi[inp2ix],forwardimlohihi[inp2ix],
               forwardimhilohi[inp2ix],forwardimlolohi[inp2ix],
               forwardimhihilo[inp2ix],forwardimlohilo[inp2ix],
               forwardimhilolo[inp2ix],forwardimlololo[inp2ix],
               forwardrehihihi[outidx],forwardrelohihi[outidx],
               forwardrehilohi[outidx],forwardrelolohi[outidx],
               forwardrehihilo[outidx],forwardrelohilo[outidx],
               forwardrehilolo[outidx],forwardrelololo[outidx],
               forwardimhihihi[outidx],forwardimlohihi[outidx],
               forwardimhilohi[outidx],forwardimlolohi[outidx],
               forwardimhihilo[outidx],forwardimlohilo[outidx],
               forwardimhilolo[outidx],forwardimlololo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_cmplx8_product(deg,
            crossrehihihi[inp1ix],crossrelohihi[inp1ix],
            crossrehilohi[inp1ix],crossrelolohi[inp1ix],
            crossrehihilo[inp1ix],crossrelohilo[inp1ix],
            crossrehilolo[inp1ix],crossrelololo[inp1ix],
            crossimhihihi[inp1ix],crossimlohihi[inp1ix],
            crossimhilohi[inp1ix],crossimlolohi[inp1ix],
            crossimhihilo[inp1ix],crossimlohilo[inp1ix],
            crossimhilolo[inp1ix],crossimlololo[inp1ix],
            inputrehihihi[inp2ix],inputrelohihi[inp2ix],
            inputrehilohi[inp2ix],inputrelolohi[inp2ix],
            inputrehihilo[inp2ix],inputrelohilo[inp2ix],
            inputrehilolo[inp2ix],inputrelololo[inp2ix],
            inputimhihihi[inp2ix],inputimlohihi[inp2ix],
            inputimhilohi[inp2ix],inputimlolohi[inp2ix],
            inputimhihilo[inp2ix],inputimlohilo[inp2ix],
            inputimhilolo[inp2ix],inputimlololo[inp2ix],
            forwardrehihihi[outidx],forwardrelohihi[outidx],
            forwardrehilohi[outidx],forwardrelolohi[outidx],
            forwardrehihilo[outidx],forwardrelohilo[outidx],
            forwardrehilolo[outidx],forwardrelololo[outidx],
            forwardimhihihi[outidx],forwardimlohihi[outidx],
            forwardimhilohi[outidx],forwardimlolohi[outidx],
            forwardimhihilo[outidx],forwardimlohilo[outidx],
            forwardimhilolo[outidx],forwardimlololo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx8_product(deg,
               inputrehihihi[inp1ix],inputrelohihi[inp1ix],
               inputrehilohi[inp1ix],inputrelolohi[inp1ix],
               inputrehihilo[inp1ix],inputrelohilo[inp1ix],
               inputrehilolo[inp1ix],inputrelololo[inp1ix],
               inputimhihihi[inp1ix],inputimlohihi[inp1ix],
               inputimhilohi[inp1ix],inputimlolohi[inp1ix],
               inputimhihilo[inp1ix],inputimlohilo[inp1ix],
               inputimhilolo[inp1ix],inputimlololo[inp1ix],
               cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
               cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
               forwardrehihihi[outidx],forwardrelohihi[outidx],
               forwardrehilohi[outidx],forwardrelolohi[outidx],
               forwardrehihilo[outidx],forwardrelohilo[outidx],
               forwardrehilolo[outidx],forwardrelololo[outidx],
               forwardimhihihi[outidx],forwardimlohihi[outidx],
               forwardimhilohi[outidx],forwardimlolohi[outidx],
               forwardimhihilo[outidx],forwardimlohilo[outidx],
               forwardimhilolo[outidx],forwardimlololo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx8_product(deg,
               forwardrehihihi[inp1ix],forwardrelohihi[inp1ix],
               forwardrehilohi[inp1ix],forwardrelolohi[inp1ix],
               forwardrehihilo[inp1ix],forwardrelohilo[inp1ix],
               forwardrehilolo[inp1ix],forwardrelololo[inp1ix],
               forwardimhihihi[inp1ix],forwardimlohihi[inp1ix],
               forwardimhilohi[inp1ix],forwardimlolohi[inp1ix],
               forwardimhihilo[inp1ix],forwardimlohilo[inp1ix],
               forwardimhilolo[inp1ix],forwardimlololo[inp1ix],
               inputrehihihi[inp2ix],inputrelohihi[inp2ix],
               inputrehilohi[inp2ix],inputrelolohi[inp2ix],
               inputrehihilo[inp2ix],inputrelohilo[inp2ix],
               inputrehilolo[inp2ix],inputrelololo[inp2ix],
               inputimhihihi[inp2ix],inputimlohihi[inp2ix],
               inputimhilohi[inp2ix],inputimlolohi[inp2ix],
               inputimhihilo[inp2ix],inputimlohilo[inp2ix],
               inputimhilolo[inp2ix],inputimlololo[inp2ix],
               forwardrehihihi[outidx],forwardrelohihi[outidx],
               forwardrehilohi[outidx],forwardrelolohi[outidx],
               forwardrehihilo[outidx],forwardrelohilo[outidx],
               forwardrehilolo[outidx],forwardrelololo[outidx],
               forwardimhihihi[outidx],forwardimlohihi[outidx],
               forwardimhilohi[outidx],forwardimlolohi[outidx],
               forwardimhihilo[outidx],forwardimlohilo[outidx],
               forwardimhilolo[outidx],forwardimlololo[outidx]);
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
            CPU_cmplx8_product(deg,
               cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
               cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
               inputrehihihi[inp2ix],inputrelohihi[inp2ix],
               inputrehilohi[inp2ix],inputrelolohi[inp2ix],
               inputrehihilo[inp2ix],inputrelohilo[inp2ix],
               inputrehilolo[inp2ix],inputrelololo[inp2ix],
               inputimhihihi[inp2ix],inputimlohihi[inp2ix],
               inputimhilohi[inp2ix],inputimlolohi[inp2ix],
               inputimhihilo[inp2ix],inputimlohilo[inp2ix],
               inputimhilolo[inp2ix],inputimlololo[inp2ix],
               backwardrehihihi[outidx],backwardrelohihi[outidx],
               backwardrehilohi[outidx],backwardrelolohi[outidx],
               backwardrehihilo[outidx],backwardrelohilo[outidx],
               backwardrehilolo[outidx],backwardrelololo[outidx],
               backwardimhihihi[outidx],backwardimlohihi[outidx],
               backwardimhilohi[outidx],backwardimlolohi[outidx],
               backwardimhihilo[outidx],backwardimlohilo[outidx],
               backwardimhilolo[outidx],backwardimlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_cmplx8_product(deg,
               cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
               cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
               backwardrehihihi[inp2ix],backwardrelohihi[inp2ix],
               backwardrehilohi[inp2ix],backwardrelolohi[inp2ix],
               backwardrehihilo[inp2ix],backwardrelohilo[inp2ix],
               backwardrehilolo[inp2ix],backwardrelololo[inp2ix],
               backwardimhihihi[inp2ix],backwardimlohihi[inp2ix],
               backwardimhilohi[inp2ix],backwardimlolohi[inp2ix],
               backwardimhihilo[inp2ix],backwardimlohilo[inp2ix],
               backwardimhilolo[inp2ix],backwardimlololo[inp2ix],
               backwardrehihihi[outidx],backwardrelohihi[outidx],
               backwardrehilohi[outidx],backwardrelolohi[outidx],
               backwardrehihilo[outidx],backwardrelohilo[outidx],
               backwardrehilolo[outidx],backwardrelololo[outidx],
               backwardimhihihi[outidx],backwardimlohihi[outidx],
               backwardimhilohi[outidx],backwardimlolohi[outidx],
               backwardimhihilo[outidx],backwardimlohilo[outidx],
               backwardimhilolo[outidx],backwardimlololo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_cmplx8_product(deg,
               inputrehihihi[inp1ix],inputrelohihi[inp1ix],
               inputrehilohi[inp1ix],inputrelolohi[inp1ix],
               inputrehihilo[inp1ix],inputrelohilo[inp1ix],
               inputrehilolo[inp1ix],inputrelololo[inp1ix],
               inputimhihihi[inp1ix],inputimlohihi[inp1ix],
               inputimhilohi[inp1ix],inputimlolohi[inp1ix],
               inputimhihilo[inp1ix],inputimlohilo[inp1ix],
               inputimhilolo[inp1ix],inputimlololo[inp1ix],
               inputrehihihi[inp2ix],inputrelohihi[inp2ix],
               inputrehilohi[inp2ix],inputrelolohi[inp2ix],
               inputrehihilo[inp2ix],inputrelohilo[inp2ix],
               inputrehilolo[inp2ix],inputrelololo[inp2ix],
               inputimhihihi[inp2ix],inputimlohihi[inp2ix],
               inputimhilohi[inp2ix],inputimlolohi[inp2ix],
               inputimhihilo[inp2ix],inputimlohilo[inp2ix],
               inputimhilolo[inp2ix],inputimlololo[inp2ix],
               backwardrehihihi[outidx],backwardrelohihi[outidx],
               backwardrehilohi[outidx],backwardrelolohi[outidx],
               backwardrehihilo[outidx],backwardrelohilo[outidx],
               backwardrehilolo[outidx],backwardrelololo[outidx],
               backwardimhihihi[outidx],backwardimlohihi[outidx],
               backwardimhilohi[outidx],backwardimlolohi[outidx],
               backwardimhihilo[outidx],backwardimlohilo[outidx],
               backwardimhilolo[outidx],backwardimlololo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_cmplx8_product(deg,
               inputrehihihi[inp1ix],inputrelohihi[inp1ix],
               inputrehilohi[inp1ix],inputrelolohi[inp1ix],
               inputrehihilo[inp1ix],inputrelohilo[inp1ix],
               inputrehilolo[inp1ix],inputrelololo[inp1ix],
               inputimhihihi[inp1ix],inputimlohihi[inp1ix],
               inputimhilohi[inp1ix],inputimlolohi[inp1ix],
               inputimhihilo[inp1ix],inputimlohilo[inp1ix],
               inputimhilolo[inp1ix],inputimlololo[inp1ix],
               backwardrehihihi[inp2ix],backwardrelohihi[inp2ix],
               backwardrehilohi[inp2ix],backwardrelolohi[inp2ix],
               backwardrehihilo[inp2ix],backwardrelohilo[inp2ix],
               backwardrehilolo[inp2ix],backwardrelololo[inp2ix],
               backwardimhihihi[inp2ix],backwardimlohihi[inp2ix],
               backwardimhilohi[inp2ix],backwardimlolohi[inp2ix],
               backwardimhihilo[inp2ix],backwardimlohilo[inp2ix],
               backwardimhilolo[inp2ix],backwardimlololo[inp2ix],
               backwardrehihihi[outidx],backwardrelohihi[outidx],
               backwardrehilohi[outidx],backwardrelolohi[outidx],
               backwardrehihilo[outidx],backwardrelohilo[outidx],
               backwardrehilolo[outidx],backwardrelololo[outidx],
               backwardimhihihi[outidx],backwardimlohihi[outidx],
               backwardimhilohi[outidx],backwardimlolohi[outidx],
               backwardimhihilo[outidx],backwardimlohilo[outidx],
               backwardimhilolo[outidx],backwardimlololo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_cmplx8_product(deg,
               backwardrehihihi[inp1ix],backwardrelohihi[inp1ix],
               backwardrehilohi[inp1ix],backwardrelolohi[inp1ix],
               backwardrehihilo[inp1ix],backwardrelohilo[inp1ix],
               backwardrehilolo[inp1ix],backwardrelololo[inp1ix],
               backwardimhihihi[inp1ix],backwardimlohihi[inp1ix],
               backwardimhilohi[inp1ix],backwardimlolohi[inp1ix],
               backwardimhihilo[inp1ix],backwardimlohilo[inp1ix],
               backwardimhilolo[inp1ix],backwardimlololo[inp1ix],
               cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
               cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
               backwardrehihihi[outidx],backwardrelohihi[outidx],
               backwardrehilohi[outidx],backwardrelolohi[outidx],
               backwardrehihilo[outidx],backwardrelohilo[outidx],
               backwardrehilolo[outidx],backwardrelololo[outidx],
               backwardimhihihi[outidx],backwardimlohihi[outidx],
               backwardimhilohi[outidx],backwardimlolohi[outidx],
               backwardimhihilo[outidx],backwardimlohilo[outidx],
               backwardimhilolo[outidx],backwardimlololo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx8_product(deg,
               backwardrehihihi[inp1ix],backwardrelohihi[inp1ix],
               backwardrehilohi[inp1ix],backwardrelolohi[inp1ix],
               backwardrehihilo[inp1ix],backwardrelohilo[inp1ix],
               backwardrehilolo[inp1ix],backwardrelololo[inp1ix],
               backwardimhihihi[inp1ix],backwardimlohihi[inp1ix],
               backwardimhilohi[inp1ix],backwardimlolohi[inp1ix],
               backwardimhihilo[inp1ix],backwardimlohilo[inp1ix],
               backwardimhilolo[inp1ix],backwardimlololo[inp1ix],
               inputrehihihi[inp2ix],inputrelohihi[inp2ix],
               inputrehilohi[inp2ix],inputrelolohi[inp2ix],
               inputrehihilo[inp2ix],inputrelohilo[inp2ix],
               inputrehilolo[inp2ix],inputrelololo[inp2ix],
               inputimhihihi[inp2ix],inputimlohihi[inp2ix],
               inputimhilohi[inp2ix],inputimlolohi[inp2ix],
               inputimhihilo[inp2ix],inputimlohilo[inp2ix],
               inputimhilolo[inp2ix],inputimlololo[inp2ix],
               backwardrehihihi[outidx],backwardrelohihi[outidx],
               backwardrehilohi[outidx],backwardrelolohi[outidx],
               backwardrehihilo[outidx],backwardrelohilo[outidx],
               backwardrehilolo[outidx],backwardrelololo[outidx],
               backwardimhihihi[outidx],backwardimlohihi[outidx],
               backwardimhilohi[outidx],backwardimlolohi[outidx],
               backwardimhihilo[outidx],backwardimlohilo[outidx],
               backwardimhilolo[outidx],backwardimlololo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_cmplx8_product(deg,
            cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
            cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
            cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
            cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
            inputrehihihi[inp2ix],inputrelohihi[inp2ix],
            inputrehilohi[inp2ix],inputrelolohi[inp2ix],
            inputrehihilo[inp2ix],inputrelohilo[inp2ix],
            inputrehilolo[inp2ix],inputrelololo[inp2ix],
            inputimhihihi[inp2ix],inputimlohihi[inp2ix],
            inputimhilohi[inp2ix],inputimlolohi[inp2ix],
            inputimhihilo[inp2ix],inputimlohilo[inp2ix],
            inputimhilolo[inp2ix],inputimlololo[inp2ix],
            crossrehihihi[outidx],crossrelohihi[outidx],
            crossrehilohi[outidx],crossrelolohi[outidx],
            crossrehihilo[outidx],crossrelohilo[outidx],
            crossrehilolo[outidx],crossrelololo[outidx],
            crossimhihihi[outidx],crossimlohihi[outidx],
            crossimhilohi[outidx],crossimlolohi[outidx],
            crossimhihilo[outidx],crossimlohilo[outidx],
            crossimhilolo[outidx],crossimlololo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx8_product(deg,
            inputrehihihi[inp1ix],inputrelohihi[inp1ix],
            inputrehilohi[inp1ix],inputrelolohi[inp1ix],
            inputrehihilo[inp1ix],inputrelohilo[inp1ix],
            inputrehilolo[inp1ix],inputrelololo[inp1ix],
            inputimhihihi[inp1ix],inputimlohihi[inp1ix],
            inputimhilohi[inp1ix],inputimlolohi[inp1ix],
            inputimhihilo[inp1ix],inputimlohilo[inp1ix],
            inputimhilolo[inp1ix],inputimlololo[inp1ix],
            forwardrehihihi[inp2ix],forwardrelohihi[inp2ix],
            forwardrehilohi[inp2ix],forwardrelolohi[inp2ix],
            forwardrehihilo[inp2ix],forwardrelohilo[inp2ix],
            forwardrehilolo[inp2ix],forwardrelololo[inp2ix],
            forwardimhihihi[inp2ix],forwardimlohihi[inp2ix],
            forwardimhilohi[inp2ix],forwardimlolohi[inp2ix],
            forwardimhihilo[inp2ix],forwardimlohilo[inp2ix],
            forwardimhilolo[inp2ix],forwardimlololo[inp2ix],
            crossrehihihi[outidx],crossrelohihi[outidx],
            crossrehilohi[outidx],crossrelolohi[outidx],
            crossrehihilo[outidx],crossrelohilo[outidx],
            crossrehilolo[outidx],crossrelololo[outidx],
            crossimhihihi[outidx],crossimlohihi[outidx],
            crossimhilohi[outidx],crossimlolohi[outidx],
            crossimhihilo[outidx],crossimlohilo[outidx],
            crossimhilolo[outidx],crossimlololo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_cmplx8_product(deg,
              forwardrehihihi[inp1ix],forwardrelohihi[inp1ix],
              forwardrehilohi[inp1ix],forwardrelolohi[inp1ix],
              forwardrehihilo[inp1ix],forwardrelohilo[inp1ix],
              forwardrehilolo[inp1ix],forwardrelololo[inp1ix],
              forwardimhihihi[inp1ix],forwardimlohihi[inp1ix],
              forwardimhilohi[inp1ix],forwardimlolohi[inp1ix],
              forwardimhihilo[inp1ix],forwardimlohilo[inp1ix],
              forwardimhilolo[inp1ix],forwardimlololo[inp1ix],
              inputrehihihi[inp2ix],inputrelohihi[inp2ix],
              inputrehilohi[inp2ix],inputrelolohi[inp2ix],
              inputrehihilo[inp2ix],inputrelohilo[inp2ix],
              inputrehilolo[inp2ix],inputrelololo[inp2ix],
              inputimhihihi[inp2ix],inputimlohihi[inp2ix],
              inputimhilohi[inp2ix],inputimlolohi[inp2ix],
              inputimhihilo[inp2ix],inputimlohilo[inp2ix],
              inputimhilolo[inp2ix],inputimlololo[inp2ix],
              crossrehihihi[outidx],crossrelohihi[outidx],
              crossrehilohi[outidx],crossrelolohi[outidx],
              crossrehihilo[outidx],crossrelohilo[outidx],
              crossrehilolo[outidx],crossrelololo[outidx],
              crossimhihihi[outidx],crossimlohihi[outidx],
              crossimhilohi[outidx],crossimlolohi[outidx],
              crossimhihilo[outidx],crossimlohilo[outidx],
              crossimhilolo[outidx],crossimlololo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_cmplx8_product(deg,
              forwardrehihihi[inp1ix],forwardrelohihi[inp1ix],
              forwardrehilohi[inp1ix],forwardrelolohi[inp1ix],
              forwardrehihilo[inp1ix],forwardrelohilo[inp1ix],
              forwardrehilolo[inp1ix],forwardrelololo[inp1ix],
              forwardimhihihi[inp1ix],forwardimlohihi[inp1ix],
              forwardimhilohi[inp1ix],forwardimlolohi[inp1ix],
              forwardimhihilo[inp1ix],forwardimlohilo[inp1ix],
              forwardimhilolo[inp1ix],forwardimlololo[inp1ix],
              backwardrehihihi[inp2ix],backwardrelohihi[inp2ix],
              backwardrehilohi[inp2ix],backwardrelolohi[inp2ix],
              backwardrehihilo[inp2ix],backwardrelohilo[inp2ix],
              backwardrehilolo[inp2ix],backwardrelololo[inp2ix],
              backwardimhihihi[inp2ix],backwardimlohihi[inp2ix],
              backwardimhilohi[inp2ix],backwardimlolohi[inp2ix],
              backwardimhihilo[inp2ix],backwardimlohilo[inp2ix],
              backwardimhilolo[inp2ix],backwardimlololo[inp2ix],
              crossrehihihi[outidx],crossrelohihi[outidx],
              crossrehilohi[outidx],crossrelolohi[outidx],
              crossrehihilo[outidx],crossrelohilo[outidx],
              crossrehilolo[outidx],crossrelololo[outidx],
              crossimhihihi[outidx],crossimlohihi[outidx],
              crossimhilohi[outidx],crossimlolohi[outidx],
              crossimhihilo[outidx],crossimlohilo[outidx],
              crossimhilolo[outidx],crossimlololo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx8_product(deg,
            backwardrehihihi[inp1ix],backwardrelohihi[inp1ix],
            backwardrehilohi[inp1ix],backwardrelolohi[inp1ix],
            backwardrehihilo[inp1ix],backwardrelohilo[inp1ix],
            backwardrehilolo[inp1ix],backwardrelololo[inp1ix],
            backwardimhihihi[inp1ix],backwardimlohihi[inp1ix],
            backwardimhilohi[inp1ix],backwardimlolohi[inp1ix],
            backwardimhihilo[inp1ix],backwardimlohilo[inp1ix],
            backwardimhilolo[inp1ix],backwardimlololo[inp1ix],
            forwardrehihihi[inp2ix],forwardrelohihi[inp2ix],
            forwardrehilohi[inp2ix],forwardrelolohi[inp2ix],
            forwardrehihilo[inp2ix],forwardrelohilo[inp2ix],
            forwardrehilolo[inp2ix],forwardrelololo[inp2ix],
            forwardimhihihi[inp2ix],forwardimlohihi[inp2ix],
            forwardimhilohi[inp2ix],forwardimlolohi[inp2ix],
            forwardimhihilo[inp2ix],forwardimlohilo[inp2ix],
            forwardimhilolo[inp2ix],forwardimlololo[inp2ix],
            crossrehihihi[outidx],crossrelohihi[outidx],
            crossrehilohi[outidx],crossrelolohi[outidx],
            crossrehihilo[outidx],crossrelohilo[outidx],
            crossrehilolo[outidx],crossrelololo[outidx],
            crossimhihihi[outidx],crossimlohihi[outidx],
            crossimhilohi[outidx],crossimlolohi[outidx],
            crossimhihilo[outidx],crossimlohilo[outidx],
            crossimhilolo[outidx],crossimlololo[outidx]);
      }
   }
}

void CPU_dbl8_add_job
 ( int deg,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double ***forwardhihihi, double ***forwardlohihi,
   double ***forwardhilohi, double ***forwardlolohi,
   double ***forwardhihilo, double ***forwardlohilo,
   double ***forwardhilolo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardlohihi,
   double ***backwardhilohi, double ***backwardlolohi, 
   double ***backwardhihilo, double ***backwardlohilo,
   double ***backwardhilolo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosslohihi,
   double ***crosshilohi, double ***crosslolohi,
   double ***crosshihilo, double ***crosslohilo,
   double ***crosshilolo, double ***crosslololo,
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
                       &forwardlohihi[updmon][updidx][i],
                       &forwardhilohi[updmon][updidx][i],
                       &forwardlolohi[updmon][updidx][i],
                       &forwardhihilo[updmon][updidx][i],
                       &forwardlohilo[updmon][updidx][i],
                       &forwardhilolo[updmon][updidx][i],
                       &forwardlololo[updmon][updidx][i],
                       csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
                       csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               odf_inc(&forwardhihihi[updmon][updidx][i],
                       &forwardlohihi[updmon][updidx][i],
                       &forwardhilohi[updmon][updidx][i],
                       &forwardlolohi[updmon][updidx][i],
                       &forwardhihilo[updmon][updidx][i],
                       &forwardlohilo[updmon][updidx][i],
                       &forwardhilolo[updmon][updidx][i],
                       &forwardlololo[updmon][updidx][i],
                       cffhihihi[incidx][i],cfflohihi[incidx][i],
                       cffhilohi[incidx][i],cfflolohi[incidx][i],
                       cffhihilo[incidx][i],cfflohilo[incidx][i],
                       cffhilolo[incidx][i],cfflololo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            odf_inc(&forwardhihihi[updmon][updidx][i],
                    &forwardlohihi[updmon][updidx][i],
                    &forwardhilohi[updmon][updidx][i],
                    &forwardlolohi[updmon][updidx][i],
                    &forwardhihilo[updmon][updidx][i],
                    &forwardlohilo[updmon][updidx][i],
                    &forwardhilolo[updmon][updidx][i],
                    &forwardlololo[updmon][updidx][i],
                    forwardhihihi[incmon][incidx][i],
                    forwardlohihi[incmon][incidx][i],
                    forwardhilohi[incmon][incidx][i],
                    forwardlolohi[incmon][incidx][i],
                    forwardhihilo[incmon][incidx][i],
                    forwardlohilo[incmon][incidx][i],
                    forwardhilolo[incmon][incidx][i],
                    forwardlololo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            odf_inc(&forwardhihihi[updmon][updidx][i],
                    &forwardlohihi[updmon][updidx][i],
                    &forwardhilohi[updmon][updidx][i],
                    &forwardlolohi[updmon][updidx][i],
                    &forwardhihilo[updmon][updidx][i],
                    &forwardlohilo[updmon][updidx][i],
                    &forwardhilolo[updmon][updidx][i],
                    &forwardlololo[updmon][updidx][i],
                    backwardhihihi[incmon][incidx][i],
                    backwardlohihi[incmon][incidx][i],
                    backwardhilohi[incmon][incidx][i],
                    backwardlolohi[incmon][incidx][i],
                    backwardhihilo[incmon][incidx][i],
                    backwardlohilo[incmon][incidx][i],
                    backwardhilolo[incmon][incidx][i],
                    backwardlololo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            odf_inc(&forwardhihihi[updmon][updidx][i],
                    &forwardlohihi[updmon][updidx][i],
                    &forwardhilohi[updmon][updidx][i],
                    &forwardlolohi[updmon][updidx][i],
                    &forwardhihilo[updmon][updidx][i],
                    &forwardlohilo[updmon][updidx][i],
                    &forwardhilolo[updmon][updidx][i],
                    &forwardlololo[updmon][updidx][i],
                    crosshihihi[incmon][incidx][i],
                    crosslohihi[incmon][incidx][i],
                    crosshilohi[incmon][incidx][i],
                    crosslolohi[incmon][incidx][i],
                    crosshihilo[incmon][incidx][i],
                    crosslohilo[incmon][incidx][i],
                    crosshilolo[incmon][incidx][i],
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
                    &backwardlohihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    cffhihihi[incidx][i],cfflohihi[incidx][i],
                    cffhilohi[incidx][i],cfflolohi[incidx][i],
                    cffhihilo[incidx][i],cfflohilo[incidx][i],
                    cffhilolo[incidx][i],cfflololo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    forwardhihihi[incmon][incidx][i],
                    forwardlohihi[incmon][incidx][i],
                    forwardhilohi[incmon][incidx][i],
                    forwardlolohi[incmon][incidx][i],
                    forwardhihilo[incmon][incidx][i],
                    forwardlohilo[incmon][incidx][i],
                    forwardhilolo[incmon][incidx][i],
                    forwardlololo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    backwardhihihi[incmon][incidx][i],
                    backwardlohihi[incmon][incidx][i],
                    backwardhilohi[incmon][incidx][i],
                    backwardlolohi[incmon][incidx][i],
                    backwardhihilo[incmon][incidx][i],
                    backwardlohilo[incmon][incidx][i],
                    backwardhilolo[incmon][incidx][i],
                    backwardlololo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            odf_inc(&backwardhihihi[updmon][updidx][i],
                    &backwardlohihi[updmon][updidx][i],
                    &backwardhilohi[updmon][updidx][i],
                    &backwardlolohi[updmon][updidx][i],
                    &backwardhihilo[updmon][updidx][i],
                    &backwardlohilo[updmon][updidx][i],
                    &backwardhilolo[updmon][updidx][i],
                    &backwardlololo[updmon][updidx][i],
                    crosshihihi[incmon][incidx][i],
                    crosslohihi[incmon][incidx][i],
                    crosshilohi[incmon][incidx][i],
                    crosslolohi[incmon][incidx][i],
                    crosshihilo[incmon][incidx][i],
                    crosslohilo[incmon][incidx][i],
                    crosshilolo[incmon][incidx][i],
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
                    &crosslohihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    cffhihihi[incidx][i],cfflohihi[incidx][i],
                    cffhilohi[incidx][i],cfflolohi[incidx][i],
                    cffhihilo[incidx][i],cfflohilo[incidx][i],
                    cffhilolo[incidx][i],cfflololo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    forwardhihihi[incmon][incidx][i],
                    forwardlohihi[incmon][incidx][i],
                    forwardhilohi[incmon][incidx][i],
                    forwardlolohi[incmon][incidx][i],
                    forwardhihilo[incmon][incidx][i],
                    forwardlohilo[incmon][incidx][i],
                    forwardhilolo[incmon][incidx][i],
                    forwardlololo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    backwardhihihi[incmon][incidx][i],
                    backwardlohihi[incmon][incidx][i],
                    backwardhilohi[incmon][incidx][i],
                    backwardlolohi[incmon][incidx][i],
                    backwardhihilo[incmon][incidx][i],
                    backwardlohilo[incmon][incidx][i],
                    backwardhilolo[incmon][incidx][i],
                    backwardlololo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            odf_inc(&crosshihihi[updmon][updidx][i],
                    &crosslohihi[updmon][updidx][i],
                    &crosshilohi[updmon][updidx][i],
                    &crosslolohi[updmon][updidx][i],
                    &crosshihilo[updmon][updidx][i],
                    &crosslohilo[updmon][updidx][i],
                    &crosshilolo[updmon][updidx][i],
                    &crosslololo[updmon][updidx][i],
                    crosshihihi[incmon][incidx][i],
                    crosslohihi[incmon][incidx][i],
                    crosshilohi[incmon][incidx][i],
                    crosslolohi[incmon][incidx][i],
                    crosshihilo[incmon][incidx][i],
                    crosslohilo[incmon][incidx][i],
                    crosshilolo[incmon][incidx][i],
                    crosslololo[incmon][incidx][i]);
      }
   }
}

void CPU_cmplx8_add_job
 ( int deg,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double ***forwardrehihihi, double ***forwardrelohihi,
   double ***forwardrehilohi, double ***forwardrelolohi,
   double ***forwardrehihilo, double ***forwardrelohilo,
   double ***forwardrehilolo, double ***forwardrelololo,
   double ***forwardimhihihi, double ***forwardimlohihi,
   double ***forwardimhilohi, double ***forwardimlolohi,
   double ***forwardimhihilo, double ***forwardimlohilo,
   double ***forwardimhilolo, double ***forwardimlololo,
   double ***backwardrehihihi, double ***backwardrelohihi,
   double ***backwardrehilohi, double ***backwardrelolohi, 
   double ***backwardrehihilo, double ***backwardrelohilo,
   double ***backwardrehilolo, double ***backwardrelololo, 
   double ***backwardimhihihi, double ***backwardimlohihi,
   double ***backwardimhilohi, double ***backwardimlolohi, 
   double ***backwardimhihilo, double ***backwardimlohilo,
   double ***backwardimhilolo, double ***backwardimlololo, 
   double ***crossrehihihi, double ***crossrelohihi,
   double ***crossrehilohi, double ***crossrelolohi,
   double ***crossrehihilo, double ***crossrelohilo,
   double ***crossrehilolo, double ***crossrelololo,
   double ***crossimhihihi, double ***crossimlohihi,
   double ***crossimhilohi, double ***crossimlolohi,
   double ***crossimhihilo, double ***crossimlohilo,
   double ***crossimhilolo, double ***crossimlololo,
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
               odf_inc(&forwardrehihihi[updmon][updidx][i],
                       &forwardrelohihi[updmon][updidx][i],
                       &forwardrehilohi[updmon][updidx][i],
                       &forwardrelolohi[updmon][updidx][i],
                       &forwardrehihilo[updmon][updidx][i],
                       &forwardrelohilo[updmon][updidx][i],
                       &forwardrehilolo[updmon][updidx][i],
                       &forwardrelololo[updmon][updidx][i],
                       cstrehihihi[i],cstrelohihi[i],
                       cstrehilohi[i],cstrelolohi[i],
                       cstrehihilo[i],cstrelohilo[i],
                       cstrehilolo[i],cstrelololo[i]);
               odf_inc(&forwardimhihihi[updmon][updidx][i],
                       &forwardimlohihi[updmon][updidx][i],
                       &forwardimhilohi[updmon][updidx][i],
                       &forwardimlolohi[updmon][updidx][i],
                       &forwardimhihilo[updmon][updidx][i],
                       &forwardimlohilo[updmon][updidx][i],
                       &forwardimhilolo[updmon][updidx][i],
                       &forwardimlololo[updmon][updidx][i],
                       cstimhihihi[i],cstimlohihi[i],
                       cstimhilohi[i],cstimlolohi[i],
                       cstimhihilo[i],cstimlohilo[i],
                       cstimhilolo[i],cstimlololo[i]);
            }
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
            {
               odf_inc(&forwardrehihihi[updmon][updidx][i],
                       &forwardrelohihi[updmon][updidx][i],
                       &forwardrehilohi[updmon][updidx][i],
                       &forwardrelolohi[updmon][updidx][i],
                       &forwardrehihilo[updmon][updidx][i],
                       &forwardrelohilo[updmon][updidx][i],
                       &forwardrehilolo[updmon][updidx][i],
                       &forwardrelololo[updmon][updidx][i],
                       cffrehihihi[incidx][i],cffrelohihi[incidx][i],
                       cffrehilohi[incidx][i],cffrelolohi[incidx][i],
                       cffrehihilo[incidx][i],cffrelohilo[incidx][i],
                       cffrehilolo[incidx][i],cffrelololo[incidx][i]);
               odf_inc(&forwardimhihihi[updmon][updidx][i],
                       &forwardimlohihi[updmon][updidx][i],
                       &forwardimhilohi[updmon][updidx][i],
                       &forwardimlolohi[updmon][updidx][i],
                       &forwardimhihilo[updmon][updidx][i],
                       &forwardimlohilo[updmon][updidx][i],
                       &forwardimhilolo[updmon][updidx][i],
                       &forwardimlololo[updmon][updidx][i],
                       cffimhihihi[incidx][i],cffimlohihi[incidx][i],
                       cffimhilohi[incidx][i],cffimlolohi[incidx][i],
                       cffimhihilo[incidx][i],cffimlohilo[incidx][i],
                       cffimhilolo[incidx][i],cffimlololo[incidx][i]);
            }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            odf_inc(&forwardrehihihi[updmon][updidx][i],
                    &forwardrelohihi[updmon][updidx][i],
                    &forwardrehilohi[updmon][updidx][i],
                    &forwardrelolohi[updmon][updidx][i],
                    &forwardrehihilo[updmon][updidx][i],
                    &forwardrelohilo[updmon][updidx][i],
                    &forwardrehilolo[updmon][updidx][i],
                    &forwardrelololo[updmon][updidx][i],
                    forwardrehihihi[incmon][incidx][i],
                    forwardrelohihi[incmon][incidx][i],
                    forwardrehilohi[incmon][incidx][i],
                    forwardrelolohi[incmon][incidx][i],
                    forwardrehihilo[incmon][incidx][i],
                    forwardrelohilo[incmon][incidx][i],
                    forwardrehilolo[incmon][incidx][i],
                    forwardrelololo[incmon][incidx][i]);
            odf_inc(&forwardimhihihi[updmon][updidx][i],
                    &forwardimlohihi[updmon][updidx][i],
                    &forwardimhilohi[updmon][updidx][i],
                    &forwardimlolohi[updmon][updidx][i],
                    &forwardimhihilo[updmon][updidx][i],
                    &forwardimlohilo[updmon][updidx][i],
                    &forwardimhilolo[updmon][updidx][i],
                    &forwardimlololo[updmon][updidx][i],
                    forwardimhihihi[incmon][incidx][i],
                    forwardimlohihi[incmon][incidx][i],
                    forwardimhilohi[incmon][incidx][i],
                    forwardimlolohi[incmon][incidx][i],
                    forwardimhihilo[incmon][incidx][i],
                    forwardimlohilo[incmon][incidx][i],
                    forwardimhilolo[incmon][incidx][i],
                    forwardimlololo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            odf_inc(&forwardrehihihi[updmon][updidx][i],
                    &forwardrelohihi[updmon][updidx][i],
                    &forwardrehilohi[updmon][updidx][i],
                    &forwardrelolohi[updmon][updidx][i],
                    &forwardrehihilo[updmon][updidx][i],
                    &forwardrelohilo[updmon][updidx][i],
                    &forwardrehilolo[updmon][updidx][i],
                    &forwardrelololo[updmon][updidx][i],
                    backwardrehihihi[incmon][incidx][i],
                    backwardrelohihi[incmon][incidx][i],
                    backwardrehilohi[incmon][incidx][i],
                    backwardrelolohi[incmon][incidx][i],
                    backwardrehihilo[incmon][incidx][i],
                    backwardrelohilo[incmon][incidx][i],
                    backwardrehilolo[incmon][incidx][i],
                    backwardrelololo[incmon][incidx][i]);
            odf_inc(&forwardimhihihi[updmon][updidx][i],
                    &forwardimlohihi[updmon][updidx][i],
                    &forwardimhilohi[updmon][updidx][i],
                    &forwardimlolohi[updmon][updidx][i],
                    &forwardimhihilo[updmon][updidx][i],
                    &forwardimlohilo[updmon][updidx][i],
                    &forwardimhilolo[updmon][updidx][i],
                    &forwardimlololo[updmon][updidx][i],
                    backwardimhihihi[incmon][incidx][i],
                    backwardimlohihi[incmon][incidx][i],
                    backwardimhilohi[incmon][incidx][i],
                    backwardimlolohi[incmon][incidx][i],
                    backwardimhihilo[incmon][incidx][i],
                    backwardimlohilo[incmon][incidx][i],
                    backwardimhilolo[incmon][incidx][i],
                    backwardimlololo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            odf_inc(&forwardrehihihi[updmon][updidx][i],
                    &forwardrelohihi[updmon][updidx][i],
                    &forwardrehilohi[updmon][updidx][i],
                    &forwardrelolohi[updmon][updidx][i],
                    &forwardrehihilo[updmon][updidx][i],
                    &forwardrelohilo[updmon][updidx][i],
                    &forwardrehilolo[updmon][updidx][i],
                    &forwardrelololo[updmon][updidx][i],
                    crossrehihihi[incmon][incidx][i],
                    crossrelohihi[incmon][incidx][i],
                    crossrehilohi[incmon][incidx][i],
                    crossrelolohi[incmon][incidx][i],
                    crossrehihilo[incmon][incidx][i],
                    crossrelohilo[incmon][incidx][i],
                    crossrehilolo[incmon][incidx][i],
                    crossrelololo[incmon][incidx][i]);
            odf_inc(&forwardimhihihi[updmon][updidx][i],
                    &forwardimlohihi[updmon][updidx][i],
                    &forwardimhilohi[updmon][updidx][i],
                    &forwardimlolohi[updmon][updidx][i],
                    &forwardimhihilo[updmon][updidx][i],
                    &forwardimlohilo[updmon][updidx][i],
                    &forwardimhilolo[updmon][updidx][i],
                    &forwardimlololo[updmon][updidx][i],
                    crossimhihihi[incmon][incidx][i],
                    crossimlohihi[incmon][incidx][i],
                    crossimhilohi[incmon][incidx][i],
                    crossimlolohi[incmon][incidx][i],
                    crossimhihilo[incmon][incidx][i],
                    crossimlohilo[incmon][incidx][i],
                    crossimhilolo[incmon][incidx][i],
                    crossimlololo[incmon][incidx][i]);
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
            odf_inc(&backwardrehihihi[updmon][updidx][i],
                    &backwardrelohihi[updmon][updidx][i],
                    &backwardrehilohi[updmon][updidx][i],
                    &backwardrelolohi[updmon][updidx][i],
                    &backwardrehihilo[updmon][updidx][i],
                    &backwardrelohilo[updmon][updidx][i],
                    &backwardrehilolo[updmon][updidx][i],
                    &backwardrelololo[updmon][updidx][i],
                    cffrehihihi[incidx][i],cffrelohihi[incidx][i],
                    cffrehilohi[incidx][i],cffrelolohi[incidx][i],
                    cffrehihilo[incidx][i],cffrelohilo[incidx][i],
                    cffrehilolo[incidx][i],cffrelololo[incidx][i]);
            odf_inc(&backwardimhihihi[updmon][updidx][i],
                    &backwardimlohihi[updmon][updidx][i],
                    &backwardimhilohi[updmon][updidx][i],
                    &backwardimlolohi[updmon][updidx][i],
                    &backwardimhihilo[updmon][updidx][i],
                    &backwardimlohilo[updmon][updidx][i],
                    &backwardimhilolo[updmon][updidx][i],
                    &backwardimlololo[updmon][updidx][i],
                    cffimhihihi[incidx][i],cffimlohihi[incidx][i],
                    cffimhilohi[incidx][i],cffimlolohi[incidx][i],
                    cffimhihilo[incidx][i],cffimlohilo[incidx][i],
                    cffimhilolo[incidx][i],cffimlololo[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            odf_inc(&backwardrehihihi[updmon][updidx][i],
                    &backwardrelohihi[updmon][updidx][i],
                    &backwardrehilohi[updmon][updidx][i],
                    &backwardrelolohi[updmon][updidx][i],
                    &backwardrehihilo[updmon][updidx][i],
                    &backwardrelohilo[updmon][updidx][i],
                    &backwardrehilolo[updmon][updidx][i],
                    &backwardrelololo[updmon][updidx][i],
                    forwardrehihihi[incmon][incidx][i],
                    forwardrelohihi[incmon][incidx][i],
                    forwardrehilohi[incmon][incidx][i],
                    forwardrelolohi[incmon][incidx][i],
                    forwardrehihilo[incmon][incidx][i],
                    forwardrelohilo[incmon][incidx][i],
                    forwardrehilolo[incmon][incidx][i],
                    forwardrelololo[incmon][incidx][i]);
            odf_inc(&backwardimhihihi[updmon][updidx][i],
                    &backwardimlohihi[updmon][updidx][i],
                    &backwardimhilohi[updmon][updidx][i],
                    &backwardimlolohi[updmon][updidx][i],
                    &backwardimhihilo[updmon][updidx][i],
                    &backwardimlohilo[updmon][updidx][i],
                    &backwardimhilolo[updmon][updidx][i],
                    &backwardimlololo[updmon][updidx][i],
                    forwardimhihihi[incmon][incidx][i],
                    forwardimlohihi[incmon][incidx][i],
                    forwardimhilohi[incmon][incidx][i],
                    forwardimlolohi[incmon][incidx][i],
                    forwardimhihilo[incmon][incidx][i],
                    forwardimlohilo[incmon][incidx][i],
                    forwardimhilolo[incmon][incidx][i],
                    forwardimlololo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            odf_inc(&backwardrehihihi[updmon][updidx][i],
                    &backwardrelohihi[updmon][updidx][i],
                    &backwardrehilohi[updmon][updidx][i],
                    &backwardrelolohi[updmon][updidx][i],
                    &backwardrehihilo[updmon][updidx][i],
                    &backwardrelohilo[updmon][updidx][i],
                    &backwardrehilolo[updmon][updidx][i],
                    &backwardrelololo[updmon][updidx][i],
                    backwardrehihihi[incmon][incidx][i],
                    backwardrelohihi[incmon][incidx][i],
                    backwardrehilohi[incmon][incidx][i],
                    backwardrelolohi[incmon][incidx][i],
                    backwardrehihilo[incmon][incidx][i],
                    backwardrelohilo[incmon][incidx][i],
                    backwardrehilolo[incmon][incidx][i],
                    backwardrelololo[incmon][incidx][i]);
            odf_inc(&backwardimhihihi[updmon][updidx][i],
                    &backwardimlohihi[updmon][updidx][i],
                    &backwardimhilohi[updmon][updidx][i],
                    &backwardimlolohi[updmon][updidx][i],
                    &backwardimhihilo[updmon][updidx][i],
                    &backwardimlohilo[updmon][updidx][i],
                    &backwardimhilolo[updmon][updidx][i],
                    &backwardimlololo[updmon][updidx][i],
                    backwardimhihihi[incmon][incidx][i],
                    backwardimlohihi[incmon][incidx][i],
                    backwardimhilohi[incmon][incidx][i],
                    backwardimlolohi[incmon][incidx][i],
                    backwardimhihilo[incmon][incidx][i],
                    backwardimlohilo[incmon][incidx][i],
                    backwardimhilolo[incmon][incidx][i],
                    backwardimlololo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            odf_inc(&backwardrehihihi[updmon][updidx][i],
                    &backwardrelohihi[updmon][updidx][i],
                    &backwardrehilohi[updmon][updidx][i],
                    &backwardrelolohi[updmon][updidx][i],
                    &backwardrehihilo[updmon][updidx][i],
                    &backwardrelohilo[updmon][updidx][i],
                    &backwardrehilolo[updmon][updidx][i],
                    &backwardrelololo[updmon][updidx][i],
                    crossrehihihi[incmon][incidx][i],
                    crossrelohihi[incmon][incidx][i],
                    crossrehilohi[incmon][incidx][i],
                    crossrelolohi[incmon][incidx][i],
                    crossrehihilo[incmon][incidx][i],
                    crossrelohilo[incmon][incidx][i],
                    crossrehilolo[incmon][incidx][i],
                    crossrelololo[incmon][incidx][i]);
            odf_inc(&backwardimhihihi[updmon][updidx][i],
                    &backwardimlohihi[updmon][updidx][i],
                    &backwardimhilohi[updmon][updidx][i],
                    &backwardimlolohi[updmon][updidx][i],
                    &backwardimhihilo[updmon][updidx][i],
                    &backwardimlohilo[updmon][updidx][i],
                    &backwardimhilolo[updmon][updidx][i],
                    &backwardimlololo[updmon][updidx][i],
                    crossimhihihi[incmon][incidx][i],
                    crossimlohihi[incmon][incidx][i],
                    crossimhilohi[incmon][incidx][i],
                    crossimlolohi[incmon][incidx][i],
                    crossimhihilo[incmon][incidx][i],
                    crossimlohilo[incmon][incidx][i],
                    crossimhilolo[incmon][incidx][i],
                    crossimlololo[incmon][incidx][i]);
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
            odf_inc(&crossrehihihi[updmon][updidx][i],
                    &crossrelohihi[updmon][updidx][i],
                    &crossrehilohi[updmon][updidx][i],
                    &crossrelolohi[updmon][updidx][i],
                    &crossrehihilo[updmon][updidx][i],
                    &crossrelohilo[updmon][updidx][i],
                    &crossrehilolo[updmon][updidx][i],
                    &crossrelololo[updmon][updidx][i],
                    cffrehihihi[incidx][i],cffrelohihi[incidx][i],
                    cffrehilohi[incidx][i],cffrelolohi[incidx][i],
                    cffrehihilo[incidx][i],cffrelohilo[incidx][i],
                    cffrehilolo[incidx][i],cffrelololo[incidx][i]);
            odf_inc(&crossimhihihi[updmon][updidx][i],
                    &crossimlohihi[updmon][updidx][i],
                    &crossimhilohi[updmon][updidx][i],
                    &crossimlolohi[updmon][updidx][i],
                    &crossimhihilo[updmon][updidx][i],
                    &crossimlohilo[updmon][updidx][i],
                    &crossimhilolo[updmon][updidx][i],
                    &crossimlololo[updmon][updidx][i],
                    cffimhihihi[incidx][i],cffimlohihi[incidx][i],
                    cffimhilohi[incidx][i],cffimlolohi[incidx][i],
                    cffimhihilo[incidx][i],cffimlohilo[incidx][i],
                    cffimhilolo[incidx][i],cffimlololo[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            odf_inc(&crossrehihihi[updmon][updidx][i],
                    &crossrelohihi[updmon][updidx][i],
                    &crossrehilohi[updmon][updidx][i],
                    &crossrelolohi[updmon][updidx][i],
                    &crossrehihilo[updmon][updidx][i],
                    &crossrelohilo[updmon][updidx][i],
                    &crossrehilolo[updmon][updidx][i],
                    &crossrelololo[updmon][updidx][i],
                    forwardrehihihi[incmon][incidx][i],
                    forwardrelohihi[incmon][incidx][i],
                    forwardrehilohi[incmon][incidx][i],
                    forwardrelolohi[incmon][incidx][i],
                    forwardrehihilo[incmon][incidx][i],
                    forwardrelohilo[incmon][incidx][i],
                    forwardrehilolo[incmon][incidx][i],
                    forwardrelololo[incmon][incidx][i]);
            odf_inc(&crossimhihihi[updmon][updidx][i],
                    &crossimlohihi[updmon][updidx][i],
                    &crossimhilohi[updmon][updidx][i],
                    &crossimlolohi[updmon][updidx][i],
                    &crossimhihilo[updmon][updidx][i],
                    &crossimlohilo[updmon][updidx][i],
                    &crossimhilolo[updmon][updidx][i],
                    &crossimlololo[updmon][updidx][i],
                    forwardimhihihi[incmon][incidx][i],
                    forwardimlohihi[incmon][incidx][i],
                    forwardimhilohi[incmon][incidx][i],
                    forwardimlolohi[incmon][incidx][i],
                    forwardimhihilo[incmon][incidx][i],
                    forwardimlohilo[incmon][incidx][i],
                    forwardimhilolo[incmon][incidx][i],
                    forwardimlololo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            odf_inc(&crossrehihihi[updmon][updidx][i],
                    &crossrelohihi[updmon][updidx][i],
                    &crossrehilohi[updmon][updidx][i],
                    &crossrelolohi[updmon][updidx][i],
                    &crossrehihilo[updmon][updidx][i],
                    &crossrelohilo[updmon][updidx][i],
                    &crossrehilolo[updmon][updidx][i],
                    &crossrelololo[updmon][updidx][i],
                    backwardrehihihi[incmon][incidx][i],
                    backwardrelohihi[incmon][incidx][i],
                    backwardrehilohi[incmon][incidx][i],
                    backwardrelolohi[incmon][incidx][i],
                    backwardrehihilo[incmon][incidx][i],
                    backwardrelohilo[incmon][incidx][i],
                    backwardrehilolo[incmon][incidx][i],
                    backwardrelololo[incmon][incidx][i]);
            odf_inc(&crossimhihihi[updmon][updidx][i],
                    &crossimlohihi[updmon][updidx][i],
                    &crossimhilohi[updmon][updidx][i],
                    &crossimlolohi[updmon][updidx][i],
                    &crossimhihilo[updmon][updidx][i],
                    &crossimlohilo[updmon][updidx][i],
                    &crossimhilolo[updmon][updidx][i],
                    &crossimlololo[updmon][updidx][i],
                    backwardimhihihi[incmon][incidx][i],
                    backwardimlohihi[incmon][incidx][i],
                    backwardimhilohi[incmon][incidx][i],
                    backwardimlolohi[incmon][incidx][i],
                    backwardimhihilo[incmon][incidx][i],
                    backwardimlohilo[incmon][incidx][i],
                    backwardimhilolo[incmon][incidx][i],
                    backwardimlololo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            odf_inc(&crossrehihihi[updmon][updidx][i],
                    &crossrelohihi[updmon][updidx][i],
                    &crossrehilohi[updmon][updidx][i],
                    &crossrelolohi[updmon][updidx][i],
                    &crossrehihilo[updmon][updidx][i],
                    &crossrelohilo[updmon][updidx][i],
                    &crossrehilolo[updmon][updidx][i],
                    &crossrelololo[updmon][updidx][i],
                    crossrehihihi[incmon][incidx][i],
                    crossrelohihi[incmon][incidx][i],
                    crossrehilohi[incmon][incidx][i],
                    crossrelolohi[incmon][incidx][i],
                    crossrehihilo[incmon][incidx][i],
                    crossrelohilo[incmon][incidx][i],
                    crossrehilolo[incmon][incidx][i],
                    crossrelololo[incmon][incidx][i]);
            odf_inc(&crossimhihihi[updmon][updidx][i],
                    &crossimlohihi[updmon][updidx][i],
                    &crossimhilohi[updmon][updidx][i],
                    &crossimlolohi[updmon][updidx][i],
                    &crossimhihilo[updmon][updidx][i],
                    &crossimlohilo[updmon][updidx][i],
                    &crossimhilolo[updmon][updidx][i],
                    &crossimlololo[updmon][updidx][i],
                    crossimhihihi[incmon][incidx][i],
                    crossimlohihi[incmon][incidx][i],
                    crossimhilohi[incmon][incidx][i],
                    crossimlolohi[incmon][incidx][i],
                    crossimhihilo[incmon][incidx][i],
                    crossimlohilo[incmon][incidx][i],
                    crossimhilolo[incmon][incidx][i],
                    crossimlololo[incmon][incidx][i]);
         }
      }
   }
}

void CPU_dbl8_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
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
   double ***forwardhihihi, double ***forwardlohihi,
   double ***forwardhilohi, double ***forwardlolohi,
   double ***forwardhihilo, double ***forwardlohilo,
   double ***forwardhilolo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardlohihi,
   double ***backwardhilohi, double ***backwardlolohi, 
   double ***backwardhihilo, double ***backwardlohilo,
   double ***backwardhilolo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosslohihi,
   double ***crosshilohi, double ***crosslolohi,
   double ***crosshihilo, double ***crosslohilo,
   double ***crosshilolo, double ***crosslololo )
{
   for(int i=0; i<=deg; i++)
   {
      outputhihihi[dim][i] = csthihihi[i];
      outputlohihi[dim][i] = cstlohihi[i];
      outputhilohi[dim][i] = csthilohi[i];
      outputlolohi[dim][i] = cstlolohi[i];
      outputhihilo[dim][i] = csthihilo[i];
      outputlohilo[dim][i] = cstlohilo[i];
      outputhilolo[dim][i] = csthilolo[i];
      outputlololo[dim][i] = cstlololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhihihi[i][j] = 0.0;
         outputlohihi[i][j] = 0.0;
         outputhilohi[i][j] = 0.0;
         outputlolohi[i][j] = 0.0;
         outputhihilo[i][j] = 0.0;
         outputlohilo[i][j] = 0.0;
         outputhilolo[i][j] = 0.0;
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
         odf_inc(&outputhihihi[dim][i],   &outputlohihi[dim][i],
                 &outputhilohi[dim][i],   &outputlolohi[dim][i],
                 &outputhihilo[dim][i],   &outputlohilo[dim][i],
                 &outputhilolo[dim][i],   &outputlololo[dim][i],
                 forwardhihihi[k][ix1][i],forwardlohihi[k][ix1][i],
                 forwardhilohi[k][ix1][i],forwardlolohi[k][ix1][i],
                 forwardhihilo[k][ix1][i],forwardlohihi[k][ix1][i],
                 forwardhilolo[k][ix1][i],forwardlololo[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            odf_inc(&outputhihihi[ix0][i],&outputlohihi[ix0][i],
                    &outputhilohi[ix0][i],&outputlolohi[ix0][i],
                    &outputhihilo[ix0][i],&outputlohilo[ix0][i],
                    &outputhilolo[ix0][i],&outputlololo[ix0][i],
                        cffhihihi[k][i],      cfflohihi[k][i],
                        cffhilohi[k][i],      cfflolohi[k][i],
                        cffhihilo[k][i],      cfflohilo[k][i],
                        cffhilolo[k][i],      cfflololo[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            odf_inc(&outputhihihi[ixn][i],   &outputlohihi[ixn][i],
                    &outputhilohi[ixn][i],   &outputlolohi[ixn][i],
                    &outputhihilo[ixn][i],   &outputlohilo[ixn][i],
                    &outputhilolo[ixn][i],   &outputlololo[ixn][i],
                    forwardhihihi[k][ix2][i],forwardlohihi[k][ix2][i],
                    forwardhilohi[k][ix2][i],forwardlolohi[k][ix2][i],
                    forwardhihilo[k][ix2][i],forwardlohilo[k][ix2][i],
                    forwardhilolo[k][ix2][i],forwardlololo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            odf_inc( &outputhihihi[ix0][i],    &outputlohihi[ix0][i],
                     &outputhilohi[ix0][i],    &outputlolohi[ix0][i],
                     &outputhihilo[ix0][i],    &outputlohilo[ix0][i],
                     &outputhilolo[ix0][i],    &outputlololo[ix0][i],
                    backwardhihihi[k][ix2][i],backwardlohihi[k][ix2][i],
                    backwardhilohi[k][ix2][i],backwardlolohi[k][ix2][i],
                    backwardhihilo[k][ix2][i],backwardlohilo[k][ix2][i],
                    backwardhilolo[k][ix2][i],backwardlololo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  odf_inc(&outputhihihi[ix0][i], &outputlohihi[ix0][i],
                          &outputhilohi[ix0][i], &outputlolohi[ix0][i],
                          &outputhihilo[ix0][i], &outputlohilo[ix0][i],
                          &outputhilolo[ix0][i], &outputlololo[ix0][i],
                            crosshihihi[k][j-1][i],crosslohihi[k][j-1][i],
                            crosshilohi[k][j-1][i],crosslolohi[k][j-1][i],
                            crosshihilo[k][j-1][i],crosslohilo[k][j-1][i],
                            crosshilolo[k][j-1][i],crosslololo[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_cmplx8_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
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
   double ***forwardrehihihi, double ***forwardrelohihi,
   double ***forwardrehilohi, double ***forwardrelolohi,
   double ***forwardrehihilo, double ***forwardrelohilo,
   double ***forwardrehilolo, double ***forwardrelololo,
   double ***forwardimhihihi, double ***forwardimlohihi,
   double ***forwardimhilohi, double ***forwardimlolohi,
   double ***forwardimhihilo, double ***forwardimlohilo,
   double ***forwardimhilolo, double ***forwardimlololo,
   double ***backwardrehihihi, double ***backwardrelohihi,
   double ***backwardrehilohi, double ***backwardrelolohi, 
   double ***backwardrehihilo, double ***backwardrelohilo,
   double ***backwardrehilolo, double ***backwardrelololo, 
   double ***backwardimhihihi, double ***backwardimlohihi,
   double ***backwardimhilohi, double ***backwardimlolohi, 
   double ***backwardimhihilo, double ***backwardimlohilo,
   double ***backwardimhilolo, double ***backwardimlololo, 
   double ***crossrehihihi, double ***crossrelohihi,
   double ***crossrehilohi, double ***crossrelolohi,
   double ***crossrehihilo, double ***crossrelohilo,
   double ***crossrehilolo, double ***crossrelololo,
   double ***crossimhihihi, double ***crossimlohihi,
   double ***crossimhilohi, double ***crossimlolohi,
   double ***crossimhihilo, double ***crossimlohilo,
   double ***crossimhilolo, double ***crossimlololo )
{
   for(int i=0; i<=deg; i++)
   {
      outputrehihihi[dim][i] = cstrehihihi[i];
      outputrelohihi[dim][i] = cstrelohihi[i];
      outputrehilohi[dim][i] = cstrehilohi[i];
      outputrelolohi[dim][i] = cstrelolohi[i];
      outputrehihilo[dim][i] = cstrehihilo[i];
      outputrelohilo[dim][i] = cstrelohilo[i];
      outputrehilolo[dim][i] = cstrehilolo[i];
      outputrelololo[dim][i] = cstrelololo[i];
      outputimhihihi[dim][i] = cstimhihihi[i];
      outputimlohihi[dim][i] = cstimlohihi[i];
      outputimhilohi[dim][i] = cstimhilohi[i];
      outputimlolohi[dim][i] = cstimlolohi[i];
      outputimhihilo[dim][i] = cstimhihilo[i];
      outputimlohilo[dim][i] = cstimlohilo[i];
      outputimhilolo[dim][i] = cstimhilolo[i];
      outputimlololo[dim][i] = cstimlololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrehihihi[i][j] = 0.0; outputrelohihi[i][j] = 0.0;
         outputrehilohi[i][j] = 0.0; outputrelolohi[i][j] = 0.0;
         outputrehihilo[i][j] = 0.0; outputrelohilo[i][j] = 0.0;
         outputrehilolo[i][j] = 0.0; outputrelololo[i][j] = 0.0;
         outputimhihihi[i][j] = 0.0; outputimlohihi[i][j] = 0.0;
         outputimhilohi[i][j] = 0.0; outputimlolohi[i][j] = 0.0;
         outputimhihilo[i][j] = 0.0; outputimlohilo[i][j] = 0.0;
         outputimhilolo[i][j] = 0.0; outputimlololo[i][j] = 0.0;
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
         odf_inc(&outputrehihihi[dim][i],&outputrelohihi[dim][i],
                 &outputrehilohi[dim][i],&outputrelolohi[dim][i],
                 &outputrehihilo[dim][i],&outputrelohilo[dim][i],
                 &outputrehilolo[dim][i],&outputrelololo[dim][i],
                 forwardrehihihi[k][ix1][i],forwardrelohihi[k][ix1][i],
                 forwardrehilohi[k][ix1][i],forwardrelolohi[k][ix1][i],
                 forwardrehihilo[k][ix1][i],forwardrelohilo[k][ix1][i],
                 forwardrehilolo[k][ix1][i],forwardrelololo[k][ix1][i]);
         odf_inc(&outputimhihihi[dim][i],&outputimlohihi[dim][i],
                 &outputimhilohi[dim][i],&outputimlolohi[dim][i],
                 &outputimhihilo[dim][i],&outputimlohilo[dim][i],
                 &outputimhilolo[dim][i],&outputimlololo[dim][i],
                 forwardimhihihi[k][ix1][i],forwardimlohihi[k][ix1][i],
                 forwardimhilohi[k][ix1][i],forwardimlolohi[k][ix1][i],
                 forwardimhihilo[k][ix1][i],forwardimlohilo[k][ix1][i],
                 forwardimhilolo[k][ix1][i],forwardimlololo[k][ix1][i]);
      }
      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
         {
            odf_inc(&outputrehihihi[ix0][i],&outputrelohihi[ix0][i],
                    &outputrehilohi[ix0][i],&outputrelolohi[ix0][i],
                    &outputrehihilo[ix0][i],&outputrelohilo[ix0][i],
                    &outputrehilolo[ix0][i],&outputrelololo[ix0][i],
                    cffrehihihi[k][i],cffrelohihi[k][i],
                    cffrehilohi[k][i],cffrelolohi[k][i],
                    cffrehihilo[k][i],cffrelohilo[k][i],
                    cffrehilolo[k][i],cffrelololo[k][i]);
            odf_inc(&outputimhihihi[ix0][i],&outputimlohihi[ix0][i],
                    &outputimhilohi[ix0][i],&outputimlolohi[ix0][i],
                    &outputimhihilo[ix0][i],&outputimlohilo[ix0][i],
                    &outputimhilolo[ix0][i],&outputimlololo[ix0][i],
                    cffimhihihi[k][i],cffimlohihi[k][i],
                    cffimhilohi[k][i],cffimlolohi[k][i],
                    cffimhihilo[k][i],cffimlohilo[k][i],
                    cffimhilolo[k][i],cffimlololo[k][i]);
         }
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            odf_inc(&outputrehihihi[ixn][i],&outputrelohihi[ixn][i],
                    &outputrehilohi[ixn][i],&outputrelolohi[ixn][i],
                    &outputrehihilo[ixn][i],&outputrelohilo[ixn][i],
                    &outputrehilolo[ixn][i],&outputrelololo[ixn][i],
                    forwardrehihihi[k][ix2][i],forwardrelohihi[k][ix2][i],
                    forwardrehilohi[k][ix2][i],forwardrelolohi[k][ix2][i],
                    forwardrehihilo[k][ix2][i],forwardrelohilo[k][ix2][i],
                    forwardrehilolo[k][ix2][i],forwardrelololo[k][ix2][i]);
            odf_inc(&outputimhihihi[ixn][i],&outputimlohihi[ixn][i],
                    &outputimhilohi[ixn][i],&outputimlolohi[ixn][i],
                    &outputimhihilo[ixn][i],&outputimlohilo[ixn][i],
                    &outputimhilolo[ixn][i],&outputimlololo[ixn][i],
                    forwardimhihihi[k][ix2][i],forwardimlohihi[k][ix2][i],
                    forwardimhilohi[k][ix2][i],forwardimlolohi[k][ix2][i],
                    forwardimhihilo[k][ix2][i],forwardimlohilo[k][ix2][i],
                    forwardimhilolo[k][ix2][i],forwardimlololo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            odf_inc(&outputrehihihi[ix0][i],&outputrelohihi[ix0][i],
                    &outputrehilohi[ix0][i],&outputrelolohi[ix0][i],
                    &outputrehihilo[ix0][i],&outputrelohilo[ix0][i],
                    &outputrehilolo[ix0][i],&outputrelololo[ix0][i],
                    backwardrehihihi[k][ix2][i],backwardrelohihi[k][ix2][i],
                    backwardrehilohi[k][ix2][i],backwardrelolohi[k][ix2][i],
                    backwardrehihilo[k][ix2][i],backwardrelohilo[k][ix2][i],
                    backwardrehilolo[k][ix2][i],backwardrelololo[k][ix2][i]);
            odf_inc(&outputimhihihi[ix0][i],&outputimlohihi[ix0][i],
                    &outputimhilohi[ix0][i],&outputimlolohi[ix0][i],
                    &outputimhihilo[ix0][i],&outputimlohilo[ix0][i],
                    &outputimhilolo[ix0][i],&outputimlololo[ix0][i],
                    backwardimhihihi[k][ix2][i],backwardimlohihi[k][ix2][i],
                    backwardimhilohi[k][ix2][i],backwardimlolohi[k][ix2][i],
                    backwardimhihilo[k][ix2][i],backwardimlohilo[k][ix2][i],
                    backwardimhilolo[k][ix2][i],backwardimlololo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
               {
                  odf_inc(&outputrehihihi[ix0][i],&outputrelohihi[ix0][i],
                          &outputrehilohi[ix0][i],&outputrelolohi[ix0][i],
                          &outputrehihilo[ix0][i],&outputrelohilo[ix0][i],
                          &outputrehilolo[ix0][i],&outputrelololo[ix0][i],
                          crossrehihihi[k][j-1][i],crossrelohihi[k][j-1][i],
                          crossrehilohi[k][j-1][i],crossrelolohi[k][j-1][i],
                          crossrehihilo[k][j-1][i],crossrelohilo[k][j-1][i],
                          crossrehilolo[k][j-1][i],crossrelololo[k][j-1][i]);
                  odf_inc(&outputimhihihi[ix0][i],&outputimlohihi[ix0][i],
                          &outputimhilohi[ix0][i],&outputimlolohi[ix0][i],
                          &outputimhihilo[ix0][i],&outputimlohilo[ix0][i],
                          &outputimhilolo[ix0][i],&outputimlololo[ix0][i],
                          crossimhihihi[k][j-1][i],crossimlohihi[k][j-1][i],
                          crossimhilohi[k][j-1][i],crossimlolohi[k][j-1][i],
                          crossimhihilo[k][j-1][i],crossimlohilo[k][j-1][i],
                          crossimhilolo[k][j-1][i],crossimlololo[k][j-1][i]);
               }
            }
         }
      }
   }
}

void CPU_dbl8_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
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
   double ***forwardhihihi, double ***forwardlohihi,
   double ***forwardhilohi, double ***forwardlolohi,
   double ***forwardhihilo, double ***forwardlohilo,
   double ***forwardhilolo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardlohihi,
   double ***backwardhilohi, double ***backwardlolohi, 
   double ***backwardhihilo, double ***backwardlohilo,
   double ***backwardhilolo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosslohihi,
   double ***crosshilohi, double ***crosslolohi,
   double ***crosshihilo, double ***crosslohilo,
   double ***crosshilolo, double ***crosslololo,
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
                 csthihihi,     cstlohihi,     csthilohi,     cstlolohi,
                 csthihilo,     cstlohilo,     csthilolo,     cstlololo,
                 cffhihihi,     cfflohihi,     cffhilohi,     cfflolohi,
                 cffhihilo,     cfflohilo,     cffhilolo,     cfflololo,
             forwardhihihi, forwardlohihi, forwardhilohi, forwardlolohi,
             forwardhihilo, forwardlohilo, forwardhilolo, forwardlololo,
            backwardhihihi,backwardlohihi,backwardhilohi,backwardlolohi,
            backwardhihilo,backwardlohilo,backwardhilolo,backwardlololo,
               crosshihihi,   crosslohihi,   crosshilohi,   crosslolohi,
               crosshihilo,   crosslohilo,   crosshilolo,   crosslololo,
            job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputhihihi[dim][i] = forwardhihihi[lastmon][lastidx][i];
      outputlohihi[dim][i] = forwardlohihi[lastmon][lastidx][i];
      outputhilohi[dim][i] = forwardhilohi[lastmon][lastidx][i];
      outputlolohi[dim][i] = forwardlolohi[lastmon][lastidx][i];
      outputhihilo[dim][i] = forwardhihilo[lastmon][lastidx][i];
      outputlohilo[dim][i] = forwardlohilo[lastmon][lastidx][i];
      outputhilolo[dim][i] = forwardhilolo[lastmon][lastidx][i];
      outputlololo[dim][i] = forwardlololo[lastmon][lastidx][i];
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
            outputhihihi[0][i] = 0.0;
            outputlohihi[0][i] = 0.0;
            outputhilohi[0][i] = 0.0;
            outputlolohi[0][i] = 0.0;
            outputhihilo[0][i] = 0.0;
            outputlohilo[0][i] = 0.0;
            outputhilolo[0][i] = 0.0;
            outputlololo[0][i] = 0.0;
         }
      }
      else
      {
         if(verbose)
            cout << "updating derivative 0 with coefficient "
                 << difidx << endl;

         for(int i=0; i<=deg; i++)
         {
            outputhihihi[0][i] = cffhihihi[difidx][i];
            outputlohihi[0][i] = cfflohihi[difidx][i];
            outputhilohi[0][i] = cffhilohi[difidx][i];
            outputlolohi[0][i] = cfflolohi[difidx][i];
            outputhihilo[0][i] = cffhihilo[difidx][i];
            outputlohilo[0][i] = cfflohilo[difidx][i];
            outputhilolo[0][i] = cffhilolo[difidx][i];
            outputlololo[0][i] = cfflololo[difidx][i];
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
         outputhihihi[0][i] = backwardhihihi[ix0][ix2][i];
         outputlohihi[0][i] = backwardlohihi[ix0][ix2][i];
         outputhilohi[0][i] = backwardhilohi[ix0][ix2][i];
         outputlolohi[0][i] = backwardlolohi[ix0][ix2][i];
         outputhihilo[0][i] = backwardhihilo[ix0][ix2][i];
         outputlohilo[0][i] = backwardlohilo[ix0][ix2][i];
         outputhilolo[0][i] = backwardhilolo[ix0][ix2][i];
         outputlololo[0][i] = backwardlololo[ix0][ix2][i];
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
               outputhihihi[k][i] = 0.0;
               outputlohihi[k][i] = 0.0;
               outputhilohi[k][i] = 0.0;
               outputlolohi[k][i] = 0.0;
               outputhihilo[k][i] = 0.0;
               outputlohilo[k][i] = 0.0;
               outputhilolo[k][i] = 0.0;
               outputlololo[k][i] = 0.0;
            }
         }
         else
         {
            if(verbose)
               cout << "updating derivative " << k 
                    << " with coefficient " << difidx << endl;

            for(int i=0; i<=deg; i++)
            {
               outputhihihi[k][i] = cffhihihi[difidx][i];
               outputlohihi[k][i] = cfflohihi[difidx][i];
               outputhilohi[k][i] = cffhilohi[difidx][i];
               outputlolohi[k][i] = cfflolohi[difidx][i];
               outputhihilo[k][i] = cffhihilo[difidx][i];
               outputlohilo[k][i] = cfflohilo[difidx][i];
               outputhilolo[k][i] = cffhilolo[difidx][i];
               outputlololo[k][i] = cfflololo[difidx][i];
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
               outputhihihi[k][i] = backwardhihihi[ix0][ix2][i];
               outputlohihi[k][i] = backwardlohihi[ix0][ix2][i];
               outputhilohi[k][i] = backwardhilohi[ix0][ix2][i];
               outputlolohi[k][i] = backwardlolohi[ix0][ix2][i];
               outputhihilo[k][i] = backwardhihilo[ix0][ix2][i];
               outputlohilo[k][i] = backwardlohilo[ix0][ix2][i];
               outputhilolo[k][i] = backwardhilolo[ix0][ix2][i];
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
               outputlohihi[k][i] = forwardlohihi[ix0][ix2][i];
               outputhilohi[k][i] = forwardhilohi[ix0][ix2][i];
               outputlolohi[k][i] = forwardlolohi[ix0][ix2][i];
               outputhihilo[k][i] = forwardhihilo[ix0][ix2][i];
               outputlohilo[k][i] = forwardlohilo[ix0][ix2][i];
               outputhilolo[k][i] = forwardhilolo[ix0][ix2][i];
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
               outputlohihi[k][i] = crosslohihi[ix0][ix2][i];
               outputhilohi[k][i] = crosshilohi[ix0][ix2][i];
               outputlolohi[k][i] = crosslolohi[ix0][ix2][i];
               outputhihilo[k][i] = crosshihilo[ix0][ix2][i];
               outputlohilo[k][i] = crosslohilo[ix0][ix2][i];
               outputhilolo[k][i] = crosshilolo[ix0][ix2][i];
               outputlololo[k][i] = crosslololo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_cmplx8_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
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
   double ***forwardrehihihi, double ***forwardrelohihi,
   double ***forwardrehilohi, double ***forwardrelolohi,
   double ***forwardrehihilo, double ***forwardrelohilo,
   double ***forwardrehilolo, double ***forwardrelololo,
   double ***forwardimhihihi, double ***forwardimlohihi,
   double ***forwardimhilohi, double ***forwardimlolohi,
   double ***forwardimhihilo, double ***forwardimlohilo,
   double ***forwardimhilolo, double ***forwardimlololo,
   double ***backwardrehihihi, double ***backwardrelohihi,
   double ***backwardrehilohi, double ***backwardrelolohi, 
   double ***backwardrehihilo, double ***backwardrelohilo,
   double ***backwardrehilolo, double ***backwardrelololo, 
   double ***backwardimhihihi, double ***backwardimlohihi,
   double ***backwardimhilohi, double ***backwardimlolohi, 
   double ***backwardimhihilo, double ***backwardimlohilo,
   double ***backwardimhilolo, double ***backwardimlololo, 
   double ***crossrehihihi, double ***crossrelohihi,
   double ***crossrehilohi, double ***crossrelolohi,
   double ***crossrehihilo, double ***crossrelohilo,
   double ***crossrehilolo, double ***crossrelololo,
   double ***crossimhihihi, double ***crossimlohihi,
   double ***crossimhilohi, double ***crossimlolohi,
   double ***crossimhihilo, double ***crossimlohilo,
   double ***crossimhilolo, double ***crossimlololo,
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

         CPU_cmplx8_add_job(deg,
            cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
            cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
            cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
            cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
            cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
            cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
            cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
            cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
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
            crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo,
            job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputrehihihi[dim][i] = forwardrehihihi[lastmon][lastidx][i];
      outputrelohihi[dim][i] = forwardrelohihi[lastmon][lastidx][i];
      outputrehilohi[dim][i] = forwardrehilohi[lastmon][lastidx][i];
      outputrelolohi[dim][i] = forwardrelolohi[lastmon][lastidx][i];
      outputrehihilo[dim][i] = forwardrehihilo[lastmon][lastidx][i];
      outputrelohilo[dim][i] = forwardrelohilo[lastmon][lastidx][i];
      outputrehilolo[dim][i] = forwardrehilolo[lastmon][lastidx][i];
      outputrelololo[dim][i] = forwardrelololo[lastmon][lastidx][i];
      outputimhihihi[dim][i] = forwardimhihihi[lastmon][lastidx][i];
      outputimlohihi[dim][i] = forwardimlohihi[lastmon][lastidx][i];
      outputimhilohi[dim][i] = forwardimhilohi[lastmon][lastidx][i];
      outputimlolohi[dim][i] = forwardimlolohi[lastmon][lastidx][i];
      outputimhihilo[dim][i] = forwardimhihilo[lastmon][lastidx][i];
      outputimlohilo[dim][i] = forwardimlohilo[lastmon][lastidx][i];
      outputimhilolo[dim][i] = forwardimhilolo[lastmon][lastidx][i];
      outputimlololo[dim][i] = forwardimlololo[lastmon][lastidx][i];
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
            outputrehihihi[0][i] = 0.0; outputrelohihi[0][i] = 0.0;
            outputrehilohi[0][i] = 0.0; outputrelolohi[0][i] = 0.0;
            outputrehihilo[0][i] = 0.0; outputrelohilo[0][i] = 0.0;
            outputrehilolo[0][i] = 0.0; outputrelololo[0][i] = 0.0;
            outputimhihihi[0][i] = 0.0; outputimlohihi[0][i] = 0.0; 
            outputimhilohi[0][i] = 0.0; outputimlolohi[0][i] = 0.0;
            outputimhihilo[0][i] = 0.0; outputimlohilo[0][i] = 0.0; 
            outputimhilolo[0][i] = 0.0; outputimlololo[0][i] = 0.0;
         }
      }
      else
      {
         if(verbose)
            cout << "updating derivative 0 with coefficient "
                 << difidx << endl;

         for(int i=0; i<=deg; i++)
         {
            outputrehihihi[0][i] = cffrehihihi[difidx][i];
            outputrelohihi[0][i] = cffrelohihi[difidx][i];
            outputrehilohi[0][i] = cffrehilohi[difidx][i];
            outputrelolohi[0][i] = cffrelolohi[difidx][i];
            outputrehihilo[0][i] = cffrehihilo[difidx][i];
            outputrelohilo[0][i] = cffrelohilo[difidx][i];
            outputrehilolo[0][i] = cffrehilolo[difidx][i];
            outputrelololo[0][i] = cffrelololo[difidx][i];
            outputimhihihi[0][i] = cffimhihihi[difidx][i];
            outputimlohihi[0][i] = cffimlohihi[difidx][i];
            outputimhilohi[0][i] = cffimhilohi[difidx][i];
            outputimlolohi[0][i] = cffimlolohi[difidx][i];
            outputimhihilo[0][i] = cffimhihilo[difidx][i];
            outputimlohilo[0][i] = cffimlohilo[difidx][i];
            outputimhilolo[0][i] = cffimhilolo[difidx][i];
            outputimlololo[0][i] = cffimlololo[difidx][i];
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
         outputrehihihi[0][i] = backwardrehihihi[ix0][ix2][i];
         outputrelohihi[0][i] = backwardrelohihi[ix0][ix2][i];
         outputrehilohi[0][i] = backwardrehilohi[ix0][ix2][i];
         outputrelolohi[0][i] = backwardrelolohi[ix0][ix2][i];
         outputrehihilo[0][i] = backwardrehihilo[ix0][ix2][i];
         outputrelohilo[0][i] = backwardrelohilo[ix0][ix2][i];
         outputrehilolo[0][i] = backwardrehilolo[ix0][ix2][i];
         outputrelololo[0][i] = backwardrelololo[ix0][ix2][i];
         outputimhihihi[0][i] = backwardimhihihi[ix0][ix2][i];
         outputimlohihi[0][i] = backwardimlohihi[ix0][ix2][i];
         outputimhilohi[0][i] = backwardimhilohi[ix0][ix2][i];
         outputimlolohi[0][i] = backwardimlolohi[ix0][ix2][i];
         outputimhihilo[0][i] = backwardimhihilo[ix0][ix2][i];
         outputimlohilo[0][i] = backwardimlohilo[ix0][ix2][i];
         outputimhilolo[0][i] = backwardimhilolo[ix0][ix2][i];
         outputimlololo[0][i] = backwardimlololo[ix0][ix2][i];
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
               outputrehihihi[k][i] = 0.0; outputrelohihi[k][i] = 0.0;
               outputrehilohi[k][i] = 0.0; outputrelolohi[k][i] = 0.0;
               outputrehihilo[k][i] = 0.0; outputrelohilo[k][i] = 0.0;
               outputrehilolo[k][i] = 0.0; outputrelololo[k][i] = 0.0;
               outputimhihihi[k][i] = 0.0; outputimlohihi[k][i] = 0.0; 
               outputimhilohi[k][i] = 0.0; outputimlolohi[k][i] = 0.0;
               outputimhihilo[k][i] = 0.0; outputimlohilo[k][i] = 0.0; 
               outputimhilolo[k][i] = 0.0; outputimlololo[k][i] = 0.0;
            }
         }
         else
         {
            if(verbose)
               cout << "updating derivative " << k 
                    << " with coefficient " << difidx << endl;

            for(int i=0; i<=deg; i++)
            {
               outputrehihihi[k][i] = cffrehihihi[difidx][i];
               outputrelohihi[k][i] = cffrelohihi[difidx][i];
               outputrehilohi[k][i] = cffrehilohi[difidx][i];
               outputrelolohi[k][i] = cffrelolohi[difidx][i];
               outputrehihilo[k][i] = cffrehihilo[difidx][i];
               outputrelohilo[k][i] = cffrelohilo[difidx][i];
               outputrehilolo[k][i] = cffrehilolo[difidx][i];
               outputrelololo[k][i] = cffrelololo[difidx][i];
               outputimhihihi[k][i] = cffimhihihi[difidx][i];
               outputimlohihi[k][i] = cffimlohihi[difidx][i];
               outputimhilohi[k][i] = cffimhilohi[difidx][i];
               outputimlolohi[k][i] = cffimlolohi[difidx][i];
               outputimhihilo[k][i] = cffimhihilo[difidx][i];
               outputimlohilo[k][i] = cffimlohilo[difidx][i];
               outputimhilolo[k][i] = cffimhilolo[difidx][i];
               outputimlololo[k][i] = cffimlololo[difidx][i];
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
               outputrehihihi[k][i] = backwardrehihihi[ix0][ix2][i];
               outputrelohihi[k][i] = backwardrelohihi[ix0][ix2][i];
               outputrehilohi[k][i] = backwardrehilohi[ix0][ix2][i];
               outputrelolohi[k][i] = backwardrelolohi[ix0][ix2][i];
               outputrehihilo[k][i] = backwardrehihilo[ix0][ix2][i];
               outputrelohilo[k][i] = backwardrelohilo[ix0][ix2][i];
               outputrehilolo[k][i] = backwardrehilolo[ix0][ix2][i];
               outputrelololo[k][i] = backwardrelololo[ix0][ix2][i];
               outputimhihihi[k][i] = backwardimhihihi[ix0][ix2][i];
               outputimlohihi[k][i] = backwardimlohihi[ix0][ix2][i];
               outputimhilohi[k][i] = backwardimhilohi[ix0][ix2][i];
               outputimlolohi[k][i] = backwardimlolohi[ix0][ix2][i];
               outputimhihilo[k][i] = backwardimhihilo[ix0][ix2][i];
               outputimlohilo[k][i] = backwardimlohilo[ix0][ix2][i];
               outputimhilolo[k][i] = backwardimhilolo[ix0][ix2][i];
               outputimlololo[k][i] = backwardimlololo[ix0][ix2][i];
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
               outputrehihihi[k][i] = forwardrehihihi[ix0][ix2][i];
               outputrelohihi[k][i] = forwardrelohihi[ix0][ix2][i];
               outputrehilohi[k][i] = forwardrehilohi[ix0][ix2][i];
               outputrelolohi[k][i] = forwardrelolohi[ix0][ix2][i];
               outputrehihilo[k][i] = forwardrehihilo[ix0][ix2][i];
               outputrelohilo[k][i] = forwardrelohilo[ix0][ix2][i];
               outputrehilolo[k][i] = forwardrehilolo[ix0][ix2][i];
               outputrelololo[k][i] = forwardrelololo[ix0][ix2][i];
               outputimhihihi[k][i] = forwardimhihihi[ix0][ix2][i];
               outputimlohihi[k][i] = forwardimlohihi[ix0][ix2][i];
               outputimhilohi[k][i] = forwardimhilohi[ix0][ix2][i];
               outputimlolohi[k][i] = forwardimlolohi[ix0][ix2][i];
               outputimhihilo[k][i] = forwardimhihilo[ix0][ix2][i];
               outputimlohilo[k][i] = forwardimlohilo[ix0][ix2][i];
               outputimhilolo[k][i] = forwardimhilolo[ix0][ix2][i];
               outputimlololo[k][i] = forwardimlololo[ix0][ix2][i];
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
               outputrehihihi[k][i] = crossrehihihi[ix0][ix2][i];
               outputrelohihi[k][i] = crossrelohihi[ix0][ix2][i];
               outputrehilohi[k][i] = crossrehilohi[ix0][ix2][i];
               outputrelolohi[k][i] = crossrelolohi[ix0][ix2][i];
               outputrehihilo[k][i] = crossrehihilo[ix0][ix2][i];
               outputrelohilo[k][i] = crossrelohilo[ix0][ix2][i];
               outputrehilolo[k][i] = crossrehilolo[ix0][ix2][i];
               outputrelololo[k][i] = crossrelololo[ix0][ix2][i];
               outputimhihihi[k][i] = crossimhihihi[ix0][ix2][i];
               outputimlohihi[k][i] = crossimlohihi[ix0][ix2][i];
               outputimhilohi[k][i] = crossimhilohi[ix0][ix2][i];
               outputimlolohi[k][i] = crossimlolohi[ix0][ix2][i];
               outputimhihilo[k][i] = crossimhihilo[ix0][ix2][i];
               outputimlohilo[k][i] = crossimlohilo[ix0][ix2][i];
               outputimhilolo[k][i] = crossimhilolo[ix0][ix2][i];
               outputimlololo[k][i] = crossimlololo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl8_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
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
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double ***forwardhihihi = new double**[nbr];
   double ***forwardlohihi = new double**[nbr];
   double ***forwardhilohi = new double**[nbr];
   double ***forwardlolohi = new double**[nbr];
   double ***forwardhihilo = new double**[nbr];
   double ***forwardlohilo = new double**[nbr];
   double ***forwardhilolo = new double**[nbr];
   double ***forwardlololo = new double**[nbr];
   double ***backwardhihihi = new double**[nbr];
   double ***backwardlohihi = new double**[nbr];
   double ***backwardhilohi = new double**[nbr];
   double ***backwardlolohi = new double**[nbr];
   double ***backwardhihilo = new double**[nbr];
   double ***backwardlohilo = new double**[nbr];
   double ***backwardhilolo = new double**[nbr];
   double ***backwardlololo = new double**[nbr];
   double ***crosshihihi = new double**[nbr];
   double ***crosslohihi = new double**[nbr];
   double ***crosshilohi = new double**[nbr];
   double ***crosslolohi = new double**[nbr];
   double ***crosshihilo = new double**[nbr];
   double ***crosslohilo = new double**[nbr];
   double ***crosshilolo = new double**[nbr];
   double ***crosslololo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardhihihi[k] = new double*[nvrk];
      forwardlohihi[k] = new double*[nvrk];
      forwardhilohi[k] = new double*[nvrk];
      forwardlolohi[k] = new double*[nvrk];
      forwardhihilo[k] = new double*[nvrk];
      forwardlohilo[k] = new double*[nvrk];
      forwardhilolo[k] = new double*[nvrk];
      forwardlololo[k] = new double*[nvrk];

      for(int i=0; i<nvrk; i++) 
      {
         forwardhihihi[k][i] = new double[deg+1];
         forwardlohihi[k][i] = new double[deg+1];
         forwardhilohi[k][i] = new double[deg+1];
         forwardlolohi[k][i] = new double[deg+1];
         forwardhihilo[k][i] = new double[deg+1];
         forwardlohilo[k][i] = new double[deg+1];
         forwardhilolo[k][i] = new double[deg+1];
         forwardlololo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardhihihi[k] = new double*[nvrk-1];
         backwardlohihi[k] = new double*[nvrk-1];
         backwardhilohi[k] = new double*[nvrk-1];
         backwardlolohi[k] = new double*[nvrk-1];
         backwardhihilo[k] = new double*[nvrk-1];
         backwardlohilo[k] = new double*[nvrk-1];
         backwardhilolo[k] = new double*[nvrk-1];
         backwardlololo[k] = new double*[nvrk-1];

         for(int i=0; i<nvrk-1; i++) 
         {
            backwardhihihi[k][i] = new double[deg+1];
            backwardlohihi[k][i] = new double[deg+1];
            backwardhilohi[k][i] = new double[deg+1];
            backwardlolohi[k][i] = new double[deg+1];
            backwardhihilo[k][i] = new double[deg+1];
            backwardlohilo[k][i] = new double[deg+1];
            backwardhilolo[k][i] = new double[deg+1];
            backwardlololo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crosshihihi[k] = new double*[nvrk-2];
         crosslohihi[k] = new double*[nvrk-2];
         crosshilohi[k] = new double*[nvrk-2];
         crosslolohi[k] = new double*[nvrk-2];
         crosshihilo[k] = new double*[nvrk-2];
         crosslohilo[k] = new double*[nvrk-2];
         crosshilolo[k] = new double*[nvrk-2];
         crosslololo[k] = new double*[nvrk-2];

         for(int i=0; i<nvrk-2; i++)
         {
            crosshihihi[k][i] = new double[deg+1];
            crosslohihi[k][i] = new double[deg+1];
            crosshilohi[k][i] = new double[deg+1];
            crosslolohi[k][i] = new double[deg+1];
            crosshihilo[k][i] = new double[deg+1];
            crosslohilo[k][i] = new double[deg+1];
            crosshilolo[k][i] = new double[deg+1];
            crosslololo[k][i] = new double[deg+1];
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

         CPU_dbl8_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffhihihi[monidx],cfflohihi[monidx],
             cffhilohi[monidx],cfflolohi[monidx],
             cffhihilo[monidx],cfflohilo[monidx],
             cffhilolo[monidx],cfflololo[monidx],
             inputhihihi,inputlohihi,inputhilohi,inputlolohi,
             inputhihilo,inputlohilo,inputhilolo,inputlololo,
              forwardhihihi[monidx], forwardlohihi[monidx],
              forwardhilohi[monidx], forwardlolohi[monidx],
              forwardhihilo[monidx], forwardlohilo[monidx],
              forwardhilolo[monidx], forwardlololo[monidx],
             backwardhihihi[monidx],backwardlohihi[monidx],
             backwardhilohi[monidx],backwardlolohi[monidx],
             backwardhihilo[monidx],backwardlohilo[monidx],
             backwardhilolo[monidx],backwardlololo[monidx],
                crosshihihi[monidx],   crosslohihi[monidx],
                crosshilohi[monidx],   crosslolohi[monidx],
                crosshihilo[monidx],   crosslohilo[monidx],
                crosshilolo[monidx],   crosslololo[monidx],job,vrb);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl8_poly_addjobs
      (dim,nbr,deg,nvr,idx,
            csthihihi,     cstlohihi,     csthilohi,     cstlolohi,
            csthihilo,     cstlohilo,     csthilolo,     cstlololo,
            cffhihihi,     cfflohihi,     cffhilohi,     cfflolohi,
            cffhihilo,     cfflohilo,     cffhilolo,     cfflololo,
          inputhihihi,   inputlohihi,   inputhilohi,   inputlolohi,
          inputhihilo,   inputlohilo,   inputhilolo,   inputlololo,
         outputhihihi,  outputlohihi,  outputhilohi,  outputlolohi,
         outputhihilo,  outputlohilo,  outputhilolo,  outputlololo,
        forwardhihihi, forwardlohihi, forwardhilohi, forwardlolohi,
        forwardhihilo, forwardlohilo, forwardhilolo, forwardlololo,
       backwardhihihi,backwardlohihi,backwardhilohi,backwardlolohi,
       backwardhihilo,backwardlohilo,backwardhilolo,backwardlololo,
          crosshihihi,   crosslohihi,   crosshilohi,   crosslolohi,
          crosshihilo,   crosslohilo,   crosshilolo,   crosslololo,
       addjobs,vrb);
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
         free(forwardhihihi[k][i]);
         free(forwardlohihi[k][i]);
         free(forwardhilohi[k][i]);
         free(forwardlolohi[k][i]);
         free(forwardhihilo[k][i]);
         free(forwardlohilo[k][i]);
         free(forwardhilolo[k][i]);
         free(forwardlololo[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardhihihi[k][i]);
                      free(backwardlohihi[k][i]);
                      free(backwardhilohi[k][i]);
                      free(backwardlolohi[k][i]);
                      free(backwardhihilo[k][i]);
                      free(backwardlohilo[k][i]);
                      free(backwardhilolo[k][i]);
                      free(backwardlololo[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crosshihihi[k][i]);
                      free(crosslohihi[k][i]);
                      free(crosshilohi[k][i]);
                      free(crosslolohi[k][i]);
                      free(crosshihilo[k][i]);
                      free(crosslohilo[k][i]);
                      free(crosshilolo[k][i]);
                      free(crosslololo[k][i]);
                   }
   }
   free(forwardhihihi); free(backwardhihihi); free(crosshihihi);
   free(forwardlohihi); free(backwardlohihi); free(crosslohihi);
   free(forwardhilohi); free(backwardhilohi); free(crosshilohi);
   free(forwardlolohi); free(backwardlolohi); free(crosslolohi);
   free(forwardhihilo); free(backwardhihilo); free(crosshihilo);
   free(forwardlohilo); free(backwardlohilo); free(crosslohilo);
   free(forwardhilolo); free(backwardhilolo); free(crosshilolo);
   free(forwardlololo); free(backwardlololo); free(crosslololo);
}

void CPU_cmplx8_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
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
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double ***forwardrehihihi = new double**[nbr];
   double ***forwardrelohihi = new double**[nbr];
   double ***forwardrehilohi = new double**[nbr];
   double ***forwardrelolohi = new double**[nbr];
   double ***forwardrehihilo = new double**[nbr];
   double ***forwardrelohilo = new double**[nbr];
   double ***forwardrehilolo = new double**[nbr];
   double ***forwardrelololo = new double**[nbr];
   double ***forwardimhihihi = new double**[nbr];
   double ***forwardimlohihi = new double**[nbr];
   double ***forwardimhilohi = new double**[nbr];
   double ***forwardimlolohi = new double**[nbr];
   double ***forwardimhihilo = new double**[nbr];
   double ***forwardimlohilo = new double**[nbr];
   double ***forwardimhilolo = new double**[nbr];
   double ***forwardimlololo = new double**[nbr];
   double ***backwardrehihihi = new double**[nbr];
   double ***backwardrelohihi = new double**[nbr];
   double ***backwardrehilohi = new double**[nbr];
   double ***backwardrelolohi = new double**[nbr];
   double ***backwardrehihilo = new double**[nbr];
   double ***backwardrelohilo = new double**[nbr];
   double ***backwardrehilolo = new double**[nbr];
   double ***backwardrelololo = new double**[nbr];
   double ***backwardimhihihi = new double**[nbr];
   double ***backwardimlohihi = new double**[nbr];
   double ***backwardimhilohi = new double**[nbr];
   double ***backwardimlolohi = new double**[nbr];
   double ***backwardimhihilo = new double**[nbr];
   double ***backwardimlohilo = new double**[nbr];
   double ***backwardimhilolo = new double**[nbr];
   double ***backwardimlololo = new double**[nbr];
   double ***crossrehihihi = new double**[nbr];
   double ***crossrelohihi = new double**[nbr];
   double ***crossrehilohi = new double**[nbr];
   double ***crossrelolohi = new double**[nbr];
   double ***crossrehihilo = new double**[nbr];
   double ***crossrelohilo = new double**[nbr];
   double ***crossrehilolo = new double**[nbr];
   double ***crossrelololo = new double**[nbr];
   double ***crossimhihihi = new double**[nbr];
   double ***crossimlohihi = new double**[nbr];
   double ***crossimhilohi = new double**[nbr];
   double ***crossimlolohi = new double**[nbr];
   double ***crossimhihilo = new double**[nbr];
   double ***crossimlohilo = new double**[nbr];
   double ***crossimhilolo = new double**[nbr];
   double ***crossimlololo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardrehihihi[k] = new double*[nvrk];
      forwardrelohihi[k] = new double*[nvrk];
      forwardrehilohi[k] = new double*[nvrk];
      forwardrelolohi[k] = new double*[nvrk];
      forwardrehihilo[k] = new double*[nvrk];
      forwardrelohilo[k] = new double*[nvrk];
      forwardrehilolo[k] = new double*[nvrk];
      forwardrelololo[k] = new double*[nvrk];
      forwardimhihihi[k] = new double*[nvrk];
      forwardimlohihi[k] = new double*[nvrk];
      forwardimhilohi[k] = new double*[nvrk];
      forwardimlolohi[k] = new double*[nvrk];
      forwardimhihilo[k] = new double*[nvrk];
      forwardimlohilo[k] = new double*[nvrk];
      forwardimhilolo[k] = new double*[nvrk];
      forwardimlololo[k] = new double*[nvrk];

      for(int i=0; i<nvrk; i++) 
      {
         forwardrehihihi[k][i] = new double[deg+1];
         forwardrelohihi[k][i] = new double[deg+1];
         forwardrehilohi[k][i] = new double[deg+1];
         forwardrelolohi[k][i] = new double[deg+1];
         forwardrehihilo[k][i] = new double[deg+1];
         forwardrelohilo[k][i] = new double[deg+1];
         forwardrehilolo[k][i] = new double[deg+1];
         forwardrelololo[k][i] = new double[deg+1];
         forwardimhihihi[k][i] = new double[deg+1];
         forwardimlohihi[k][i] = new double[deg+1];
         forwardimhilohi[k][i] = new double[deg+1];
         forwardimlolohi[k][i] = new double[deg+1];
         forwardimhihilo[k][i] = new double[deg+1];
         forwardimlohilo[k][i] = new double[deg+1];
         forwardimhilolo[k][i] = new double[deg+1];
         forwardimlololo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardrehihihi[k] = new double*[nvrk-1];
         backwardrelohihi[k] = new double*[nvrk-1];
         backwardrehilohi[k] = new double*[nvrk-1];
         backwardrelolohi[k] = new double*[nvrk-1];
         backwardrehihilo[k] = new double*[nvrk-1];
         backwardrelohilo[k] = new double*[nvrk-1];
         backwardrehilolo[k] = new double*[nvrk-1];
         backwardrelololo[k] = new double*[nvrk-1];
         backwardimhihihi[k] = new double*[nvrk-1];
         backwardimlohihi[k] = new double*[nvrk-1];
         backwardimhilohi[k] = new double*[nvrk-1];
         backwardimlolohi[k] = new double*[nvrk-1];
         backwardimhihilo[k] = new double*[nvrk-1];
         backwardimlohilo[k] = new double*[nvrk-1];
         backwardimhilolo[k] = new double*[nvrk-1];
         backwardimlololo[k] = new double*[nvrk-1];

         for(int i=0; i<nvrk-1; i++) 
         {
            backwardrehihihi[k][i] = new double[deg+1];
            backwardrelohihi[k][i] = new double[deg+1];
            backwardrehilohi[k][i] = new double[deg+1];
            backwardrelolohi[k][i] = new double[deg+1];
            backwardrehihilo[k][i] = new double[deg+1];
            backwardrelohilo[k][i] = new double[deg+1];
            backwardrehilolo[k][i] = new double[deg+1];
            backwardrelololo[k][i] = new double[deg+1];
            backwardimhihihi[k][i] = new double[deg+1];
            backwardimlohihi[k][i] = new double[deg+1];
            backwardimhilohi[k][i] = new double[deg+1];
            backwardimlolohi[k][i] = new double[deg+1];
            backwardimhihilo[k][i] = new double[deg+1];
            backwardimlohilo[k][i] = new double[deg+1];
            backwardimhilolo[k][i] = new double[deg+1];
            backwardimlololo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crossrehihihi[k] = new double*[nvrk-2];
         crossrelohihi[k] = new double*[nvrk-2];
         crossrehilohi[k] = new double*[nvrk-2];
         crossrelolohi[k] = new double*[nvrk-2];
         crossrehihilo[k] = new double*[nvrk-2];
         crossrelohilo[k] = new double*[nvrk-2];
         crossrehilolo[k] = new double*[nvrk-2];
         crossrelololo[k] = new double*[nvrk-2];
         crossimhihihi[k] = new double*[nvrk-2];
         crossimlohihi[k] = new double*[nvrk-2];
         crossimhilohi[k] = new double*[nvrk-2];
         crossimlolohi[k] = new double*[nvrk-2];
         crossimhihilo[k] = new double*[nvrk-2];
         crossimlohilo[k] = new double*[nvrk-2];
         crossimhilolo[k] = new double*[nvrk-2];
         crossimlololo[k] = new double*[nvrk-2];

         for(int i=0; i<nvrk-2; i++)
         {
            crossrehihihi[k][i] = new double[deg+1];
            crossrelohihi[k][i] = new double[deg+1];
            crossrehilohi[k][i] = new double[deg+1];
            crossrelolohi[k][i] = new double[deg+1];
            crossrehihilo[k][i] = new double[deg+1];
            crossrelohilo[k][i] = new double[deg+1];
            crossrehilolo[k][i] = new double[deg+1];
            crossrelololo[k][i] = new double[deg+1];
            crossimhihihi[k][i] = new double[deg+1];
            crossimlohihi[k][i] = new double[deg+1];
            crossimhilohi[k][i] = new double[deg+1];
            crossimlolohi[k][i] = new double[deg+1];
            crossimhihilo[k][i] = new double[deg+1];
            crossimlohilo[k][i] = new double[deg+1];
            crossimhilolo[k][i] = new double[deg+1];
            crossimlololo[k][i] = new double[deg+1];
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

         CPU_cmplx8_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffrehihihi[monidx],cffrelohihi[monidx],
             cffrehilohi[monidx],cffrelolohi[monidx],
             cffrehihilo[monidx],cffrelohilo[monidx],
             cffrehilolo[monidx],cffrelololo[monidx],
             cffimhihihi[monidx],cffimlohihi[monidx],
             cffimhilohi[monidx],cffimlolohi[monidx],
             cffimhihilo[monidx],cffimlohilo[monidx],
             cffimhilolo[monidx],cffimlololo[monidx],
             inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
             inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
             inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
             inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
             forwardrehihihi[monidx],forwardrelohihi[monidx],
             forwardrehilohi[monidx],forwardrelolohi[monidx],
             forwardrehihilo[monidx],forwardrelohilo[monidx],
             forwardrehilolo[monidx],forwardrelololo[monidx],
             forwardimhihihi[monidx],forwardimlohihi[monidx],
             forwardimhilohi[monidx],forwardimlolohi[monidx],
             forwardimhihilo[monidx],forwardimlohilo[monidx],
             forwardimhilolo[monidx],forwardimlololo[monidx],
             backwardrehihihi[monidx],backwardrelohihi[monidx],
             backwardrehilohi[monidx],backwardrelolohi[monidx],
             backwardrehihilo[monidx],backwardrelohilo[monidx],
             backwardrehilolo[monidx],backwardrelololo[monidx],
             backwardimhihihi[monidx],backwardimlohihi[monidx],
             backwardimhilohi[monidx],backwardimlolohi[monidx],
             backwardimhihilo[monidx],backwardimlohilo[monidx],
             backwardimhilolo[monidx],backwardimlololo[monidx],
             crossrehihihi[monidx],crossrelohihi[monidx],
             crossrehilohi[monidx],crossrelolohi[monidx],
             crossrehihilo[monidx],crossrelohilo[monidx],
             crossrehilolo[monidx],crossrelololo[monidx],
             crossimhihihi[monidx],crossimlohihi[monidx],
             crossimhilohi[monidx],crossimlolohi[monidx],
             crossimhihilo[monidx],crossimlohilo[monidx],
             crossimhilolo[monidx],crossimlololo[monidx],job,vrb);
      }
   }
 /*
   CPU_cmplx8_poly_updates
      (dim,nbr,deg,nvr,idx,
       cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
       cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
       cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
       cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
       inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
       inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
       inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
       outputrehihihi,outputrelohihi,outputrehilohi,outputrelolohi,
       outputrehihilo,outputrelohilo,outputrehilolo,outputrelololo,
       outputimhihihi,outputimlohihi,outputimhilohi,outputimlolohi,
       outputimhihilo,outputimlohilo,outputimhilolo,outputimlololo,
       forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
       forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
       forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
       forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
       backwardrehihihi,backwardrelohihi,backwardrehilohi,backwardrelolohi,
       backwardrehihilo,backwardrelohilo,backwardrehilolo,backwardrelololo,
       backwardimhihihi,backwardimlohihi,backwardimhilohi,backwardimlolohi,
       backwardimhihilo,backwardimlohilo,backwardimhilolo,backwardimlololo,
       crossrehihihi,crossrelohihi,crossrehilohi,crossrelolohi,
       crossrehihilo,crossrelohilo,crossrehilolo,crossrelololo,
       crossimhihihi,crossimlohihi,crossimhilohi,crossimlolohi,
       crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo);
  */
   CPU_cmplx8_poly_addjobs
      (dim,nbr,deg,nvr,idx,
       cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
       cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
       cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
       cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
       inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
       inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
       inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
       outputrehihihi,outputrelohihi,outputrehilohi,outputrelolohi,
       outputrehihilo,outputrelohilo,outputrehilolo,outputrelololo,
       outputimhihihi,outputimlohihi,outputimhilohi,outputimlolohi,
       outputimhihilo,outputimlohilo,outputimhilolo,outputimlololo,
       forwardrehihihi,forwardrelohihi,forwardrehilohi,forwardrelolohi,
       forwardrehihilo,forwardrelohilo,forwardrehilolo,forwardrelololo,
       forwardimhihihi,forwardimlohihi,forwardimhilohi,forwardimlolohi,
       forwardimhihilo,forwardimlohilo,forwardimhilolo,forwardimlololo,
       backwardrehihihi,backwardrelohihi,backwardrehilohi,backwardrelolohi,
       backwardrehihilo,backwardrelohilo,backwardrehilolo,backwardrelololo,
       backwardimhihihi,backwardimlohihi,backwardimhilohi,backwardimlolohi,
       backwardimhihilo,backwardimlohilo,backwardimhilolo,backwardimlololo,
       crossrehihihi,crossrelohihi,crossrehilohi,crossrelolohi,
       crossrehihilo,crossrelohilo,crossrehilolo,crossrelololo,
       crossimhihihi,crossimlohihi,crossimhilohi,crossimlolohi,
       crossimhihilo,crossimlohilo,crossimhilolo,crossimlololo,
       addjobs,vrb);
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
         free(forwardrehihihi[k][i]); free(forwardrelohihi[k][i]);
         free(forwardrehilohi[k][i]); free(forwardrelolohi[k][i]);
         free(forwardrehihilo[k][i]); free(forwardrelohilo[k][i]);
         free(forwardrehilolo[k][i]); free(forwardrelololo[k][i]);
         free(forwardimhihihi[k][i]); free(forwardimlohihi[k][i]);
         free(forwardimhilohi[k][i]); free(forwardimlolohi[k][i]);
         free(forwardimhihilo[k][i]); free(forwardimlohilo[k][i]);
         free(forwardimhilolo[k][i]); free(forwardimlololo[k][i]);
      }
      if(nvrk > 1)
         for(int i=0; i<nvrk-1; i++)
         {
            free(backwardrehihihi[k][i]); free(backwardrelohihi[k][i]);
            free(backwardrehilohi[k][i]); free(backwardrelolohi[k][i]);
            free(backwardrehihilo[k][i]); free(backwardrelohilo[k][i]);
            free(backwardrehilolo[k][i]); free(backwardrelololo[k][i]);
            free(backwardimhihihi[k][i]); free(backwardimlohihi[k][i]);
            free(backwardimhilohi[k][i]); free(backwardimlolohi[k][i]);
            free(backwardimhihilo[k][i]); free(backwardimlohilo[k][i]);
            free(backwardimhilolo[k][i]); free(backwardimlololo[k][i]);
         }
      if(nvrk > 2)
         for(int i=0; i<nvrk-2; i++)
         {
            free(crossrehihihi[k][i]); free(crossrelohihi[k][i]);
            free(crossrehilohi[k][i]); free(crossrelolohi[k][i]);
            free(crossrehihilo[k][i]); free(crossrelohilo[k][i]);
            free(crossrehilolo[k][i]); free(crossrelololo[k][i]);
            free(crossimhihihi[k][i]); free(crossimlohihi[k][i]);
            free(crossimhilohi[k][i]); free(crossimlolohi[k][i]);
            free(crossimhihilo[k][i]); free(crossimlohilo[k][i]);
            free(crossimhilolo[k][i]); free(crossimlololo[k][i]);
         }
   }
   free(forwardrehihihi); free(backwardrehihihi); free(crossrehihihi);
   free(forwardrelohihi); free(backwardrelohihi); free(crossrelohihi);
   free(forwardrehilohi); free(backwardrehilohi); free(crossrehilohi);
   free(forwardrelolohi); free(backwardrelolohi); free(crossrelolohi);
   free(forwardrehihilo); free(backwardrehihilo); free(crossrehihilo);
   free(forwardrelohilo); free(backwardrelohilo); free(crossrelohilo);
   free(forwardrehilolo); free(backwardrehilolo); free(crossrehilolo);
   free(forwardrelololo); free(backwardrelololo); free(crossrelololo);
   free(forwardimhihihi); free(backwardimhihihi); free(crossimhihihi);
   free(forwardimlohihi); free(backwardimlohihi); free(crossimlohihi);
   free(forwardimhilohi); free(backwardimhilohi); free(crossimhilohi);
   free(forwardimlolohi); free(backwardimlolohi); free(crossimlolohi);
   free(forwardimhihilo); free(backwardimhihilo); free(crossimhihilo);
   free(forwardimlohilo); free(backwardimlohilo); free(crossimlohilo);
   free(forwardimhilolo); free(backwardimhilolo); free(crossimhilolo);
   free(forwardimlololo); free(backwardimlololo); free(crossimlololo);
}
