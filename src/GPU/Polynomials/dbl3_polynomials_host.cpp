/* The file dbl3_polynomials_host.cpp defines functions specified
 * in dbl3_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "triple_double_functions.h"
#include "dbl3_convolutions_host.h"
#include "dbl3_monomials_host.h"
#include "dbl3_polynomials_host.h"

void CPU_dbl3_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   double **forwardhi, double **forwardmi, double **forwardlo,
   double **backwardhi, double **backwardmi, double **backwardlo,
   double **crosshi, double **crossmi, double **crosslo, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl3_product(deg,inputhi[ix1],inputmi[ix1],inputlo[ix1],
                              cffhi[i],cffmi[i],cfflo[i],
                              forwardhi[0],forwardmi[0],forwardlo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            tdf_inc(&outputhi[dim][j],&outputmi[dim][j],&outputlo[dim][j],
                    forwardhi[0][j],forwardmi[0][j],forwardlo[0][j]);
            // output[ix1][j] += cff[i][j];
            tdf_inc(&outputhi[ix1][j],&outputmi[ix1][j],&outputlo[ix1][j],
                    cffhi[i][j],cffmi[i][j],cfflo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl3_product(deg,cffhi[i],cffmi[i],cfflo[i],
                              inputhi[ix1],inputmi[ix1],inputlo[ix1],
                              forwardhi[0],forwardmi[0],forwardlo[0]);
         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
            tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                    forwardhi[0][j],forwardmi[0][j],forwardlo[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl3_product(deg,cffhi[i],cffmi[i],cfflo[i],
                              inputhi[ix2],inputmi[ix2],inputlo[ix2],
                              backwardhi[0],backwardmi[0],backwardlo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
            tdf_inc(&outputhi[ix1][j],&outputmi[ix1][j],&outputlo[ix1][j],
                    backwardhi[0][j],backwardmi[0][j],backwardlo[0][j]);

         CPU_dbl3_product(deg,forwardhi[0],forwardmi[0],forwardlo[0],
                              inputhi[ix2],inputmi[ix2],inputlo[ix2],
                              forwardhi[1],forwardmi[1],forwardlo[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
            tdf_inc(&outputhi[dim][j],&outputmi[dim][j],&outputlo[dim][j],
                    forwardhi[1][j],forwardmi[1][j],forwardlo[1][j]);
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl3_speel
            (nvr[i],deg,idx[i],cffhi[i],cffmi[i],cfflo[i],
             inputhi,inputmi,inputlo,forwardhi,forwardmi,forwardlo,
             backwardhi,backwardmi,backwardlo,crosshi,crossmi,crosslo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
            tdf_inc(&outputhi[dim][j],&outputmi[dim][j],&outputlo[dim][j],
                    forwardhi[ix1][j],forwardmi[ix1][j],forwardlo[ix1][j]);

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
            tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                    forwardhi[ix1][j],forwardmi[ix1][j],forwardlo[ix1][j]);

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
            tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                    backwardhi[ix1][j],backwardmi[ix1][j],backwardlo[ix1][j]);

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
               tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                       crosshi[k-1][j],crossmi[k-1][j],crosslo[k-1][j]);
         }
      }
   }
}

void CPU_cmplx3_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
   double **forwardrehi, double **forwardremi, double **forwardrelo,
   double **forwardimhi, double **forwardimmi, double **forwardimlo,
   double **backwardrehi, double **backwardremi, double **backwardrelo,
   double **backwardimhi, double **backwardimmi, double **backwardimlo,
   double **crossrehi, double **crossremi, double **crossrelo,
   double **crossimhi, double **crossimmi, double **crossimlo, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_cmplx3_product(deg,
            inputrehi[ix1],inputremi[ix1],inputrelo[ix1],
            inputimhi[ix1],inputimmi[ix1],inputimlo[ix1],
            cffrehi[i],cffremi[i],cffrelo[i],cffimhi[i],cffimmi[i],cffimlo[i],
            forwardrehi[0],forwardremi[0],forwardrelo[0],
            forwardimhi[0],forwardimmi[0],forwardimlo[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            tdf_inc
               (&outputrehi[dim][j],&outputremi[dim][j],&outputrelo[dim][j],
                forwardrehi[0][j],forwardremi[0][j],forwardrelo[0][j]);
            tdf_inc
               (&outputimhi[dim][j],&outputimmi[dim][j],&outputimlo[dim][j],
                forwardimhi[0][j],forwardimmi[0][j],forwardimlo[0][j]);
            // output[ix1][j] += cff[i][j];
            tdf_inc
               (&outputrehi[ix1][j],&outputremi[ix1][j],&outputrelo[ix1][j],
                cffrehi[i][j],cffremi[i][j],cffrelo[i][j]);
            tdf_inc
               (&outputimhi[ix1][j],&outputimmi[ix1][j],&outputimlo[ix1][j],
                cffimhi[i][j],cffimmi[i][j],cffimlo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_cmplx3_product(deg,
            cffrehi[i],cffremi[i],cffrelo[i],cffimhi[i],cffimmi[i],cffimlo[i],
            inputrehi[ix1],inputremi[ix1],inputrelo[ix1],
            inputimhi[ix1],inputimmi[ix1],inputimlo[ix1],
            forwardrehi[0],forwardremi[0],forwardrelo[0],
            forwardimhi[0],forwardimmi[0],forwardimlo[0]);

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
         {
            tdf_inc
               (&outputrehi[ix2][j],&outputremi[ix2][j],&outputrelo[ix2][j],
                forwardrehi[0][j],forwardremi[0][j],forwardrelo[0][j]);
            tdf_inc
               (&outputimhi[ix2][j],&outputimmi[ix2][j],&outputimlo[ix2][j],
                forwardimhi[0][j],forwardimmi[0][j],forwardimlo[0][j]);
         }
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_cmplx3_product(deg,
            cffrehi[i],cffremi[i],cffrelo[i],cffimhi[i],cffimmi[i],cffimlo[i],
            inputrehi[ix2],inputremi[ix2],inputrelo[ix2],
            inputimhi[ix2],inputimmi[ix2],inputimlo[ix2],
            backwardrehi[0],backwardremi[0],backwardrelo[0],
            backwardimhi[0],backwardimmi[0],backwardimlo[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
         {
            tdf_inc
               (&outputrehi[ix1][j],&outputremi[ix1][j],&outputrelo[ix1][j],
                backwardrehi[0][j],backwardremi[0][j],backwardrelo[0][j]);
            tdf_inc
               (&outputimhi[ix1][j],&outputimmi[ix1][j],&outputimlo[ix1][j],
                backwardimhi[0][j],backwardimmi[0][j],backwardimlo[0][j]);
         }
         CPU_cmplx3_product(deg,
            forwardrehi[0],forwardremi[0],forwardrelo[0],
            forwardimhi[0],forwardimmi[0],forwardimlo[0],
            inputrehi[ix2],inputremi[ix2],inputrelo[ix2],
            inputimhi[ix2],inputimmi[ix2],inputimlo[ix2],
            forwardrehi[1],forwardremi[1],forwardrelo[1],
            forwardimhi[1],forwardimmi[1],forwardimlo[1]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
         {
            tdf_inc
               (&outputrehi[dim][j],&outputremi[dim][j],&outputrelo[dim][j],
                forwardrehi[1][j],forwardremi[1][j],forwardrelo[1][j]);
            tdf_inc
               (&outputimhi[dim][j],&outputimmi[dim][j],&outputimlo[dim][j],
                forwardimhi[1][j],forwardimmi[1][j],forwardimlo[1][j]);
         }
      }
      else if(nvr[i] > 2)
      {
         CPU_cmplx3_speel
            (nvr[i],deg,idx[i],cffrehi[i],cffremi[i],cffrelo[i],
             cffimhi[i],cffimmi[i],cffimlo[i],
             inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
             forwardrehi,forwardremi,forwardrelo,
             forwardimhi,forwardimmi,forwardimlo,
             backwardrehi,backwardremi,backwardrelo,
             backwardimhi,backwardimmi,backwardimlo,
             crossrehi,crossremi,crossrelo,crossimhi,crossimmi,crossimlo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
         {
            tdf_inc
               (&outputrehi[dim][j],&outputremi[dim][j],&outputrelo[dim][j],
                forwardrehi[ix1][j],forwardremi[ix1][j],forwardrelo[ix1][j]);
            tdf_inc
               (&outputimhi[dim][j],&outputimmi[dim][j],&outputimlo[dim][j],
                forwardimhi[ix1][j],forwardimmi[ix1][j],forwardimlo[ix1][j]);
         }
         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
         {
            tdf_inc
               (&outputrehi[ix2][j],&outputremi[ix2][j],&outputrelo[ix2][j],
                forwardrehi[ix1][j],forwardremi[ix1][j],forwardrelo[ix1][j]);
            tdf_inc
               (&outputimhi[ix2][j],&outputimmi[ix2][j],&outputimlo[ix2][j],
                forwardimhi[ix1][j],forwardimmi[ix1][j],forwardimlo[ix1][j]);
         }
         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
         {
            tdf_inc
               (&outputrehi[ix2][j],&outputremi[ix2][j],&outputrelo[ix2][j],
                backwardrehi[ix1][j],backwardremi[ix1][j],
                backwardrelo[ix1][j]);
            tdf_inc
               (&outputimhi[ix2][j],&outputimmi[ix2][j],&outputimlo[ix2][j],
                backwardimhi[ix1][j],backwardimmi[ix1][j],
                backwardimlo[ix1][j]);
         }
         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
            {
               tdf_inc
                  (&outputrehi[ix2][j],&outputremi[ix2][j],&outputrelo[ix2][j],
                   crossrehi[k-1][j],crossremi[k-1][j],crossrelo[k-1][j]);
               tdf_inc
                  (&outputimhi[ix2][j],&outputimmi[ix2][j],&outputimlo[ix2][j],
                   crossimhi[k-1][j],crossimmi[k-1][j],crossimlo[k-1][j]);
            }
         }
      }
   }
}

void CPU_dbl3_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo, 
   double **outputhi, double **outputmi, double **outputlo,
   double *elapsedsec, bool verbose )
{
   double **forwardhi = new double*[dim];
   double **forwardmi = new double*[dim];
   double **forwardlo = new double*[dim];
   double **backwardhi = new double*[dim-1]; // in case dim = 2
   double **backwardmi = new double*[dim-1];
   double **backwardlo = new double*[dim-1];
   double **crosshi = new double*[dim-1];    // in case dim = 2
   double **crossmi = new double*[dim-1];
   double **crosslo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardhi[i] = new double[deg+1]; forwardmi[i] = new double[deg+1];
      forwardlo[i] = new double[deg+1];
      backwardhi[i] = new double[deg+1]; backwardmi[i] = new double[deg+1];
      backwardlo[i] = new double[deg+1];
      crosshi[i] = new double[deg+1]; crossmi[i] = new double[deg+1];
      crosslo[i] = new double[deg+1];
   }
   forwardhi[dim-1] = new double[deg+1];
   forwardmi[dim-1] = new double[deg+1];
   forwardlo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputhi[dim][i] = csthi[i];
      outputmi[dim][i] = cstmi[i];
      outputlo[dim][i] = cstlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhi[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputlo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl3_poly_speel
      (dim,nbr,deg,nvr,idx,cffhi,cffmi,cfflo,inputhi,inputmi,inputlo,
       outputhi,outputmi,outputlo,forwardhi,forwardmi,forwardlo,
       backwardhi,backwardmi,backwardlo,crosshi,crossmi,crosslo,verbose);
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
      free(forwardhi[i]); free(backwardhi[i]); free(crosshi[i]);
      free(forwardmi[i]); free(backwardmi[i]); free(crossmi[i]);
      free(forwardlo[i]); free(backwardlo[i]); free(crosslo[i]);
   }
   free(forwardhi[dim-1]); free(forwardmi[dim-1]); free(forwardlo[dim-1]);
   free(forwardhi); free(backwardhi); free(crosshi);
   free(forwardmi); free(backwardmi); free(crossmi);
   free(forwardlo); free(backwardlo); free(crosslo);
}

void CPU_cmplx3_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo, 
   double **inputimhi, double **inputimmi, double **inputimlo, 
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
   double *elapsedsec, bool verbose )
{
   double **forwardrehi = new double*[dim];
   double **forwardremi = new double*[dim];
   double **forwardrelo = new double*[dim];
   double **forwardimhi = new double*[dim];
   double **forwardimmi = new double*[dim];
   double **forwardimlo = new double*[dim];
   double **backwardrehi = new double*[dim-1]; // in case dim = 2
   double **backwardremi = new double*[dim-1];
   double **backwardrelo = new double*[dim-1];
   double **backwardimhi = new double*[dim-1];
   double **backwardimmi = new double*[dim-1];
   double **backwardimlo = new double*[dim-1];
   double **crossrehi = new double*[dim-1];    // in case dim = 2
   double **crossremi = new double*[dim-1];
   double **crossrelo = new double*[dim-1];
   double **crossimhi = new double*[dim-1];
   double **crossimmi = new double*[dim-1];
   double **crossimlo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardrehi[i] = new double[deg+1];
      forwardremi[i] = new double[deg+1];
      forwardrelo[i] = new double[deg+1];
      forwardimhi[i] = new double[deg+1];
      forwardimmi[i] = new double[deg+1];
      forwardimlo[i] = new double[deg+1];
      backwardrehi[i] = new double[deg+1];
      backwardremi[i] = new double[deg+1];
      backwardrelo[i] = new double[deg+1];
      backwardimhi[i] = new double[deg+1];
      backwardimmi[i] = new double[deg+1];
      backwardimlo[i] = new double[deg+1];
      crossrehi[i] = new double[deg+1];
      crossremi[i] = new double[deg+1];
      crossrelo[i] = new double[deg+1];
      crossimhi[i] = new double[deg+1];
      crossimmi[i] = new double[deg+1];
      crossimlo[i] = new double[deg+1];
   }
   forwardrehi[dim-1] = new double[deg+1];
   forwardremi[dim-1] = new double[deg+1];
   forwardrelo[dim-1] = new double[deg+1];
   forwardimhi[dim-1] = new double[deg+1];
   forwardimmi[dim-1] = new double[deg+1];
   forwardimlo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputrehi[dim][i] = cstrehi[i];
      outputremi[dim][i] = cstremi[i];
      outputrelo[dim][i] = cstrelo[i];
      outputimhi[dim][i] = cstimhi[i];
      outputimmi[dim][i] = cstimmi[i];
      outputimlo[dim][i] = cstimlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrehi[i][j] = 0.0;
         outputremi[i][j] = 0.0;
         outputrelo[i][j] = 0.0;
         outputimhi[i][j] = 0.0;
         outputimmi[i][j] = 0.0;
         outputimlo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_cmplx3_poly_speel
      (dim,nbr,deg,nvr,idx,
       cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
       inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
       outputrehi,outputremi,outputrelo,outputimhi,outputimmi,outputimlo,
       forwardrehi,forwardremi,forwardrelo,
       forwardimhi,forwardimmi,forwardimlo,
       backwardrehi,backwardremi,backwardrelo,
       backwardimhi,backwardimmi,backwardimlo,
       crossrehi,crossremi,crossrelo,crossimhi,crossimmi,crossimlo,verbose);
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
      free(forwardrehi[i]); free(backwardrehi[i]); free(crossrehi[i]);
      free(forwardremi[i]); free(backwardremi[i]); free(crossremi[i]);
      free(forwardrelo[i]); free(backwardrelo[i]); free(crossrelo[i]);
      free(forwardimhi[i]); free(backwardimhi[i]); free(crossimhi[i]);
      free(forwardimmi[i]); free(backwardimmi[i]); free(crossimmi[i]);
      free(forwardimlo[i]); free(backwardimlo[i]); free(crossimlo[i]);
   }
   free(forwardrehi[dim-1]); free(forwardremi[dim-1]);
   free(forwardrelo[dim-1]);
   free(forwardimhi[dim-1]); free(forwardimmi[dim-1]);
   free(forwardimlo[dim-1]);
   free(forwardrehi); free(backwardrehi); free(crossrehi);
   free(forwardremi); free(backwardremi); free(crossremi);
   free(forwardrelo); free(backwardrelo); free(crossrelo);
   free(forwardimhi); free(backwardimhi); free(crossimhi);
   free(forwardimmi); free(backwardimmi); free(crossimmi);
   free(forwardimlo); free(backwardimlo); free(crossimlo);
}

void CPU_dbl3_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhi, double *cffmi, double *cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **forwardhi, double **forwardmi, double **forwardlo,
   double **backwardhi, double **backwardmi, double **backwardlo,
   double **crosshi, double **crossmi, double **crosslo,
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
         CPU_dbl3_product(deg,
            cffhi,cffmi,cfflo,
            inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
            forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               cffhi,cffmi,cfflo,
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               forwardhi[inp2ix],forwardmi[inp2ix],forwardlo[inp2ix],
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            crosshi[inp1ix],crossmi[inp1ix],crosslo[inp1ix],
            inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
            forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               cffhi,cffmi,cfflo,
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               forwardhi[inp1ix],forwardmi[inp1ix],forwardlo[inp1ix],
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
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
            CPU_dbl3_product(deg,
               cffhi,cffmi,cfflo,
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               cffhi,cffmi,cfflo,
               backwardhi[inp2ix],backwardmi[inp2ix],backwardlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               backwardhi[inp2ix],backwardmi[inp2ix],backwardlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl3_product(deg,
               backwardhi[inp1ix],backwardmi[inp1ix],backwardlo[inp1ix],
               cffhi,cffmi,cfflo,
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               backwardhi[inp1ix],backwardmi[inp1ix],backwardlo[inp1ix],
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            cffhi,cffmi,cfflo,
            inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
            crosshi[outidx],crossmi[outidx],crosslo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
            forwardhi[inp2ix],forwardmi[inp2ix],forwardlo[inp2ix],
            crosshi[outidx],crossmi[outidx],crosslo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl3_product(deg,
              forwardhi[inp1ix],forwardmi[inp1ix],forwardlo[inp1ix],
              inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
              crosshi[outidx],crossmi[outidx],crosslo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl3_product(deg,
              forwardhi[inp1ix],forwardmi[inp1ix],forwardlo[inp1ix],
              backwardhi[inp2ix],backwardmi[inp2ix],backwardlo[inp2ix],
              crosshi[outidx],crossmi[outidx],crosslo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            backwardhi[inp1ix],backwardmi[inp1ix],backwardlo[inp1ix],
            forwardhi[inp2ix],forwardmi[inp2ix],forwardlo[inp2ix],
            crosshi[outidx],crossmi[outidx],crosslo[outidx]);
      }
   }
}

void CPU_cmplx3_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **forwardrehi, double **forwardremi, double **forwardrelo,
   double **forwardimhi, double **forwardimmi, double **forwardimlo,
   double **backwardrehi, double **backwardremi, double **backwardrelo,
   double **backwardimhi, double **backwardimmi, double **backwardimlo,
   double **crossrehi, double **crossremi, double **crossrelo,
   double **crossimhi, double **crossimmi, double **crossimlo,
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
         CPU_cmplx3_product(deg,
            cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
            inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
            inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
            forwardrehi[outidx],forwardremi[outidx],forwardrelo[outidx],
            forwardimhi[outidx],forwardimmi[outidx],forwardimlo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx3_product(deg,
               inputrehi[inp1ix],inputremi[inp1ix],inputrelo[inp1ix],
               inputimhi[inp1ix],inputimmi[inp1ix],inputimlo[inp1ix],
               cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
               forwardrehi[outidx],forwardremi[outidx],forwardrelo[outidx],
               forwardimhi[outidx],forwardimmi[outidx],forwardimlo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_cmplx3_product(deg,
               inputrehi[inp1ix],inputremi[inp1ix],inputrelo[inp1ix],
               inputimhi[inp1ix],inputimmi[inp1ix],inputimlo[inp1ix],
               forwardrehi[inp2ix],forwardremi[inp2ix],forwardrelo[inp2ix],
               forwardimhi[inp2ix],forwardimmi[inp2ix],forwardimlo[inp2ix],
               forwardrehi[outidx],forwardremi[outidx],forwardrelo[outidx],
               forwardimhi[outidx],forwardimmi[outidx],forwardimlo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_cmplx3_product(deg,
            crossrehi[inp1ix],crossremi[inp1ix],crossrelo[inp1ix],
            crossimhi[inp1ix],crossimmi[inp1ix],crossimlo[inp1ix],
            inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
            inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
            forwardrehi[outidx],forwardremi[outidx],forwardrelo[outidx],
            forwardimhi[outidx],forwardimmi[outidx],forwardimlo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx3_product(deg,
               inputrehi[inp1ix],inputremi[inp1ix],inputrelo[inp1ix],
               inputimhi[inp1ix],inputimmi[inp1ix],inputimlo[inp1ix],
               cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
               forwardrehi[outidx],forwardremi[outidx],forwardrelo[outidx],
               forwardimhi[outidx],forwardimmi[outidx],forwardimlo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx3_product(deg,
               forwardrehi[inp1ix],forwardremi[inp1ix],forwardrelo[inp1ix],
               forwardimhi[inp1ix],forwardimmi[inp1ix],forwardimlo[inp1ix],
               inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
               inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
               forwardrehi[outidx],forwardremi[outidx],forwardrelo[outidx],
               forwardimhi[outidx],forwardimmi[outidx],forwardimlo[outidx]);
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
            CPU_cmplx3_product(deg,
               cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
               inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
               inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
               backwardrehi[outidx],backwardremi[outidx],backwardrelo[outidx],
               backwardimhi[outidx],backwardimmi[outidx],backwardimlo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_cmplx3_product(deg,
               cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
               backwardrehi[inp2ix],backwardremi[inp2ix],backwardrelo[inp2ix],
               backwardimhi[inp2ix],backwardimmi[inp2ix],backwardimlo[inp2ix],
               backwardrehi[outidx],backwardremi[outidx],backwardrelo[outidx],
               backwardimhi[outidx],backwardimmi[outidx],backwardimlo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_cmplx3_product(deg,
               inputrehi[inp1ix],inputremi[inp1ix],inputrelo[inp1ix],
               inputimhi[inp1ix],inputimmi[inp1ix],inputimlo[inp1ix],
               inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
               inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
               backwardrehi[outidx],backwardremi[outidx],backwardrelo[outidx],
               backwardimhi[outidx],backwardimmi[outidx],backwardimlo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_cmplx3_product(deg,
               inputrehi[inp1ix],inputremi[inp1ix],inputrelo[inp1ix],
               inputimhi[inp1ix],inputimmi[inp1ix],inputimlo[inp1ix],
               backwardrehi[inp2ix],backwardremi[inp2ix],backwardrelo[inp2ix],
               backwardimhi[inp2ix],backwardimmi[inp2ix],backwardimlo[inp2ix],
               backwardrehi[outidx],backwardremi[outidx],backwardrelo[outidx],
               backwardimhi[outidx],backwardimmi[outidx],backwardimlo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_cmplx3_product(deg,
               backwardrehi[inp1ix],backwardremi[inp1ix],backwardrelo[inp1ix],
               backwardimhi[inp1ix],backwardimmi[inp1ix],backwardimlo[inp1ix],
               cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
               backwardrehi[outidx],backwardremi[outidx],backwardrelo[outidx],
               backwardimhi[outidx],backwardimmi[outidx],backwardimlo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx3_product(deg,
               backwardrehi[inp1ix],backwardremi[inp1ix],backwardrelo[inp1ix],
               backwardimhi[inp1ix],backwardimmi[inp1ix],backwardimlo[inp1ix],
               inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
               inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
               backwardrehi[outidx],backwardremi[outidx],backwardrelo[outidx],
               backwardimhi[outidx],backwardimmi[outidx],backwardimlo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_cmplx3_product(deg,cffrehi,cffremi,cffrelo,
                                cffimhi,cffimmi,cffimlo,
            inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
            inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
            crossrehi[outidx],crossremi[outidx],crossrelo[outidx],
            crossimhi[outidx],crossimmi[outidx],crossimlo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx3_product(deg,
            inputrehi[inp1ix],inputremi[inp1ix],inputrelo[inp1ix],
            inputimhi[inp1ix],inputimmi[inp1ix],inputimlo[inp1ix],
            forwardrehi[inp2ix],forwardremi[inp2ix],forwardrelo[inp2ix],
            forwardimhi[inp2ix],forwardimmi[inp2ix],forwardimlo[inp2ix],
            crossrehi[outidx],crossremi[outidx],crossrelo[outidx],
            crossimhi[outidx],crossimmi[outidx],crossimlo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_cmplx3_product(deg,
              forwardrehi[inp1ix],forwardremi[inp1ix],forwardrelo[inp1ix],
              forwardimhi[inp1ix],forwardimmi[inp1ix],forwardimlo[inp1ix],
              inputrehi[inp2ix],inputremi[inp2ix],inputrelo[inp2ix],
              inputimhi[inp2ix],inputimmi[inp2ix],inputimlo[inp2ix],
              crossrehi[outidx],crossremi[outidx],crossrelo[outidx],
              crossimhi[outidx],crossimmi[outidx],crossimlo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_cmplx3_product(deg,
              forwardrehi[inp1ix],forwardremi[inp1ix],forwardrelo[inp1ix],
              forwardimhi[inp1ix],forwardimmi[inp1ix],forwardimlo[inp1ix],
              backwardrehi[inp2ix],backwardremi[inp2ix],backwardrelo[inp2ix],
              backwardimhi[inp2ix],backwardimmi[inp2ix],backwardimlo[inp2ix],
              crossrehi[outidx],crossremi[outidx],crossrelo[outidx],
              crossimhi[outidx],crossimmi[outidx],crossimlo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx3_product(deg,
            backwardrehi[inp1ix],backwardremi[inp1ix],backwardrelo[inp1ix],
            backwardimhi[inp1ix],backwardimmi[inp1ix],backwardimlo[inp1ix],
            forwardrehi[inp2ix],forwardremi[inp2ix],forwardrelo[inp2ix],
            forwardimhi[inp2ix],forwardimmi[inp2ix],forwardimlo[inp2ix],
            crossrehi[outidx],crossremi[outidx],crossrelo[outidx],
            crossimhi[outidx],crossimmi[outidx],crossimlo[outidx]);
      }
   }
}

void CPU_dbl3_add_job
 ( int deg, double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double ***forwardhi, double ***forwardmi, double ***forwardlo,
   double ***backwardhi, double ***backwardmi,  double ***backwardlo, 
   double ***crosshi, double ***crossmi, double ***crosslo,
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
               tdf_inc(&forwardhi[updmon][updidx][i],
                       &forwardmi[updmon][updidx][i],
                       &forwardlo[updmon][updidx][i],
                       csthi[i],cstmi[i],cstlo[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               tdf_inc(&forwardhi[updmon][updidx][i],
                       &forwardmi[updmon][updidx][i],
                       &forwardlo[updmon][updidx][i],
                       cffhi[incidx][i],cffmi[incidx][i],cfflo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            tdf_inc(&forwardhi[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardlo[updmon][updidx][i],
                    forwardhi[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardlo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            tdf_inc(&forwardhi[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardlo[updmon][updidx][i],
                    backwardhi[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardlo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            tdf_inc(&forwardhi[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardlo[updmon][updidx][i],
                    crosshi[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crosslo[incmon][incidx][i]);
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    cffhi[incidx][i],cffmi[incidx][i],cfflo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    forwardhi[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardlo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    backwardhi[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardlo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    crosshi[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crosslo[incmon][incidx][i]);
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    cffhi[incidx][i],cffmi[incidx][i],cfflo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    forwardhi[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardlo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    backwardhi[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardlo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    crosshi[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crosslo[incmon][incidx][i]);
      }
   }
}

void CPU_cmplx3_add_job
 ( int deg, double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double ***forwardrehi, double ***forwardremi, double ***forwardrelo,
   double ***forwardimhi, double ***forwardimmi, double ***forwardimlo,
   double ***backwardrehi, double ***backwardremi, double ***backwardrelo, 
   double ***backwardimhi, double ***backwardimmi, double ***backwardimlo, 
   double ***crossrehi, double ***crossremi, double ***crossrelo,
   double ***crossimhi, double ***crossimmi, double ***crossimlo,
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
               tdf_inc(&forwardrehi[updmon][updidx][i],
                       &forwardremi[updmon][updidx][i],
                       &forwardrelo[updmon][updidx][i],
                       cstrehi[i],cstremi[i],cstrelo[i]);
               tdf_inc(&forwardimhi[updmon][updidx][i],
                       &forwardimmi[updmon][updidx][i],
                       &forwardimlo[updmon][updidx][i],
                       cstimhi[i],cstimmi[i],cstimlo[i]);
            }
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
            {
               tdf_inc(&forwardrehi[updmon][updidx][i],
                       &forwardremi[updmon][updidx][i],
                       &forwardrelo[updmon][updidx][i],
                       cffrehi[incidx][i],cffremi[incidx][i],
                       cffrelo[incidx][i]);
               tdf_inc(&forwardimhi[updmon][updidx][i],
                       &forwardimmi[updmon][updidx][i],
                       &forwardimlo[updmon][updidx][i],
                       cffimhi[incidx][i],cffimmi[incidx][i],
                       cffimlo[incidx][i]);
            }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            tdf_inc(&forwardrehi[updmon][updidx][i],
                    &forwardremi[updmon][updidx][i],
                    &forwardrelo[updmon][updidx][i],
                    forwardrehi[incmon][incidx][i],
                    forwardremi[incmon][incidx][i],
                    forwardrelo[incmon][incidx][i]);
            tdf_inc(&forwardimhi[updmon][updidx][i],
                    &forwardimmi[updmon][updidx][i],
                    &forwardimlo[updmon][updidx][i],
                    forwardimhi[incmon][incidx][i],
                    forwardimmi[incmon][incidx][i],
                    forwardimlo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            tdf_inc(&forwardrehi[updmon][updidx][i],
                    &forwardremi[updmon][updidx][i],
                    &forwardrelo[updmon][updidx][i],
                    backwardrehi[incmon][incidx][i],
                    backwardremi[incmon][incidx][i],
                    backwardrelo[incmon][incidx][i]);
            tdf_inc(&forwardimhi[updmon][updidx][i],
                    &forwardimmi[updmon][updidx][i],
                    &forwardimlo[updmon][updidx][i],
                    backwardimhi[incmon][incidx][i],
                    backwardimmi[incmon][incidx][i],
                    backwardimlo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            tdf_inc(&forwardrehi[updmon][updidx][i],
                    &forwardremi[updmon][updidx][i],
                    &forwardrelo[updmon][updidx][i],
                    crossrehi[incmon][incidx][i],
                    crossremi[incmon][incidx][i],
                    crossrelo[incmon][incidx][i]);
            tdf_inc(&forwardimhi[updmon][updidx][i],
                    &forwardimmi[updmon][updidx][i],
                    &forwardimlo[updmon][updidx][i],
                    crossimhi[incmon][incidx][i],
                    crossimmi[incmon][incidx][i],
                    crossimlo[incmon][incidx][i]);
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
            tdf_inc(&backwardrehi[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrelo[updmon][updidx][i],
                    cffrehi[incidx][i],cffremi[incidx][i],cffrelo[incidx][i]);
            tdf_inc(&backwardimhi[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimlo[updmon][updidx][i],
                    cffimhi[incidx][i],cffimmi[incidx][i],cffimlo[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            tdf_inc(&backwardrehi[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrelo[updmon][updidx][i],
                    forwardrehi[incmon][incidx][i],
                    forwardremi[incmon][incidx][i],
                    forwardrelo[incmon][incidx][i]);
            tdf_inc(&backwardimhi[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimlo[updmon][updidx][i],
                    forwardimhi[incmon][incidx][i],
                    forwardimmi[incmon][incidx][i],
                    forwardimlo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            tdf_inc(&backwardrehi[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrelo[updmon][updidx][i],
                    backwardrehi[incmon][incidx][i],
                    backwardremi[incmon][incidx][i],
                    backwardrelo[incmon][incidx][i]);
            tdf_inc(&backwardimhi[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimlo[updmon][updidx][i],
                    backwardimhi[incmon][incidx][i],
                    backwardimmi[incmon][incidx][i],
                    backwardimlo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            tdf_inc(&backwardrehi[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrelo[updmon][updidx][i],
                    crossrehi[incmon][incidx][i],
                    crossremi[incmon][incidx][i],
                    crossrelo[incmon][incidx][i]);
            tdf_inc(&backwardimhi[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimlo[updmon][updidx][i],
                    crossimhi[incmon][incidx][i],
                    crossimmi[incmon][incidx][i],
                    crossimlo[incmon][incidx][i]);
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
            tdf_inc(&crossrehi[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrelo[updmon][updidx][i],
                    cffrehi[incidx][i],cffremi[incidx][i],cffrelo[incidx][i]);
            tdf_inc(&crossimhi[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimlo[updmon][updidx][i],
                    cffimhi[incidx][i],cffimmi[incidx][i],cffimlo[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            tdf_inc(&crossrehi[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrelo[updmon][updidx][i],
                    forwardrehi[incmon][incidx][i],
                    forwardremi[incmon][incidx][i],
                    forwardrelo[incmon][incidx][i]);
            tdf_inc(&crossimhi[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimlo[updmon][updidx][i],
                    forwardimhi[incmon][incidx][i],
                    forwardimmi[incmon][incidx][i],
                    forwardimlo[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            tdf_inc(&crossrehi[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrelo[updmon][updidx][i],
                    backwardrehi[incmon][incidx][i],
                    backwardremi[incmon][incidx][i],
                    backwardrelo[incmon][incidx][i]);
            tdf_inc(&crossimhi[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimlo[updmon][updidx][i],
                    backwardimhi[incmon][incidx][i],
                    backwardimmi[incmon][incidx][i],
                    backwardimlo[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            tdf_inc(&crossrehi[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrelo[updmon][updidx][i],
                    crossrehi[incmon][incidx][i],
                    crossremi[incmon][incidx][i],
                    crossrelo[incmon][incidx][i]);
            tdf_inc(&crossimhi[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimlo[updmon][updidx][i],
                    crossimhi[incmon][incidx][i],
                    crossimmi[incmon][incidx][i],
                    crossimlo[incmon][incidx][i]);
         }
      }
   }
}


void CPU_dbl3_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   double ***forwardhi, double ***forwardmi, double ***forwardlo,
   double ***backwardhi, double ***backwardmi, double ***backwardlo,
   double ***crosshi, double ***crossmi, double ***crosslo )
{
   for(int i=0; i<=deg; i++)
   {
      outputhi[dim][i] = csthi[i];
      outputmi[dim][i] = cstmi[i];
      outputlo[dim][i] = cstlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhi[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputlo[i][j] = 0.0;
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
         tdf_inc(&outputhi[dim][i],&outputmi[dim][i],&outputlo[dim][i],
                 forwardhi[k][ix1][i],forwardmi[k][ix1][i],
                 forwardlo[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    cffhi[k][i],cffmi[k][i],cfflo[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            tdf_inc(&outputhi[ixn][i],&outputmi[ixn][i],&outputlo[ixn][i],
                    forwardhi[k][ix2][i],forwardmi[k][ix2][i],
                    forwardlo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    backwardhi[k][ix2][i],backwardmi[k][ix2][i],
                    backwardlo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],
                          &outputlo[ix0][i],crosshi[k][j-1][i],
                          crossmi[k][j-1][i],crosslo[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_cmplx3_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
   double ***forwardrehi, double ***forwardremi, double ***forwardrelo,
   double ***forwardimhi, double ***forwardimmi, double ***forwardimlo,
   double ***backwardrehi, double ***backwardremi, double ***backwardrelo,
   double ***backwardimhi, double ***backwardimmi, double ***backwardimlo,
   double ***crossrehi, double ***crossremi, double ***crossrelo,
   double ***crossimhi, double ***crossimmi, double ***crossimlo )
{
   for(int i=0; i<=deg; i++)
   {
      outputrehi[dim][i] = cstrehi[i];
      outputremi[dim][i] = cstremi[i];
      outputrelo[dim][i] = cstrelo[i];
      outputimhi[dim][i] = cstimhi[i];
      outputimmi[dim][i] = cstimmi[i];
      outputimlo[dim][i] = cstimlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrehi[i][j] = 0.0; outputremi[i][j] = 0.0;
         outputrelo[i][j] = 0.0;
         outputimhi[i][j] = 0.0; outputimmi[i][j] = 0.0;
         outputimlo[i][j] = 0.0;
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
         tdf_inc(&outputrehi[dim][i],&outputremi[dim][i],&outputrelo[dim][i],
                 forwardrehi[k][ix1][i],forwardremi[k][ix1][i],
                 forwardrelo[k][ix1][i]);
         tdf_inc(&outputimhi[dim][i],&outputimmi[dim][i],&outputimlo[dim][i],
                 forwardimhi[k][ix1][i],forwardimmi[k][ix1][i],
                 forwardimlo[k][ix1][i]);
      }
      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
         {
            tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                    &outputrelo[ix0][i],
                    cffrehi[k][i],cffremi[k][i],cffrelo[k][i]);
            tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                    &outputimlo[ix0][i],
                    cffimhi[k][i],cffimmi[k][i],cffimlo[k][i]);
         }
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            tdf_inc(&outputrehi[ixn][i],&outputremi[ixn][i],
                    &outputrelo[ixn][i],
                    forwardrehi[k][ix2][i],forwardremi[k][ix2][i],
                    forwardrelo[k][ix2][i]);
            tdf_inc(&outputimhi[ixn][i],&outputimmi[ixn][i],
                    &outputimlo[ixn][i],
                    forwardimhi[k][ix2][i],forwardimmi[k][ix2][i],
                    forwardimlo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                    &outputrelo[ix0][i],
                    backwardrehi[k][ix2][i],backwardremi[k][ix2][i],
                    backwardrelo[k][ix2][i]);
            tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                    &outputimlo[ix0][i],
                    backwardimhi[k][ix2][i],backwardimmi[k][ix2][i],
                    backwardimlo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
               {
                  tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                          &outputrelo[ix0][i],
                          crossrehi[k][j-1][i],crossremi[k][j-1][i],
                          crossrelo[k][j-1][i]);
                  tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                          &outputimlo[ix0][i],
                          crossimhi[k][j-1][i],crossimmi[k][j-1][i],
                          crossimlo[k][j-1][i]);
               }
            }
         }
      }
   }
}

void CPU_dbl3_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   double ***forwardhi, double ***forwardmi, double ***forwardlo,
   double ***backwardhi, double ***backwardmi, double ***backwardlo,
   double ***crosshi, double ***crossmi, double ***crosslo,
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

         CPU_dbl3_add_job(deg,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
            forwardhi,forwardmi,forwardlo,backwardhi,backwardmi,backwardlo,
            crosshi,crossmi,crosslo,job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputhi[dim][i] = forwardhi[lastmon][lastidx][i];
      outputmi[dim][i] = forwardmi[lastmon][lastidx][i];
      outputlo[dim][i] = forwardlo[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputhi[0][i] = 0.0;
         outputmi[0][i] = 0.0;
         outputlo[0][i] = 0.0;
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
         outputhi[0][i] = backwardhi[ix0][ix2][i];
         outputmi[0][i] = backwardmi[ix0][ix2][i];
         outputlo[0][i] = backwardlo[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputhi[k][i] = 0.0;
            outputmi[k][i] = 0.0;
            outputlo[k][i] = 0.0;
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
               outputhi[k][i] = backwardhi[ix0][ix2][i];
               outputmi[k][i] = backwardmi[ix0][ix2][i];
               outputlo[k][i] = backwardlo[ix0][ix2][i];
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
               outputhi[k][i] = forwardhi[ix0][ix2][i];
               outputmi[k][i] = forwardmi[ix0][ix2][i];
               outputlo[k][i] = forwardlo[ix0][ix2][i];
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
               outputhi[k][i] = crosshi[ix0][ix2][i];
               outputmi[k][i] = crossmi[ix0][ix2][i];
               outputlo[k][i] = crosslo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_cmplx3_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
   double ***forwardrehi, double ***forwardremi, double ***forwardrelo,
   double ***forwardimhi, double ***forwardimmi, double ***forwardimlo,
   double ***backwardrehi, double ***backwardremi, double ***backwardrelo,
   double ***backwardimhi, double ***backwardimmi, double ***backwardimlo,
   double ***crossrehi, double ***crossremi, double ***crossrelo,
   double ***crossimhi, double ***crossimmi, double ***crossimlo,
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

         CPU_cmplx3_add_job(deg,
            cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
            cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
            forwardrehi,forwardremi,forwardrelo,
            forwardimhi,forwardimmi,forwardimlo,
            backwardrehi,backwardremi,backwardrelo,
            backwardimhi,backwardimmi,backwardimlo,
            crossrehi,crossremi,crossrelo,
            crossimhi,crossimmi,crossimlo,job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputrehi[dim][i] = forwardrehi[lastmon][lastidx][i];
      outputremi[dim][i] = forwardremi[lastmon][lastidx][i];
      outputrelo[dim][i] = forwardrelo[lastmon][lastidx][i];
      outputimhi[dim][i] = forwardimhi[lastmon][lastidx][i];
      outputimmi[dim][i] = forwardimmi[lastmon][lastidx][i];
      outputimlo[dim][i] = forwardimlo[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputrehi[0][i] = 0.0; outputremi[0][i] = 0.0;
         outputrelo[0][i] = 0.0;
         outputimhi[0][i] = 0.0; outputimmi[0][i] = 0.0; 
         outputimlo[0][i] = 0.0;
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
         outputrehi[0][i] = backwardrehi[ix0][ix2][i];
         outputremi[0][i] = backwardremi[ix0][ix2][i];
         outputrelo[0][i] = backwardrelo[ix0][ix2][i];
         outputimhi[0][i] = backwardimhi[ix0][ix2][i];
         outputimmi[0][i] = backwardimmi[ix0][ix2][i];
         outputimlo[0][i] = backwardimlo[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputrehi[k][i] = 0.0; outputremi[k][i] = 0.0;
            outputrelo[k][i] = 0.0;
            outputimhi[k][i] = 0.0; outputimmi[k][i] = 0.0; 
            outputimlo[k][i] = 0.0;
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
               outputrehi[k][i] = backwardrehi[ix0][ix2][i];
               outputremi[k][i] = backwardremi[ix0][ix2][i];
               outputrelo[k][i] = backwardrelo[ix0][ix2][i];
               outputimhi[k][i] = backwardimhi[ix0][ix2][i];
               outputimmi[k][i] = backwardimmi[ix0][ix2][i];
               outputimlo[k][i] = backwardimlo[ix0][ix2][i];
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
               outputrehi[k][i] = forwardrehi[ix0][ix2][i];
               outputremi[k][i] = forwardremi[ix0][ix2][i];
               outputrelo[k][i] = forwardrelo[ix0][ix2][i];
               outputimhi[k][i] = forwardimhi[ix0][ix2][i];
               outputimmi[k][i] = forwardimmi[ix0][ix2][i];
               outputimlo[k][i] = forwardimlo[ix0][ix2][i];
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
               outputrehi[k][i] = crossrehi[ix0][ix2][i];
               outputremi[k][i] = crossremi[ix0][ix2][i];
               outputrelo[k][i] = crossrelo[ix0][ix2][i];
               outputimhi[k][i] = crossimhi[ix0][ix2][i];
               outputimmi[k][i] = crossimmi[ix0][ix2][i];
               outputimlo[k][i] = crossimlo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl3_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose )
{
   double ***forwardhi = new double**[nbr];
   double ***forwardmi = new double**[nbr];
   double ***forwardlo = new double**[nbr];
   double ***backwardhi = new double**[nbr];
   double ***backwardmi = new double**[nbr];
   double ***backwardlo = new double**[nbr];
   double ***crosshi = new double**[nbr];
   double ***crossmi = new double**[nbr];
   double ***crosslo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardhi[k] = new double*[nvrk];
      forwardmi[k] = new double*[nvrk];
      forwardlo[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardhi[k][i] = new double[deg+1];
         forwardmi[k][i] = new double[deg+1];
         forwardlo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardhi[k] = new double*[nvrk-1];
         backwardmi[k] = new double*[nvrk-1];
         backwardlo[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardhi[k][i] = new double[deg+1];
            backwardmi[k][i] = new double[deg+1];
            backwardlo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crosshi[k] = new double*[nvrk-2];
         crossmi[k] = new double*[nvrk-2];
         crosslo[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crosshi[k][i] = new double[deg+1];
            crossmi[k][i] = new double[deg+1];
            crosslo[k][i] = new double[deg+1];
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

         CPU_dbl3_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffhi[monidx],cffmi[monidx],cfflo[monidx],
             inputhi,inputmi,inputlo,
             forwardhi[monidx],forwardmi[monidx],forwardlo[monidx],
             backwardhi[monidx],backwardmi[monidx],backwardlo[monidx],
             crosshi[monidx],crossmi[monidx],crosslo[monidx],job,verbose);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl3_poly_addjobs
      (dim,nbr,deg,nvr,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
       inputhi,inputmi,inputlo,outputhi,outputmi,outputlo,
       forwardhi,forwardmi,forwardlo,backwardhi,backwardmi,backwardlo,
       crosshi,crossmi,crosslo,addjobs,verbose);
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
         free(forwardhi[k][i]);
         free(forwardmi[k][i]);
         free(forwardlo[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardhi[k][i]);
                      free(backwardmi[k][i]);
                      free(backwardlo[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crosshi[k][i]);
                      free(crossmi[k][i]);
                      free(crosslo[k][i]);
                   }
   }
   free(forwardhi); free(backwardhi); free(crosshi);
   free(forwardmi); free(backwardmi); free(crossmi);
   free(forwardlo); free(backwardlo); free(crosslo);
}

void CPU_cmplx3_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose )
{
   double ***forwardrehi = new double**[nbr];
   double ***forwardremi = new double**[nbr];
   double ***forwardrelo = new double**[nbr];
   double ***forwardimhi = new double**[nbr];
   double ***forwardimmi = new double**[nbr];
   double ***forwardimlo = new double**[nbr];
   double ***backwardrehi = new double**[nbr];
   double ***backwardremi = new double**[nbr];
   double ***backwardrelo = new double**[nbr];
   double ***backwardimhi = new double**[nbr];
   double ***backwardimmi = new double**[nbr];
   double ***backwardimlo = new double**[nbr];
   double ***crossrehi = new double**[nbr];
   double ***crossremi = new double**[nbr];
   double ***crossrelo = new double**[nbr];
   double ***crossimhi = new double**[nbr];
   double ***crossimmi = new double**[nbr];
   double ***crossimlo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardrehi[k] = new double*[nvrk];
      forwardremi[k] = new double*[nvrk];
      forwardrelo[k] = new double*[nvrk];
      forwardimhi[k] = new double*[nvrk];
      forwardimmi[k] = new double*[nvrk];
      forwardimlo[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardrehi[k][i] = new double[deg+1];
         forwardremi[k][i] = new double[deg+1];
         forwardrelo[k][i] = new double[deg+1];
         forwardimhi[k][i] = new double[deg+1];
         forwardimmi[k][i] = new double[deg+1];
         forwardimlo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardrehi[k] = new double*[nvrk-1];
         backwardremi[k] = new double*[nvrk-1];
         backwardrelo[k] = new double*[nvrk-1];
         backwardimhi[k] = new double*[nvrk-1];
         backwardimmi[k] = new double*[nvrk-1];
         backwardimlo[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardrehi[k][i] = new double[deg+1];
            backwardremi[k][i] = new double[deg+1];
            backwardrelo[k][i] = new double[deg+1];
            backwardimhi[k][i] = new double[deg+1];
            backwardimmi[k][i] = new double[deg+1];
            backwardimlo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crossrehi[k] = new double*[nvrk-2];
         crossremi[k] = new double*[nvrk-2];
         crossrelo[k] = new double*[nvrk-2];
         crossimhi[k] = new double*[nvrk-2];
         crossimmi[k] = new double*[nvrk-2];
         crossimlo[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crossrehi[k][i] = new double[deg+1];
            crossremi[k][i] = new double[deg+1];
            crossrelo[k][i] = new double[deg+1];
            crossimhi[k][i] = new double[deg+1];
            crossimmi[k][i] = new double[deg+1];
            crossimlo[k][i] = new double[deg+1];
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

         CPU_cmplx3_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffrehi[monidx],cffremi[monidx],cffrelo[monidx],
             cffimhi[monidx],cffimmi[monidx],cffimlo[monidx],
             inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
             forwardrehi[monidx],forwardremi[monidx],forwardrelo[monidx],
             forwardimhi[monidx],forwardimmi[monidx],forwardimlo[monidx],
             backwardrehi[monidx],backwardremi[monidx],backwardrelo[monidx],
             backwardimhi[monidx],backwardimmi[monidx],backwardimlo[monidx],
             crossrehi[monidx],crossremi[monidx],crossrelo[monidx],
             crossimhi[monidx],crossimmi[monidx],crossimlo[monidx],
             job,verbose);
      }
   }
   //CPU_cmplx2_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cffrehi,cffrelo,cffimhi,cffimlo,
   //    inputrehi,inputrelo,outputrehi,outputrelo,outputimhi,outputimlo,
   //    forwardrehi,forwardrelo,forwardimhi,forwardimlo,
   //    backwardrehi,backwardrelo,backwardimhi,backwardimlo,
   //    crossrehi,crossrelo,crossimhi,crossimlo);
   CPU_cmplx3_poly_addjobs
      (dim,nbr,deg,nvr,idx,
       cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
       cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
       inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
       outputrehi,outputremi,outputrelo,outputimhi,outputimmi,outputimlo,
       forwardrehi,forwardremi,forwardrelo,
       forwardimhi,forwardimmi,forwardimlo,
       backwardrehi,backwardremi,backwardrelo,
       backwardimhi,backwardimmi,backwardimlo,
       crossrehi,crossremi,crossrelo,crossimhi,crossimmi,crossimlo,
       addjobs,verbose);

   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(verbose)
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
         free(forwardrehi[k][i]); free(forwardremi[k][i]);
         free(forwardrelo[k][i]);
         free(forwardimhi[k][i]); free(forwardimmi[k][i]);
         free(forwardimlo[k][i]);
      }
      if(nvrk > 1)
         for(int i=0; i<nvrk-1; i++)
         {
            free(backwardrehi[k][i]); free(backwardremi[k][i]);
            free(backwardrelo[k][i]);
            free(backwardimhi[k][i]); free(backwardimmi[k][i]);
            free(backwardimlo[k][i]);
         }
      if(nvrk > 2)
         for(int i=0; i<nvrk-2; i++)
         {
            free(crossrehi[k][i]); free(crossremi[k][i]);
            free(crossrelo[k][i]);
            free(crossimhi[k][i]); free(crossimmi[k][i]);
            free(crossimlo[k][i]);
         }
   }
   free(forwardrehi); free(backwardrehi); free(crossrehi);
   free(forwardremi); free(backwardremi); free(crossremi);
   free(forwardrelo); free(backwardrelo); free(crossrelo);
   free(forwardimhi); free(backwardimhi); free(crossimhi);
   free(forwardimmi); free(backwardimmi); free(crossimmi);
   free(forwardimlo); free(backwardimlo); free(crossimlo);
}
