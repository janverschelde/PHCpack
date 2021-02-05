/* The file dbl5_polynomials_host.cpp defines functions specified
 * in dbl5_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "penta_double_functions.h"
#include "dbl5_convolutions_host.h"
#include "dbl5_monomials_host.h"
#include "dbl5_polynomials_host.h"

void CPU_dbl5_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double **forwardtb, double **forwardix, double **forwardmi,
   double **forwardrg, double **forwardpk,
   double **backwardtb, double **backwardix, double **backwardmi,
   double **backwardrg, double **backwardpk,
   double **crosstb, double **crossix, double **crossmi,
   double **crossrg, double **crosspk, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl5_product(deg,inputtb[ix1],inputix[ix1],inputmi[ix1],
                              inputrg[ix1],inputpk[ix1],
                                cfftb[i],    cffix[i],    cffmi[i],
                                cffrg[i],     cffpk[i],
                            forwardtb[0],forwardix[0],forwardmi[0],
                            forwardrg[0],forwardpk[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            pdf_inc(&outputtb[dim][j], &outputix[dim][j], &outputmi[dim][j],
                    &outputrg[dim][j],&outputpk[dim][j],
                    forwardtb[0][j],  forwardix[0][j],    forwardmi[0][j],
                    forwardrg[0][j],  forwardpk[0][j]);
            // output[ix1][j] += cff[i][j];
            pdf_inc(&outputtb[ix1][j],&outputix[ix1][j],&outputmi[ix1][j],
                    &outputrg[ix1][j],&outputpk[ix1][j],
                        cfftb[i][j],      cffix[i][j],      cffmi[i][j],
                        cffrg[i][j],      cffpk[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl5_product(deg, cfftb[i],    cffix[i],    cffmi[i],
                               cffrg[i],    cffpk[i],
                             inputtb[ix1],inputix[ix1],inputmi[ix1],
                             inputrg[ix1],inputpk[ix1],
                           forwardtb[0],forwardix[0],forwardmi[0],
                           forwardrg[0],forwardpk[0]);
         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
            pdf_inc(&outputtb[ix2][j],&outputix[ix2][j],&outputmi[ix2][j],
                    &outputrg[ix2][j],&outputpk[ix2][j],
                    forwardtb[0][j],   forwardix[0][j], forwardmi[0][j],
                    forwardrg[0][j],   forwardpk[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl5_product(deg, cfftb[i],     cffix[i],     cffmi[i],    
                               cffrg[i],     cffpk[i],
                             inputtb[ix2], inputix[ix2], inputmi[ix2], 
                             inputrg[ix2], inputpk[ix2],
                          backwardtb[0],backwardix[0],backwardmi[0],
                          backwardrg[0],backwardpk[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
            pdf_inc( &outputtb[ix1][j],&outputix[ix1][j],&outputmi[ix1][j],
                     &outputrg[ix1][j],&outputpk[ix1][j],
                    backwardtb[0][j], backwardix[0][j], backwardmi[0][j],
                    backwardrg[0][j], backwardpk[0][j]);

         CPU_dbl5_product(deg,forwardtb[0],forwardix[0],forwardmi[0],
                              forwardrg[0],forwardpk[0],
                                inputtb[ix2],inputix[ix2],inputmi[ix2],
                                inputrg[ix2],inputpk[ix2],
                              forwardtb[1],forwardix[1],forwardmi[1],
                              forwardrg[1],forwardpk[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
            pdf_inc(&outputtb[dim][j],&outputix[dim][j],&outputmi[dim][j],
                    &outputrg[dim][j],&outputpk[dim][j],
                    forwardtb[1][j],  forwardix[1][j],  forwardmi[1][j],
                    forwardrg[1][j],  forwardpk[1][j]);
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl5_speel(nvr[i],deg,idx[i],
                 cfftb[i],  cffix[i],  cffmi[i],  cffrg[i],  cffpk[i],
               inputtb,   inputix,   inputmi,   inputrg,   inputpk,
             forwardtb, forwardix, forwardmi, forwardrg, forwardpk,
            backwardtb,backwardix,backwardmi,backwardrg,backwardpk,
               crosstb,   crossix,   crossmi,   crossrg,   crosspk);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
            pdf_inc(&outputtb[dim][j],&outputix[dim][j],&outputmi[dim][j],
                    &outputrg[dim][j],&outputpk[dim][j],
                    forwardtb[ix1][j],forwardix[ix1][j],forwardmi[ix1][j],
                    forwardrg[ix1][j],forwardpk[ix1][j]);

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
            pdf_inc(&outputtb[ix2][j],&outputix[ix2][j],&outputmi[ix2][j],
                    &outputrg[ix2][j],&outputpk[ix2][j],
                    forwardtb[ix1][j],forwardix[ix1][j],forwardmi[ix1][j],
                    forwardrg[ix1][j],forwardpk[ix1][j]);

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
            pdf_inc( &outputtb[ix2][j], &outputix[ix2][j],  &outputmi[ix2][j],
                     &outputrg[ix2][j], &outputpk[ix2][j],
                    backwardtb[ix1][j],backwardix[ix1][j], backwardmi[ix1][j],
                    backwardrg[ix1][j],backwardpk[ix1][j]);

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
               pdf_inc(&outputtb[ix2][j],&outputix[ix2][j],&outputmi[ix2][j],
                       &outputrg[ix2][j],&outputpk[ix2][j],
                         crosstb[k-1][j],  crossix[k-1][j],  crossmi[k-1][j],
                         crossrg[k-1][j],  crosspk[k-1][j]);
         }
      }
   }
}

void CPU_cmplx5_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffretb, double **cffreix, double **cffremi, 
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk,
   double **forwardretb, double **forwardreix, double **forwardremi,
   double **forwardrerg, double **forwardrepk,
   double **forwardimtb, double **forwardimix, double **forwardimmi,
   double **forwardimrg, double **forwardimpk,
   double **backwardretb, double **backwardreix, double **backwardremi,
   double **backwardrerg, double **backwardrepk,
   double **backwardimtb, double **backwardimix, double **backwardimmi,
   double **backwardimrg, double **backwardimpk,
   double **crossretb, double **crossreix, double **crossremi,
   double **crossrerg, double **crossrepk,
   double **crossimtb, double **crossimix, double **crossimmi,
   double **crossimrg, double **crossimpk, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_cmplx5_product(deg,
            inputretb[ix1],inputreix[ix1],inputremi[ix1],
            inputrerg[ix1],inputrepk[ix1],
            inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
            inputimrg[ix1],inputimpk[ix1],
            cffretb[i],cffreix[i],cffremi[i],cffrerg[i],cffrepk[i],
            cffimtb[i],cffimix[i],cffimmi[i],cffimrg[i],cffimpk[i],
            forwardretb[0],forwardreix[0],forwardremi[0],
            forwardrerg[0],forwardrepk[0],
            forwardimtb[0],forwardimix[0],forwardimmi[0],
            forwardimrg[0],forwardimpk[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            pdf_inc
               (&outputretb[dim][j],&outputreix[dim][j],&outputremi[dim][j],
                &outputrerg[dim][j],&outputrepk[dim][j],
                forwardretb[0][j],forwardreix[0][j],forwardremi[0][j],
                forwardrerg[0][j],forwardrepk[0][j]);
            pdf_inc
               (&outputimtb[dim][j],&outputimix[dim][j],&outputimmi[dim][j],
                &outputimrg[dim][j],&outputimpk[dim][j],
                forwardimtb[0][j],forwardimix[0][j],forwardimmi[0][j],
                forwardimrg[0][j],forwardimpk[0][j]);
            // output[ix1][j] += cff[i][j];
            pdf_inc
               (&outputretb[ix1][j],&outputreix[ix1][j],&outputremi[ix1][j],
                &outputrerg[ix1][j],&outputrepk[ix1][j],
                cffretb[i][j],cffreix[i][j],cffremi[i][j],
                cffrerg[i][j],cffrepk[i][j]);
            pdf_inc
               (&outputimtb[ix1][j],&outputimix[ix1][j],&outputimmi[ix1][j],
                &outputimrg[ix1][j],&outputimpk[ix1][j],
                cffimtb[i][j],cffimix[i][j],cffimmi[i][j],
                cffimrg[i][j],cffimpk[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_cmplx5_product(deg,
            cffretb[i],cffreix[i],cffremi[i],cffrerg[i],cffrepk[i],
            cffimtb[i],cffimix[i],cffimmi[i],cffimrg[i],cffimpk[i],
            inputretb[ix1],inputreix[ix1],inputremi[ix1],
            inputrerg[ix1],inputrepk[ix1],
            inputimtb[ix1],inputimix[ix1],inputimmi[ix1],
            inputimrg[ix1],inputimpk[ix1],
            forwardretb[0],forwardreix[0],forwardremi[0],
            forwardrerg[0],forwardrepk[0],
            forwardimtb[0],forwardimix[0],forwardimmi[0],
            forwardimrg[0],forwardimpk[0]);

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
         {
            pdf_inc
               (&outputretb[ix2][j],&outputreix[ix2][j],&outputremi[ix2][j],
                &outputrerg[ix2][j],&outputrepk[ix2][j],
                forwardretb[0][j],forwardreix[0][j],forwardremi[0][j],
                forwardrerg[0][j],forwardrepk[0][j]);
            pdf_inc
               (&outputimtb[ix2][j],&outputimix[ix2][j],&outputimmi[ix2][j],
                &outputimrg[ix2][j],&outputimpk[ix2][j],
                forwardimtb[0][j],forwardimix[0][j],forwardimmi[0][j],
                forwardimrg[0][j],forwardimpk[0][j]);
         }
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_cmplx5_product(deg,
            cffretb[i],cffreix[i],cffremi[i],cffrerg[i],cffrepk[i],
            cffimtb[i],cffimix[i],cffimmi[i],cffimrg[i],cffimpk[i],
            inputretb[ix2],inputreix[ix2],inputremi[ix2],
            inputrerg[ix2],inputrepk[ix2],
            inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
            inputimrg[ix2],inputimpk[ix2],
            backwardretb[0],backwardreix[0],backwardremi[0],
            backwardrerg[0],backwardrepk[0],
            backwardimtb[0],backwardimix[0],backwardimmi[0],
            backwardimrg[0],backwardimpk[0]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
         {
            pdf_inc
               (&outputretb[ix1][j],&outputreix[ix1][j],&outputremi[ix1][j],
                &outputrerg[ix1][j],&outputrepk[ix1][j],
                backwardretb[0][j],backwardreix[0][j],backwardremi[0][j],
                backwardrerg[0][j],backwardrepk[0][j]);
            pdf_inc
               (&outputimtb[ix1][j],&outputimix[ix1][j],&outputimmi[ix1][j],
                &outputimrg[ix1][j],&outputimpk[ix1][j],
                backwardimtb[0][j],backwardimix[0][j],backwardimmi[0][j],
                backwardimrg[0][j],backwardimpk[0][j]);
         }
         CPU_cmplx5_product(deg,
            forwardretb[0],forwardreix[0],forwardremi[0],
            forwardrerg[0],forwardrepk[0],
            forwardimtb[0],forwardimix[0],forwardimmi[0],
            forwardimrg[0],forwardimpk[0],
            inputretb[ix2],inputreix[ix2],inputremi[ix2],
            inputrerg[ix2],inputrepk[ix2],
            inputimtb[ix2],inputimix[ix2],inputimmi[ix2],
            inputimrg[ix2],inputimpk[ix2],
            forwardretb[1],forwardreix[1],forwardremi[1],
            forwardrerg[1],forwardrepk[1],
            forwardimtb[1],forwardimix[1],forwardimmi[1],
            forwardimrg[1],forwardimpk[1]);

         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
         {
            pdf_inc
               (&outputretb[dim][j],&outputreix[dim][j],&outputremi[dim][j],
                &outputrerg[dim][j],&outputrepk[dim][j],
                forwardretb[1][j],forwardreix[1][j],forwardremi[1][j],
                forwardrerg[1][j],forwardrepk[1][j]);
            pdf_inc
               (&outputimtb[dim][j],&outputimix[dim][j],&outputimmi[dim][j],
                &outputimrg[dim][j],&outputimpk[dim][j],
                forwardimtb[1][j],forwardimix[1][j],forwardimmi[1][j],
                forwardimrg[1][j],forwardimpk[1][j]);
         }
      }
      else if(nvr[i] > 2)
      {
         CPU_cmplx5_speel
            (nvr[i],deg,idx[i],
             cffretb[i],cffreix[i],cffremi[i],cffrerg[i],cffrepk[i],
             cffimtb[i],cffimix[i],cffimmi[i],cffimrg[i],cffimpk[i],
             inputretb,inputreix,inputremi,inputrerg,inputrepk,
             inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
             forwardretb,forwardreix,forwardremi,forwardrerg,forwardrepk,
             forwardimtb,forwardimix,forwardimmi,forwardimrg,forwardimpk,
             backwardretb,backwardreix,backwardremi,backwardrerg,backwardrepk,
             backwardimtb,backwardimix,backwardimmi,backwardimrg,backwardimpk,
             crossretb,crossreix,crossremi,crossrerg,crossrepk,
             crossimtb,crossimix,crossimmi,crossimrg,crossimpk);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
         {
            pdf_inc
               (&outputretb[dim][j],&outputreix[dim][j],&outputremi[dim][j],
                &outputrerg[dim][j],&outputrepk[dim][j],
                forwardretb[ix1][j],forwardreix[ix1][j],forwardremi[ix1][j],
                forwardrerg[ix1][j],forwardrepk[ix1][j]);
            pdf_inc
               (&outputimtb[dim][j],&outputimix[dim][j],&outputimmi[dim][j],
                &outputimrg[dim][j],&outputimpk[dim][j],
                forwardimtb[ix1][j],forwardimix[ix1][j],forwardimmi[ix1][j],
                forwardimrg[ix1][j],forwardimpk[ix1][j]);
         }
         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
         {
            pdf_inc
               (&outputretb[ix2][j],&outputreix[ix2][j],&outputremi[ix2][j],
                &outputrerg[ix2][j],&outputrepk[ix2][j],
                forwardretb[ix1][j],forwardremi[ix1][j],forwardreix[ix1][j],
                forwardrerg[ix1][j],forwardrepk[ix1][j]);
            pdf_inc
               (&outputimtb[ix2][j],&outputimix[ix2][j],&outputimmi[ix2][j],
                &outputimrg[ix2][j],&outputimpk[ix2][j],
                forwardimtb[ix1][j],forwardimix[ix1][j],forwardimmi[ix1][j],
                forwardimrg[ix1][j],forwardimpk[ix1][j]);
         }
         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
         {
            pdf_inc
               (&outputretb[ix2][j],&outputreix[ix2][j],&outputremi[ix2][j],
                &outputrerg[ix2][j],&outputrepk[ix2][j],
                backwardretb[ix1][j],backwardreix[ix1][j],
                backwardremi[ix1][j],backwardrerg[ix1][j],
                backwardrepk[ix1][j]);
            pdf_inc
               (&outputimtb[ix2][j],&outputimix[ix2][j],&outputimmi[ix2][j],
                &outputimrg[ix2][j],&outputimpk[ix2][j],
                backwardimtb[ix1][j],backwardimix[ix1][j],
                backwardimmi[ix1][j],backwardimrg[ix1][j],
                backwardimpk[ix1][j]);
         }
         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
            {
               pdf_inc
                  (&outputretb[ix2][j],&outputreix[ix2][j],
                   &outputremi[ix2][j],&outputrerg[ix2][j],
                   &outputrepk[ix2][j],
                   crossretb[k-1][j],crossreix[k-1][j],crossremi[k-1][j],
                   crossrerg[k-1][j],crossrepk[k-1][j]);
               pdf_inc
                  (&outputimtb[ix2][j],&outputimix[ix2][j],
                   &outputimmi[ix2][j],&outputimrg[ix2][j],
                   &outputimpk[ix2][j],
                   crossimtb[k-1][j],crossimix[k-1][j],crossimmi[k-1][j],
                   crossimrg[k-1][j],crossimpk[k-1][j]);
            }
         }
      }
   }
}

void CPU_dbl5_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double *elapsedsec, bool verbose )
{
   double **forwardtb = new double*[dim];
   double **forwardix = new double*[dim];
   double **forwardmi = new double*[dim];
   double **forwardrg = new double*[dim];
   double **forwardpk = new double*[dim];
   double **backwardtb = new double*[dim-1]; // in case dim = 2
   double **backwardix = new double*[dim-1];
   double **backwardmi = new double*[dim-1];
   double **backwardrg = new double*[dim-1];
   double **backwardpk = new double*[dim-1];
   double **crosstb = new double*[dim-1];    // in case dim = 2
   double **crossix = new double*[dim-1];
   double **crossmi = new double*[dim-1];
   double **crossrg = new double*[dim-1];
   double **crosspk = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardtb[i] = new double[deg+1];
      forwardix[i] = new double[deg+1];
      forwardmi[i] = new double[deg+1];
      forwardrg[i] = new double[deg+1];
      forwardpk[i] = new double[deg+1];
      backwardtb[i] = new double[deg+1];
      backwardix[i] = new double[deg+1];
      backwardmi[i] = new double[deg+1];
      backwardrg[i] = new double[deg+1];
      backwardpk[i] = new double[deg+1];
      crosstb[i] = new double[deg+1];
      crossix[i] = new double[deg+1];
      crossmi[i] = new double[deg+1];
      crossrg[i] = new double[deg+1];
      crosspk[i] = new double[deg+1];
   }
   forwardtb[dim-1] = new double[deg+1];
   forwardix[dim-1] = new double[deg+1];
   forwardmi[dim-1] = new double[deg+1];
   forwardrg[dim-1] = new double[deg+1];
   forwardpk[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputtb[dim][i] = csttb[i];
      outputix[dim][i] = cstix[i];
      outputmi[dim][i] = cstmi[i];
      outputrg[dim][i] = cstrg[i];
      outputpk[dim][i] = cstpk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputtb[i][j] = 0.0;
         outputix[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputrg[i][j] = 0.0;
         outputpk[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl5_poly_speel
      (dim,nbr,deg,nvr,idx,
            cfftb,     cffix,     cffmi,     cffrg,     cffpk,
          inputtb,   inputix,   inputmi,   inputrg,   inputpk,
         outputtb,  outputix,  outputmi,  outputrg,  outputpk,
        forwardtb, forwardix, forwardmi, forwardrg, forwardpk,
       backwardtb,backwardix,backwardmi,backwardrg,backwardpk,
          crosstb,   crossix,   crossmi,   crossrg,   crosspk,verbose);
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
      free(forwardtb[i]); free(backwardtb[i]); free(crosstb[i]);
      free(forwardix[i]); free(backwardix[i]); free(crossix[i]);
      free(forwardmi[i]); free(backwardmi[i]); free(crossmi[i]);
      free(forwardrg[i]); free(backwardrg[i]); free(crossrg[i]);
      free(forwardpk[i]); free(backwardpk[i]); free(crosspk[i]);
   }
   free(forwardtb[dim-1]); free(forwardix[dim-1]); free(forwardmi[dim-1]);
   free(forwardrg[dim-1]); free(forwardpk[dim-1]);
   free(forwardtb); free(backwardtb); free(crosstb);
   free(forwardix); free(backwardix); free(crossix);
   free(forwardmi); free(backwardmi); free(crossmi);
   free(forwardrg); free(backwardrg); free(crossrg);
   free(forwardpk); free(backwardpk); free(crosspk);
}

void CPU_cmplx5_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk, 
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk, 
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk,
   double *elapsedsec, bool verbose )
{
   double **forwardretb = new double*[dim];
   double **forwardreix = new double*[dim];
   double **forwardremi = new double*[dim];
   double **forwardrerg = new double*[dim];
   double **forwardrepk = new double*[dim];
   double **forwardimtb = new double*[dim];
   double **forwardimix = new double*[dim];
   double **forwardimmi = new double*[dim];
   double **forwardimrg = new double*[dim];
   double **forwardimpk = new double*[dim];
   double **backwardretb = new double*[dim-1]; // in case dim = 2
   double **backwardreix = new double*[dim-1];
   double **backwardremi = new double*[dim-1];
   double **backwardrerg = new double*[dim-1];
   double **backwardrepk = new double*[dim-1];
   double **backwardimtb = new double*[dim-1];
   double **backwardimix = new double*[dim-1];
   double **backwardimmi = new double*[dim-1];
   double **backwardimrg = new double*[dim-1];
   double **backwardimpk = new double*[dim-1];
   double **crossretb = new double*[dim-1];    // in case dim = 2
   double **crossreix = new double*[dim-1];
   double **crossremi = new double*[dim-1];
   double **crossrerg = new double*[dim-1];
   double **crossrepk = new double*[dim-1];
   double **crossimtb = new double*[dim-1];
   double **crossimix = new double*[dim-1];
   double **crossimmi = new double*[dim-1];
   double **crossimrg = new double*[dim-1];
   double **crossimpk = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardretb[i] = new double[deg+1];
      forwardreix[i] = new double[deg+1];
      forwardremi[i] = new double[deg+1];
      forwardrerg[i] = new double[deg+1];
      forwardrepk[i] = new double[deg+1];
      forwardimtb[i] = new double[deg+1];
      forwardimix[i] = new double[deg+1];
      forwardimmi[i] = new double[deg+1];
      forwardimrg[i] = new double[deg+1];
      forwardimpk[i] = new double[deg+1];
      backwardretb[i] = new double[deg+1];
      backwardreix[i] = new double[deg+1];
      backwardremi[i] = new double[deg+1];
      backwardrerg[i] = new double[deg+1];
      backwardrepk[i] = new double[deg+1];
      backwardimtb[i] = new double[deg+1];
      backwardimix[i] = new double[deg+1];
      backwardimmi[i] = new double[deg+1];
      backwardimrg[i] = new double[deg+1];
      backwardimpk[i] = new double[deg+1];
      crossretb[i] = new double[deg+1];
      crossreix[i] = new double[deg+1];
      crossremi[i] = new double[deg+1];
      crossrerg[i] = new double[deg+1];
      crossrepk[i] = new double[deg+1];
      crossimtb[i] = new double[deg+1];
      crossimix[i] = new double[deg+1];
      crossimmi[i] = new double[deg+1];
      crossimrg[i] = new double[deg+1];
      crossimpk[i] = new double[deg+1];
   }
   forwardretb[dim-1] = new double[deg+1];
   forwardreix[dim-1] = new double[deg+1];
   forwardremi[dim-1] = new double[deg+1];
   forwardrerg[dim-1] = new double[deg+1];
   forwardrepk[dim-1] = new double[deg+1];
   forwardimtb[dim-1] = new double[deg+1];
   forwardimix[dim-1] = new double[deg+1];
   forwardimmi[dim-1] = new double[deg+1];
   forwardimrg[dim-1] = new double[deg+1];
   forwardimpk[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputretb[dim][i] = cstretb[i];
      outputreix[dim][i] = cstreix[i];
      outputremi[dim][i] = cstremi[i];
      outputrerg[dim][i] = cstrerg[i];
      outputrepk[dim][i] = cstrepk[i];
      outputimtb[dim][i] = cstimtb[i];
      outputimix[dim][i] = cstimix[i];
      outputimmi[dim][i] = cstimmi[i];
      outputimrg[dim][i] = cstimrg[i];
      outputimpk[dim][i] = cstimpk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputretb[i][j] = 0.0;
         outputreix[i][j] = 0.0;
         outputremi[i][j] = 0.0;
         outputrerg[i][j] = 0.0;
         outputrepk[i][j] = 0.0;
         outputimtb[i][j] = 0.0;
         outputimix[i][j] = 0.0;
         outputimmi[i][j] = 0.0;
         outputimrg[i][j] = 0.0;
         outputimpk[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_cmplx5_poly_speel
      (dim,nbr,deg,nvr,idx,
       cffretb,cffreix,cffremi,cffrerg,cffrepk,
       cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
       inputretb,inputreix,inputremi,inputrerg,inputrepk,
       inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
       outputretb,outputreix,outputremi,outputrerg,outputrepk,
       outputimtb,outputimix,outputimmi,outputimrg,outputimpk,
       forwardretb,forwardreix,forwardremi,forwardrerg,forwardrepk,
       forwardimtb,forwardimix,forwardimmi,forwardimrg,forwardimpk,
       backwardretb,backwardreix,backwardremi,backwardrerg,backwardrepk,
       backwardimtb,backwardimix,backwardimmi,backwardimrg,backwardimpk,
       crossretb,crossreix,crossremi,crossrerg,crossrepk,
       crossimtb,crossimix,crossimmi,crossimrg,crossimpk,verbose);
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
      free(forwardretb[i]); free(backwardretb[i]); free(crossretb[i]);
      free(forwardreix[i]); free(backwardreix[i]); free(crossreix[i]);
      free(forwardremi[i]); free(backwardremi[i]); free(crossremi[i]);
      free(forwardrerg[i]); free(backwardrerg[i]); free(crossrerg[i]);
      free(forwardrepk[i]); free(backwardrepk[i]); free(crossrepk[i]);
      free(forwardimtb[i]); free(backwardimtb[i]); free(crossimtb[i]);
      free(forwardimix[i]); free(backwardimix[i]); free(crossimix[i]);
      free(forwardimmi[i]); free(backwardimmi[i]); free(crossimmi[i]);
      free(forwardimrg[i]); free(backwardimrg[i]); free(crossimrg[i]);
      free(forwardimpk[i]); free(backwardimpk[i]); free(crossimpk[i]);
   }
   free(forwardretb[dim-1]); free(forwardreix[dim-1]);
   free(forwardremi[dim-1]); free(forwardrerg[dim-1]);
   free(forwardrepk[dim-1]);
   free(forwardimtb[dim-1]); free(forwardimix[dim-1]);
   free(forwardimmi[dim-1]); free(forwardimrg[dim-1]);
   free(forwardimpk[dim-1]);
   free(forwardretb); free(backwardretb); free(crossretb);
   free(forwardreix); free(backwardreix); free(crossreix);
   free(forwardremi); free(backwardremi); free(crossremi);
   free(forwardrerg); free(backwardrerg); free(crossrerg);
   free(forwardrepk); free(backwardrepk); free(crossrepk);
   free(forwardimtb); free(backwardimtb); free(crossimtb);
   free(forwardimix); free(backwardimix); free(crossimix);
   free(forwardimmi); free(backwardimmi); free(crossimmi);
   free(forwardimrg); free(backwardimrg); free(crossimrg);
   free(forwardimpk); free(backwardimpk); free(crossimpk);
}

void CPU_dbl5_conv_job
 ( int deg, int nvr, int *idx,
   double *cfftb, double *cffix, double *cffmi,
   double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **forwardtb, double **forwardix, double **forwardmi,
   double **forwardrg, double **forwardpk,
   double **backwardtb, double **backwardix, double **backwardmi,
   double **backwardrg, double **backwardpk,
   double **crosstb, double **crossix, double **crossmi,
   double **crossrg, double **crosspk,
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
         CPU_dbl5_product(deg,
                cfftb,cffix,cffmi,cffrg,cffpk,
              inputtb[inp2ix],  inputix[inp2ix],  inputmi[inp2ix],
              inputrg[inp2ix],  inputpk[inp2ix],
            forwardtb[outidx],forwardix[outidx],forwardmi[outidx],
            forwardrg[outidx],forwardpk[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl5_product(deg,
               inputtb[inp1ix],inputix[inp1ix],inputmi[inp1ix],
               inputrg[inp1ix],inputpk[inp1ix],
               cfftb,cffix,cffmi,cffrg,cffpk,
               forwardtb[outidx],forwardix[outidx],forwardmi[outidx],
               forwardrg[outidx],forwardpk[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl5_product(deg,
                 inputtb[inp1ix],  inputix[inp1ix],  inputmi[inp1ix],
                 inputrg[inp1ix],  inputpk[inp1ix],
               forwardtb[inp2ix],forwardix[inp2ix],forwardmi[inp2ix],
               forwardrg[inp2ix],forwardpk[inp2ix],
               forwardtb[outidx],forwardix[outidx],forwardmi[outidx],
               forwardrg[outidx],forwardpk[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl5_product(deg,
              crosstb[inp1ix],  crossix[inp1ix],  crossmi[inp1ix],
              crossrg[inp1ix],  crosspk[inp1ix],
              inputtb[inp2ix],  inputix[inp2ix],  inputmi[inp2ix],
              inputrg[inp2ix],  inputpk[inp2ix],
            forwardtb[outidx],forwardix[outidx],forwardmi[outidx],
            forwardrg[outidx],forwardpk[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl5_product(deg,
               inputtb[inp1ix],inputix[inp1ix],inputmi[inp1ix],
               inputrg[inp1ix],inputpk[inp1ix],
               cfftb,cffix,cffmi,cffrg,cffpk,
               forwardtb[outidx],forwardix[outidx],forwardmi[outidx],
               forwardrg[outidx],forwardpk[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl5_product(deg,
               forwardtb[inp1ix],forwardix[inp1ix],forwardmi[inp1ix],
               forwardrg[inp1ix],forwardpk[inp1ix],
                 inputtb[inp2ix],  inputix[inp2ix],  inputmi[inp2ix],
                 inputrg[inp2ix],  inputpk[inp2ix],
               forwardtb[outidx],forwardix[outidx],forwardmi[outidx],
               forwardrg[outidx],forwardpk[outidx]);
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
            CPU_dbl5_product(deg,
                    cfftb,cffix,cffmi,cffrg,cffpk,
                  inputtb[inp2ix],   inputix[inp2ix],   inputmi[inp2ix],
                  inputrg[inp2ix],   inputpk[inp2ix],
               backwardtb[outidx],backwardix[outidx],backwardmi[outidx],
               backwardrg[outidx],backwardpk[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl5_product(deg,
               cfftb,cffix,cffmi,cffrg,cffpk,
               backwardtb[inp2ix],backwardix[inp2ix],backwardmi[inp2ix],
               backwardrg[inp2ix],backwardpk[inp2ix],
               backwardtb[outidx],backwardix[outidx],backwardmi[outidx],
               backwardrg[outidx],backwardpk[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl5_product(deg,
                  inputtb[inp1ix],   inputix[inp1ix],   inputmi[inp1ix],
                  inputrg[inp1ix],   inputpk[inp1ix],
                  inputtb[inp2ix],   inputix[inp2ix],   inputmi[inp2ix],
                  inputrg[inp2ix],   inputpk[inp2ix],
               backwardtb[outidx],backwardix[outidx],backwardmi[outidx],
               backwardrg[outidx],backwardpk[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl5_product(deg,
                  inputtb[inp1ix],   inputix[inp1ix],   inputmi[inp1ix],
                  inputrg[inp1ix],   inputpk[inp1ix],
               backwardtb[inp2ix],backwardix[inp2ix],backwardmi[inp2ix],
               backwardrg[inp2ix],backwardpk[inp2ix],
               backwardtb[outidx],backwardix[outidx],backwardmi[outidx],
               backwardrg[outidx],backwardpk[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl5_product(deg,
               backwardtb[inp1ix],backwardix[inp1ix],backwardmi[inp1ix],
               backwardrg[inp1ix],backwardpk[inp1ix],
               cfftb,cffix,cffmi,cffrg,cffpk,
               backwardtb[outidx],backwardix[outidx],backwardmi[outidx],
               backwardrg[outidx],backwardpk[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl5_product(deg,
               backwardtb[inp1ix],backwardix[inp1ix],backwardmi[inp1ix],
               backwardrg[inp1ix],backwardpk[inp1ix],
                  inputtb[inp2ix],   inputix[inp2ix],   inputmi[inp2ix],
                  inputrg[inp2ix],   inputpk[inp2ix],
               backwardtb[outidx],backwardix[outidx],backwardmi[outidx],
               backwardrg[outidx],backwardpk[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl5_product(deg,
            cfftb,cffix,cffmi,cffrg,cffpk,
            inputtb[inp2ix],inputix[inp2ix],inputmi[inp2ix],
            inputrg[inp2ix],inputpk[inp2ix],
            crosstb[outidx],crossix[outidx],crossmi[outidx],
            crossrg[outidx],crosspk[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl5_product(deg,
              inputtb[inp1ix],  inputix[inp1ix],  inputmi[inp1ix],
              inputrg[inp1ix],  inputpk[inp1ix],
            forwardtb[inp2ix],forwardix[inp2ix],forwardmi[inp2ix],
            forwardrg[inp2ix],forwardpk[inp2ix],
              crosstb[outidx],  crossix[outidx],  crossmi[outidx],
              crossrg[outidx],  crosspk[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl5_product(deg,
              forwardtb[inp1ix],forwardix[inp1ix],forwardmi[inp1ix],
              forwardrg[inp1ix],forwardpk[inp1ix],
                inputtb[inp2ix],  inputix[inp2ix],  inputmi[inp2ix],
                inputrg[inp2ix],  inputpk[inp2ix],
                crosstb[outidx],  crossix[outidx],  crossmi[outidx],
                crossrg[outidx],  crosspk[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl5_product(deg,
               forwardtb[inp1ix], forwardix[inp1ix], forwardmi[inp1ix],
               forwardrg[inp1ix], forwardpk[inp1ix],
              backwardtb[inp2ix],backwardix[inp2ix],backwardmi[inp2ix],
              backwardrg[inp2ix],backwardpk[inp2ix],
                 crosstb[outidx],   crossix[outidx],   crossmi[outidx],
                 crossrg[outidx],   crosspk[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl5_product(deg,
            backwardtb[inp1ix],backwardix[inp1ix],backwardmi[inp1ix],
            backwardrg[inp1ix],backwardpk[inp1ix],
             forwardtb[inp2ix], forwardix[inp2ix], forwardmi[inp2ix],
             forwardrg[inp2ix], forwardpk[inp2ix],
               crosstb[outidx],   crossix[outidx],   crossmi[outidx],
               crossrg[outidx],   crosspk[outidx]);
      }
   }
}

void CPU_cmplx5_conv_job
 ( int deg, int nvr, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **forwardretb, double **forwardreix, double **forwardremi,
   double **forwardrerg, double **forwardrepk,
   double **forwardimtb, double **forwardimix, double **forwardimmi,
   double **forwardimrg, double **forwardimpk,
   double **backwardretb, double **backwardreix, double **backwardremi,
   double **backwardrerg, double **backwardrepk,
   double **backwardimtb, double **backwardimix, double **backwardimmi,
   double **backwardimrg, double **backwardimpk,
   double **crossretb, double **crossreix, double **crossremi,
   double **crossrerg, double **crossrepk,
   double **crossimtb, double **crossimix, double **crossimmi,
   double **crossimrg, double **crossimpk,
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
         CPU_cmplx5_product(deg,
            cffretb,cffreix,cffremi,cffrerg,cffrepk,
            cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
            inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
            inputrerg[inp2ix],inputrepk[inp2ix],
            inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
            inputimrg[inp2ix],inputimpk[inp2ix],
            forwardretb[outidx],forwardreix[outidx],forwardremi[outidx],
            forwardrerg[outidx],forwardrepk[outidx],
            forwardimtb[outidx],forwardimix[outidx],forwardimmi[outidx],
            forwardimrg[outidx],forwardimpk[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx5_product(deg,
               inputretb[inp1ix],inputreix[inp1ix],inputremi[inp1ix],
               inputrerg[inp1ix],inputrepk[inp1ix],
               inputimtb[inp1ix],inputimix[inp1ix],inputimmi[inp1ix],
               inputimrg[inp1ix],inputimpk[inp1ix],
               cffretb,cffreix,cffremi,cffrerg,cffrepk,
               cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
               forwardretb[outidx],forwardreix[outidx],forwardremi[outidx],
               forwardrerg[outidx],forwardrepk[outidx],
               forwardimtb[outidx],forwardimix[outidx],forwardimmi[outidx],
               forwardimrg[outidx],forwardimpk[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_cmplx5_product(deg,
               inputretb[inp1ix],inputreix[inp1ix],inputremi[inp1ix],
               inputrerg[inp1ix],inputrepk[inp1ix],
               inputimtb[inp1ix],inputimix[inp1ix],inputimmi[inp1ix],
               inputimrg[inp1ix],inputimpk[inp1ix],
               forwardretb[inp2ix],forwardreix[inp2ix],forwardremi[inp2ix],
               forwardrerg[inp2ix],forwardrepk[inp2ix],
               forwardimtb[inp2ix],forwardimix[inp2ix],forwardimmi[inp2ix],
               forwardimrg[inp2ix],forwardimpk[inp2ix],
               forwardretb[outidx],forwardreix[outidx],forwardremi[outidx],
               forwardrerg[outidx],forwardrepk[outidx],
               forwardimtb[outidx],forwardimix[outidx],forwardimmi[outidx],
               forwardimrg[outidx],forwardimpk[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_cmplx5_product(deg,
            crossretb[inp1ix],crossreix[inp1ix],crossremi[inp1ix],
            crossrerg[inp1ix],crossrepk[inp1ix],
            crossimtb[inp1ix],crossimix[inp1ix],crossimmi[inp1ix],
            crossimrg[inp1ix],crossimpk[inp1ix],
            inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
            inputrerg[inp2ix],inputrepk[inp2ix],
            inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
            inputimrg[inp2ix],inputimpk[inp2ix],
            forwardretb[outidx],forwardreix[outidx],forwardremi[outidx],
            forwardrerg[outidx],forwardrepk[outidx],
            forwardimtb[outidx],forwardimix[outidx],forwardimmi[outidx],
            forwardimrg[outidx],forwardimpk[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx5_product(deg,
               inputretb[inp1ix],inputreix[inp1ix],inputremi[inp1ix],
               inputrerg[inp1ix],inputrepk[inp1ix],
               inputimtb[inp1ix],inputimix[inp1ix],inputimmi[inp1ix],
               inputimrg[inp1ix],inputimpk[inp1ix],
               cffretb,cffreix,cffremi,cffrerg,cffrepk,
               cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
               forwardretb[outidx],forwardreix[outidx],forwardremi[outidx],
               forwardrerg[outidx],forwardrepk[outidx],
               forwardimtb[outidx],forwardimix[outidx],forwardimmi[outidx],
               forwardimrg[outidx],forwardimpk[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx5_product(deg,
               forwardretb[inp1ix],forwardreix[inp1ix],forwardremi[inp1ix],
               forwardrerg[inp1ix],forwardrepk[inp1ix],
               forwardimtb[inp1ix],forwardimix[inp1ix],forwardimmi[inp1ix],
               forwardimrg[inp1ix],forwardimpk[inp1ix],
               inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
               inputrerg[inp2ix],inputrepk[inp2ix],
               inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
               inputimrg[inp2ix],inputimpk[inp2ix],
               forwardretb[outidx],forwardreix[outidx],forwardremi[outidx],
               forwardrerg[outidx],forwardrepk[outidx],
               forwardimtb[outidx],forwardimix[outidx],forwardimmi[outidx],
               forwardimrg[outidx],forwardimpk[outidx]);
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
            CPU_cmplx5_product(deg,
               cffretb,cffreix,cffremi,cffrerg,cffrepk,
               cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
               inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
               inputrerg[inp2ix],inputrepk[inp2ix],
               inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
               inputimrg[inp2ix],inputimpk[inp2ix],
               backwardretb[outidx],backwardreix[outidx],backwardremi[outidx],
               backwardrerg[outidx],backwardrepk[outidx],
               backwardimtb[outidx],backwardimix[outidx],backwardimmi[outidx],
               backwardimrg[outidx],backwardimpk[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_cmplx5_product(deg,
               cffretb,cffreix,cffremi,cffrerg,cffrepk,
               cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
               backwardretb[inp2ix],backwardreix[inp2ix],backwardremi[inp2ix],
               backwardrerg[inp2ix],backwardrepk[inp2ix],
               backwardimtb[inp2ix],backwardimix[inp2ix],backwardimmi[inp2ix],
               backwardimrg[inp2ix],backwardimpk[inp2ix],
               backwardretb[outidx],backwardreix[outidx],backwardremi[outidx],
               backwardrerg[outidx],backwardrepk[outidx],
               backwardimtb[outidx],backwardimix[outidx],backwardimmi[outidx],
               backwardimrg[outidx],backwardimpk[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_cmplx5_product(deg,
               inputretb[inp1ix],inputreix[inp1ix],inputremi[inp1ix],
               inputrerg[inp1ix],inputrepk[inp1ix],
               inputimtb[inp1ix],inputimix[inp1ix],inputimmi[inp1ix],
               inputimrg[inp1ix],inputimpk[inp1ix],
               inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
               inputrerg[inp2ix],inputrepk[inp2ix],
               inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
               inputimrg[inp2ix],inputimpk[inp2ix],
               backwardretb[outidx],backwardreix[outidx],backwardremi[outidx],
               backwardrerg[outidx],backwardrepk[outidx],
               backwardimtb[outidx],backwardimix[outidx],backwardimmi[outidx],
               backwardimrg[outidx],backwardimpk[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_cmplx5_product(deg,
               inputretb[inp1ix],inputreix[inp1ix],inputremi[inp1ix],
               inputrerg[inp1ix],inputrepk[inp1ix],
               inputimtb[inp1ix],inputimix[inp1ix],inputimmi[inp1ix],
               inputimrg[inp1ix],inputimpk[inp1ix],
               backwardretb[inp2ix],backwardreix[inp2ix],backwardremi[inp2ix],
               backwardrerg[inp2ix],backwardrepk[inp2ix],
               backwardimtb[inp2ix],backwardimix[inp2ix],backwardimmi[inp2ix],
               backwardimrg[inp2ix],backwardimpk[inp2ix],
               backwardretb[outidx],backwardreix[outidx],backwardremi[outidx],
               backwardrerg[outidx],backwardrepk[outidx],
               backwardimtb[outidx],backwardimix[outidx],backwardimmi[outidx],
               backwardimrg[outidx],backwardimpk[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_cmplx5_product(deg,
               backwardretb[inp1ix],backwardreix[inp1ix],backwardremi[inp1ix],
               backwardrerg[inp1ix],backwardrepk[inp1ix],
               backwardimtb[inp1ix],backwardimix[inp1ix],backwardimmi[inp1ix],
               backwardimrg[inp1ix],backwardimpk[inp1ix],
               cffretb,cffreix,cffremi,cffrerg,cffrepk,
               cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
               backwardretb[outidx],backwardreix[outidx],backwardremi[outidx],
               backwardrerg[outidx],backwardrepk[outidx],
               backwardimtb[outidx],backwardimix[outidx],backwardimmi[outidx],
               backwardimrg[outidx],backwardimpk[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx5_product(deg,
               backwardretb[inp1ix],backwardreix[inp1ix],backwardremi[inp1ix],
               backwardrerg[inp1ix],backwardrepk[inp1ix],
               backwardimtb[inp1ix],backwardimix[inp1ix],backwardimmi[inp1ix],
               backwardimrg[inp1ix],backwardimpk[inp1ix],
               inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
               inputrerg[inp2ix],inputrepk[inp2ix],
               inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
               inputimrg[inp2ix],inputimpk[inp2ix],
               backwardretb[outidx],backwardreix[outidx],backwardremi[outidx],
               backwardrerg[outidx],backwardrepk[outidx],
               backwardimtb[outidx],backwardimix[outidx],backwardimmi[outidx],
               backwardimrg[outidx],backwardimpk[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_cmplx5_product(deg,
            cffretb,cffreix,cffremi,cffrerg,cffrepk,
            cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
            inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
            inputrerg[inp2ix],inputrepk[inp2ix],
            inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
            inputimrg[inp2ix],inputimpk[inp2ix],
            crossretb[outidx],crossreix[outidx],crossremi[outidx],
            crossrerg[outidx],crossrepk[outidx],
            crossimtb[outidx],crossimix[outidx],crossimmi[outidx],
            crossimrg[outidx],crossimpk[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx5_product(deg,
            inputretb[inp1ix],inputreix[inp1ix],inputremi[inp1ix],
            inputrerg[inp1ix],inputrepk[inp1ix],
            inputimtb[inp1ix],inputimix[inp1ix],inputimmi[inp1ix],
            inputimrg[inp1ix],inputimpk[inp1ix],
            forwardretb[inp2ix],forwardreix[inp2ix],forwardremi[inp2ix],
            forwardrerg[inp2ix],forwardrepk[inp2ix],
            forwardimtb[inp2ix],forwardimix[inp2ix],forwardimmi[inp2ix],
            forwardimrg[inp2ix],forwardimpk[inp2ix],
            crossretb[outidx],crossreix[outidx],crossremi[outidx],
            crossrerg[outidx],crossrepk[outidx],
            crossimtb[outidx],crossimix[outidx],crossimmi[outidx],
            crossimrg[outidx],crossimpk[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_cmplx5_product(deg,
              forwardretb[inp1ix],forwardreix[inp1ix],forwardremi[inp1ix],
              forwardrerg[inp1ix],forwardrepk[inp1ix],
              forwardimtb[inp1ix],forwardimix[inp1ix],forwardimmi[inp1ix],
              forwardimrg[inp1ix],forwardimpk[inp1ix],
              inputretb[inp2ix],inputreix[inp2ix],inputremi[inp2ix],
              inputrerg[inp2ix],inputrepk[inp2ix],
              inputimtb[inp2ix],inputimix[inp2ix],inputimmi[inp2ix],
              inputimrg[inp2ix],inputimpk[inp2ix],
              crossretb[outidx],crossreix[outidx],crossremi[outidx],
              crossrerg[outidx],crossrepk[outidx],
              crossimtb[outidx],crossimix[outidx],crossimmi[outidx],
              crossimrg[outidx],crossimpk[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_cmplx5_product(deg,
              forwardretb[inp1ix],forwardreix[inp1ix],forwardremi[inp1ix],
              forwardrerg[inp1ix],forwardrepk[inp1ix],
              forwardimtb[inp1ix],forwardimix[inp1ix],forwardimmi[inp1ix],
              forwardimrg[inp1ix],forwardimpk[inp1ix],
              backwardretb[inp2ix],backwardreix[inp2ix],backwardremi[inp2ix],
              backwardrerg[inp2ix],backwardrepk[inp2ix],
              backwardimtb[inp2ix],backwardimix[inp2ix],backwardimmi[inp2ix],
              backwardimrg[inp2ix],backwardimpk[inp2ix],
              crossretb[outidx],crossreix[outidx],crossremi[outidx],
              crossrerg[outidx],crossrepk[outidx],
              crossimtb[outidx],crossimix[outidx],crossimmi[outidx],
              crossimrg[outidx],crossimpk[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx5_product(deg,
            backwardretb[inp1ix],backwardreix[inp1ix],backwardremi[inp1ix],
            backwardrerg[inp1ix],backwardrepk[inp1ix],
            backwardimtb[inp1ix],backwardimix[inp1ix],backwardimmi[inp1ix],
            backwardimrg[inp1ix],backwardimpk[inp1ix],
            forwardretb[inp2ix],forwardreix[inp2ix],forwardremi[inp2ix],
            forwardrerg[inp2ix],forwardrepk[inp2ix],
            forwardimtb[inp2ix],forwardimix[inp2ix],forwardimmi[inp2ix],
            forwardimrg[inp2ix],forwardimpk[inp2ix],
            crossretb[outidx],crossreix[outidx],crossremi[outidx],
            crossrerg[outidx],crossrepk[outidx],
            crossimtb[outidx],crossimix[outidx],crossimmi[outidx],
            crossimrg[outidx],crossimpk[outidx]);
      }
   }
}

void CPU_dbl5_add_job
 ( int deg,
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double ***forwardtb, double ***forwardix, double ***forwardmi,
   double ***forwardrg, double ***forwardpk,
   double ***backwardtb, double ***backwardix, double ***backwardmi,
   double ***backwardrg, double ***backwardpk, 
   double ***crosstb, double ***crossix, double ***crossmi,
   double ***crossrg, double ***crosspk,
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
               pdf_inc(&forwardtb[updmon][updidx][i],
                       &forwardix[updmon][updidx][i],
                       &forwardmi[updmon][updidx][i],
                       &forwardrg[updmon][updidx][i],
                       &forwardpk[updmon][updidx][i],
                       csttb[i],cstix[i],cstmi[i],cstrg[i],cstpk[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               pdf_inc(&forwardtb[updmon][updidx][i],
                       &forwardix[updmon][updidx][i],
                       &forwardmi[updmon][updidx][i],
                       &forwardrg[updmon][updidx][i],
                       &forwardpk[updmon][updidx][i],
                       cfftb[incidx][i],cffix[incidx][i],cffmi[incidx][i],
                       cffrg[incidx][i],cffpk[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            pdf_inc(&forwardtb[updmon][updidx][i],
                    &forwardix[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardrg[updmon][updidx][i],
                    &forwardpk[updmon][updidx][i],
                    forwardtb[incmon][incidx][i],
                    forwardix[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardrg[incmon][incidx][i],
                    forwardpk[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            pdf_inc(&forwardtb[updmon][updidx][i],
                    &forwardix[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardrg[updmon][updidx][i],
                    &forwardpk[updmon][updidx][i],
                    backwardtb[incmon][incidx][i],
                    backwardix[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardrg[incmon][incidx][i],
                    backwardpk[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            pdf_inc(&forwardtb[updmon][updidx][i],
                    &forwardix[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardrg[updmon][updidx][i],
                    &forwardpk[updmon][updidx][i],
                    crosstb[incmon][incidx][i],
                    crossix[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crossrg[incmon][incidx][i],
                    crosspk[incmon][incidx][i]);
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
            pdf_inc(&backwardtb[updmon][updidx][i],
                    &backwardix[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardrg[updmon][updidx][i],
                    &backwardpk[updmon][updidx][i],
                    cfftb[incidx][i],cffix[incidx][i],cffmi[incidx][i],
                    cffrg[incidx][i],cffpk[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            pdf_inc(&backwardtb[updmon][updidx][i],
                    &backwardix[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardrg[updmon][updidx][i],
                    &backwardpk[updmon][updidx][i],
                    forwardtb[incmon][incidx][i],
                    forwardix[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardrg[incmon][incidx][i],
                    forwardpk[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            pdf_inc(&backwardtb[updmon][updidx][i],
                    &backwardix[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardrg[updmon][updidx][i],
                    &backwardpk[updmon][updidx][i],
                    backwardtb[incmon][incidx][i],
                    backwardix[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardrg[incmon][incidx][i],
                    backwardpk[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            pdf_inc(&backwardtb[updmon][updidx][i],
                    &backwardix[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardrg[updmon][updidx][i],
                    &backwardpk[updmon][updidx][i],
                    crosstb[incmon][incidx][i],
                    crossix[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crossrg[incmon][incidx][i],
                    crosspk[incmon][incidx][i]);
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
            pdf_inc(&crosstb[updmon][updidx][i],
                    &crossix[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crossrg[updmon][updidx][i],
                    &crosspk[updmon][updidx][i],
                    cfftb[incidx][i],cffix[incidx][i],cffmi[incidx][i],
                    cffrg[incidx][i],cffpk[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            pdf_inc(&crosstb[updmon][updidx][i],
                    &crossix[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crossrg[updmon][updidx][i],
                    &crosspk[updmon][updidx][i],
                    forwardtb[incmon][incidx][i],
                    forwardix[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardrg[incmon][incidx][i],
                    forwardpk[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            pdf_inc(&crosstb[updmon][updidx][i],
                    &crossix[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crossrg[updmon][updidx][i],
                    &crosspk[updmon][updidx][i],
                    backwardtb[incmon][incidx][i],
                    backwardix[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardrg[incmon][incidx][i],
                    backwardpk[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            pdf_inc(&crosstb[updmon][updidx][i],
                    &crossix[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crossrg[updmon][updidx][i],
                    &crosspk[updmon][updidx][i],
                    crosstb[incmon][incidx][i],
                    crossix[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crossrg[incmon][incidx][i],
                    crosspk[incmon][incidx][i]);
      }
   }
}

void CPU_cmplx5_add_job
 ( int deg, double *cstretb, double *cstreix, double *cstremi,
            double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk,
   double ***forwardretb, double ***forwardreix, double ***forwardremi,
   double ***forwardrerg, double ***forwardrepk,
   double ***forwardimtb, double ***forwardimix, double ***forwardimmi,
   double ***forwardimrg, double ***forwardimpk,
   double ***backwardretb, double ***backwardreix, double ***backwardremi,
   double ***backwardrerg, double ***backwardrepk, 
   double ***backwardimtb, double ***backwardimix, double ***backwardimmi,
   double ***backwardimrg, double ***backwardimpk, 
   double ***crossretb, double ***crossreix, double ***crossremi,
   double ***crossrerg, double ***crossrepk,
   double ***crossimtb, double ***crossimix, double ***crossimmi,
   double ***crossimrg, double ***crossimpk,
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
               pdf_inc(&forwardretb[updmon][updidx][i],
                       &forwardreix[updmon][updidx][i],
                       &forwardremi[updmon][updidx][i],
                       &forwardrerg[updmon][updidx][i],
                       &forwardrepk[updmon][updidx][i],
                       cstretb[i],cstreix[i],cstremi[i],
                       cstrerg[i],cstrepk[i]);
               pdf_inc(&forwardimtb[updmon][updidx][i],
                       &forwardimix[updmon][updidx][i],
                       &forwardimmi[updmon][updidx][i],
                       &forwardimrg[updmon][updidx][i],
                       &forwardimpk[updmon][updidx][i],
                       cstimtb[i],cstimix[i],cstimmi[i],
                       cstimrg[i],cstimpk[i]);
            }
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
            {
               pdf_inc(&forwardretb[updmon][updidx][i],
                       &forwardreix[updmon][updidx][i],
                       &forwardremi[updmon][updidx][i],
                       &forwardrerg[updmon][updidx][i],
                       &forwardrepk[updmon][updidx][i],
                       cffretb[incidx][i],cffreix[incidx][i],
                       cffremi[incidx][i],cffrerg[incidx][i],
                       cffrepk[incidx][i]);
               pdf_inc(&forwardimtb[updmon][updidx][i],
                       &forwardimix[updmon][updidx][i],
                       &forwardimmi[updmon][updidx][i],
                       &forwardimrg[updmon][updidx][i],
                       &forwardimpk[updmon][updidx][i],
                       cffimtb[incidx][i],cffimix[incidx][i],
                       cffimmi[incidx][i],cffimrg[incidx][i],
                       cffimpk[incidx][i]);
            }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            pdf_inc(&forwardretb[updmon][updidx][i],
                    &forwardreix[updmon][updidx][i],
                    &forwardremi[updmon][updidx][i],
                    &forwardrerg[updmon][updidx][i],
                    &forwardrepk[updmon][updidx][i],
                    forwardretb[incmon][incidx][i],
                    forwardreix[incmon][incidx][i],
                    forwardremi[incmon][incidx][i],
                    forwardrerg[incmon][incidx][i],
                    forwardrepk[incmon][incidx][i]);
            pdf_inc(&forwardimtb[updmon][updidx][i],
                    &forwardimix[updmon][updidx][i],
                    &forwardimmi[updmon][updidx][i],
                    &forwardimrg[updmon][updidx][i],
                    &forwardimpk[updmon][updidx][i],
                    forwardimtb[incmon][incidx][i],
                    forwardimix[incmon][incidx][i],
                    forwardimmi[incmon][incidx][i],
                    forwardimrg[incmon][incidx][i],
                    forwardimpk[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            pdf_inc(&forwardretb[updmon][updidx][i],
                    &forwardreix[updmon][updidx][i],
                    &forwardremi[updmon][updidx][i],
                    &forwardrerg[updmon][updidx][i],
                    &forwardrepk[updmon][updidx][i],
                    backwardretb[incmon][incidx][i],
                    backwardreix[incmon][incidx][i],
                    backwardremi[incmon][incidx][i],
                    backwardrerg[incmon][incidx][i],
                    backwardrepk[incmon][incidx][i]);
            pdf_inc(&forwardimtb[updmon][updidx][i],
                    &forwardimix[updmon][updidx][i],
                    &forwardimmi[updmon][updidx][i],
                    &forwardimrg[updmon][updidx][i],
                    &forwardimpk[updmon][updidx][i],
                    backwardimtb[incmon][incidx][i],
                    backwardimix[incmon][incidx][i],
                    backwardimmi[incmon][incidx][i],
                    backwardimrg[incmon][incidx][i],
                    backwardimpk[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            pdf_inc(&forwardretb[updmon][updidx][i],
                    &forwardreix[updmon][updidx][i],
                    &forwardremi[updmon][updidx][i],
                    &forwardrerg[updmon][updidx][i],
                    &forwardrepk[updmon][updidx][i],
                    crossretb[incmon][incidx][i],
                    crossreix[incmon][incidx][i],
                    crossremi[incmon][incidx][i],
                    crossrerg[incmon][incidx][i],
                    crossrepk[incmon][incidx][i]);
            pdf_inc(&forwardimtb[updmon][updidx][i],
                    &forwardimix[updmon][updidx][i],
                    &forwardimmi[updmon][updidx][i],
                    &forwardimrg[updmon][updidx][i],
                    &forwardimpk[updmon][updidx][i],
                    crossimtb[incmon][incidx][i],
                    crossimix[incmon][incidx][i],
                    crossimmi[incmon][incidx][i],
                    crossimrg[incmon][incidx][i],
                    crossimpk[incmon][incidx][i]);
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
            pdf_inc(&backwardretb[updmon][updidx][i],
                    &backwardreix[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrerg[updmon][updidx][i],
                    &backwardrepk[updmon][updidx][i],
                    cffretb[incidx][i],cffreix[incidx][i],
                    cffremi[incidx][i],cffrerg[incidx][i],
                    cffrepk[incidx][i]);
            pdf_inc(&backwardimtb[updmon][updidx][i],
                    &backwardimix[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimrg[updmon][updidx][i],
                    &backwardimpk[updmon][updidx][i],
                    cffimtb[incidx][i],cffimix[incidx][i],
                    cffimmi[incidx][i],cffimrg[incidx][i],
                    cffimpk[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            pdf_inc(&backwardretb[updmon][updidx][i],
                    &backwardreix[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrerg[updmon][updidx][i],
                    &backwardrepk[updmon][updidx][i],
                    forwardretb[incmon][incidx][i],
                    forwardreix[incmon][incidx][i],
                    forwardremi[incmon][incidx][i],
                    forwardrerg[incmon][incidx][i],
                    forwardrepk[incmon][incidx][i]);
            pdf_inc(&backwardimtb[updmon][updidx][i],
                    &backwardimix[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimrg[updmon][updidx][i],
                    &backwardimpk[updmon][updidx][i],
                    forwardimtb[incmon][incidx][i],
                    forwardimix[incmon][incidx][i],
                    forwardimmi[incmon][incidx][i],
                    forwardimrg[incmon][incidx][i],
                    forwardimpk[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            pdf_inc(&backwardretb[updmon][updidx][i],
                    &backwardreix[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrerg[updmon][updidx][i],
                    &backwardrepk[updmon][updidx][i],
                    backwardretb[incmon][incidx][i],
                    backwardreix[incmon][incidx][i],
                    backwardremi[incmon][incidx][i],
                    backwardrerg[incmon][incidx][i],
                    backwardrepk[incmon][incidx][i]);
            pdf_inc(&backwardimtb[updmon][updidx][i],
                    &backwardimix[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimrg[updmon][updidx][i],
                    &backwardimpk[updmon][updidx][i],
                    backwardimtb[incmon][incidx][i],
                    backwardimix[incmon][incidx][i],
                    backwardimmi[incmon][incidx][i],
                    backwardimrg[incmon][incidx][i],
                    backwardimpk[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            pdf_inc(&backwardretb[updmon][updidx][i],
                    &backwardreix[updmon][updidx][i],
                    &backwardremi[updmon][updidx][i],
                    &backwardrerg[updmon][updidx][i],
                    &backwardrepk[updmon][updidx][i],
                    crossretb[incmon][incidx][i],
                    crossreix[incmon][incidx][i],
                    crossremi[incmon][incidx][i],
                    crossrerg[incmon][incidx][i],
                    crossrepk[incmon][incidx][i]);
            pdf_inc(&backwardimtb[updmon][updidx][i],
                    &backwardimix[updmon][updidx][i],
                    &backwardimmi[updmon][updidx][i],
                    &backwardimrg[updmon][updidx][i],
                    &backwardimpk[updmon][updidx][i],
                    crossimtb[incmon][incidx][i],
                    crossimix[incmon][incidx][i],
                    crossimmi[incmon][incidx][i],
                    crossimrg[incmon][incidx][i],
                    crossimpk[incmon][incidx][i]);
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
            pdf_inc(&crossretb[updmon][updidx][i],
                    &crossreix[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrerg[updmon][updidx][i],
                    &crossrepk[updmon][updidx][i],
                    cffretb[incidx][i],cffreix[incidx][i],
                    cffremi[incidx][i],cffrerg[incidx][i],
                    cffrepk[incidx][i]);
            pdf_inc(&crossimtb[updmon][updidx][i],
                    &crossimix[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimrg[updmon][updidx][i],
                    &crossimpk[updmon][updidx][i],
                    cffimtb[incidx][i],cffimix[incidx][i],
                    cffimmi[incidx][i],cffimrg[incidx][i],
                    cffimpk[incidx][i]);
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
         {
            pdf_inc(&crossretb[updmon][updidx][i],
                    &crossreix[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrerg[updmon][updidx][i],
                    &crossrepk[updmon][updidx][i],
                    forwardretb[incmon][incidx][i],
                    forwardreix[incmon][incidx][i],
                    forwardremi[incmon][incidx][i],
                    forwardrerg[incmon][incidx][i],
                    forwardrepk[incmon][incidx][i]);
            pdf_inc(&crossimtb[updmon][updidx][i],
                    &crossimix[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimrg[updmon][updidx][i],
                    &crossimpk[updmon][updidx][i],
                    forwardimtb[incmon][incidx][i],
                    forwardimix[incmon][incidx][i],
                    forwardimmi[incmon][incidx][i],
                    forwardimrg[incmon][incidx][i],
                    forwardimpk[incmon][incidx][i]);
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
         {
            pdf_inc(&crossretb[updmon][updidx][i],
                    &crossreix[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrerg[updmon][updidx][i],
                    &crossrepk[updmon][updidx][i],
                    backwardretb[incmon][incidx][i],
                    backwardreix[incmon][incidx][i],
                    backwardremi[incmon][incidx][i],
                    backwardrerg[incmon][incidx][i],
                    backwardrepk[incmon][incidx][i]);
            pdf_inc(&crossimtb[updmon][updidx][i],
                    &crossimix[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimrg[updmon][updidx][i],
                    &crossimpk[updmon][updidx][i],
                    backwardimtb[incmon][incidx][i],
                    backwardimix[incmon][incidx][i],
                    backwardimmi[incmon][incidx][i],
                    backwardimrg[incmon][incidx][i],
                    backwardimpk[incmon][incidx][i]);
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
         {
            pdf_inc(&crossretb[updmon][updidx][i],
                    &crossreix[updmon][updidx][i],
                    &crossremi[updmon][updidx][i],
                    &crossrerg[updmon][updidx][i],
                    &crossrepk[updmon][updidx][i],
                    crossretb[incmon][incidx][i],
                    crossreix[incmon][incidx][i],
                    crossremi[incmon][incidx][i],
                    crossrerg[incmon][incidx][i],
                    crossrepk[incmon][incidx][i]);
            pdf_inc(&crossimtb[updmon][updidx][i],
                    &crossimix[updmon][updidx][i],
                    &crossimmi[updmon][updidx][i],
                    &crossimrg[updmon][updidx][i],
                    &crossimpk[updmon][updidx][i],
                    crossimtb[incmon][incidx][i],
                    crossimix[incmon][incidx][i],
                    crossimmi[incmon][incidx][i],
                    crossimrg[incmon][incidx][i],
                    crossimpk[incmon][incidx][i]);
         }
      }
   }
}

void CPU_dbl5_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double ***forwardtb, double ***forwardix, double ***forwardmi,
   double ***forwardrg, double ***forwardpk,
   double ***backwardtb, double ***backwardix, double ***backwardmi,
   double ***backwardrg, double ***backwardpk, 
   double ***crosstb, double ***crossix, double ***crossmi,
   double ***crossrg, double ***crosspk )
{
   for(int i=0; i<=deg; i++)
   {
      outputtb[dim][i] = csttb[i];
      outputix[dim][i] = cstix[i];
      outputmi[dim][i] = cstmi[i];
      outputrg[dim][i] = cstrg[i];
      outputpk[dim][i] = cstpk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputtb[i][j] = 0.0;
         outputix[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputrg[i][j] = 0.0;
         outputpk[i][j] = 0.0;
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
         pdf_inc(&outputtb[dim][i],   &outputix[dim][i],
                 &outputmi[dim][i],
                 &outputrg[dim][i],   &outputpk[dim][i],
                 forwardtb[k][ix1][i],forwardix[k][ix1][i],
                 forwardmi[k][ix1][i],
                 forwardrg[k][ix1][i],forwardpk[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            pdf_inc(&outputtb[ix0][i],&outputix[ix0][i],&outputmi[ix0][i],
                    &outputrg[ix0][i],&outputpk[ix0][i],
                        cfftb[k][i],      cffix[k][i],      cffmi[k][i],
                        cffrg[k][i],      cffpk[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            pdf_inc(&outputtb[ixn][i],   &outputix[ixn][i],
                    &outputmi[ixn][i],
                    &outputrg[ixn][i],   &outputpk[ixn][i],
                    forwardtb[k][ix2][i],forwardix[k][ix2][i],
                    forwardmi[k][ix2][i],
                    forwardrg[k][ix2][i],forwardpk[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            pdf_inc( &outputtb[ix0][i],    &outputix[ix0][i],
                     &outputmi[ix0][i],
                     &outputrg[ix0][i],    &outputpk[ix0][i],
                    backwardtb[k][ix2][i],backwardix[k][ix2][i],
                    backwardmi[k][ix2][i],
                    backwardrg[k][ix2][i],backwardpk[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  pdf_inc(&outputtb[ix0][i], &outputix[ix0][i],
                          &outputmi[ix0][i],
                          &outputrg[ix0][i], &outputpk[ix0][i],
                            crosstb[k][j-1][i],crossix[k][j-1][i],
                            crossmi[k][j-1][i],
                            crossrg[k][j-1][i],crosspk[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_dbl5_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double ***forwardtb, double ***forwardix, double ***forwardmi,
   double ***forwardrg, double ***forwardpk,
   double ***backwardtb, double ***backwardix, double ***backwardmi,
   double ***backwardrg, double ***backwardpk, 
   double ***crosstb, double ***crossix, double ***crossmi,
   double ***crossrg, double ***crosspk,
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

         CPU_dbl5_add_job(deg,
                 csttb,     cstix,     cstmi,     cstrg,     cstpk,
                 cfftb,     cffix,     cffmi,     cffrg,     cffpk,
             forwardtb, forwardix, forwardmi, forwardrg, forwardpk,
            backwardtb,backwardix,backwardmi,backwardrg,backwardpk,
               crosstb,   crossix,   crossmi,   crossrg,   crosspk,
            job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputtb[dim][i] = forwardtb[lastmon][lastidx][i];
      outputix[dim][i] = forwardix[lastmon][lastidx][i];
      outputmi[dim][i] = forwardmi[lastmon][lastidx][i];
      outputrg[dim][i] = forwardrg[lastmon][lastidx][i];
      outputpk[dim][i] = forwardpk[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputtb[0][i] = 0.0;
         outputix[0][i] = 0.0;
         outputmi[0][i] = 0.0;
         outputrg[0][i] = 0.0;
         outputpk[0][i] = 0.0;
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
         outputtb[0][i] = backwardtb[ix0][ix2][i];
         outputix[0][i] = backwardix[ix0][ix2][i];
         outputmi[0][i] = backwardmi[ix0][ix2][i];
         outputrg[0][i] = backwardrg[ix0][ix2][i];
         outputpk[0][i] = backwardpk[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputtb[k][i] = 0.0;
            outputix[k][i] = 0.0;
            outputmi[k][i] = 0.0;
            outputrg[k][i] = 0.0;
            outputpk[k][i] = 0.0;
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
               outputtb[k][i] = backwardtb[ix0][ix2][i];
               outputix[k][i] = backwardix[ix0][ix2][i];
               outputmi[k][i] = backwardmi[ix0][ix2][i];
               outputrg[k][i] = backwardrg[ix0][ix2][i];
               outputpk[k][i] = backwardpk[ix0][ix2][i];
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
               outputtb[k][i] = forwardtb[ix0][ix2][i];
               outputix[k][i] = forwardix[ix0][ix2][i];
               outputmi[k][i] = forwardmi[ix0][ix2][i];
               outputrg[k][i] = forwardrg[ix0][ix2][i];
               outputpk[k][i] = forwardpk[ix0][ix2][i];
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
               outputtb[k][i] = crosstb[ix0][ix2][i];
               outputix[k][i] = crossix[ix0][ix2][i];
               outputmi[k][i] = crossmi[ix0][ix2][i];
               outputrg[k][i] = crossrg[ix0][ix2][i];
               outputpk[k][i] = crosspk[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl5_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose )
{
   double ***forwardtb = new double**[nbr];
   double ***forwardix = new double**[nbr];
   double ***forwardmi = new double**[nbr];
   double ***forwardrg = new double**[nbr];
   double ***forwardpk = new double**[nbr];
   double ***backwardtb = new double**[nbr];
   double ***backwardix = new double**[nbr];
   double ***backwardmi = new double**[nbr];
   double ***backwardrg = new double**[nbr];
   double ***backwardpk = new double**[nbr];
   double ***crosstb = new double**[nbr];
   double ***crossix = new double**[nbr];
   double ***crossmi = new double**[nbr];
   double ***crossrg = new double**[nbr];
   double ***crosspk = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardtb[k] = new double*[nvrk];
      forwardix[k] = new double*[nvrk];
      forwardmi[k] = new double*[nvrk];
      forwardrg[k] = new double*[nvrk];
      forwardpk[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardtb[k][i] = new double[deg+1];
         forwardix[k][i] = new double[deg+1];
         forwardmi[k][i] = new double[deg+1];
         forwardrg[k][i] = new double[deg+1];
         forwardpk[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardtb[k] = new double*[nvrk-1];
         backwardix[k] = new double*[nvrk-1];
         backwardmi[k] = new double*[nvrk-1];
         backwardrg[k] = new double*[nvrk-1];
         backwardpk[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardtb[k][i] = new double[deg+1];
            backwardix[k][i] = new double[deg+1];
            backwardmi[k][i] = new double[deg+1];
            backwardrg[k][i] = new double[deg+1];
            backwardpk[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crosstb[k] = new double*[nvrk-2];
         crossix[k] = new double*[nvrk-2];
         crossmi[k] = new double*[nvrk-2];
         crossrg[k] = new double*[nvrk-2];
         crosspk[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crosstb[k][i] = new double[deg+1];
            crossix[k][i] = new double[deg+1];
            crossmi[k][i] = new double[deg+1];
            crossrg[k][i] = new double[deg+1];
            crosspk[k][i] = new double[deg+1];
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

         CPU_dbl5_conv_job
            (deg,nvr[monidx],idx[monidx],
             cfftb[monidx],cffix[monidx],cffmi[monidx],
             cffrg[monidx],cffpk[monidx],
             inputtb,inputix,inputmi,inputrg,inputpk,
              forwardtb[monidx], forwardix[monidx], forwardmi[monidx],
              forwardrg[monidx], forwardpk[monidx],
             backwardtb[monidx],backwardix[monidx],backwardmi[monidx],
             backwardrg[monidx],backwardpk[monidx],
                crosstb[monidx],   crossix[monidx],   crossmi[monidx],
                crossrg[monidx],   crosspk[monidx],job,verbose);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl5_poly_addjobs
      (dim,nbr,deg,nvr,idx,csttb,cstix,cstmi,cstrg,cstpk,
            cfftb,     cffix,     cffmi,     cffrg,     cffpk,
          inputtb,   inputix,   inputmi,   inputrg,   inputpk,
         outputtb,  outputix,  outputmi,  outputrg,  outputpk,
        forwardtb, forwardix, forwardmi, forwardrg, forwardpk,
       backwardtb,backwardix,backwardmi,backwardrg,backwardpk,
          crosstb,   crossix,   crossmi,   crossrg,   crosspk,addjobs,verbose);
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
         free(forwardtb[k][i]);
         free(forwardix[k][i]);
         free(forwardmi[k][i]);
         free(forwardrg[k][i]);
         free(forwardpk[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardtb[k][i]);
                      free(backwardix[k][i]);
                      free(backwardmi[k][i]);
                      free(backwardrg[k][i]);
                      free(backwardpk[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crosstb[k][i]);
                      free(crossix[k][i]);
                      free(crossmi[k][i]);
                      free(crossrg[k][i]);
                      free(crosspk[k][i]);
                   }
   }
   free(forwardtb); free(backwardtb); free(crosstb);
   free(forwardix); free(backwardix); free(crossix);
   free(forwardmi); free(backwardmi); free(crossmi);
   free(forwardrg); free(backwardrg); free(crossrg);
   free(forwardpk); free(backwardpk); free(crosspk);
}
