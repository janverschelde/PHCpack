/* The file dbl5_polynomials_host.cpp defines functions specified
 * in dbl5_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
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

void CPU_dbl5_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk, bool verbose )
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

   CPU_dbl5_poly_speel
      (dim,nbr,deg,nvr,idx,
            cfftb,     cffix,     cffmi,     cffrg,     cffpk,
          inputtb,   inputix,   inputmi,   inputrg,   inputpk,
         outputtb,  outputix,  outputmi,  outputrg,  outputpk,
        forwardtb, forwardix, forwardmi, forwardrg, forwardpk,
       backwardtb,backwardix,backwardmi,backwardrg,backwardpk,
          crosstb,   crossix,   crossmi,   crossrg,   crosspk,verbose);

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
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, bool verbose )
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
          crosstb,   crossix,   crossmi,   crossrg,   crosspk,addjobs,true);

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
