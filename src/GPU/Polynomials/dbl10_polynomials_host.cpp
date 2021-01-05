/* The file dbl10_polynomials_host.cpp defines functions specified
 * in dbl10_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "deca_double_functions.h"
#include "dbl10_convolutions_host.h"
#include "dbl10_monomials_host.h"
#include "dbl10_polynomials_host.h"

void CPU_dbl10_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   double **forwardrtb, double **forwardrix, double **forwardrmi,
   double **forwardrrg, double **forwardrpk,
   double **forwardltb, double **forwardlix, double **forwardlmi,
   double **forwardlrg, double **forwardlpk,
   double **backwardrtb, double **backwardrix, double **backwardrmi,
   double **backwardrrg, double **backwardrpk,
   double **backwardltb, double **backwardlix, double **backwardlmi,
   double **backwardlrg, double **backwardlpk,
   double **crossrtb, double **crossrix, double **crossrmi,
   double **crossrrg, double **crossrpk,
   double **crossltb, double **crosslix, double **crosslmi,
   double **crosslrg, double **crosslpk, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl10_product(deg,inputrtb[ix1],inputrix[ix1],inputrmi[ix1],
                               inputrrg[ix1],inputrpk[ix1],
                               inputltb[ix1],inputlix[ix1],inputlmi[ix1],
                               inputlrg[ix1],inputlpk[ix1],
                                 cffrtb[i],    cffrix[i],    cffrmi[i],
                                 cffrrg[i],    cffrpk[i],
                                 cffltb[i],    cfflix[i],    cfflmi[i],
                                 cfflrg[i],    cfflpk[i],
                             forwardrtb[0],forwardrix[0],forwardrmi[0],
                             forwardrrg[0],forwardrpk[0],
                             forwardltb[0],forwardlix[0],forwardlmi[0],
                             forwardlrg[0],forwardlpk[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            daf_inc(&outputrtb[dim][j],&outputrix[dim][j],&outputrmi[dim][j],
                    &outputrrg[dim][j],&outputrpk[dim][j],
                    &outputltb[dim][j],&outputlix[dim][j],&outputlmi[dim][j],
                    &outputlrg[dim][j],&outputlpk[dim][j],
                    forwardrtb[0][j],  forwardrix[0][j],  forwardrmi[0][j],
                    forwardrrg[0][j],  forwardrpk[0][j],
                    forwardltb[0][j],  forwardlix[0][j],  forwardlmi[0][j],
                    forwardlrg[0][j],  forwardlpk[0][j]);
            // output[ix1][j] += cff[i][j];
            daf_inc(&outputrtb[ix1][j],&outputrix[ix1][j],&outputrmi[ix1][j],
                    &outputrrg[ix1][j],&outputrpk[ix1][j],
                    &outputltb[ix1][j],&outputlix[ix1][j],&outputlmi[ix1][j],
                    &outputlrg[ix1][j],&outputlpk[ix1][j],
                        cffrtb[i][j],      cffrix[i][j],      cffrmi[i][j],
                        cffrrg[i][j],      cffrpk[i][j],
                        cffltb[i][j],      cfflix[i][j],      cfflmi[i][j],
                        cfflrg[i][j],      cfflpk[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl10_product(deg,cffrtb[i],    cffrix[i],    cffrmi[i],
                               cffrrg[i],    cffrpk[i],
                               cffltb[i],    cfflix[i],    cfflmi[i],
                               cfflrg[i],    cfflpk[i],
                             inputrtb[ix1],inputrix[ix1],inputrmi[ix1],
                             inputrrg[ix1],inputrpk[ix1],
                             inputltb[ix1],inputlix[ix1],inputlmi[ix1],
                             inputlrg[ix1],inputlpk[ix1],
                           forwardrtb[0],forwardrix[0],forwardrmi[0],
                           forwardrrg[0],forwardrpk[0],
                           forwardltb[0],forwardlix[0],forwardlmi[0],
                           forwardlrg[0],forwardlpk[0]);
         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
            daf_inc
               (&outputrtb[ix2][j],&outputrix[ix2][j],&outputrmi[ix2][j],
                &outputrrg[ix2][j],&outputrpk[ix2][j],
                &outputltb[ix2][j],&outputlix[ix2][j],&outputlmi[ix2][j],
                &outputlrg[ix2][j],&outputlpk[ix2][j],
                forwardrtb[0][j],  forwardrix[0][j],  forwardrmi[0][j],
                forwardrrg[0][j],  forwardrpk[0][j],
                forwardltb[0][j],  forwardlix[0][j],  forwardlmi[0][j],
                forwardlrg[0][j],  forwardlpk[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl10_product(deg,cffrtb[i],     cffrix[i],     cffrmi[i],
                               cffrrg[i],     cffrpk[i],
                               cffltb[i],     cfflix[i],     cfflmi[i],
                               cfflrg[i],     cfflpk[i],
                             inputrtb[ix2], inputrix[ix2], inputrmi[ix2],
                             inputrrg[ix2], inputrpk[ix2],
                             inputltb[ix2], inputlix[ix2], inputlmi[ix2],
                             inputlrg[ix2], inputlpk[ix2],
                          backwardrtb[0],backwardrix[0],backwardrmi[0],
                          backwardrrg[0],backwardrpk[0],
                          backwardltb[0],backwardlix[0],backwardlmi[0],
                          backwardlrg[0],backwardlpk[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
            daf_inc
               ( &outputrtb[ix1][j],&outputrix[ix1][j],&outputrmi[ix1][j],
                 &outputrrg[ix1][j],&outputrpk[ix1][j],
                 &outputltb[ix1][j],&outputlix[ix1][j],&outputlmi[ix1][j],
                 &outputlrg[ix1][j],&outputlpk[ix1][j],
                backwardrtb[0][j], backwardrix[0][j], backwardrmi[0][j],
                backwardrrg[0][j], backwardrpk[0][j],
                backwardltb[0][j], backwardlix[0][j], backwardlmi[0][j],
                backwardlrg[0][j], backwardlpk[0][j]);

         CPU_dbl10_product(deg,forwardrtb[0],forwardrix[0],forwardrmi[0],
                               forwardrrg[0],forwardrpk[0],
                               forwardltb[0],forwardlix[0],forwardlmi[0],
                               forwardlrg[0],forwardlpk[0],
                                 inputrtb[ix2],inputrix[ix2],inputrmi[ix2],
                                 inputrrg[ix2],inputrpk[ix2],
                                 inputltb[ix2],inputlix[ix2],inputlmi[ix2],
                                 inputlrg[ix2],inputlpk[ix2],
                               forwardrtb[1],forwardrix[1],forwardrmi[1],
                               forwardrrg[1],forwardrpk[1],
                               forwardltb[1],forwardlix[1],forwardlmi[1],
                               forwardlrg[1],forwardlpk[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
            daf_inc
               (&outputrtb[dim][j],&outputrix[dim][j],&outputrmi[dim][j],
                &outputrrg[dim][j],&outputrpk[dim][j],
                &outputltb[dim][j],&outputlix[dim][j],&outputlmi[dim][j],
                &outputlrg[dim][j],&outputlpk[dim][j],
                forwardrtb[1][j],  forwardrix[1][j],  forwardrmi[1][j],
                forwardrrg[1][j],  forwardrpk[1][j],
                forwardltb[1][j],  forwardlix[1][j],  forwardlmi[1][j],
                forwardlrg[1][j],  forwardlpk[1][j]);
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl10_speel(nvr[i],deg,idx[i],
                 cffrtb[i],  cffrix[i],  cffrmi[i],  cffrrg[i],  cffrpk[i],
                 cffltb[i],  cfflix[i],  cfflmi[i],  cfflrg[i],  cfflpk[i],
               inputrtb,   inputrix,   inputrmi,   inputrrg,   inputrpk,
               inputltb,   inputlix,   inputlmi,   inputlrg,   inputlpk,
             forwardrtb, forwardrix, forwardrmi, forwardrrg, forwardrpk,
             forwardltb, forwardlix, forwardlmi, forwardlrg, forwardlpk,
            backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
            backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,
               crossrtb,   crossrix,   crossrmi,   crossrrg,   crossrpk,
               crossltb,   crosslix,   crosslmi,   crosslrg,   crosslpk);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
            daf_inc
               (&outputrtb[dim][j],&outputrix[dim][j],&outputrmi[dim][j],
                &outputrrg[dim][j],&outputrpk[dim][j],
                &outputltb[dim][j],&outputlix[dim][j],&outputlmi[dim][j],
                &outputlrg[dim][j],&outputlpk[dim][j],
                forwardrtb[ix1][j],forwardrix[ix1][j],forwardrmi[ix1][j],
                forwardrrg[ix1][j],forwardrpk[ix1][j],
                forwardltb[ix1][j],forwardlix[ix1][j],forwardlmi[ix1][j],
                forwardlrg[ix1][j],forwardlpk[ix1][j]);

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
            daf_inc
               (&outputrtb[ix2][j],&outputrix[ix2][j],&outputrmi[ix2][j],
                &outputrrg[ix2][j],&outputrpk[ix2][j],
                &outputltb[ix2][j],&outputlix[ix2][j],&outputlmi[ix2][j],
                &outputlrg[ix2][j],&outputlpk[ix2][j],
                forwardrtb[ix1][j],forwardrix[ix1][j],forwardrmi[ix1][j],
                forwardrrg[ix1][j],forwardrpk[ix1][j],
                forwardltb[ix1][j],forwardlix[ix1][j],forwardlmi[ix1][j],
                forwardlrg[ix1][j],forwardlpk[ix1][j]);

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
            daf_inc
               ( &outputrtb[ix2][j], &outputrix[ix2][j],  &outputrmi[ix2][j],
                 &outputrrg[ix2][j], &outputrpk[ix2][j],
                 &outputltb[ix2][j], &outputlix[ix2][j],  &outputlmi[ix2][j],
                 &outputlrg[ix2][j], &outputlpk[ix2][j],
                backwardrtb[ix1][j],backwardrix[ix1][j], backwardrmi[ix1][j],
                backwardrrg[ix1][j],backwardrpk[ix1][j],
                backwardltb[ix1][j],backwardlix[ix1][j], backwardlmi[ix1][j],
                backwardlrg[ix1][j],backwardlpk[ix1][j]);

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
               daf_inc
                  (&outputrtb[ix2][j],&outputrix[ix2][j],&outputrmi[ix2][j],
                   &outputrrg[ix2][j],&outputrpk[ix2][j],
                   &outputltb[ix2][j],&outputlix[ix2][j],&outputlmi[ix2][j],
                   &outputlrg[ix2][j],&outputlpk[ix2][j],
                     crossrtb[k-1][j],  crossrix[k-1][j],  crossrmi[k-1][j],
                     crossrrg[k-1][j],  crossrpk[k-1][j],
                     crossltb[k-1][j],  crosslix[k-1][j],  crosslmi[k-1][j],
                     crosslrg[k-1][j],  crosslpk[k-1][j]);
         }
      }
   }
}

void CPU_dbl10_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk, 
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk, 
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   double *elapsedsec, bool verbose )
{
   double **forwardrtb = new double*[dim];
   double **forwardrix = new double*[dim];
   double **forwardrmi = new double*[dim];
   double **forwardrrg = new double*[dim];
   double **forwardrpk = new double*[dim];
   double **forwardltb = new double*[dim];
   double **forwardlix = new double*[dim];
   double **forwardlmi = new double*[dim];
   double **forwardlrg = new double*[dim];
   double **forwardlpk = new double*[dim];
   double **backwardrtb = new double*[dim-1]; // in case dim = 2
   double **backwardrix = new double*[dim-1];
   double **backwardrmi = new double*[dim-1];
   double **backwardrrg = new double*[dim-1];
   double **backwardrpk = new double*[dim-1];
   double **backwardltb = new double*[dim-1];
   double **backwardlix = new double*[dim-1];
   double **backwardlmi = new double*[dim-1];
   double **backwardlrg = new double*[dim-1];
   double **backwardlpk = new double*[dim-1];
   double **crossrtb = new double*[dim-1];    // in case dim = 2
   double **crossrix = new double*[dim-1];
   double **crossrmi = new double*[dim-1];
   double **crossrrg = new double*[dim-1];
   double **crossrpk = new double*[dim-1];
   double **crossltb = new double*[dim-1];
   double **crosslix = new double*[dim-1];
   double **crosslmi = new double*[dim-1];
   double **crosslrg = new double*[dim-1];
   double **crosslpk = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardrtb[i] = new double[deg+1];
      forwardrix[i] = new double[deg+1];
      forwardrmi[i] = new double[deg+1];
      forwardrrg[i] = new double[deg+1];
      forwardrpk[i] = new double[deg+1];
      forwardltb[i] = new double[deg+1];
      forwardlix[i] = new double[deg+1];
      forwardlmi[i] = new double[deg+1];
      forwardlrg[i] = new double[deg+1];
      forwardlpk[i] = new double[deg+1];
      backwardrtb[i] = new double[deg+1];
      backwardrix[i] = new double[deg+1];
      backwardrmi[i] = new double[deg+1];
      backwardrrg[i] = new double[deg+1];
      backwardrpk[i] = new double[deg+1];
      backwardltb[i] = new double[deg+1];
      backwardlix[i] = new double[deg+1];
      backwardlmi[i] = new double[deg+1];
      backwardlrg[i] = new double[deg+1];
      backwardlpk[i] = new double[deg+1];
      crossrtb[i] = new double[deg+1];
      crossrix[i] = new double[deg+1];
      crossrmi[i] = new double[deg+1];
      crossrrg[i] = new double[deg+1];
      crossrpk[i] = new double[deg+1];
      crossltb[i] = new double[deg+1];
      crosslix[i] = new double[deg+1];
      crosslmi[i] = new double[deg+1];
      crosslrg[i] = new double[deg+1];
      crosslpk[i] = new double[deg+1];
   }
   forwardrtb[dim-1] = new double[deg+1];
   forwardrix[dim-1] = new double[deg+1];
   forwardrmi[dim-1] = new double[deg+1];
   forwardrrg[dim-1] = new double[deg+1];
   forwardrpk[dim-1] = new double[deg+1];
   forwardltb[dim-1] = new double[deg+1];
   forwardlix[dim-1] = new double[deg+1];
   forwardlmi[dim-1] = new double[deg+1];
   forwardlrg[dim-1] = new double[deg+1];
   forwardlpk[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputrtb[dim][i] = cstrtb[i];
      outputrix[dim][i] = cstrix[i];
      outputrmi[dim][i] = cstrmi[i];
      outputrrg[dim][i] = cstrrg[i];
      outputrpk[dim][i] = cstrpk[i];
      outputltb[dim][i] = cstltb[i];
      outputlix[dim][i] = cstlix[i];
      outputlmi[dim][i] = cstlmi[i];
      outputlrg[dim][i] = cstlrg[i];
      outputlpk[dim][i] = cstlpk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrtb[i][j] = 0.0;
         outputrix[i][j] = 0.0;
         outputrmi[i][j] = 0.0;
         outputrrg[i][j] = 0.0;
         outputrpk[i][j] = 0.0;
         outputltb[i][j] = 0.0;
         outputlix[i][j] = 0.0;
         outputlmi[i][j] = 0.0;
         outputlrg[i][j] = 0.0;
         outputlpk[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl10_poly_speel
      (dim,nbr,deg,nvr,idx,
            cffrtb,     cffrix,     cffrmi,     cffrrg,     cffrpk,
            cffltb,     cfflix,     cfflmi,     cfflrg,     cfflpk,
          inputrtb,   inputrix,   inputrmi,   inputrrg,   inputrpk,
          inputltb,   inputlix,   inputlmi,   inputlrg,   inputlpk,
         outputrtb,  outputrix,  outputrmi,  outputrrg,  outputrpk,
         outputltb,  outputlix,  outputlmi,  outputlrg,  outputlpk,
        forwardrtb, forwardrix, forwardrmi, forwardrrg, forwardrpk,
        forwardltb, forwardlix, forwardlmi, forwardlrg, forwardlpk,
       backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
       backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,
          crossrtb,   crossrix,   crossrmi,   crossrrg,   crossrpk,
          crossltb,   crosslix,   crosslmi,   crosslrg,   crosslpk,verbose);
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
      free(forwardrtb[i]); free(backwardrtb[i]); free(crossrtb[i]);
      free(forwardrix[i]); free(backwardrix[i]); free(crossrix[i]);
      free(forwardrmi[i]); free(backwardrmi[i]); free(crossrmi[i]);
      free(forwardrrg[i]); free(backwardrrg[i]); free(crossrrg[i]);
      free(forwardrpk[i]); free(backwardrpk[i]); free(crossrpk[i]);
      free(forwardltb[i]); free(backwardltb[i]); free(crossltb[i]);
      free(forwardlix[i]); free(backwardlix[i]); free(crosslix[i]);
      free(forwardlmi[i]); free(backwardlmi[i]); free(crosslmi[i]);
      free(forwardlrg[i]); free(backwardlrg[i]); free(crosslrg[i]);
      free(forwardlpk[i]); free(backwardlpk[i]); free(crosslpk[i]);
   }
   free(forwardrtb[dim-1]); free(forwardrix[dim-1]); free(forwardrmi[dim-1]);
   free(forwardrrg[dim-1]); free(forwardrpk[dim-1]);
   free(forwardltb[dim-1]); free(forwardlix[dim-1]); free(forwardlmi[dim-1]);
   free(forwardlrg[dim-1]); free(forwardlpk[dim-1]);
   free(forwardrtb); free(backwardrtb); free(crossrtb);
   free(forwardrix); free(backwardrix); free(crossrix);
   free(forwardrmi); free(backwardrmi); free(crossrmi);
   free(forwardrrg); free(backwardrrg); free(crossrrg);
   free(forwardrpk); free(backwardrpk); free(crossrpk);
   free(forwardltb); free(backwardltb); free(crossltb);
   free(forwardlix); free(backwardlix); free(crosslix);
   free(forwardlmi); free(backwardlmi); free(crosslmi);
   free(forwardlrg); free(backwardlrg); free(crosslrg);
   free(forwardlpk); free(backwardlpk); free(crosslpk);
}

void CPU_dbl10_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi,
   double *cffrrg, double *cffrpk,
   double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **forwardrtb, double **forwardrix, double **forwardrmi,
   double **forwardrrg, double **forwardrpk,
   double **forwardltb, double **forwardlix, double **forwardlmi,
   double **forwardlrg, double **forwardlpk,
   double **backwardrtb, double **backwardrix, double **backwardrmi,
   double **backwardrrg, double **backwardrpk,
   double **backwardltb, double **backwardlix, double **backwardlmi,
   double **backwardlrg, double **backwardlpk,
   double **crossrtb, double **crossrix, double **crossrmi,
   double **crossrrg, double **crossrpk,
   double **crossltb, double **crosslix, double **crosslmi,
   double **crosslrg, double **crosslpk,
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
         CPU_dbl10_product(deg,
                cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
                cffltb,cfflix,cfflmi,cfflrg,cfflpk,
              inputrtb[inp2ix],  inputrix[inp2ix],  inputrmi[inp2ix],
              inputrrg[inp2ix],  inputrpk[inp2ix],
              inputltb[inp2ix],  inputlix[inp2ix],  inputlmi[inp2ix],
              inputlrg[inp2ix],  inputlpk[inp2ix],
            forwardrtb[outidx],forwardrix[outidx],forwardrmi[outidx],
            forwardrrg[outidx],forwardrpk[outidx],
            forwardltb[outidx],forwardlix[outidx],forwardlmi[outidx],
            forwardlrg[outidx],forwardlpk[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl10_product(deg,
               inputrtb[inp1ix],inputrix[inp1ix],inputrmi[inp1ix],
               inputrrg[inp1ix],inputrpk[inp1ix],
               inputltb[inp1ix],inputlix[inp1ix],inputlmi[inp1ix],
               inputlrg[inp1ix],inputlpk[inp1ix],
               cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
               cffltb,cfflix,cfflmi,cfflrg,cfflpk,
               forwardrtb[outidx],forwardrix[outidx],forwardrmi[outidx],
               forwardrrg[outidx],forwardrpk[outidx],
               forwardltb[outidx],forwardlix[outidx],forwardlmi[outidx],
               forwardlrg[outidx],forwardlpk[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl10_product(deg,
                 inputrtb[inp1ix],  inputrix[inp1ix],  inputrmi[inp1ix],
                 inputrrg[inp1ix],  inputrpk[inp1ix],
                 inputltb[inp1ix],  inputlix[inp1ix],  inputlmi[inp1ix],
                 inputlrg[inp1ix],  inputlpk[inp1ix],
               forwardrtb[inp2ix],forwardrix[inp2ix],forwardrmi[inp2ix],
               forwardrrg[inp2ix],forwardrpk[inp2ix],
               forwardltb[inp2ix],forwardlix[inp2ix],forwardlmi[inp2ix],
               forwardlrg[inp2ix],forwardlpk[inp2ix],
               forwardrtb[outidx],forwardrix[outidx],forwardrmi[outidx],
               forwardrrg[outidx],forwardrpk[outidx],
               forwardltb[outidx],forwardlix[outidx],forwardlmi[outidx],
               forwardlrg[outidx],forwardlpk[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl10_product(deg,
              crossrtb[inp1ix],  crossrix[inp1ix],  crossrmi[inp1ix],
              crossrrg[inp1ix],  crossrpk[inp1ix],
              crossltb[inp1ix],  crosslix[inp1ix],  crosslmi[inp1ix],
              crosslrg[inp1ix],  crosslpk[inp1ix],
              inputrtb[inp2ix],  inputrix[inp2ix],  inputrmi[inp2ix],
              inputrrg[inp2ix],  inputrpk[inp2ix],
              inputltb[inp2ix],  inputlix[inp2ix],  inputlmi[inp2ix],
              inputlrg[inp2ix],  inputlpk[inp2ix],
            forwardrtb[outidx],forwardrix[outidx],forwardrmi[outidx],
            forwardrrg[outidx],forwardrpk[outidx],
            forwardltb[outidx],forwardlix[outidx],forwardlmi[outidx],
            forwardlrg[outidx],forwardlpk[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl10_product(deg,
               inputrtb[inp1ix],inputrix[inp1ix],inputrmi[inp1ix],
               inputrrg[inp1ix],inputrpk[inp1ix],
               inputltb[inp1ix],inputlix[inp1ix],inputlmi[inp1ix],
               inputlrg[inp1ix],inputlpk[inp1ix],
               cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
               cffltb,cfflix,cfflmi,cfflrg,cfflpk,
               forwardrtb[outidx],forwardrix[outidx],forwardrmi[outidx],
               forwardrrg[outidx],forwardrpk[outidx],
               forwardltb[outidx],forwardlix[outidx],forwardlmi[outidx],
               forwardlrg[outidx],forwardlpk[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl10_product(deg,
               forwardrtb[inp1ix],forwardrix[inp1ix],forwardrmi[inp1ix],
               forwardrrg[inp1ix],forwardrpk[inp1ix],
               forwardltb[inp1ix],forwardlix[inp1ix],forwardlmi[inp1ix],
               forwardlrg[inp1ix],forwardlpk[inp1ix],
                 inputrtb[inp2ix],  inputrix[inp2ix],  inputrmi[inp2ix],
                 inputrrg[inp2ix],  inputrpk[inp2ix],
                 inputltb[inp2ix],  inputlix[inp2ix],  inputlmi[inp2ix],
                 inputlrg[inp2ix],  inputlpk[inp2ix],
               forwardrtb[outidx],forwardrix[outidx],forwardrmi[outidx],
               forwardrrg[outidx],forwardrpk[outidx],
               forwardltb[outidx],forwardlix[outidx],forwardlmi[outidx],
               forwardlrg[outidx],forwardlpk[outidx]);
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
            CPU_dbl10_product(deg,
                    cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
                    cffltb,cfflix,cfflmi,cfflrg,cfflpk,
                  inputrtb[inp2ix],   inputrix[inp2ix],   inputrmi[inp2ix],
                  inputrrg[inp2ix],   inputrpk[inp2ix],
                  inputltb[inp2ix],   inputlix[inp2ix],   inputlmi[inp2ix],
                  inputlrg[inp2ix],   inputlpk[inp2ix],
               backwardrtb[outidx],backwardrix[outidx],backwardrmi[outidx],
               backwardrrg[outidx],backwardrpk[outidx],
               backwardltb[outidx],backwardlix[outidx],backwardlmi[outidx],
               backwardlrg[outidx],backwardlpk[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl10_product(deg,
               cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
               cffltb,cfflix,cfflmi,cfflrg,cfflpk,
               backwardrtb[inp2ix],backwardrix[inp2ix],backwardrmi[inp2ix],
               backwardrrg[inp2ix],backwardrpk[inp2ix],
               backwardltb[inp2ix],backwardlix[inp2ix],backwardlmi[inp2ix],
               backwardlrg[inp2ix],backwardlpk[inp2ix],
               backwardrtb[outidx],backwardrix[outidx],backwardrmi[outidx],
               backwardrrg[outidx],backwardrpk[outidx],
               backwardltb[outidx],backwardlix[outidx],backwardlmi[outidx],
               backwardlrg[outidx],backwardlpk[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl10_product(deg,
                  inputrtb[inp1ix],   inputrix[inp1ix],   inputrmi[inp1ix],
                  inputrrg[inp1ix],   inputrpk[inp1ix],
                  inputltb[inp1ix],   inputlix[inp1ix],   inputlmi[inp1ix],
                  inputlrg[inp1ix],   inputlpk[inp1ix],
                  inputrtb[inp2ix],   inputrix[inp2ix],   inputrmi[inp2ix],
                  inputrrg[inp2ix],   inputrpk[inp2ix],
                  inputltb[inp2ix],   inputlix[inp2ix],   inputlmi[inp2ix],
                  inputlrg[inp2ix],   inputlpk[inp2ix],
               backwardrtb[outidx],backwardrix[outidx],backwardrmi[outidx],
               backwardrrg[outidx],backwardrpk[outidx],
               backwardltb[outidx],backwardlix[outidx],backwardlmi[outidx],
               backwardlrg[outidx],backwardlpk[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl10_product(deg,
                  inputrtb[inp1ix],   inputrix[inp1ix],   inputrmi[inp1ix],
                  inputrrg[inp1ix],   inputrpk[inp1ix],
                  inputltb[inp1ix],   inputlix[inp1ix],   inputlmi[inp1ix],
                  inputlrg[inp1ix],   inputlpk[inp1ix],
               backwardrtb[inp2ix],backwardrix[inp2ix],backwardrmi[inp2ix],
               backwardrrg[inp2ix],backwardrpk[inp2ix],
               backwardltb[inp2ix],backwardlix[inp2ix],backwardlmi[inp2ix],
               backwardlrg[inp2ix],backwardlpk[inp2ix],
               backwardrtb[outidx],backwardrix[outidx],backwardrmi[outidx],
               backwardrrg[outidx],backwardrpk[outidx],
               backwardltb[outidx],backwardlix[outidx],backwardlmi[outidx],
               backwardlrg[outidx],backwardlpk[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl10_product(deg,
               backwardrtb[inp1ix],backwardrix[inp1ix],backwardrmi[inp1ix],
               backwardrrg[inp1ix],backwardrpk[inp1ix],
               backwardltb[inp1ix],backwardlix[inp1ix],backwardlmi[inp1ix],
               backwardlrg[inp1ix],backwardlpk[inp1ix],
               cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
               cffltb,cfflix,cfflmi,cfflrg,cfflpk,
               backwardrtb[outidx],backwardrix[outidx],backwardrmi[outidx],
               backwardrrg[outidx],backwardrpk[outidx],
               backwardltb[outidx],backwardlix[outidx],backwardlmi[outidx],
               backwardlrg[outidx],backwardlpk[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl10_product(deg,
               backwardrtb[inp1ix],backwardrix[inp1ix],backwardrmi[inp1ix],
               backwardrrg[inp1ix],backwardrpk[inp1ix],
               backwardltb[inp1ix],backwardlix[inp1ix],backwardlmi[inp1ix],
               backwardlrg[inp1ix],backwardlpk[inp1ix],
                  inputrtb[inp2ix],   inputrix[inp2ix],   inputrmi[inp2ix],
                  inputrrg[inp2ix],   inputrpk[inp2ix],
                  inputltb[inp2ix],   inputlix[inp2ix],   inputlmi[inp2ix],
                  inputlrg[inp2ix],   inputlpk[inp2ix],
               backwardrtb[outidx],backwardrix[outidx],backwardrmi[outidx],
               backwardrrg[outidx],backwardrpk[outidx],
               backwardltb[outidx],backwardlix[outidx],backwardlmi[outidx],
               backwardlrg[outidx],backwardlpk[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl10_product(deg,
            cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
            cffltb,cfflix,cfflmi,cfflrg,cfflpk,
            inputrtb[inp2ix],inputrix[inp2ix],inputrmi[inp2ix],
            inputrrg[inp2ix],inputrpk[inp2ix],
            inputltb[inp2ix],inputlix[inp2ix],inputlmi[inp2ix],
            inputlrg[inp2ix],inputlpk[inp2ix],
            crossrtb[outidx],crossrix[outidx],crossrmi[outidx],
            crossrrg[outidx],crossrpk[outidx],
            crossltb[outidx],crosslix[outidx],crosslmi[outidx],
            crosslrg[outidx],crosslpk[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl10_product(deg,
              inputrtb[inp1ix],  inputrix[inp1ix],  inputrmi[inp1ix],
              inputrrg[inp1ix],  inputrpk[inp1ix],
              inputltb[inp1ix],  inputlix[inp1ix],  inputlmi[inp1ix],
              inputlrg[inp1ix],  inputlpk[inp1ix],
            forwardrtb[inp2ix],forwardrix[inp2ix],forwardrmi[inp2ix],
            forwardrrg[inp2ix],forwardrpk[inp2ix],
            forwardltb[inp2ix],forwardlix[inp2ix],forwardlmi[inp2ix],
            forwardlrg[inp2ix],forwardlpk[inp2ix],
              crossrtb[outidx],  crossrix[outidx],  crossrmi[outidx],
              crossrrg[outidx],  crossrpk[outidx],
              crossltb[outidx],  crosslix[outidx],  crosslmi[outidx],
              crosslrg[outidx],  crosslpk[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl10_product(deg,
              forwardrtb[inp1ix],forwardrix[inp1ix],forwardrmi[inp1ix],
              forwardrrg[inp1ix],forwardrpk[inp1ix],
              forwardltb[inp1ix],forwardlix[inp1ix],forwardlmi[inp1ix],
              forwardlrg[inp1ix],forwardlpk[inp1ix],
                inputrtb[inp2ix],  inputrix[inp2ix],  inputrmi[inp2ix],
                inputrrg[inp2ix],  inputrpk[inp2ix],
                inputltb[inp2ix],  inputlix[inp2ix],  inputlmi[inp2ix],
                inputlrg[inp2ix],  inputlpk[inp2ix],
                crossrtb[outidx],  crossrix[outidx],  crossrmi[outidx],
                crossrrg[outidx],  crossrpk[outidx],
                crossltb[outidx],  crosslix[outidx],  crosslmi[outidx],
                crosslrg[outidx],  crosslpk[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl10_product(deg,
               forwardrtb[inp1ix], forwardrix[inp1ix], forwardrmi[inp1ix],
               forwardrrg[inp1ix], forwardrpk[inp1ix],
               forwardltb[inp1ix], forwardlix[inp1ix], forwardlmi[inp1ix],
               forwardlrg[inp1ix], forwardlpk[inp1ix],
              backwardrtb[inp2ix],backwardrix[inp2ix],backwardrmi[inp2ix],
              backwardrrg[inp2ix],backwardrpk[inp2ix],
              backwardltb[inp2ix],backwardlix[inp2ix],backwardlmi[inp2ix],
              backwardlrg[inp2ix],backwardlpk[inp2ix],
                 crossrtb[outidx],   crossrix[outidx],   crossrmi[outidx],
                 crossrrg[outidx],   crossrpk[outidx],
                 crossltb[outidx],   crosslix[outidx],   crosslmi[outidx],
                 crosslrg[outidx],   crosslpk[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl10_product(deg,
            backwardrtb[inp1ix],backwardrix[inp1ix],backwardrmi[inp1ix],
            backwardrrg[inp1ix],backwardrpk[inp1ix],
            backwardltb[inp1ix],backwardlix[inp1ix],backwardlmi[inp1ix],
            backwardlrg[inp1ix],backwardlpk[inp1ix],
             forwardrtb[inp2ix], forwardrix[inp2ix], forwardrmi[inp2ix],
             forwardrrg[inp2ix], forwardrpk[inp2ix],
             forwardltb[inp2ix], forwardlix[inp2ix], forwardlmi[inp2ix],
             forwardlrg[inp2ix], forwardlpk[inp2ix],
               crossrtb[outidx],   crossrix[outidx],   crossrmi[outidx],
               crossrrg[outidx],   crossrpk[outidx],
               crossltb[outidx],   crosslix[outidx],   crosslmi[outidx],
               crosslrg[outidx],   crosslpk[outidx]);
      }
   }
}

void CPU_dbl10_add_job
 ( int deg,
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double ***forwardrtb, double ***forwardrix, double ***forwardrmi,
   double ***forwardrrg, double ***forwardrpk,
   double ***forwardltb, double ***forwardlix, double ***forwardlmi,
   double ***forwardlrg, double ***forwardlpk,
   double ***backwardrtb, double ***backwardrix, double ***backwardrmi,
   double ***backwardrrg, double ***backwardrpk, 
   double ***backwardltb, double ***backwardlix, double ***backwardlmi,
   double ***backwardlrg, double ***backwardlpk, 
   double ***crossrtb, double ***crossrix, double ***crossrmi,
   double ***crossrrg, double ***crossrpk,
   double ***crossltb, double ***crosslix, double ***crosslmi,
   double ***crosslrg, double ***crosslpk,
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
               daf_inc(&forwardrtb[updmon][updidx][i],
                       &forwardrix[updmon][updidx][i],
                       &forwardrmi[updmon][updidx][i],
                       &forwardrrg[updmon][updidx][i],
                       &forwardrpk[updmon][updidx][i],
                       &forwardltb[updmon][updidx][i],
                       &forwardlix[updmon][updidx][i],
                       &forwardlmi[updmon][updidx][i],
                       &forwardlrg[updmon][updidx][i],
                       &forwardlpk[updmon][updidx][i],
                       cstrtb[i],cstrix[i],cstrmi[i],cstrrg[i],cstrpk[i],
                       cstltb[i],cstlix[i],cstlmi[i],cstlrg[i],cstlpk[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               daf_inc(&forwardrtb[updmon][updidx][i],
                       &forwardrix[updmon][updidx][i],
                       &forwardrmi[updmon][updidx][i],
                       &forwardrrg[updmon][updidx][i],
                       &forwardrpk[updmon][updidx][i],
                       &forwardltb[updmon][updidx][i],
                       &forwardlix[updmon][updidx][i],
                       &forwardlmi[updmon][updidx][i],
                       &forwardlrg[updmon][updidx][i],
                       &forwardlpk[updmon][updidx][i],
                       cffrtb[incidx][i],cffrix[incidx][i],cffrmi[incidx][i],
                       cffrrg[incidx][i],cffrpk[incidx][i],
                       cffltb[incidx][i],cfflix[incidx][i],cfflmi[incidx][i],
                       cfflrg[incidx][i],cfflpk[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            daf_inc(&forwardrtb[updmon][updidx][i],
                    &forwardrix[updmon][updidx][i],
                    &forwardrmi[updmon][updidx][i],
                    &forwardrrg[updmon][updidx][i],
                    &forwardrpk[updmon][updidx][i],
                    &forwardltb[updmon][updidx][i],
                    &forwardlix[updmon][updidx][i],
                    &forwardlmi[updmon][updidx][i],
                    &forwardlrg[updmon][updidx][i],
                    &forwardlpk[updmon][updidx][i],
                    forwardrtb[incmon][incidx][i],
                    forwardrix[incmon][incidx][i],
                    forwardrmi[incmon][incidx][i],
                    forwardrrg[incmon][incidx][i],
                    forwardrpk[incmon][incidx][i],
                    forwardltb[incmon][incidx][i],
                    forwardlix[incmon][incidx][i],
                    forwardlmi[incmon][incidx][i],
                    forwardlrg[incmon][incidx][i],
                    forwardlpk[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            daf_inc(&forwardrtb[updmon][updidx][i],
                    &forwardrix[updmon][updidx][i],
                    &forwardrmi[updmon][updidx][i],
                    &forwardrrg[updmon][updidx][i],
                    &forwardrpk[updmon][updidx][i],
                    &forwardltb[updmon][updidx][i],
                    &forwardlix[updmon][updidx][i],
                    &forwardlmi[updmon][updidx][i],
                    &forwardlrg[updmon][updidx][i],
                    &forwardlpk[updmon][updidx][i],
                    backwardrtb[incmon][incidx][i],
                    backwardrix[incmon][incidx][i],
                    backwardrmi[incmon][incidx][i],
                    backwardrrg[incmon][incidx][i],
                    backwardrpk[incmon][incidx][i],
                    backwardltb[incmon][incidx][i],
                    backwardlix[incmon][incidx][i],
                    backwardlmi[incmon][incidx][i],
                    backwardlrg[incmon][incidx][i],
                    backwardlpk[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            daf_inc(&forwardrtb[updmon][updidx][i],
                    &forwardrix[updmon][updidx][i],
                    &forwardrmi[updmon][updidx][i],
                    &forwardrrg[updmon][updidx][i],
                    &forwardrpk[updmon][updidx][i],
                    &forwardltb[updmon][updidx][i],
                    &forwardlix[updmon][updidx][i],
                    &forwardlmi[updmon][updidx][i],
                    &forwardlrg[updmon][updidx][i],
                    &forwardlpk[updmon][updidx][i],
                    crossrtb[incmon][incidx][i],
                    crossrix[incmon][incidx][i],
                    crossrmi[incmon][incidx][i],
                    crossrrg[incmon][incidx][i],
                    crossrpk[incmon][incidx][i],
                    crossltb[incmon][incidx][i],
                    crosslix[incmon][incidx][i],
                    crosslmi[incmon][incidx][i],
                    crosslrg[incmon][incidx][i],
                    crosslpk[incmon][incidx][i]);
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
            daf_inc(&backwardrtb[updmon][updidx][i],
                    &backwardrix[updmon][updidx][i],
                    &backwardrmi[updmon][updidx][i],
                    &backwardrrg[updmon][updidx][i],
                    &backwardrpk[updmon][updidx][i],
                    &backwardltb[updmon][updidx][i],
                    &backwardlix[updmon][updidx][i],
                    &backwardlmi[updmon][updidx][i],
                    &backwardlrg[updmon][updidx][i],
                    &backwardlpk[updmon][updidx][i],
                    cffrtb[incidx][i],cffrix[incidx][i],cffrmi[incidx][i],
                    cffrrg[incidx][i],cffrpk[incidx][i],
                    cffltb[incidx][i],cfflix[incidx][i],cfflmi[incidx][i],
                    cfflrg[incidx][i],cfflpk[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            daf_inc(&backwardrtb[updmon][updidx][i],
                    &backwardrix[updmon][updidx][i],
                    &backwardrmi[updmon][updidx][i],
                    &backwardrrg[updmon][updidx][i],
                    &backwardrpk[updmon][updidx][i],
                    &backwardltb[updmon][updidx][i],
                    &backwardlix[updmon][updidx][i],
                    &backwardlmi[updmon][updidx][i],
                    &backwardlrg[updmon][updidx][i],
                    &backwardlpk[updmon][updidx][i],
                    forwardrtb[incmon][incidx][i],
                    forwardrix[incmon][incidx][i],
                    forwardrmi[incmon][incidx][i],
                    forwardrrg[incmon][incidx][i],
                    forwardrpk[incmon][incidx][i],
                    forwardltb[incmon][incidx][i],
                    forwardlix[incmon][incidx][i],
                    forwardlmi[incmon][incidx][i],
                    forwardlrg[incmon][incidx][i],
                    forwardlpk[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            daf_inc(&backwardrtb[updmon][updidx][i],
                    &backwardrix[updmon][updidx][i],
                    &backwardrmi[updmon][updidx][i],
                    &backwardrrg[updmon][updidx][i],
                    &backwardrpk[updmon][updidx][i],
                    &backwardltb[updmon][updidx][i],
                    &backwardlix[updmon][updidx][i],
                    &backwardlmi[updmon][updidx][i],
                    &backwardlrg[updmon][updidx][i],
                    &backwardlpk[updmon][updidx][i],
                    backwardrtb[incmon][incidx][i],
                    backwardrix[incmon][incidx][i],
                    backwardrmi[incmon][incidx][i],
                    backwardrrg[incmon][incidx][i],
                    backwardrpk[incmon][incidx][i],
                    backwardltb[incmon][incidx][i],
                    backwardlix[incmon][incidx][i],
                    backwardlmi[incmon][incidx][i],
                    backwardlrg[incmon][incidx][i],
                    backwardlpk[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            daf_inc(&backwardrtb[updmon][updidx][i],
                    &backwardrix[updmon][updidx][i],
                    &backwardrmi[updmon][updidx][i],
                    &backwardrrg[updmon][updidx][i],
                    &backwardrpk[updmon][updidx][i],
                    &backwardltb[updmon][updidx][i],
                    &backwardlix[updmon][updidx][i],
                    &backwardlmi[updmon][updidx][i],
                    &backwardlrg[updmon][updidx][i],
                    &backwardlpk[updmon][updidx][i],
                    crossrtb[incmon][incidx][i],
                    crossrix[incmon][incidx][i],
                    crossrmi[incmon][incidx][i],
                    crossrrg[incmon][incidx][i],
                    crossrpk[incmon][incidx][i],
                    crossltb[incmon][incidx][i],
                    crosslix[incmon][incidx][i],
                    crosslmi[incmon][incidx][i],
                    crosslrg[incmon][incidx][i],
                    crosslpk[incmon][incidx][i]);
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
            daf_inc(&crossrtb[updmon][updidx][i],
                    &crossrix[updmon][updidx][i],
                    &crossrmi[updmon][updidx][i],
                    &crossrrg[updmon][updidx][i],
                    &crossrpk[updmon][updidx][i],
                    &crossltb[updmon][updidx][i],
                    &crosslix[updmon][updidx][i],
                    &crosslmi[updmon][updidx][i],
                    &crosslrg[updmon][updidx][i],
                    &crosslpk[updmon][updidx][i],
                    cffrtb[incidx][i],cffrix[incidx][i],cffrmi[incidx][i],
                    cffrrg[incidx][i],cffrpk[incidx][i],
                    cffltb[incidx][i],cfflix[incidx][i],cfflmi[incidx][i],
                    cfflrg[incidx][i],cfflpk[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            daf_inc(&crossrtb[updmon][updidx][i],
                    &crossrix[updmon][updidx][i],
                    &crossrmi[updmon][updidx][i],
                    &crossrrg[updmon][updidx][i],
                    &crossrpk[updmon][updidx][i],
                    &crossltb[updmon][updidx][i],
                    &crosslix[updmon][updidx][i],
                    &crosslmi[updmon][updidx][i],
                    &crosslrg[updmon][updidx][i],
                    &crosslpk[updmon][updidx][i],
                    forwardrtb[incmon][incidx][i],
                    forwardrix[incmon][incidx][i],
                    forwardrmi[incmon][incidx][i],
                    forwardrrg[incmon][incidx][i],
                    forwardrpk[incmon][incidx][i],
                    forwardltb[incmon][incidx][i],
                    forwardlix[incmon][incidx][i],
                    forwardlmi[incmon][incidx][i],
                    forwardlrg[incmon][incidx][i],
                    forwardlpk[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            daf_inc(&crossrtb[updmon][updidx][i],
                    &crossrix[updmon][updidx][i],
                    &crossrmi[updmon][updidx][i],
                    &crossrrg[updmon][updidx][i],
                    &crossrpk[updmon][updidx][i],
                    &crossltb[updmon][updidx][i],
                    &crosslix[updmon][updidx][i],
                    &crosslmi[updmon][updidx][i],
                    &crosslrg[updmon][updidx][i],
                    &crosslpk[updmon][updidx][i],
                    backwardrtb[incmon][incidx][i],
                    backwardrix[incmon][incidx][i],
                    backwardrmi[incmon][incidx][i],
                    backwardrrg[incmon][incidx][i],
                    backwardrpk[incmon][incidx][i],
                    backwardltb[incmon][incidx][i],
                    backwardlix[incmon][incidx][i],
                    backwardlmi[incmon][incidx][i],
                    backwardlrg[incmon][incidx][i],
                    backwardlpk[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            daf_inc(&crossrtb[updmon][updidx][i],
                    &crossrix[updmon][updidx][i],
                    &crossrmi[updmon][updidx][i],
                    &crossrrg[updmon][updidx][i],
                    &crossrpk[updmon][updidx][i],
                    &crossltb[updmon][updidx][i],
                    &crosslix[updmon][updidx][i],
                    &crosslmi[updmon][updidx][i],
                    &crosslrg[updmon][updidx][i],
                    &crosslpk[updmon][updidx][i],
                    crossrtb[incmon][incidx][i],
                    crossrix[incmon][incidx][i],
                    crossrmi[incmon][incidx][i],
                    crossrrg[incmon][incidx][i],
                    crossrpk[incmon][incidx][i],
                    crossltb[incmon][incidx][i],
                    crosslix[incmon][incidx][i],
                    crosslmi[incmon][incidx][i],
                    crosslrg[incmon][incidx][i],
                    crosslpk[incmon][incidx][i]);
      }
   }
}

void CPU_dbl10_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk, 
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk, 
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   double ***forwardrtb, double ***forwardrix, double ***forwardrmi,
   double ***forwardrrg, double ***forwardrpk,
   double ***forwardltb, double ***forwardlix, double ***forwardlmi,
   double ***forwardlrg, double ***forwardlpk,
   double ***backwardrtb, double ***backwardrix, double ***backwardrmi,
   double ***backwardrrg, double ***backwardrpk, 
   double ***backwardltb, double ***backwardlix, double ***backwardlmi,
   double ***backwardlrg, double ***backwardlpk, 
   double ***crossrtb, double ***crossrix, double ***crossrmi,
   double ***crossrrg, double ***crossrpk,
   double ***crossltb, double ***crosslix, double ***crosslmi,
   double ***crosslrg, double ***crosslpk )
{
   for(int i=0; i<=deg; i++)
   {
      outputrtb[dim][i] = cstrtb[i];
      outputrix[dim][i] = cstrix[i];
      outputrmi[dim][i] = cstrmi[i];
      outputrrg[dim][i] = cstrrg[i];
      outputrpk[dim][i] = cstrpk[i];
      outputltb[dim][i] = cstltb[i];
      outputlix[dim][i] = cstlix[i];
      outputlmi[dim][i] = cstlmi[i];
      outputlrg[dim][i] = cstlrg[i];
      outputlpk[dim][i] = cstlpk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputrtb[i][j] = 0.0;
         outputrix[i][j] = 0.0;
         outputrmi[i][j] = 0.0;
         outputrrg[i][j] = 0.0;
         outputrpk[i][j] = 0.0;
         outputltb[i][j] = 0.0;
         outputlix[i][j] = 0.0;
         outputlmi[i][j] = 0.0;
         outputlrg[i][j] = 0.0;
         outputlpk[i][j] = 0.0;
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
         daf_inc(&outputrtb[dim][i],   &outputrix[dim][i],
                 &outputrmi[dim][i],
                 &outputrrg[dim][i],   &outputrpk[dim][i],
                 &outputltb[dim][i],   &outputlix[dim][i],
                 &outputlmi[dim][i],
                 &outputlrg[dim][i],   &outputlpk[dim][i],
                 forwardrtb[k][ix1][i],forwardrix[k][ix1][i],
                 forwardrmi[k][ix1][i],
                 forwardrrg[k][ix1][i],forwardrpk[k][ix1][i],
                 forwardltb[k][ix1][i],forwardlix[k][ix1][i],
                 forwardlmi[k][ix1][i],
                 forwardlrg[k][ix1][i],forwardlpk[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            daf_inc(&outputrtb[ix0][i],&outputrix[ix0][i],&outputrmi[ix0][i],
                    &outputrrg[ix0][i],&outputrpk[ix0][i],
                    &outputltb[ix0][i],&outputlix[ix0][i],&outputlmi[ix0][i],
                    &outputlrg[ix0][i],&outputlpk[ix0][i],
                        cffrtb[k][i],      cffrix[k][i],      cffrmi[k][i],
                        cffrrg[k][i],      cffrpk[k][i],
                        cffltb[k][i],      cfflix[k][i],      cfflmi[k][i],
                        cfflrg[k][i],      cfflpk[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            daf_inc(&outputrtb[ixn][i],   &outputrix[ixn][i],
                    &outputrmi[ixn][i],
                    &outputrrg[ixn][i],   &outputrpk[ixn][i],
                    &outputltb[ixn][i],   &outputlix[ixn][i],
                    &outputlmi[ixn][i],
                    &outputlrg[ixn][i],   &outputlpk[ixn][i],
                    forwardrtb[k][ix2][i],forwardrix[k][ix2][i],
                    forwardrmi[k][ix2][i],
                    forwardrrg[k][ix2][i],forwardrpk[k][ix2][i],
                    forwardltb[k][ix2][i],forwardlix[k][ix2][i],
                    forwardlmi[k][ix2][i],
                    forwardlrg[k][ix2][i],forwardlpk[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            daf_inc( &outputrtb[ix0][i],    &outputrix[ix0][i],
                     &outputrmi[ix0][i],
                     &outputrrg[ix0][i],    &outputrpk[ix0][i],
                     &outputltb[ix0][i],    &outputlix[ix0][i],
                     &outputlmi[ix0][i],
                     &outputlrg[ix0][i],    &outputlpk[ix0][i],
                    backwardrtb[k][ix2][i],backwardrix[k][ix2][i],
                    backwardrmi[k][ix2][i],
                    backwardrrg[k][ix2][i],backwardrpk[k][ix2][i],
                    backwardltb[k][ix2][i],backwardlix[k][ix2][i],
                    backwardlmi[k][ix2][i],
                    backwardlrg[k][ix2][i],backwardlpk[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  daf_inc(&outputrtb[ix0][i], &outputrix[ix0][i],
                          &outputrmi[ix0][i],
                          &outputrrg[ix0][i], &outputrpk[ix0][i],
                          &outputltb[ix0][i], &outputlix[ix0][i],
                          &outputlmi[ix0][i],
                          &outputlrg[ix0][i], &outputlpk[ix0][i],
                            crossrtb[k][j-1][i],crossrix[k][j-1][i],
                            crossrmi[k][j-1][i],
                            crossrrg[k][j-1][i],crossrpk[k][j-1][i],
                            crossltb[k][j-1][i],crosslix[k][j-1][i],
                            crosslmi[k][j-1][i],
                            crosslrg[k][j-1][i],crosslpk[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_dbl10_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk, 
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk, 
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   double ***forwardrtb, double ***forwardrix, double ***forwardrmi,
   double ***forwardrrg, double ***forwardrpk,
   double ***forwardltb, double ***forwardlix, double ***forwardlmi,
   double ***forwardlrg, double ***forwardlpk,
   double ***backwardrtb, double ***backwardrix, double ***backwardrmi,
   double ***backwardrrg, double ***backwardrpk, 
   double ***backwardltb, double ***backwardlix, double ***backwardlmi,
   double ***backwardlrg, double ***backwardlpk, 
   double ***crossrtb, double ***crossrix, double ***crossrmi,
   double ***crossrrg, double ***crossrpk,
   double ***crossltb, double ***crosslix, double ***crosslmi,
   double ***crosslrg, double ***crosslpk,
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

         CPU_dbl10_add_job(deg,
                 cstrtb,     cstrix,     cstrmi,     cstrrg,     cstrpk,
                 cstltb,     cstlix,     cstlmi,     cstlrg,     cstlpk,
                 cffrtb,     cffrix,     cffrmi,     cffrrg,     cffrpk,
                 cffltb,     cfflix,     cfflmi,     cfflrg,     cfflpk,
             forwardrtb, forwardrix, forwardrmi, forwardrrg, forwardrpk,
             forwardltb, forwardlix, forwardlmi, forwardlrg, forwardlpk,
            backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
            backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,
               crossrtb,   crossrix,   crossrmi,   crossrrg,   crossrpk,
               crossltb,   crosslix,   crosslmi,   crosslrg,   crosslpk,
            job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputrtb[dim][i] = forwardrtb[lastmon][lastidx][i];
      outputrix[dim][i] = forwardrix[lastmon][lastidx][i];
      outputrmi[dim][i] = forwardrmi[lastmon][lastidx][i];
      outputrrg[dim][i] = forwardrrg[lastmon][lastidx][i];
      outputrpk[dim][i] = forwardrpk[lastmon][lastidx][i];
      outputltb[dim][i] = forwardltb[lastmon][lastidx][i];
      outputlix[dim][i] = forwardlix[lastmon][lastidx][i];
      outputlmi[dim][i] = forwardlmi[lastmon][lastidx][i];
      outputlrg[dim][i] = forwardlrg[lastmon][lastidx][i];
      outputlpk[dim][i] = forwardlpk[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputrtb[0][i] = 0.0;
         outputrix[0][i] = 0.0;
         outputrmi[0][i] = 0.0;
         outputrrg[0][i] = 0.0;
         outputrpk[0][i] = 0.0;
         outputltb[0][i] = 0.0;
         outputlix[0][i] = 0.0;
         outputlmi[0][i] = 0.0;
         outputlrg[0][i] = 0.0;
         outputlpk[0][i] = 0.0;
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
         outputrtb[0][i] = backwardrtb[ix0][ix2][i];
         outputrix[0][i] = backwardrix[ix0][ix2][i];
         outputrmi[0][i] = backwardrmi[ix0][ix2][i];
         outputrrg[0][i] = backwardrrg[ix0][ix2][i];
         outputrpk[0][i] = backwardrpk[ix0][ix2][i];
         outputltb[0][i] = backwardltb[ix0][ix2][i];
         outputlix[0][i] = backwardlix[ix0][ix2][i];
         outputlmi[0][i] = backwardlmi[ix0][ix2][i];
         outputlrg[0][i] = backwardlrg[ix0][ix2][i];
         outputlpk[0][i] = backwardlpk[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputrtb[k][i] = 0.0;
            outputrix[k][i] = 0.0;
            outputrmi[k][i] = 0.0;
            outputrrg[k][i] = 0.0;
            outputrpk[k][i] = 0.0;
            outputltb[k][i] = 0.0;
            outputlix[k][i] = 0.0;
            outputlmi[k][i] = 0.0;
            outputlrg[k][i] = 0.0;
            outputlpk[k][i] = 0.0;
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
               outputrtb[k][i] = backwardrtb[ix0][ix2][i];
               outputrix[k][i] = backwardrix[ix0][ix2][i];
               outputrmi[k][i] = backwardrmi[ix0][ix2][i];
               outputrrg[k][i] = backwardrrg[ix0][ix2][i];
               outputrpk[k][i] = backwardrpk[ix0][ix2][i];
               outputltb[k][i] = backwardltb[ix0][ix2][i];
               outputlix[k][i] = backwardlix[ix0][ix2][i];
               outputlmi[k][i] = backwardlmi[ix0][ix2][i];
               outputlrg[k][i] = backwardlrg[ix0][ix2][i];
               outputlpk[k][i] = backwardlpk[ix0][ix2][i];
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
               outputrtb[k][i] = forwardrtb[ix0][ix2][i];
               outputrix[k][i] = forwardrix[ix0][ix2][i];
               outputrmi[k][i] = forwardrmi[ix0][ix2][i];
               outputrrg[k][i] = forwardrrg[ix0][ix2][i];
               outputrpk[k][i] = forwardrpk[ix0][ix2][i];
               outputltb[k][i] = forwardltb[ix0][ix2][i];
               outputlix[k][i] = forwardlix[ix0][ix2][i];
               outputlmi[k][i] = forwardlmi[ix0][ix2][i];
               outputlrg[k][i] = forwardlrg[ix0][ix2][i];
               outputlpk[k][i] = forwardlpk[ix0][ix2][i];
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
               outputrtb[k][i] = crossrtb[ix0][ix2][i];
               outputrix[k][i] = crossrix[ix0][ix2][i];
               outputrmi[k][i] = crossrmi[ix0][ix2][i];
               outputrrg[k][i] = crossrrg[ix0][ix2][i];
               outputrpk[k][i] = crossrpk[ix0][ix2][i];
               outputltb[k][i] = crossltb[ix0][ix2][i];
               outputlix[k][i] = crosslix[ix0][ix2][i];
               outputlmi[k][i] = crosslmi[ix0][ix2][i];
               outputlrg[k][i] = crosslrg[ix0][ix2][i];
               outputlpk[k][i] = crosslpk[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl10_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk, 
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk, 
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose )
{
   double ***forwardrtb = new double**[nbr];
   double ***forwardrix = new double**[nbr];
   double ***forwardrmi = new double**[nbr];
   double ***forwardrrg = new double**[nbr];
   double ***forwardrpk = new double**[nbr];
   double ***forwardltb = new double**[nbr];
   double ***forwardlix = new double**[nbr];
   double ***forwardlmi = new double**[nbr];
   double ***forwardlrg = new double**[nbr];
   double ***forwardlpk = new double**[nbr];
   double ***backwardrtb = new double**[nbr];
   double ***backwardrix = new double**[nbr];
   double ***backwardrmi = new double**[nbr];
   double ***backwardrrg = new double**[nbr];
   double ***backwardrpk = new double**[nbr];
   double ***backwardltb = new double**[nbr];
   double ***backwardlix = new double**[nbr];
   double ***backwardlmi = new double**[nbr];
   double ***backwardlrg = new double**[nbr];
   double ***backwardlpk = new double**[nbr];
   double ***crossrtb = new double**[nbr];
   double ***crossrix = new double**[nbr];
   double ***crossrmi = new double**[nbr];
   double ***crossrrg = new double**[nbr];
   double ***crossrpk = new double**[nbr];
   double ***crossltb = new double**[nbr];
   double ***crosslix = new double**[nbr];
   double ***crosslmi = new double**[nbr];
   double ***crosslrg = new double**[nbr];
   double ***crosslpk = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardrtb[k] = new double*[nvrk];
      forwardrix[k] = new double*[nvrk];
      forwardrmi[k] = new double*[nvrk];
      forwardrrg[k] = new double*[nvrk];
      forwardrpk[k] = new double*[nvrk];
      forwardltb[k] = new double*[nvrk];
      forwardlix[k] = new double*[nvrk];
      forwardlmi[k] = new double*[nvrk];
      forwardlrg[k] = new double*[nvrk];
      forwardlpk[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardrtb[k][i] = new double[deg+1];
         forwardrix[k][i] = new double[deg+1];
         forwardrmi[k][i] = new double[deg+1];
         forwardrrg[k][i] = new double[deg+1];
         forwardrpk[k][i] = new double[deg+1];
         forwardltb[k][i] = new double[deg+1];
         forwardlix[k][i] = new double[deg+1];
         forwardlmi[k][i] = new double[deg+1];
         forwardlrg[k][i] = new double[deg+1];
         forwardlpk[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardrtb[k] = new double*[nvrk-1];
         backwardrix[k] = new double*[nvrk-1];
         backwardrmi[k] = new double*[nvrk-1];
         backwardrrg[k] = new double*[nvrk-1];
         backwardrpk[k] = new double*[nvrk-1];
         backwardltb[k] = new double*[nvrk-1];
         backwardlix[k] = new double*[nvrk-1];
         backwardlmi[k] = new double*[nvrk-1];
         backwardlrg[k] = new double*[nvrk-1];
         backwardlpk[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardrtb[k][i] = new double[deg+1];
            backwardrix[k][i] = new double[deg+1];
            backwardrmi[k][i] = new double[deg+1];
            backwardrrg[k][i] = new double[deg+1];
            backwardrpk[k][i] = new double[deg+1];
            backwardltb[k][i] = new double[deg+1];
            backwardlix[k][i] = new double[deg+1];
            backwardlmi[k][i] = new double[deg+1];
            backwardlrg[k][i] = new double[deg+1];
            backwardlpk[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crossrtb[k] = new double*[nvrk-2];
         crossrix[k] = new double*[nvrk-2];
         crossrmi[k] = new double*[nvrk-2];
         crossrrg[k] = new double*[nvrk-2];
         crossrpk[k] = new double*[nvrk-2];
         crossltb[k] = new double*[nvrk-2];
         crosslix[k] = new double*[nvrk-2];
         crosslmi[k] = new double*[nvrk-2];
         crosslrg[k] = new double*[nvrk-2];
         crosslpk[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crossrtb[k][i] = new double[deg+1];
            crossrix[k][i] = new double[deg+1];
            crossrmi[k][i] = new double[deg+1];
            crossrrg[k][i] = new double[deg+1];
            crossrpk[k][i] = new double[deg+1];
            crossltb[k][i] = new double[deg+1];
            crosslix[k][i] = new double[deg+1];
            crosslmi[k][i] = new double[deg+1];
            crosslrg[k][i] = new double[deg+1];
            crosslpk[k][i] = new double[deg+1];
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

         CPU_dbl10_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffrtb[monidx],cffrix[monidx],cffrmi[monidx],
             cffrrg[monidx],cffrpk[monidx],
             cffltb[monidx],cfflix[monidx],cfflmi[monidx],
             cfflrg[monidx],cfflpk[monidx],
             inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
             inputltb,inputlix,inputlmi,inputlrg,inputlpk,
              forwardrtb[monidx], forwardrix[monidx], forwardrmi[monidx],
              forwardrrg[monidx], forwardrpk[monidx],
              forwardltb[monidx], forwardlix[monidx], forwardlmi[monidx],
              forwardlrg[monidx], forwardlpk[monidx],
             backwardrtb[monidx],backwardrix[monidx],backwardrmi[monidx],
             backwardrrg[monidx],backwardrpk[monidx],
             backwardltb[monidx],backwardlix[monidx],backwardlmi[monidx],
             backwardlrg[monidx],backwardlpk[monidx],
                crossrtb[monidx],   crossrix[monidx],   crossrmi[monidx],
                crossrrg[monidx],   crossrpk[monidx],
                crossltb[monidx],   crosslix[monidx],   crosslmi[monidx],
                crosslrg[monidx],   crosslpk[monidx],job,verbose);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl10_poly_addjobs
      (dim,nbr,deg,nvr,idx,
            cstrtb,     cstrix,     cstrmi,     cstrrg,     cstrpk,
            cstltb,     cstlix,     cstlmi,     cstlrg,     cstlpk,
            cffrtb,     cffrix,     cffrmi,     cffrrg,     cffrpk,
            cffltb,     cfflix,     cfflmi,     cfflrg,     cfflpk,
          inputrtb,   inputrix,   inputrmi,   inputrrg,   inputrpk,
          inputltb,   inputlix,   inputlmi,   inputlrg,   inputlpk,
         outputrtb,  outputrix,  outputrmi,  outputrrg,  outputrpk,
         outputltb,  outputlix,  outputlmi,  outputlrg,  outputlpk,
        forwardrtb, forwardrix, forwardrmi, forwardrrg, forwardrpk,
        forwardltb, forwardlix, forwardlmi, forwardlrg, forwardlpk,
       backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
       backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,
          crossrtb,   crossrix,   crossrmi,   crossrrg,   crossrpk,
          crossltb,   crosslix,   crosslmi,   crosslrg,   crosslpk,
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
         free(forwardrtb[k][i]);
         free(forwardrix[k][i]);
         free(forwardrmi[k][i]);
         free(forwardrrg[k][i]);
         free(forwardrpk[k][i]);
         free(forwardltb[k][i]);
         free(forwardlix[k][i]);
         free(forwardlmi[k][i]);
         free(forwardlrg[k][i]);
         free(forwardlpk[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardrtb[k][i]);
                      free(backwardrix[k][i]);
                      free(backwardrmi[k][i]);
                      free(backwardrrg[k][i]);
                      free(backwardrpk[k][i]);
                      free(backwardltb[k][i]);
                      free(backwardlix[k][i]);
                      free(backwardlmi[k][i]);
                      free(backwardlrg[k][i]);
                      free(backwardlpk[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crossrtb[k][i]);
                      free(crossrix[k][i]);
                      free(crossrmi[k][i]);
                      free(crossrrg[k][i]);
                      free(crossrpk[k][i]);
                      free(crossltb[k][i]);
                      free(crosslix[k][i]);
                      free(crosslmi[k][i]);
                      free(crosslrg[k][i]);
                      free(crosslpk[k][i]);
                   }
   }
   free(forwardrtb); free(backwardrtb); free(crossrtb);
   free(forwardrix); free(backwardrix); free(crossrix);
   free(forwardrmi); free(backwardrmi); free(crossrmi);
   free(forwardrrg); free(backwardrrg); free(crossrrg);
   free(forwardrpk); free(backwardrpk); free(crossrpk);
   free(forwardltb); free(backwardltb); free(crossltb);
   free(forwardlix); free(backwardlix); free(crosslix);
   free(forwardlmi); free(backwardlmi); free(crosslmi);
   free(forwardlrg); free(backwardlrg); free(crosslrg);
   free(forwardlpk); free(backwardlpk); free(crosslpk);
}
