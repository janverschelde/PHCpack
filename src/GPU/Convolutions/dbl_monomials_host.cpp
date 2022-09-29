/* The file dbl_monomials_host.cpp defines functions specified
 * in dbl_monomials_host.h. */

/* The algorithm to compute forward, backward, and cross products
 * (denoted respectively by arrays f, b, and c)
 * for a monomial cff*x[0]*x[1]* .. *x[n-1] goes as follows:
 *
 * f[0] := cff*x[0]
 * for i from 1 to n-1 do f[i] := f[i-1]*x[i]
 * if n > 2 then
 *    b[0] := x[n-1]*x[n-2]
 *    for i from 1 to n-3 do b[i] := b[i-1]*x[n-2-i]
 *    b[n-3] := b[n-3]*cff
 *    if n = 3 then
 *       c[0] = f[0]*x[2]
 *    else
 *       for i from 0 to n-4 do c[i] := f[i]*b[n-4-i]
 *       c[n-3] := f[n-3]*x[n-1]
 *
 * Compared to the evaluation and differentiation of a product of variables,
 * (without coefficient cff), two extra multiplications must be done,
 * but this is better than n+1 multiplications with cff afterwards. */

#include <cstdlib>
#include <iostream>
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"

using namespace std;

void CPU_dbl_speel
 ( int nvr, int deg, int *idx, double *cff, double **input,
   double **forward, double **backward, double **cross,
   bool verbose, int monidx )
{
   int ix1 = idx[0];
   int ix2;

   CPU_dbl_product(deg,cff,input[ix1],forward[0]); // f[0] = cff*x[0] 
   if(verbose) cout << "monomial " << monidx << " : ";
   if(verbose) cout << "cff * input[" << ix1 << "] to f[0]" << endl;

   for(int i=1; i<nvr; i++)
   {                                               // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_dbl_product(deg,forward[i-1],input[ix2],forward[i]);
      if(verbose) cout << "monomial " << monidx << " : ";
      if(verbose) cout << "f[" << i-1 << "] * "
                       << "input[" << ix2 << "] to f[" << i << "]" << endl;
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1];                            // b[0] = x[n-1]*x[n-2]
      ix2 = idx[nvr-2];
      CPU_dbl_product(deg,input[ix1],input[ix2],backward[0]);
      if(verbose) cout << "monomial " << monidx << " : ";
      if(verbose) cout << "input[" << ix1 << "] * "
                       << "input[" << ix2 << "] to b[0]" << endl;
      for(int i=1; i<nvr-2; i++)
      {                                            // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];
         CPU_dbl_product(deg,backward[i-1],input[ix2],backward[i]);
         if(verbose) cout << "monomial " << monidx << " : ";
         if(verbose) cout << "b[" << i-1 << "] * "
                          << "input[" << ix2 << "] to b[" << i << "]" << endl;
      }
                                                   // b[n-3] = cff*b[n-3]
      CPU_dbl_product(deg,backward[nvr-3],cff,cross[0]);
      if(verbose) cout << "monomial " << monidx << " : ";
      if(verbose) cout << "b[" << nvr-3 << "] * cff to c[0]" << endl;
      // cross[0] is work space, cannot write into backward[nvr-3]
      for(int i=0; i<=deg; i++) backward[nvr-3][i] = cross[0][i];

      if(nvr == 3)
      {                                            // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_dbl_product(deg,forward[0],input[ix2],cross[0]);
         if(verbose) cout << "monomial " << monidx << " : ";
         if(verbose) cout << "f[0] * input[" << ix2 << "] to c[0]" << endl;
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                         // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_dbl_product(deg,forward[i],backward[ix2],cross[i]);
            if(verbose) cout << "monomial " << monidx << " : ";
            if(verbose) cout << "f[" << i << "] * b[" << ix2
                             << "] to c[" << i << "]" << endl;
         }
                                                   // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_dbl_product(deg,forward[nvr-3],input[ix2],cross[nvr-3]);
         if(verbose) cout << "monomial " << monidx << " : ";
         if(verbose) cout << "f[" << nvr-3 << "] * input[" << ix2
                          << "] to c[" << nvr-3 << "]" << endl;
      }
   }
}

void CPU_cmplx_speel
 ( int nvr, int deg, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim, double **forwardre,
   double **forwardim, double **backwardre, double **backwardim,
   double **crossre, double **crossim )
{
   int ix1 = idx[0];
   int ix2;
                                                           // f[0] = cff*x[0]
   CPU_cmplx_product(deg,cffre,cffim,inputre[ix1],inputim[ix1],
                                     forwardre[0],forwardim[0]);
   for(int i=1; i<nvr; i++)
   {                                                    // f[i] = f[i-1]*x[i]
      ix2 = idx[i];
      CPU_cmplx_product(deg,forwardre[i-1],forwardim[i-1],inputre[ix2],
                        inputim[ix2],forwardre[i],forwardim[i]);
   }
   if(nvr > 2)
   {                                                  // b[0] = x[n-1]*x[n-2]
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];
      CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],inputre[ix2],
                            inputim[ix2],backwardre[0],backwardim[0]);
      for(int i=1; i<nvr-2; i++)
      {                                             // b[i] = b[i-1]*x[x-2-i]
         ix2 = idx[nvr-2-i];
         CPU_cmplx_product(deg,backwardre[i-1],backwardim[i-1],inputre[ix2],
                               inputim[ix2],backwardre[i],backwardim[i]);
      }
                                                       // b[n-3] = b[n-3]*cff
      CPU_cmplx_product(deg,backwardre[nvr-3],backwardim[nvr-3],cffre,cffim,
                            crossre[0],crossim[0]);    // cross is work space
      for(int i=0; i<=deg; i++)
      {
         backwardre[nvr-3][i] = crossre[0][i];
         backwardim[nvr-3][i] = crossim[0][i];
      }
      if(nvr == 3)
      {                                                   // c[0] = f[0]*x[2]
         ix2 = idx[2];
         CPU_cmplx_product(deg,forwardre[0],forwardim[0],inputre[ix2],
                               inputim[ix2],crossre[0],crossim[0]);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                            // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;
            CPU_cmplx_product(deg,forwardre[i],forwardim[i],backwardre[ix2],
                              backwardim[ix2],crossre[i],crossim[i]);
         }
                                                    // c[n-3] = f[n-3]*x[n-1]
         ix2 = idx[nvr-1];
         CPU_cmplx_product(deg,forwardre[nvr-3],forwardim[nvr-3],inputre[ix2],
                           inputim[ix2],crossre[nvr-3],crossim[nvr-3]);
      }
   }
}

void CPU_dbl_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cff,
   double **input, double **output )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_dbl_product(deg,input[ix],cff,output[dim]);

      for(int i=0; i<=deg; i++) output[ix][i] = cff[i];
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_dbl_product(deg,input[ix1],input[ix2],output[dim]);
      // CPU_dbl_product(deg,output[dim],cff,output[dim]); // wrong!
      double *acc = new double[deg+1];
      CPU_dbl_product(deg,output[dim],cff,acc);
      for(int i=0; i<=deg; i++) output[dim][i] = acc[i];
      free(acc);

      CPU_dbl_product(deg,cff,input[ix1],output[ix2]);
      CPU_dbl_product(deg,cff,input[ix2],output[ix1]);
   }
   else
   {
      double **forward = new double*[nvr];
      double **backward = new double*[nvr-2];
      double **cross = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forward[i] = new double[deg+1];
         backward[i] = new double[deg+1];
         cross[i] = new double[deg+1];
      }
      forward[nvr-2] = new double[deg+1];
      forward[nvr-1] = new double[deg+1];

      CPU_dbl_speel(nvr,deg,idx,cff,input,forward,backward,cross);

      // assign the value of the monomial

      for(int i=0; i<deg+1; i++) output[dim][i] = forward[nvr-1][i];

      if(nvr > 2)
      {
         int ix = idx[nvr-1];       // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++) output[ix][i] = forward[nvr-2][i];

         ix = idx[0];               // derivative with respect to x[0]

         for(int i=0; i<deg+1; i++) output[ix][i] = backward[nvr-3][i];

         for(int k=1; k<nvr-1; k++)
         {                          // derivative with respect to x[k]
            ix = idx[k];
            for(int i=0; i<deg+1; i++) output[ix][i] = cross[k-1][i];
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forward[i]); free(backward[i]); free(cross[i]);
      }
      free(forward[nvr-2]); free(forward[nvr-1]);
      free(forward); free(backward); free(cross);
   }
}

void CPU_cmplx_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim, double **outputre, double **outputim )
{
   if(nvr == 1)
   {
      int ix = idx[0];

      CPU_cmplx_product(deg,inputre[ix],inputim[ix],cffre,cffim,
                            outputre[dim],outputim[dim]);

      for(int i=0; i<=deg; i++) 
      {
         outputre[ix][i] = cffre[i]; outputim[ix][i] = cffim[i];
      }
   }
   else if(nvr == 2)
   {
      int ix1 = idx[0]; int ix2 = idx[1];

      CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],
                            inputre[ix2],inputim[ix2],
                            outputre[dim],outputim[dim]);
      // CPU_cmplx_product(deg,outputre[dim],outputim[dim],cffre,cffim,
      //                       outputre[dim],outputim[dim]);
      double *accre = new double[deg+1];
      double *accim = new double[deg+1];
      CPU_cmplx_product
         (deg,outputre[dim],outputim[dim],cffre,cffim,accre,accim);
      for(int i=0; i<=deg; i++)
      {
         outputre[dim][i] = accre[i];
         outputim[dim][i] = accim[i];
      }
      free(accre);
      free(accim);

      CPU_cmplx_product(deg,cffre,cffim,inputre[ix1],inputim[ix1],
                                        outputre[ix2],outputim[ix2]);
      CPU_cmplx_product(deg,cffre,cffim,inputre[ix2],inputim[ix2],
                                        outputre[ix1],outputim[ix1]);
   }
   else
   {
      double **forwardre = new double*[nvr];
      double **forwardim = new double*[nvr];
      double **backwardre = new double*[nvr-2];
      double **backwardim = new double*[nvr-2];
      double **crossre = new double*[nvr-2];
      double **crossim = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forwardre[i] = new double[deg+1]; forwardim[i] = new double[deg+1];
         backwardre[i] = new double[deg+1]; backwardim[i] = new double[deg+1];
         crossre[i] = new double[deg+1]; crossim[i] = new double[deg+1];
      }
      forwardre[nvr-2] = new double[deg+1];
      forwardim[nvr-2] = new double[deg+1];
      forwardre[nvr-1] = new double[deg+1];
      forwardim[nvr-1] = new double[deg+1];

      CPU_cmplx_speel(nvr,deg,idx,cffre,cffim,inputre,inputim,forwardre,
                      forwardim,backwardre,backwardim,crossre,crossim);

      for(int i=0; i<deg+1; i++)          // assign value of the monomial
      {
         outputre[dim][i] = forwardre[nvr-1][i];
         outputim[dim][i] = forwardim[nvr-1][i];
      }

      if(nvr > 2)
      {
         int ix = idx[nvr-1];        // derivative with respect to x[n-1]

         for(int i=0; i<deg+1; i++)
         {
            outputre[ix][i] = forwardre[nvr-2][i];
            outputim[ix][i] = forwardim[nvr-2][i];
         }

         ix = idx[0];                  // derivative with respect to x[0]
         for(int i=0; i<deg+1; i++)
         {
            outputre[ix][i] = backwardre[nvr-3][i];
            outputim[ix][i] = backwardim[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];               // derivative with respect to x[k]
            for(int i=0; i<deg+1; i++)
            {
               outputre[ix][i] = crossre[k-1][i];
               outputim[ix][i] = crossim[k-1][i];
            }
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forwardre[i]); free(backwardre[i]); free(crossre[i]);
         free(forwardim[i]); free(backwardim[i]); free(crossim[i]);
      }
      free(forwardre[nvr-2]); free(forwardre[nvr-1]); free(forwardre);
      free(forwardim[nvr-2]); free(forwardim[nvr-1]); free(forwardim);
      free(backwardre); free(crossre);
      free(backwardim); free(crossim);
   }
}
