/* The file dbl_monomials_host.cpp defines functions specified
 * in dbl_monomials_host.h. */

#include <cstdlib>
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"

void CPU_dbl_speel
 ( int nvr, int deg, int *idx, double **input,
   double **forward, double **backward, double **cross )
{
   int ix1 = idx[0];
   int ix2 = idx[1];

   CPU_dbl_product(deg,input[ix1],input[ix2],forward[0]);
   for(int i=2; i<nvr; i++)
   {
      ix2 = idx[i];
      CPU_dbl_product(deg,forward[i-2],input[ix2],forward[i-1]);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1];
      ix2 = idx[nvr-2];
      CPU_dbl_product(deg,input[ix1],input[ix2],backward[0]);
      for(int i=1; i<nvr-2; i++)
      {
         ix2 = idx[nvr-2-i];
         CPU_dbl_product(deg,backward[i-1],input[ix2],backward[i]);
      }
      if(nvr == 3)
      {
         ix1 = idx[0];
         ix2 = idx[2];
         CPU_dbl_product(deg,input[ix1],input[ix2],cross[0]);
      }
      else
      {
         ix1 = idx[0];
         ix2 = nvr-3;
         CPU_dbl_product(deg,input[ix1],backward[ix2],cross[0]);
         for(int i=1; i<nvr-3; i++)
         {
            ix2 = nvr-3-i;
            CPU_dbl_product(deg,forward[i-1],backward[ix2],cross[i]);
         }
         ix2 = idx[nvr-1];
         CPU_dbl_product(deg,forward[nvr-4],input[ix2],cross[nvr-3]);
      }
   }
}

void CPU_cmplx_speel
 ( int nvr, int deg, int *idx, double **inputre, double **inputim,
   double **forwardre, double **forwardim, double **backwardre,
   double **backwardim, double **crossre, double **crossim )
{
   int ix1 = idx[0];
   int ix2 = idx[1];

   CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],
                         inputre[ix2],inputim[ix2],
                         forwardre[0],forwardim[0]);
   for(int i=2; i<nvr; i++)
   {
      ix2 = idx[i];
      CPU_cmplx_product(deg,forwardre[i-2],forwardim[i-2],
                            inputre[ix2],inputim[ix2],
                            forwardre[i-1],forwardim[i-1]);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];

      CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],
                            inputre[ix2],inputim[ix2],
                            backwardre[0],backwardim[0]);
      for(int i=1; i<nvr-2; i++)
      {
         ix2 = idx[nvr-2-i];

         CPU_cmplx_product(deg,backwardre[i-1],backwardim[i-1],
                               inputre[ix2],inputim[ix2],
                               backwardre[i],backwardim[i]);
      }
      if(nvr == 3)
      {
         ix1 = idx[0]; ix2 = idx[2];

         CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],
                               inputre[ix2],inputim[ix2],
                               crossre[0],crossim[0]);
      }
      else
      {
         ix1 = idx[0]; ix2 = nvr-3;

         CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],
                               backwardre[ix2],backwardim[ix2],
                               crossre[0],crossim[0]);
         for(int i=1; i<nvr-3; i++)
         {
            ix2 = nvr-3-i;

            CPU_cmplx_product(deg,forwardre[i-1],forwardim[i-1],
                                  backwardre[ix2],backwardim[ix2],
                                  crossre[i],crossim[i]);
         }
         ix2 = idx[nvr-1];

         CPU_cmplx_product(deg,forwardre[nvr-4],forwardim[nvr-4],
                               inputre[ix2],inputim[ix2],
                               crossre[nvr-3],crossim[nvr-3]);
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
      CPU_dbl_product(deg,output[dim],cff,output[dim]);

      CPU_dbl_product(deg,cff,input[ix1],output[ix2]);
      CPU_dbl_product(deg,cff,input[ix2],output[ix1]);
   }
   else
   {
      double **forward = new double*[nvr-1];
      double **backward = new double*[nvr-2];
      double **cross = new double*[nvr-2];

      for(int i=0; i<nvr-2; i++)
      {
         forward[i] = new double[deg+1];
         backward[i] = new double[deg+1];
         cross[i] = new double[deg+1];
      }
      forward[nvr-2] = new double[deg+1];

      CPU_dbl_speel(nvr,deg,idx,input,forward,backward,cross);
      for(int i=0; i<deg+1; i++) output[dim][i] = forward[nvr-2][i];

      // CPU_dbl_product(deg,cff,forward[nvr-2],output[dim]);

      if(nvr > 2)
      {
         int ix = idx[nvr-1];

         // CPU_dbl_product(deg,cff,forward[nvr-3],output[ix]);
         for(int i=0; i<deg+1; i++) output[ix][i] = forward[nvr-3][i];
         ix = idx[0];
         // CPU_dbl_product(deg,cff,backward[nvr-3],output[ix]);
         for(int i=0; i<deg+1; i++) output[ix][i] = backward[nvr-3][i];

         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];
            // CPU_dbl_product(deg,cff,cross[k-1],output[ix]);
            for(int i=0; i<deg+1; i++) output[ix][i] = cross[k-1][i];
         }
      }
      for(int i=0; i<nvr-2; i++)
      {
         free(forward[i]); free(backward[i]); free(cross[i]);
      }
      free(forward[nvr-2]); free(forward); free(backward); free(cross);
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
      CPU_cmplx_product(deg,outputre[dim],outputim[dim],cffre,cffim,
                            outputre[dim],outputim[dim]);

      CPU_cmplx_product(deg,cffre,cffim,inputre[ix1],inputim[ix1],
                                        outputre[ix2],outputim[ix2]);
      CPU_cmplx_product(deg,cffre,cffim,inputre[ix2],inputim[ix1],
                                        outputre[ix1],outputim[ix1]);
   }
   else
   {
      double **forwardre = new double*[nvr-1];
      double **forwardim = new double*[nvr-1];
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

      CPU_cmplx_speel(nvr,deg,idx,inputre,inputim,forwardre,forwardim,
                      backwardre,backwardim,crossre,crossim);

      for(int i=0; i<deg+1; i++)
      {
         outputre[dim][i] = forwardre[nvr-2][i];
         outputim[dim][i] = forwardim[nvr-2][i];
      }

      // CPU_dbl_product(deg,cff,forward[nvr-2],output[dim]);

      if(nvr > 2)
      {
         int ix = idx[nvr-1];

         // CPU_dbl_product(deg,cff,forward[nvr-3],output[ix]);
         for(int i=0; i<deg+1; i++) 
         {
            outputre[ix][i] = forwardre[nvr-3][i];
            outputim[ix][i] = forwardim[nvr-3][i];
         }

         ix = idx[0];
         // CPU_dbl_product(deg,cff,backward[nvr-3],output[ix]);
         for(int i=0; i<deg+1; i++)
         {
            outputre[ix][i] = backwardre[nvr-3][i];
            outputim[ix][i] = backwardim[nvr-3][i];
         }
         for(int k=1; k<nvr-1; k++)
         {
            ix = idx[k];
            // CPU_dbl_product(deg,cff,cross[k-1],output[ix]);
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
      free(forwardre[nvr-2]); free(forwardre);
      free(forwardim[nvr-2]); free(forwardim);
      free(backwardre); free(crossre);
      free(backwardim); free(crossim);
   }
}
