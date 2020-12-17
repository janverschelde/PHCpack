/* The file dbl_polynomials_host.cpp defines functions specified
 * in dbl_polynomials_host.h. */

#include <cstdlib>
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"
#include "dbl_polynomials_host.h"

void CPU_dbl_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cff, double **input, double **output,
   double **forward, double **backward, double **cross )
{
   int ix1,ix2;

   for(int i=0; i<nbr-1; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl_product(deg,input[ix1],cff[i],forward[0]);
         for(int j=0; j<=deg; j++)
         {
            output[dim][j] += forward[0][j];
            output[ix1][j] += cff[i][j];
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl_product(deg,input[ix1],input[ix2],forward[0]);
         CPU_dbl_product(deg,forward[0],cff[i],forward[0]);

         CPU_dbl_product(deg,cff[i],input[ix1],cross[0]);
         CPU_dbl_product(deg,cff[i],input[ix2],backward[0]);

         for(int j=0; j<=deg; j++)
         {
            output[dim][j] += forward[0][j];
            output[ix1][j] += backward[i][j];
            output[ix2][j] += cross[i][j];
         }
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl_speel(nvr[i],deg,idx[i],cff[i],input,forward,backward,cross);

         ix1 = nvr[i]-1;
         for(int j=0; j<=deg; j++)     // update the value of the polynomial
            output[dim][j] += forward[ix1][j];

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) output[ix2][j] += forward[ix1][j];

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) output[ix2][j] += backward[ix1][j];

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) output[ix2][j] += cross[k-1][j];
         }
      }
   }
}

void CPU_dbl_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output )
{
   if(nbr == 1)
      CPU_dbl_evaldiff(dim,nvr[0],deg,idx[0],cff[0],input,output);
   else if(nbr > 1)
   {
      double **forward = new double*[dim];
      double **backward = new double*[dim-2];
      double **cross = new double*[dim-2];

      for(int i=0; i<dim-2; i++)
      {
         forward[i] = new double[deg+1];
         backward[i] = new double[deg+1];
         cross[i] = new double[deg+1];
      }
      forward[dim-2] = new double[deg+1];
      forward[dim-1] = new double[deg+1];

      for(int i=0; i<=deg; i++) output[dim][i] = cst[i];
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++) output[i][j] = 0.0;

      CPU_dbl_poly_speel
         (dim,nbr,deg,nvr,idx,cff,input,output,forward,backward,cross);

      for(int i=0; i<dim-2; i++)
      {
         free(forward[i]); free(backward[i]); free(cross[i]);
      }
      free(forward[dim-2]); free(forward[dim-1]);
      free(forward); free(backward); free(cross);
   }
}
