// The file dbl_comfacmon_kernels.cu defines the functions with prototypes in
// the file dbl_comfacmon_kernels.h.

#include <iostream>
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"

using namespace std;

void GPU_dbl_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac,
   int **expfac, double **cff, double *acc, double **input,
   double ***output, int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_dbl_product(deg,input[idxvar],cff[i],acc);
                  for(int L=0; L<=deg; L++) cff[i][L] = acc[L];
               }
            }
         }
      }
   }
   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << cff[i][j] << endl;
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl_evaldiff(dim,nvr[i],deg,idx[i],cff[i],input,output[i]);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << output[i][dim][j] << endl;
      }
   }

   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
                  output[i][idxvar][k] = factor*output[i][idxvar][k];
            }
         }
      }
   }
}
