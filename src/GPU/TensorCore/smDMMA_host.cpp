/* The file smDMMA_host.cpp contains the definitions of the functions with
 * prototypes in smDMMA_host.h. */

#include "smDMMA_dims.h"
#include "smDMMA_host.h"

void init_host_matrices ( double *a, double *b, double *c, int nbrange )
{
   for(int i=0; i<K_GLOBAL; i++) // #rows: M_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<i; j++) a[i*K_GLOBAL+j] = 0.0;
      for(int j=i; j<K_GLOBAL; j++) a[i*K_GLOBAL+j] = 1.0;
   }
   for(int i=K_GLOBAL; i<M_GLOBAL; i++) // row major with M_GLOBAL rows
      for(int j=0; j<K_GLOBAL; j++) a[i*K_GLOBAL+j] = 0.0;

   for(int i=0; i<K_GLOBAL; i++) // #columns: N_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<=i; j++) b[i*K_GLOBAL+j] = 1.0;
      for(int j=i+1; j<K_GLOBAL; j++) b[i*K_GLOBAL+j] = 0.0;
   }
   for(int i=K_GLOBAL; i<N_GLOBAL; i++) // #columns: N_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<K_GLOBAL; j++) 
      {
         if(j % nbrange == 0)
            b[i*K_GLOBAL+j] = 1.0;
         else
            b[i*K_GLOBAL+j] = 2.0*b[i*K_GLOBAL+j-1];
      }
   }

   for(int t=0; t<M_GLOBAL*N_GLOBAL; t++) c[t] = 0.0;
}

void matMultiplyOnHost
 ( double *A, double *B, double *C, float alpha, float beta,
   int numARows, int numAColumns, int numBRows, int numBColumns,
   int numCRows, int numCColumns )
{
   for(int i = 0; i < numCRows; i++)
   {
      for(int j = 0; j < numCColumns; j++)
      {
         double temp = 0.0;

         for(int k = 0; k < numAColumns; k++)
         {
            // B matrix is column major. A matrix is row major.
            temp += A[i*numAColumns+k] * B[j*numBRows+k];
         }
         C[i*numCColumns+j] = temp*alpha + beta*C[i*numCColumns+j];
        }
    }
}
