/* Tests the setup of random double double matrices as single indexed
 * matrices for input to tensor core matrix multiplication. */

#include <iostream>
#include <cstdlib>
#include "smDMMA_dims.h"
#include "smDMMA_host.h"

using namespace std;

int main ( void )
{
   double *A = NULL;
   double *B = NULL;
   double *C = NULL;

   A = (double*) malloc(sizeof(double) * M_GLOBAL * K_GLOBAL);
   B = (double*) malloc(sizeof(double) * K_GLOBAL * N_GLOBAL);
   C = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);

   random_double_double_matrices
      (A,B,C,M_GLOBAL,K_GLOBAL,K_GLOBAL,N_GLOBAL,M_GLOBAL,N_GLOBAL,1);

   return 0;
}
