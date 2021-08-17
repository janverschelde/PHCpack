/* The file dbl_tabs_flopcounts.cpp defines the functions with prototypes in
 * the file dbl_tabs_flopcounts.h. */

#include "dbl_tabs_flopcounts.h"

void flopcount_dbl_invert_tiles
 ( int nbt, int szt, long long int *add, long long int *mul,
   long long int *div )
{
   const int nbthreads = nbt*szt;
   const int dim = szt;

   *add += nbthreads*(dim-1)*(dim-2)/2;
   *mul += nbthreads*(dim-1)*(dim-2)/2;
   *div += nbthreads*(dim-1);
}

void flopcount_cmplx_invert_tiles
 ( int nbt, int szt, long long int *add, long long int *mul,
   long long int *div )
{
   const int nbthreads = nbt*szt;
   const int dim = szt;

   *add += 4*nbthreads*(dim-1)*(dim-2)/2;
   *mul += 4*nbthreads*(dim-1)*(dim-2)/2;
   *div += 2*nbthreads*(dim-1);
}

void flopcount_dbl_multiply_inverse
 ( int szt, long long int *add, long long int *mul )
{
   *add += szt*szt;
   *mul += szt*szt;
}

void flopcount_cmplx_multiply_inverse
 ( int szt, long long int *add, long long int *mul )
{
   *add += szt*szt;
   *mul += szt*szt;
}

void flopcount_dbl_back_substitute
 ( int nblocks, int szt, long long int *add, long long int *mul )
{
   const int nbthreads = nblocks*szt;

   *add += nbthreads*(szt + 1); // dim equals szt, one extract subtraction
   *mul += nbthreads*szt;
}

void flopcount_cmplx_back_substitute
 ( int nblocks, int szt, long long int *add, long long int *mul )
{
   const int nbthreads = nblocks*szt;

   *add += 4*nbthreads*szt + 2*nbthreads;
   *mul += 4*nbthreads*szt;
}
