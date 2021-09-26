/* The file dbl_tabs_flopcounts.cpp defines the functions with prototypes in
 * the file dbl_tabs_flopcounts.h. */

#include <iostream>
using namespace std;

#include <climits>
#include "dbl_tabs_flopcounts.h"

void update_counters
 ( long long int *cnt, double *cntover, int inc )
{
   if(inc > 0)
   {
      long long int prevcnt = *cnt;
      if(*cnt + inc > 0) *cnt += inc;
      if(*cnt > prevcnt) 
         return;
      else
      {
         cout << "Overflow at cnt = " << *cnt
              << " and inc = " << inc << endl;
         // long long int overinc = *cnt + inc;
         // overinc = overinc + LONG_MAX;
         // *cnt = (overinc < 0) ? -overinc : overinc;
         *cnt = 0;
         *cntover += ((double) LONG_MAX);
      }
   }
}

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

void overflopcount_dbl_invert_tiles
 ( int nbt, int szt,
   long long int *add, double *addover,
   long long int *mul, double *mulover,
   long long int *div, double *divover )
{
   const int nbthreads = nbt*szt;
   const int dim = szt;
   int inc;

   // *add += nbthreads*(dim-1)*(dim-2)/2;
   inc = nbthreads*(dim-1)*(dim-2)/2;
   update_counters(add,addover,inc);
   // *mul += nbthreads*(dim-1)*(dim-2)/2;
   inc = nbthreads*(dim-1)*(dim-2)/2;
   update_counters(mul,mulover,inc);
   // *div += nbthreads*(dim-1);
   inc = nbthreads*(dim-1);
   update_counters(div,divover,inc);
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

void overflopcount_cmplx_invert_tiles
 ( int nbt, int szt,
   long long int *add, double *addover,
   long long int *mul, double *mulover,
   long long int *div, double *divover )
{
   const int nbthreads = nbt*szt;
   const int dim = szt;
   int inc;

   // *add += 4*nbthreads*(dim-1)*(dim-2)/2;
   inc = 4*nbthreads*(dim-1)*(dim-2)/2;
   update_counters(add,addover,inc);
   // *mul += 4*nbthreads*(dim-1)*(dim-2)/2;
   inc = 4*nbthreads*(dim-1)*(dim-2)/2;
   update_counters(mul,mulover,inc);
   // *div += 2*nbthreads*(dim-1);
   inc = 2*nbthreads*(dim-1);
   update_counters(div,divover,inc);
}

void flopcount_dbl_multiply_inverse
 ( int szt, long long int *add, long long int *mul )
{
   *add += szt*szt;
   *mul += szt*szt;
}

void overflopcount_dbl_multiply_inverse
 ( int szt, long long int *add, double *addover,
   long long int *mul, double *mulover )
{
   int inc;

   // *add += szt*szt;
   inc = szt*szt;
   update_counters(add,addover,inc);
   // *mul += szt*szt;
   inc = szt*szt;
   update_counters(mul,mulover,inc);
}

void flopcount_cmplx_multiply_inverse
 ( int szt, long long int *add, long long int *mul )
{
   *add += szt*szt;
   *mul += szt*szt;
}

void overflopcount_cmplx_multiply_inverse
 ( int szt, long long int *add, double *addover,
   long long int *mul, double *mulover )
{
   int inc;

   // *add += szt*szt;
   inc = szt*szt;
   update_counters(add,addover,inc);
   // *mul += szt*szt;
   inc = szt*szt;
   update_counters(mul,mulover,inc);
}

void flopcount_dbl_back_substitute
 ( int nblocks, int szt, long long int *add, long long int *mul )
{
   const int nbthreads = nblocks*szt;

   *add += nbthreads*(szt + 1); // dim equals szt, one extract subtraction
   *mul += nbthreads*szt;
}

void overflopcount_dbl_back_substitute
 ( int nblocks, int szt,
   long long int *add, double *addover,
   long long int *mul, double *mulover )
{
   const int nbthreads = nblocks*szt;
   int inc;

   // *add += nbthreads*(szt + 1); // dim equals szt, one extract subtraction
   inc = nbthreads*(szt + 1);
   update_counters(add,addover,inc);
   // *mul += nbthreads*szt;
   inc = nbthreads*szt;
   update_counters(mul,mulover,inc);
}

void flopcount_cmplx_back_substitute
 ( int nblocks, int szt, long long int *add, long long int *mul )
{
   const int nbthreads = nblocks*szt;

   *add += 4*nbthreads*szt + 2*nbthreads;
   *mul += 4*nbthreads*szt;
}

void overflopcount_cmplx_back_substitute
 ( int nblocks, int szt,
   long long int *add, double *addover,
   long long int *mul, double *mulover )
{
   const int nbthreads = nblocks*szt;
   int inc;

   // *add += 4*nbthreads*szt + 2*nbthreads;
   inc = 4*nbthreads*szt + 2*nbthreads;
   update_counters(add,addover,inc);
   // *mul += 4*nbthreads*szt;
   inc = 4*nbthreads*szt;
   update_counters(mul,mulover,inc);
}
