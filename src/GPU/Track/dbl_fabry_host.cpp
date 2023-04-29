// The file dbl_fabry_host.cpp defines the functions specified in
// the file dbl_fabry_host.h.

#include <iostream>
#include "dbl_fabry_host.h"

using namespace std;

int dbl_fabry_ratio
 ( int deg, double *input, double *ratio, int vrblvl )
{
   double den = input[deg];
   double num = input[deg-1];

   *ratio = 0.0;

   if(den + 1.0 != 1.0)
   {
      *ratio = num/den;
      if(*ratio < 0.0) *ratio = -(*ratio);
   }
   else
      if(vrblvl > 0)
         cout << "Denominator " << den << " is too small!" << endl;

   return 0;
}

int dbl_fabry_vector
 ( int dim, int deg, double **input, double *ratios, int vrblvl )
{
   int fail;

   for(int i=0; i<dim; i++)
   {
      fail = dbl_fabry_ratio(deg,input[i],&ratios[i],vrblvl);
      if(vrblvl > 0)
         cout << "ratio[" << i << "] : " << ratios[i] << endl;
   }
   return fail;
}

int dbl_fabry_smallest
 ( int dim, double *ratios, double *step, int vrblvl )
{
   *step = ratios[0];

   for(int i=1; i<dim; i++)
      if(ratios[i] < *step) *step = ratios[i];

   if(vrblvl > 0) cout << "step : " << *step << endl;

   return 0;
}

int dbl_fabry_step
 ( int dim, int deg, double **input, double *ratios, double *step,
   int vrblvl )
{
   int fail;

   fail = dbl_fabry_vector(dim,deg,input,ratios,vrblvl);

   fail = dbl_fabry_smallest(dim,ratios,step,vrblvl);

   return fail;
}
