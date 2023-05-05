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

int dbl_fabry_predict1
 ( int deg, double *input, double step, double *output, int vrblvl )
{
   double fail = 0;            // we do not expect failure ...
   double den = input[deg-1];
   double num = input[deg];
   double y = 0.0;

   if(den + 1.0 == 1.0)
   {
      fail = 1;                // report failure

      if(vrblvl > 0)
         cout << "Denominator " << den << " is too small!" << endl;
      // then just evaluate the series ...
      for(int i=deg; i>0; i--) y = y + input[i]*step;
      *output = y + input[0];
   }
   else
   {
      double b = -num/den;  // coefficient of t in denominator
      double a;             // coefficient of numerator

      for(int i=deg-1; i>0; i--)
      {
         a = input[i] + b*input[i-1];
         y = y + a*step;
      }
      y = y + input[0];
      *output = y/(1.0 + b*step);
   }
   return fail;
}

int dbl_fabry_predictor
 ( int dim, int deg, double **input, double step, double *output, int vrblvl )
{
   int fail = 0;

   for(int i=0; i<dim; i++)
   {
      fail += dbl_fabry_predict1(deg,input[i],step,&output[i],vrblvl);
      if(vrblvl > 1)
         cout << "prediction[" << i << "] : " << output[i] << endl;
   }
   if(vrblvl > 1) cout << "#failed predictions : " << fail << endl; 

   return fail;
}
