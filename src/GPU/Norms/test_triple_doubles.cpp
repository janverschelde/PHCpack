// Test on triple double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "triple_double_functions.h"
#include "dbl3_sqrt_kernels.h"

using namespace std;

int my_sqrt ( double *hi, double *mi, double *lo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the triple double
 *   given by high, middle and low parts respectively in hi, mi, and lo,
 *   in as many interations as the value of max_steps.
 *   Returns in hi, mi, and lo the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 7;
   double twohi_h = 2.0;
   double twomi_h = 0.0;
   double twolo_h = 0.0;
   double twohi_d = 2.0;
   double twomi_d = 0.0;
   double twolo_d = 0.0;

   int fail = my_sqrt(&twohi_h,&twomi_h,&twolo_h,max);

   cout << "Test on triple doubles ";
   if(fail != 0)
      cout << "failed!" << endl;
   else
   {
      cout << "passed." << endl;

      GPU_dbl3_sqrt(&twohi_d,&twomi_d,&twolo_d,max);

      cout << "GPU computed sqrt :" << endl;
      tdf_write_doubles(twohi_d,twomi_d,twolo_d);

      double err = abs(twohi_h - twohi_d)
                 + abs(twomi_h - twomi_d)
                 + abs(twolo_h - twolo_d);
      cout << "  error : " << err << endl; 

      if(err < 1.0e-44)
         cout << "GPU test on triple doubles passed." << endl;
      else
         cout << "GPU test on triple doubles failed!" << endl;
   }
   return 0;
}

int my_sqrt ( double *hi, double *mi, double *lo, int max_steps )
{
   const double tol = 1.0e-44;

   double x_hi,x_mi,x_lo,z_hi,z_mi,z_lo,y_hi,y_mi,y_lo,e_hi,e_mi,e_lo;
   double a_hi,a_mi,a_lo;
   int i;

   x_hi = *hi; x_mi = *mi; x_lo = *lo;

   cout << "\nRunning Newton's method for ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      tdf_copy(x_hi,x_mi,x_lo,&y_hi,&y_mi,&y_lo); // copy for comparison
      tdf_sqr(x_hi,x_mi,x_lo,&z_hi,&z_mi,&z_lo);  // z = x*x
      cout << "x*x : " << endl;
      tdf_write_doubles(z_hi,z_mi,z_lo); cout << endl;
      tdf_inc(&z_hi,&z_mi,&z_lo,*hi,*mi,*lo);     // z += input number
      tdf_div(z_hi,z_mi,z_lo,x_hi,x_mi,x_lo,&z_hi,&z_mi,&z_lo); // z = z/x
      tdf_mul_td_d(z_hi,z_mi,z_lo,0.5,&z_hi,&z_mi,&z_lo);       // z *= 0.5
      tdf_copy(z_hi,z_mi,z_lo,&x_hi,&x_mi,&x_lo);
      cout << "after step " << i << " : " << endl;
      tdf_write_doubles(x_hi,x_mi,x_lo);
      tdf_sub(x_hi,x_mi,x_lo,y_hi,y_mi,y_lo,&e_hi,&e_mi,&e_lo); // error
      tdf_abs(e_hi,e_mi,e_lo,&a_hi,&a_mi,&a_lo);
      cout << "  error : "<< a_hi << endl;
   }

   *hi = x_hi; *mi = x_mi; *lo = x_lo;

   return int(a_hi + a_mi + a_lo > tol);
}
