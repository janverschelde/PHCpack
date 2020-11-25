// Test on double double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "double_double.h"
#include "double_double_functions.h"
#include "dbl2_sqrt_kernels.h"

using namespace std;

int my_sqrt ( double *hi, double *lo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root
 *   of the double double with high part in hi and low part in lo,
 *   in as many iterations as the value of max_steps.
 *   Returns in hi and lo the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 7;
   double twohi_h = 2.0;
   double twolo_h = 0.0;
   double twohi_d = 2.0;
   double twolo_d = 0.0;

   cout << scientific << setprecision(3);

   int fail = my_sqrt(&twohi_h,&twolo_h,max);

   if(fail != 0)
      cout << "Test failed!" << endl;
   else
   {
      cout << "Test passed." << endl;

      GPU_dbl2_sqrt(&twohi_d,&twolo_d,max);

      cout << "GPU sqrt : ";

      double result[2];
      result[0] = twohi_d;
      result[1] = twolo_d;
      dd_write(result,32);

      double err = abs(twohi_h - twohi_d)
                 + abs(twolo_h - twolo_d);
      cout << "  error : " << err << endl; 

      if(err < 1.0e-28)
         cout << "GPU test on double doubles passed." << endl;
      else
         cout << "GPU test on double doubles failed!" << endl;
   }
   return 0;
}

int my_sqrt ( double *hi, double *lo, int max_steps )
{
   const double tol = 1.0e-28;

   double n[2],x[2],y[2],z[2],e[2],a[2];
   double x_hi,x_lo,z_hi,z_lo,y_hi,y_lo,e_hi,e_lo;

   int i;

   x_hi = *hi; x_lo = *lo; n[0] = x_hi; n[1] = x_lo; // dd_copy(n,x);

   cout << "\nRunning Newton's method for sqrt ...\n";

   for(i=1; i <= max_steps; i++)
   {
      y_hi = x_hi; y_lo = x_lo;
      ddf_sqr(x_hi,x_lo,&z_hi,&z_lo);           // z = x*x
      ddf_inc(&z_hi,&z_lo,n[0],n[1]);           // z += n
      ddf_div(z_hi,z_lo,x_hi,x_lo,&z_hi,&z_lo); // z = z/x
      ddf_mlt_d(&z_hi,&z_lo,0.5);               // z *= 0.5
      z[0] = z_hi; z[1] = z_lo;
      cout << "  step " << i << " : "; dd_write(z,32);
      x_hi = z_hi; x_lo = z_lo;                 // dd_copy(z,x);
      ddf_sub(x_hi,x_lo,y_hi,y_lo,&e_hi,&e_lo); // dd_sub(x,y,e);
      ddf_abs(e_hi,e_lo,&a[0],&a[1]);
      cout << "  error : " << a[0] << endl;
   }
   *hi = x_hi; *lo = x_lo;

   return int(a[0] + a[1] > tol);
}
